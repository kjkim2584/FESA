#include "analysis/solid_dynamic_reaction.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>

#include "adapters/openmp/thread_runtime.h"
#include "analysis/linear_dynamic_element.h"
#include "analysis/numerical_support.h"
#include "model/model_storage.h"

namespace fesa {
namespace {
using numerical_support::Buffers;
struct ReactionTask {
  const HhtSolidForceInput* hht = nullptr;  // Borrowed through join.
  Vector* origins = nullptr;  // Borrowed disjoint per-element output ranges.
  int origin_offset = 0;
  const int* origin_offsets = nullptr; ///< Borrowed prefix through joined work.
  Buffers* buffers = nullptr;  // Borrowed exclusive scratch through join.
  const Vector* displacement = nullptr;  // Borrowed global physical motion.
  const Vector* velocity = nullptr;
  const Vector* acceleration = nullptr;
  MassMode mode = MassMode::kConsistent;
};
StatusCode ComputeContribution(int item, int worker, void* context,
                               Diagnostic* d) {
  auto* task = static_cast<ReactionTask*>(context);
  auto* b = task->buffers;
  auto* local = b->workers[worker];
  const auto* element = b->elements[item];
  const int size = element->dof_count();
  std::unique_ptr<numerical_support::WorkerBuffers> exact;
  if (local->u->Size() != size) {
    numerical_support::WorkerBuffers* raw = nullptr;
    const auto status = numerical_support::WorkerBuffers::Create(
        local->element, &raw, d, size);
    exact.reset(raw);
    if (status != StatusCode::kOk) return status;
    local = raw;
  }
  MassMode physical;
  auto status = ResolveDynamicMassMode(element, task->mode, &physical, d);
  if (status != StatusCode::kOk) return status;
  d->element_id = element->Id();
  const int* indices = b->global_indices + b->element_dofs * item;
  status = local->u->GatherFrom(task->displacement, indices, size, 0, d);
  if (status != StatusCode::kOk) return status;
  status =
      element->ComputeInternalForce(local->u, local->element, local->force, d);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<Vector> velocity, acceleration, inertia, damping;
  for (auto* owner : {&velocity, &acceleration, &inertia, &damping}) {
    Vector* raw = nullptr;
    status = Vector::Create(size, &raw, d);
    owner->reset(raw);
    if (status != StatusCode::kOk) return status;
  }
  status = velocity->GatherFrom(task->velocity, indices, size, 0, d);
  if (status != StatusCode::kOk) return status;
  status = acceleration->GatherFrom(task->acceleration, indices, size, 0, d);
  if (status != StatusCode::kOk) return status;
  ElementMotionInput motion;
  motion.linear_displacement = local->u;
  motion.velocity = velocity.get();
  motion.acceleration = acceleration.get();
  ElementDynamicResult result;
  result.inertia_force = inertia.get();
  result.damping_force = damping.get();
  TrialEvaluationFailure failure = TrialEvaluationFailure::kNone;
  status = element->EvaluateDynamic(&motion, physical, local->element,
                                    DynamicEvaluationRequest::kForces, &result,
                                    &failure, d);
  if (status != StatusCode::kOk) return status;
  int identity[Buffers::kMaxElementDofs] = {};
  for (int i = 0; i < size; ++i) identity[i] = i;
  if (task->hht) {
    const int offset = task->origin_offsets[task->origin_offset + item];
    status = task->origins->GatherFrom(local->force, identity, size, offset, d);
    if (status != StatusCode::kOk) return status;
    status =
        task->origins->GatherFrom(inertia.get(), identity, size, offset + size, d);
    if (status != StatusCode::kOk) return status;
    if (HasMaterialRayleighDamping(element)) {
      status = task->origins->GatherFrom(damping.get(), identity, size,
                                          offset + 2 * size, d);
      if (status != StatusCode::kOk) return status;
      status = local->force->AddScaled(1, damping.get(), d);
      if (status != StatusCode::kOk) return status;
    }
    status = local->force->Scale(-(1 + task->hht->alpha), d);
    if (status != StatusCode::kOk) return status;
    status = local->force->AddScaled(-1, inertia.get(), d);
    if (status != StatusCode::kOk) return status;
    status = local->u->GatherFrom(task->hht->previous_displacement, indices, size,
                                  0, d);
    if (status != StatusCode::kOk) return status;
    // Reuse scratch only after its physical contribution has been consumed.
    status = element->ComputeInternalForce(local->u, local->element,
                                           damping.get(), d);
    if (status != StatusCode::kOk) return status;
    status = local->force->AddScaled(task->hht->alpha, damping.get(), d);
    if (status != StatusCode::kOk) return status;
    if (HasMaterialRayleighDamping(element) && task->hht->alpha != 0) {
      status = velocity->GatherFrom(task->hht->previous_velocity, indices,
                                     size, 0, d);
      if (status == StatusCode::kOk)
        status = element->EvaluateDynamic(&motion, physical, local->element,
            DynamicEvaluationRequest::kForces, &result, &failure, d);
      if (status == StatusCode::kOk)
        status = local->force->AddScaled(task->hht->alpha, damping.get(), d);
      if (status != StatusCode::kOk) return status;
    }
  } else {
    status = local->force->AddScaled(1, inertia.get(), d);
    if (status != StatusCode::kOk) return status;
    status = local->force->AddScaled(1, damping.get(), d);
    if (status != StatusCode::kOk) return status;
  }
  return b->vector_values->GatherFrom(local->force, identity, size, b->vector_offsets[item], d);
}

StatusCode ComputeSolidDynamicForce(
    const AnalysisModel* model, const DofManager* dofs, MassMode mode,
    const Vector* displacement, const Vector* velocity,
    const Vector* acceleration, const Vector* external,
    ElementWorkspace** workers, int worker_count, Vector* reaction,
    const HhtSolidForceInput* hht, Vector** out_origins, Diagnostic* d) {
  if (!model || !dofs || worker_count <= 0 || !reaction ||
      (mode != MassMode::kConsistent && mode != MassMode::kRowSum &&
       mode != MassMode::kDefault) ||
      dofs->NodeCount() <= 0 || model->Elements()->empty() ||
      model->Nodes()->size() != static_cast<std::size_t>(dofs->NodeCount()) ||
      reaction->Size() != dofs->GlobalCount())
    return model_storage::Invalid("Solid dynamic reaction inputs", d);
  bool shell = false;
  auto status = ValidateLinearDynamicModel(model, mode, &shell, d);
  if (status == StatusCode::kOk)
    status = ValidateLinearDynamicWorkers(shell, workers, worker_count, d);
  if (status != StatusCode::kOk) return status;
  int maximum = 0;
  status = ValidateAssemblyCounts(model, dofs, d, &maximum, nullptr);
  if (status != StatusCode::kOk) return status;
  for (const auto* input : {displacement, velocity, acceleration, external}) {
    if (!input || input == reaction || input->Size() != dofs->GlobalCount())
      return model_storage::Invalid("Solid dynamic reaction vector", d);
    for (int i = 0; i < input->Size(); ++i)
      if (!std::isfinite(input->Data()[i]))
        return model_storage::Invalid("Solid dynamic reaction finite input", d);
  }
  int origin_count = 0;
  std::unique_ptr<int[]> origin_offsets;
  if (hht) {
    if (!std::isfinite(hht->alpha) || hht->alpha < -1. / 3 || hht->alpha > 0 ||
        !hht->previous_displacement ||
        hht->previous_displacement->Size() != dofs->GlobalCount())
      return model_storage::Invalid("HHT force history", d);
    for (const auto* input :
         {hht->previous_displacement, hht->previous_velocity,
          hht->external_origins}) {
      if (!input) continue;
      if (input == reaction)
        return model_storage::Invalid("HHT output alias", d);
      for (int i = 0; i < input->Size(); ++i)
        if (!std::isfinite(input->Data()[i]))
          return model_storage::Invalid("HHT finite history/origins", d);
    }
    const int external_count =
        hht->external_origins ? hht->external_origins->Size() : 0;
    origin_offsets.reset(new (std::nothrow) int[model->Elements()->size() + 1]);
    if (!origin_offsets)
      return model_storage::AllocationFailure("HHT origin prefix", d);
    origin_offsets[0] = 0;
    int index = 0;
    for (const auto& entry : *model->Elements()) {
      const bool damped = HasMaterialRayleighDamping(entry.second);
      if (damped && hht->alpha != 0 &&
          (!hht->previous_velocity ||
           hht->previous_velocity->Size() != dofs->GlobalCount()))
        return model_storage::Invalid("HHT damping previous velocity", d);
      const int width = entry.second->dof_count() * (damped ? 3 : 2);
      if (origin_offsets[index] >
          std::numeric_limits<int>::max() - external_count - width)
        return model_storage::Invalid("HHT flux contribution count", d);
      origin_offsets[index + 1] = origin_offsets[index] + width;
      ++index;
    }
    origin_count = origin_offsets[index] + external_count;
  }
  for (const auto& entry : *model->Nodes())
    if (dofs->NodeDofCount(entry.first) != (shell ? 6 : 3))
      return model_storage::Invalid("Solid dynamic reaction 3N mapping", d);
  for (int i = 0; i < dofs->GlobalCount(); ++i)
    if (dofs->Equation(i) == -1 &&
        (velocity->Data()[i] != 0 || acceleration->Data()[i] != 0))
      return model_storage::Invalid("Solid dynamic reaction fixed motion", d);
  status = ValidateAssemblyCounts(model, dofs, d);
  if (status != StatusCode::kOk) return status;
  Vector* raw_force = nullptr;
  status = Vector::Create(dofs->GlobalCount(), &raw_force, d);
  std::unique_ptr<Vector> force(raw_force);
  if (status != StatusCode::kOk) return status;
  status = force->CopyFrom(external, d);
  if (status != StatusCode::kOk) return status;
  status = force->Scale(hht ? 1 : -1, d);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<Vector> origins;
  if (hht) {
    Vector* raw_origins = nullptr;
    status = Vector::Create(origin_count, &raw_origins, d);
    origins.reset(raw_origins);
    if (status != StatusCode::kOk) return status;
    if (hht->external_origins) {
      const int offset = origin_offsets[model->Elements()->size()];
      for (int i = 0; i < hht->external_origins->Size();) {
        const int count = std::min(24, hht->external_origins->Size() - i);
        int indices[24] = {};
        for (int j = 0; j < count; ++j) indices[j] = i + j;
        status = origins->GatherFrom(hht->external_origins, indices, count,
                                     offset + i, d);
        if (status != StatusCode::kOk) return status;
        i += count;
      }
    }
  }
  const int capacity = static_cast<int>(std::min(
      model->Elements()->size(), static_cast<std::size_t>(Buffers::kElements)));
  Buffers* raw = nullptr;
  status = Buffers::Create(false, capacity, &raw, d, workers, worker_count, maximum,
                           false);
  std::unique_ptr<Buffers> buffers(raw);
  if (status != StatusCode::kOk) return status;
  ReactionTask task;
  task.hht = hht;
  task.origins = origins.get();
  task.origin_offsets = origin_offsets.get();
  task.buffers = raw;
  task.displacement = displacement;
  task.velocity = velocity;
  task.acceleration = acceleration;
  task.mode = mode;
  auto next = model->Elements()->begin();
  while (next != model->Elements()->end()) {
    int count = 0;
    raw->vector_offsets[0] = 0;
    while (count < capacity && next != model->Elements()->end()) {
      raw->elements[count] = next->second;
      const int size = next->second->dof_count();
      int* indices = raw->global_indices + raw->element_dofs * count;
      status =
          numerical_support::ElementIndices(next->second, dofs, indices, d);
      if (status != StatusCode::kOk) return status;
      const int offset = raw->vector_offsets[count];
      for (int i = 0; i < size; ++i)
        raw->vector_indices[offset + i] = indices[i];
      raw->vector_offsets[count + 1] = offset + size;
      ++count;
      ++next;
    }
    status = ExecuteElementTasks(count, worker_count, ComputeContribution,
                                 &task, raw->diagnostics, nullptr, d);
    if (status != StatusCode::kOk) return status;
    status = force->ScatterAdd(raw->vector_values, raw->vector_indices,
                               raw->vector_offsets[count], d);
    if (status != StatusCode::kOk) return status;
    if (hht) task.origin_offset += count;
  }
  for (int i = 0; !hht && i < dofs->GlobalCount(); ++i) {
    if (dofs->Equation(i) == -1) continue;
    status = force->Set(i, 0, d);
    if (status != StatusCode::kOk) return status;
  }
  status = reaction->CopyFrom(force.get(), d);
  if (status == StatusCode::kOk && hht) *out_origins = origins.release();
  return status;
}
}  // namespace

StatusCode ComputeSolidDynamicReaction(
    const AnalysisModel* model, const DofManager* dofs, MassMode mode,
    const Vector* displacement, const Vector* velocity,
    const Vector* acceleration, const Vector* external,
    ElementWorkspace** workers, int worker_count, Vector* reaction,
    Diagnostic* d) {
  return ComputeSolidDynamicForce(model, dofs, mode, displacement, velocity,
                                  acceleration, external, workers, worker_count,
                                  reaction, nullptr, nullptr, d);
}
StatusCode ComputeSolidHhtResidual(
    const AnalysisModel* model, const DofManager* dofs, MassMode mode,
    const Vector* displacement, const Vector* velocity,
    const Vector* acceleration, const Vector* weighted_external,
    const HhtSolidForceInput* hht, ElementWorkspace** workers, int worker_count,
    Vector* residual, Vector** origins, Diagnostic* d) {
  if (!origins) return model_storage::Invalid("HHT origins output", d);
  *origins = nullptr;
  if (!hht) return model_storage::Invalid("HHT force input", d);
  return ComputeSolidDynamicForce(model, dofs, mode, displacement, velocity,
                                  acceleration, weighted_external, workers,
                                  worker_count, residual, hht, origins, d);
}
}  // namespace fesa
