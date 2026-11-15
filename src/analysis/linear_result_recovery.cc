#include "analysis/linear_result_recovery.h"

#include <algorithm>
#include <memory>

#include "adapters/openmp/thread_runtime.h"
#include "analysis/numerical_support.h"
#include "model/model_storage.h"

namespace fesa {
namespace {
struct RecoveryTaskContext {
  numerical_support::Buffers* buffers = nullptr;
  const Vector* candidate = nullptr;
  const BeamLineLoadState* line_loads = nullptr;
};
StatusCode RecoverElement(int item, int worker, void* context, Diagnostic* d) {
  auto* task = static_cast<RecoveryTaskContext*>(context);
  auto* buffers = task->buffers;
  auto* local = buffers->workers[worker];
  Element* element = buffers->elements[item];
  d->element_id = element->Id();
  const int size = element->dof_count();
  std::unique_ptr<numerical_support::WorkerBuffers> owned;
  if (local->u->Size() != size) {
    numerical_support::WorkerBuffers* raw = nullptr;
    const auto status =
        numerical_support::WorkerBuffers::Create(local->element, &raw, d, size);
    owned.reset(raw);
    if (status != StatusCode::kOk) return status;
    local = raw;
  }
  auto status = local->u->GatherFrom(
      task->candidate, buffers->global_indices + buffers->element_dofs * item,
      size, 0, d);
  if (status != StatusCode::kOk) return status;
  if (task->line_loads) {
    status = numerical_support::EvaluateLineLoad(task->line_loads, element,
                                                 local, d);
    if (status != StatusCode::kOk) return status;
    status = element->UpdateLoadedResults(local->u, local->line_local,
                                          local->element, d);
  } else {
    status = element->UpdateResults(local->u, local->element, d);
  }
  return status;
}
}  // namespace
StatusCode RecoverLinearResults(const AnalysisModel* model,
                                const DofManager* dofs, AnalysisState* state,
                                int step, int increment,
                                ElementWorkspace** workers, int worker_count,
                                Diagnostic* d,
                                void (*observer)(AnalysisStage)) {
  if (!model || !dofs || !state || !model->Conditions() || step <= 0 ||
      increment < 0 || worker_count <= 0 || !state->HasValidCandidate() ||
      !state->CandidateDisplacement() ||
      (increment == 0 &&
       (!state->Dynamics() || !state->Dynamics()->CandidateVelocity())))
    return model_storage::Invalid("Linear recovery inputs/candidate", d);
  // Legacy linear states have no stored node IDs. The caller validates the
  // model mapping; recovery retains their existing physical-vector contract.
  if (state->GlobalCount() != dofs->GlobalCount())
    return model_storage::Invalid("Linear recovery physical shape", d);
  Vector* raw = nullptr;
  auto status = Vector::Create(12, &raw, d);
  std::unique_ptr<Vector> values(raw);
  if (status != StatusCode::kOk) return status;
  bool started = false;
  for (const auto& entry : *model->Nodes()) {
    if (d) {
      d->node_id = entry.first;
      d->element_id = 0;
    }
    const int components = dofs->NodeDofCount(entry.first);
    if (values->Size() != 2 * components) {
      Vector* resized = nullptr;
      status = Vector::Create(2 * components, &resized, d);
      if (status != StatusCode::kOk) return status;
      values.reset(resized);
    }
    int indices[6] = {};
    for (int i = 0; i < components; ++i)
      indices[i] = dofs->GlobalDof(entry.first, i + 1);
    status = values->GatherFrom(state->CandidateDisplacement(), indices,
                                components, 0, d);
    if (status != StatusCode::kOk) return status;
    status = values->GatherFrom(state->CandidateReaction(), indices, components,
                                components, d);
    if (status != StatusCode::kOk) return status;
    if (!started) {
      state->MarkResultsUnconfirmed(step, increment);
      started = true;
    }
    if (observer) observer(AnalysisStage::kNodeUpdate);
    status = entry.second->UpdateResults(values.get(), d);
    if (status != StatusCode::kOk) return status;
  }
  if (model->Elements()->empty()) return StatusCode::kOk;
  using numerical_support::Buffers;
  const int capacity = static_cast<int>(std::min(
      model->Elements()->size(), static_cast<std::size_t>(Buffers::kElements)));
  Buffers* buffers_raw = nullptr;
  int maximum_dofs = 12;
  bool beam_only = true;
  for (const auto& entry : *model->Elements()) {
    maximum_dofs = std::max(maximum_dofs, entry.second->dof_count());
    beam_only = beam_only && entry.second->family() == ElementFamily::kBeam;
  }
  status = Buffers::Create(false, capacity, &buffers_raw, d, workers,
                           worker_count, maximum_dofs, beam_only);
  std::unique_ptr<Buffers> buffers(buffers_raw);
  if (status != StatusCode::kOk) return status;
  const auto* line_loads = model->Conditions()->LineLoadState();
  status = numerical_support::PrepareLineLoads(line_loads, buffers_raw, d);
  if (status != StatusCode::kOk) return status;
  RecoveryTaskContext task;
  task.buffers = buffers_raw;
  task.candidate = state->CandidateDisplacement();
  task.line_loads = line_loads;
  auto next = model->Elements()->begin();
  while (next != model->Elements()->end()) {
    int count = 0;
    while (count < capacity && next != model->Elements()->end()) {
      buffers->elements[count] = next->second;
      status = numerical_support::ElementIndices(
          next->second, dofs,
          buffers->global_indices + buffers->element_dofs * count, d);
      if (status != StatusCode::kOk) {
        if (d) d->element_id = next->first;
        return status;
      }
      ++count;
      ++next;
    }
    status = ExecuteElementTasks(count, worker_count, RecoverElement, &task,
                                 buffers->diagnostics, nullptr, d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}
}  // namespace fesa
