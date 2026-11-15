#include "analysis/dynamic_initialization.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <new>
#include <utility>

#include "adapters/openmp/thread_runtime.h"
#include "analysis/finite_shell_initialization.h"
#include "analysis/linear_dynamic_element.h"
#include "analysis/linear_shell_nullspace.h"
#include "analysis/numerical_support.h"
#include "model/model_storage.h"

namespace fesa {
StatusCode ValidateDynamicConstraintTransition(
    const EffectiveConditions* previous, const EffectiveConditions* current,
    Diagnostic* d) {
  if (!current) return model_storage::Invalid("Dynamic constraint context", d);
  if (!previous)
    return current->StepNumber() == 1
               ? StatusCode::kOk
               : model_storage::Invalid("Missing previous dynamic conditions",
                                        d);
  if (previous->SourceDomain() != current->SourceDomain() ||
      previous->StepNumber() != current->StepNumber() - 1)
    return model_storage::Invalid("Dynamic condition history", d);
  std::map<std::pair<int, int>, double> old;
  for (const auto* boundary : *previous->Boundaries())
    old.emplace(std::make_pair(boundary->node_id, boundary->definition->Dof()),
                boundary->definition->Value());
  bool same = old.size() == current->Boundaries()->size();
  for (const auto* boundary : *current->Boundaries()) {
    const auto found = old.find(
        std::make_pair(boundary->node_id, boundary->definition->Dof()));
    same = same && found != old.end() &&
           found->second == boundary->definition->Value();
  }
  return same ? StatusCode::kOk
              : ReportFailure(StatusCode::kUnsupported,
                              DiagnosticStage::kValidation,
                              "Dynamic constraints must remain unchanged", d);
}

namespace {
using numerical_support::Buffers;

StatusCode ValidateModel(const AnalysisModel* model, const DofManager* dofs,
                         MassMode mode, LinearSolver* solver,
                         ElementWorkspace** workers, int worker_count,
                         Diagnostic* d) {
  if (!model || !dofs || !solver || worker_count <= 0 ||
      (mode != MassMode::kConsistent && mode != MassMode::kRowSum &&
       mode != MassMode::kDefault) ||
      dofs->NodeCount() <= 0 || model->Elements()->empty() ||
      model->Nodes()->size() != static_cast<std::size_t>(dofs->NodeCount()))
    return model_storage::Invalid("Initial acceleration inputs", d);
  for (const auto& entry : *model->Nodes())
    if (dofs->NodeDofCount(entry.first) != 3)
      return model_storage::Invalid("Initial acceleration 3N mapping", d);
  for (const auto& entry : *model->Elements()) {
    const auto* element = entry.second;
    if (element->family() != ElementFamily::kSolid ||
        !element->SupportsDynamicEvaluation(DynamicKinematics::kLinear)) {
      const auto status = ReportFailure(
          StatusCode::kUnsupported, DiagnosticStage::kValidation,
          "Initial acceleration requires supported linear Solid", d);
      if (d) d->element_id = entry.first;
      return status;
    }
    MassMode physical;
    const auto resolved = ResolveDynamicMassMode(element, mode, &physical, d);
    if (resolved != StatusCode::kOk) return resolved;
  }
  auto status = ValidateAssemblyCounts(model, dofs, d);
  if (status != StatusCode::kOk || dofs->FreeCount() == 0) return status;
  if (!workers) return model_storage::Invalid("Initial mass workers", d);
  for (int i = 0; i < worker_count; ++i) {
    if (!workers[i]) return model_storage::Invalid("Initial mass worker", d);
    for (int j = 0; j < i; ++j)
      if (workers[i] == workers[j])
        return model_storage::Invalid("Aliased initial mass workers", d);
  }
  return StatusCode::kOk;
}

StatusCode Validate(const AnalysisModel* model, const DofManager* dofs,
                    MassMode mode, const Vector* force, LinearSolver* solver,
                    ElementWorkspace** workers, int worker_count,
                    const Vector* acceleration, Diagnostic* d) {
  auto status =
      ValidateModel(model, dofs, mode, solver, workers, worker_count, d);
  if (status != StatusCode::kOk) return status;
  if (!force || !acceleration || force == acceleration ||
      force->Size() != dofs->GlobalCount() ||
      acceleration->Size() != dofs->GlobalCount())
    return model_storage::Invalid("Initial acceleration vectors", d);
  for (int i = 0; i < force->Size(); ++i)
    if (!std::isfinite(force->Data()[i]))
      return model_storage::Invalid("Initial acceleration physical force", d);
  return StatusCode::kOk;
}

class InitialCandidateGuard {
 public:
  explicit InitialCandidateGuard(AnalysisState* state) : state_(state) {}
  InitialCandidateGuard(const InitialCandidateGuard&) = delete;
  InitialCandidateGuard& operator=(const InitialCandidateGuard&) = delete;
  InitialCandidateGuard(InitialCandidateGuard&&) = delete;
  InitialCandidateGuard& operator=(InitialCandidateGuard&&) = delete;
  ~InitialCandidateGuard() {
    if (state_ && !ready_) state_->DiscardCandidate();
  }
  void KeepCandidate() { ready_ = true; }

 private:
  bool ready_ = false;
  AnalysisState* state_;  // Borrowed for the initialization call only.
};

struct InitialForceTask {
  Buffers* buffers = nullptr;  // Borrowed exclusive scratch through join.
  const Vector* displacement = nullptr;  // Borrowed physical initial u.
  const Vector* velocity = nullptr;  // Borrowed initial physical velocity.
  MassMode mode = MassMode::kConsistent;
};

StatusCode ComputeInitialForce(int item, int worker, void* context,
                               Diagnostic* d) {
  auto* task = static_cast<InitialForceTask*>(context);
  auto* buffers = task->buffers;
  auto* local = buffers->workers[worker];
  const auto* element = buffers->elements[item];
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
  status = local->u->GatherFrom(
      task->displacement, buffers->global_indices + buffers->element_dofs * item, size, 0, d);
  if (status != StatusCode::kOk) return status;
  status =
      element->ComputeInternalForce(local->u, local->element, local->force, d);
  if (status != StatusCode::kOk) return status;
  if (HasMaterialRayleighDamping(element)) {
    std::unique_ptr<Vector> fields[4];
    for (auto& field : fields) {
      Vector* raw = nullptr;
      status = Vector::Create(size, &raw, d);
      field.reset(raw);
      if (status != StatusCode::kOk) return status;
    }
    status = fields[0]->GatherFrom(task->velocity,
        buffers->global_indices + buffers->element_dofs * item, size, 0, d);
    ElementMotionInput input;
    input.linear_displacement = local->u;
    input.velocity = fields[0].get();
    input.acceleration = fields[1].get();
    ElementDynamicResult result;
    result.inertia_force = fields[2].get();
    result.damping_force = fields[3].get();
    TrialEvaluationFailure failure = TrialEvaluationFailure::kNone;
    if (status == StatusCode::kOk)
      status = element->EvaluateDynamic(&input, physical, local->element,
          DynamicEvaluationRequest::kForces, &result, &failure, d);
    if (status == StatusCode::kOk)
      status = local->force->AddScaled(1, fields[3].get(), d);
    if (status != StatusCode::kOk) return status;
  }
  status = local->force->Scale(-1, d);
  if (status != StatusCode::kOk) return status;
  int indices[Buffers::kMaxElementDofs] = {};
  for (int i = 0; i < size; ++i) indices[i] = i;
  return buffers->vector_values->GatherFrom(local->force, indices, size,
                                            buffers->vector_offsets[item], d);
}

StatusCode InitialNetForce(const AnalysisModel* model, const DofManager* dofs,
                           const Vector* displacement, const Vector* external,
                           ElementWorkspace** workers, int worker_count,
                           Vector* net_force, Diagnostic* d,
                           const Vector* velocity, MassMode mode) {
  int maximum = 0;
  auto status = ValidateAssemblyCounts(model, dofs, d, &maximum, nullptr);
  if (status != StatusCode::kOk) return status;
  status = net_force->CopyFrom(external, d);
  if (status != StatusCode::kOk) return status;
  const int capacity = static_cast<int>(std::min(
      model->Elements()->size(), static_cast<std::size_t>(Buffers::kElements)));
  Buffers* raw = nullptr;
  status = Buffers::Create(false, capacity, &raw, d, workers, worker_count, maximum,
                           false);
  std::unique_ptr<Buffers> buffers(raw);
  if (status != StatusCode::kOk) return status;
  InitialForceTask task;
  task.buffers = raw;
  task.displacement = displacement;
  task.velocity = velocity;
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
    status = ExecuteElementTasks(count, worker_count, ComputeInitialForce,
                                 &task, raw->diagnostics, nullptr, d);
    if (status != StatusCode::kOk) return status;
    status = net_force->ScatterAdd(raw->vector_values, raw->vector_indices,
                                   raw->vector_offsets[count], d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}

// Borrowed only during finite initial mass assembly. net_force is caller-owned
// scratch, never the physical external input or an AnalysisState field.
struct FiniteMassInput {
  const ElementTrialState* pose = nullptr;
  const Vector* velocity = nullptr;
  Vector* net_force = nullptr;
};

StatusCode ComputeFiniteInitialMass(const Element* element, const int* global,
                                    const FiniteMassInput* input, MassMode mode,
                                    Matrix* mass, Vector* negative_force,
                                    Diagnostic* d) {
  Matrix *raw_q = nullptr, *raw_committed = nullptr;
  auto status = Matrix::Create(4, 4, &raw_q, d);
  std::unique_ptr<Matrix> q(raw_q);
  if (status == StatusCode::kOk)
    status = Matrix::Create(4, 4, &raw_committed, d);
  std::unique_ptr<Matrix> committed(raw_committed);
  std::unique_ptr<Vector> fields[6];
  const int sizes[] = {12, 12, 24, 24, 24, 24};
  for (int i = 0; i < 6 && status == StatusCode::kOk; ++i) {
    Vector* raw = nullptr;
    status = Vector::Create(sizes[i], &raw, d);
    fields[i].reset(raw);
  }
  // The common Element API requires all five distinct derivative outputs.
  std::unique_ptr<Matrix> derivatives[4];
  for (int i = 0; i < 4 && status == StatusCode::kOk; ++i) {
    Matrix* raw = nullptr;
    status = Matrix::Create(24, 24, &raw, d);
    derivatives[i].reset(raw);
  }
  if (status != StatusCode::kOk) return status;
  int translations[12] = {};
  for (int i = 0; i < 4; ++i) {
    const int row = global[6 * i] / 6;
    for (int c = 0; c < 3; ++c) translations[3 * i + c] = 3 * row + c;
    status = q->CopyBlockFrom(input->pose->orientations, row, 0, 1, 4, i, 0, d);
    if (status == StatusCode::kOk)
      status = committed->CopyBlockFrom(input->pose->committed_orientations,
                                        row, 0, 1, 4, i, 0, d);
    if (status != StatusCode::kOk) return status;
  }
  status =
      fields[0]->GatherFrom(input->pose->translations, translations, 12, 0, d);
  if (status == StatusCode::kOk)
    status =
        fields[1]->GatherFrom(input->pose->body_spin, translations, 12, 0, d);
  if (status == StatusCode::kOk)
    status = fields[2]->GatherFrom(input->velocity, global, 24, 0, d);
  if (status != StatusCode::kOk) return status;
  ElementTrialState local;
  local.translations = fields[0].get();
  local.body_spin = fields[1].get();
  local.orientations = q.get();
  local.committed_orientations = committed.get();
  ElementMotionInput motion;
  motion.kinematics = DynamicKinematics::kFiniteRotation;
  motion.finite_state = &local;
  motion.velocity = fields[2].get();
  motion.acceleration = fields[3].get();  // Zero: inertia force is gyro only.
  ElementDynamicResult result;
  result.inertia_force = fields[4].get();
  result.damping_force = fields[5].get();
  result.inertia_acceleration = mass;
  result.inertia_displacement = derivatives[0].get();
  result.inertia_velocity = derivatives[1].get();
  result.damping_displacement = derivatives[2].get();
  result.damping_velocity = derivatives[3].get();
  TrialEvaluationFailure failure = TrialEvaluationFailure::kNone;
  status = element->EvaluateDynamic(
      &motion, mode, nullptr, DynamicEvaluationRequest::kForcesAndLinearization,
      &result, &failure, d);
  if (status == StatusCode::kOk)
    status = element->EvaluateTrial(&local, nullptr,
                                    ElementEvaluationRequest::kInternalForce,
                                    negative_force, nullptr, &failure, d);
  if (status == StatusCode::kOk)
    status = negative_force->AddScaled(1, fields[4].get(), d);
  if (status == StatusCode::kOk && HasMaterialRayleighDamping(element))
    status = negative_force->AddScaled(1, fields[5].get(), d);
  return status == StatusCode::kOk ? negative_force->Scale(-1, d) : status;
}

struct MassTask {
  const DofManager* dofs = nullptr;  // Borrowed through joined dispatch.
  Buffers* buffers = nullptr;        // Exclusive scratch owned by coordinator.
  MassMode mode = MassMode::kConsistent;
  const FiniteMassInput* finite = nullptr;  // Optional, borrowed through join.
};

StatusCode PrepareMass(Element* element, const DofManager* dofs, Buffers* b,
                       int item, Diagnostic* d) {
  b->elements[item] = element;
  int* global = b->global_indices + b->element_dofs * item;
  auto status = numerical_support::ElementIndices(element, dofs, global, d);
  if (status != StatusCode::kOk) return status;
  int count = b->matrix_offsets[item];
  for (int i = 0; i < element->dof_count(); ++i) {
    const int row = dofs->Equation(global[i]);
    if (row < 0) continue;
    for (int j = 0; j < element->dof_count(); ++j) {
      const int column = dofs->Equation(global[j]);
      if (column < row) continue;
      b->rows[count] = row;
      b->columns[count++] = column;
    }
  }
  b->matrix_offsets[item + 1] = count;
  return StatusCode::kOk;
}

StatusCode ComputeMass(int item, int worker, void* context, Diagnostic* d) {
  auto* task = static_cast<MassTask*>(context);
  auto* b = task->buffers;
  const auto* element = b->elements[item];
  auto* local = b->workers[worker];
  const int local_size = element->dof_count();
  std::unique_ptr<numerical_support::WorkerBuffers> exact;
  if (local->u->Size() != local_size) {
    numerical_support::WorkerBuffers* raw = nullptr;
    const auto status = numerical_support::WorkerBuffers::Create(
        local->element, &raw, d, local_size);
    exact.reset(raw);
    if (status != StatusCode::kOk) return status;
    local = raw;
  }
  MassMode physical;
  const auto resolved = ResolveDynamicMassMode(element, task->mode, &physical, d);
  if (resolved != StatusCode::kOk) return resolved;
  d->element_id = element->Id();
  // This exclusive square scratch holds M, not K, for the entire dispatch.
  const int* global = b->global_indices + b->element_dofs * item;
  auto status = task->finite ? ComputeFiniteInitialMass(
                                   element, global, task->finite, physical,
                                   local->stiffness, local->force, d)
                             : element->ComputeMass(physical, local->element,
                                                    local->stiffness, d);
  if (status != StatusCode::kOk) return status;
  if (task->finite) {
    int indices[24] = {};
    for (int i = 0; i < 24; ++i) indices[i] = i;
    status =
        b->vector_values->GatherFrom(local->force, indices, 24, 24 * item, d);
    if (status != StatusCode::kOk) return status;
  }
  int entries[Buffers::kMaxUpperEntries] = {};
  int count = 0;
  const int size = element->dof_count();
  for (int i = 0; i < size; ++i) {
    const int row = task->dofs->Equation(global[i]);
    if (row < 0) continue;
    for (int j = 0; j < size; ++j)
      if (task->dofs->Equation(global[j]) >= row)
        entries[count++] = size * i + j;
  }
  return local->stiffness->GatherTo(entries, count, b->matrix_values,
                                    b->matrix_offsets[item], d);
}

StatusCode AssembleMass(const AnalysisModel* model, const DofManager* dofs,
                        MassMode mode, ElementWorkspace** workers,
                        int worker_count, int size, SparseMatrix** out,
                        Diagnostic* d,
                        const FiniteMassInput* finite = nullptr) {
  *out = nullptr;
  int maximum = 0;
  auto status = ValidateAssemblyCounts(model, dofs, d, &maximum, nullptr);
  if (status != StatusCode::kOk) return status;
  SparseMatrix* matrix = nullptr;
  status =
      SparseMatrix::Create(size, 0, nullptr, nullptr, nullptr, &matrix, d);
  std::unique_ptr<SparseMatrix> owner(matrix);
  if (status != StatusCode::kOk) return status;
  const int capacity = static_cast<int>(std::min(
      model->Elements()->size(), static_cast<std::size_t>(Buffers::kElements)));
  Buffers* raw = nullptr;
  status = Buffers::Create(true, capacity, &raw, d, workers, worker_count, maximum,
                           false);
  std::unique_ptr<Buffers> buffers(raw);
  if (status != StatusCode::kOk) return status;
  MassTask task;
  task.dofs = dofs;
  task.buffers = raw;
  task.mode = mode;
  task.finite = finite;
  auto next = model->Elements()->begin();
  while (next != model->Elements()->end()) {
    int count = 0;
    raw->matrix_offsets[0] = 0;
    while (count < capacity && next != model->Elements()->end()) {
      status = PrepareMass(next->second, dofs, raw, count, d);
      if (status != StatusCode::kOk) return status;
      if (finite)
        for (int i = 0; i < 24; ++i)
          raw->vector_indices[24 * count + i] =
              raw->global_indices[24 * count + i];
      ++count;
      ++next;
    }
    status = ExecuteElementTasks(count, worker_count, ComputeMass, &task,
                                 raw->diagnostics, nullptr, d);
    if (status != StatusCode::kOk) return status;
    if (finite) {
      status = finite->net_force->ScatterAdd(
          raw->vector_values, raw->vector_indices, 24 * count, d);
      if (status != StatusCode::kOk) return status;
    }
    SparseMatrix* chunk = nullptr;
    status = SparseMatrix::Create(size, raw->matrix_offsets[count], raw->rows,
                                  raw->columns, raw->matrix_values, &chunk, d);
    std::unique_ptr<SparseMatrix> chunk_owner(chunk);
    if (status != StatusCode::kOk) return status;
    SparseMatrix* merged = nullptr;
    status = SparseMatrix::Add(owner.get(), chunk, &merged, d);
    if (status != StatusCode::kOk) return status;
    owner.reset(merged);
  }
  *out = owner.release();
  return StatusCode::kOk;
}
}  // namespace

namespace {
StatusCode PrepareInitialDynamicState(
    bool shell, const Vector* external_rate, SolidDynamicStart start,
    const AnalysisModel* model, const DofManager* dofs,
    const EffectiveConditions* conditions, MassMode mode,
    const Vector* external, const Vector* initial_velocity,
    LinearSolver* solver, ElementWorkspace** workers, int worker_count,
    double total_time, AnalysisState* state, Diagnostic* d) {
  InitialCandidateGuard guard(state);
  if (!state || !conditions || !model || model->Conditions() != conditions ||
      state->HasValidCandidate() || !std::isfinite(total_time) ||
      total_time < 0 ||
      state->StepNumber() == std::numeric_limits<int>::max() ||
      conditions->StepNumber() != state->StepNumber() + 1)
    return model_storage::Invalid("Initial dynamic state inputs/step", d);
  const auto* dynamics = state->Dynamics();
  const bool has_motion = dynamics && dynamics->CommittedVelocity();
  if ((start != SolidDynamicStart::kDirect &&
       start != SolidDynamicStart::kStaticPreload &&
       start != SolidDynamicStart::kContinuous &&
       start != SolidDynamicStart::kLoadJump) ||
      (start == SolidDynamicStart::kDirect &&
       (state->StepNumber() != 0 || total_time != 0 || has_motion)) ||
      (start == SolidDynamicStart::kStaticPreload &&
       (state->StepNumber() <= 0 || has_motion || total_time <= 0)) ||
      ((start == SolidDynamicStart::kContinuous ||
        start == SolidDynamicStart::kLoadJump) &&
       (!has_motion || total_time != dynamics->TotalTime())) ||
      (start != SolidDynamicStart::kDirect && initial_velocity))
    return model_storage::Invalid("Initial dynamic start mode", d);
  auto status = state->ValidateMapping(dofs, d);
  if (status != StatusCode::kOk) return status;
  LinearShellNullspace* raw_space = nullptr;
  if (shell) {
    if (!solver || worker_count <= 0 ||
        (mode != MassMode::kConsistent && mode != MassMode::kRowSum))
      return model_storage::Invalid("Initial shell model inputs", d);
    status = LinearShellNullspace::Create(model, dofs, &raw_space, d);
  } else {
    status = ValidateModel(model, dofs, mode, solver, workers, worker_count, d);
  }
  std::unique_ptr<LinearShellNullspace> space(raw_space);
  if (status != StatusCode::kOk) return status;
  if (start != SolidDynamicStart::kContinuous) {
    if (!shell && !workers)
      return model_storage::Invalid("Initial structural workers", d);
    for (int i = 0; !shell && i < worker_count; ++i) {
      if (!workers[i])
        return model_storage::Invalid("Initial structural worker", d);
      for (int j = 0; j < i; ++j)
        if (workers[i] == workers[j])
          return model_storage::Invalid("Aliased initial structural workers",
                                        d);
    }
    if (!external || (shell && !external_rate))
      return model_storage::Invalid("Missing initial external force", d);
    for (const auto* vector :
         {external, initial_velocity, shell ? external_rate : nullptr}) {
      if (!vector) continue;
      if (vector->Size() != state->GlobalCount())
        return model_storage::Invalid("Initial physical input shape", d);
      for (int i = 0; i < vector->Size(); ++i)
        if (!std::isfinite(vector->Data()[i]))
          return model_storage::Invalid("Initial physical input finite values",
                                        d);
    }
  }
  for (const auto* boundary : *conditions->Boundaries()) {
    const int index =
        dofs->GlobalDof(boundary->node_id, boundary->definition->Dof());
    if (index < 0 || dofs->Equation(index) != -1)
      return model_storage::Invalid("Initial boundary mapping", d);
    if ((start == SolidDynamicStart::kContinuous ||
         start == SolidDynamicStart::kLoadJump) &&
        (dynamics->CommittedVelocity()->Data()[index] != 0 ||
         dynamics->CommittedAcceleration()->Data()[index] != 0))
      return ReportFailure(StatusCode::kUnsupported,
                           DiagnosticStage::kValidation,
                           "Nonzero committed motion on fixed DOF", d);
    if (initial_velocity && initial_velocity->Data()[index] != 0) {
      status =
          ReportFailure(StatusCode::kInputError, DiagnosticStage::kValidation,
                        "Nonzero initial velocity on fixed DOF", d);
      if (d) {
        d->node_id = boundary->node_id;
        d->condition_dof = boundary->definition->Dof();
      }
      return status;
    }
    if (start != SolidDynamicStart::kDirect &&
        state->CommittedDisplacement()->Data()[index] !=
            boundary->definition->Value())
      return ReportFailure(
          StatusCode::kUnsupported, DiagnosticStage::kValidation,
          "Changed prescribed displacement at dynamic start", d);
  }
  status = shell ? state->EnableLinearDynamics(d)
                 : state->EnableTranslationalDynamics(d);
  if (status != StatusCode::kOk) return status;
  status = state->PrepareCandidate(d);
  if (status != StatusCode::kOk) return status;
  for (const auto* boundary : *conditions->Boundaries()) {
    const int index =
        dofs->GlobalDof(boundary->node_id, boundary->definition->Dof());
    status = state->CandidateDisplacement()->Set(
        index, boundary->definition->Value(), d);
    if (status != StatusCode::kOk) return status;
  }
  if (start == SolidDynamicStart::kContinuous) {
    status = state->StageDynamicKinematics(dynamics->CommittedVelocity(),
                                           dynamics->CommittedAcceleration(), 0,
                                           total_time, d);
  } else {
    std::unique_ptr<Vector> velocity, acceleration, net_force;
    for (auto* owner : {&velocity, &acceleration, &net_force}) {
      Vector* raw = nullptr;
      status = Vector::Create(state->GlobalCount(), &raw, d);
      owner->reset(raw);
      if (status != StatusCode::kOk) return status;
    }
    if (initial_velocity) {
      status = velocity->CopyFrom(initial_velocity, d);
      if (status != StatusCode::kOk) return status;
    }
    if (start == SolidDynamicStart::kLoadJump) {
      status = velocity->CopyFrom(dynamics->CommittedVelocity(), d);
      if (status != StatusCode::kOk) return status;
    }
    if (shell) {
      status = space->CheckCompatibility(state->CandidateDisplacement(),
                                         external, false, d);
      if (status == StatusCode::kOk)
        status =
            space->CheckCompatibility(velocity.get(), external_rate, true, d);
      if (status != StatusCode::kOk) return status;
    }
    status =
        InitialNetForce(model, dofs, state->CandidateDisplacement(), external,
                        workers, worker_count, net_force.get(), d,
                        velocity.get(), mode);
    if (status != StatusCode::kOk) return status;
    if (shell) {
      Vector* raw_gamma = nullptr;
      if (space->Count() > 0)
        status = Vector::Create(space->Count(), &raw_gamma, d);
      std::unique_ptr<Vector> gamma(raw_gamma);
      // N4: constant basis, Kzr=0 and right p_ddot=0 for approved CLOAD.
      if (status == StatusCode::kOk)
        status = ComputeInitialShellAcceleration(
            model, dofs, mode, net_force.get(), raw_gamma, solver, worker_count,
            acceleration.get(), d);
    } else {
      status = ComputeInitialSolidAcceleration(
          model, dofs, mode, net_force.get(), solver, workers, worker_count,
          acceleration.get(), d);
    }
    if (status != StatusCode::kOk) return status;
    status = state->StageDynamicKinematics(velocity.get(), acceleration.get(),
                                           0, total_time, d);
  }
  if (status == StatusCode::kOk) guard.KeepCandidate();
  return status;
}

}  // namespace
StatusCode PrepareInitialSolidDynamicState(
    SolidDynamicStart start, const AnalysisModel* model, const DofManager* dofs,
    const EffectiveConditions* conditions, MassMode mode,
    const Vector* external, const Vector* initial_velocity,
    LinearSolver* solver, ElementWorkspace** workers, int worker_count,
    double total_time, AnalysisState* state, Diagnostic* d) {
  return PrepareInitialDynamicState(
      false, nullptr, start, model, dofs, conditions, mode, external,
      initial_velocity, solver, workers, worker_count, total_time, state, d);
}
StatusCode PrepareInitialShellDynamicState(
    SolidDynamicStart start, const AnalysisModel* model, const DofManager* dofs,
    const EffectiveConditions* conditions, MassMode mode,
    const Vector* external, const Vector* external_rate,
    const Vector* initial_velocity, LinearSolver* solver, int worker_count,
    double total_time, AnalysisState* state, Diagnostic* d) {
  return PrepareInitialDynamicState(
      true, external_rate, start, model, dofs, conditions, mode, external,
      initial_velocity, solver, nullptr, worker_count, total_time, state, d);
}

StatusCode ComputeInitialSolidAcceleration(
    const AnalysisModel* model, const DofManager* dofs, MassMode mode,
    const Vector* net_force, LinearSolver* solver, ElementWorkspace** workers,
    int worker_count, Vector* acceleration, Diagnostic* d) {
  auto status = Validate(model, dofs, mode, net_force, solver, workers,
                         worker_count, acceleration, d);
  if (status != StatusCode::kOk) return status;
  Vector* full = nullptr;
  status = Vector::Create(dofs->GlobalCount(), &full, d);
  std::unique_ptr<Vector> full_owner(full);
  if (status != StatusCode::kOk) return status;
  const int size = dofs->FreeCount();
  if (size == 0) return acceleration->CopyFrom(full, d);
  SparseMatrix* mass = nullptr;
  status =
      AssembleMass(model, dofs, mode, workers, worker_count, size, &mass, d);
  std::unique_ptr<SparseMatrix> mass_owner(mass);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<int[]> indices(new (std::nothrow) int[size]);
  if (!indices)
    return model_storage::AllocationFailure("Initial free indices", d);
  for (int i = 0; i < size; ++i) indices[i] = dofs->GlobalDofForEquation(i);
  Vector* rhs = nullptr;
  status = Vector::Create(size, &rhs, d);
  std::unique_ptr<Vector> rhs_owner(rhs);
  if (status != StatusCode::kOk) return status;
  status = rhs->GatherFrom(net_force, indices.get(), size, 0, d);
  if (status != StatusCode::kOk) return status;
  EquationSystem* system = nullptr;
  status = EquationSystem::Create(
      mass, rhs, MatrixType::kSymmetricPositiveDefinite, &system, d);
  std::unique_ptr<EquationSystem> system_owner(system);
  if (status != StatusCode::kOk) return status;
  mass_owner.release();
  rhs_owner.release();
  Vector* free = nullptr;
  status = Vector::Create(size, &free, d);
  std::unique_ptr<Vector> free_owner(free);
  if (status != StatusCode::kOk) return status;
  status = solver->Solve(system, free, d);
  if (status != StatusCode::kOk) return status;
  status = full->ScatterAdd(free, indices.get(), size, d);
  if (status != StatusCode::kOk) return status;
  return acceleration->CopyFrom(full, d);
}
namespace {
// Common N4/N8 sparse mass range solver; space fixes the chosen basis.
// finite optionally assembles the physical current force into net_force.
StatusCode SolveInitialShellMass(
    const AnalysisModel* model, const DofManager* dofs, MassMode mode,
    const LinearShellNullspace* space, const Vector* net_force,
    const Vector* null_acceleration, LinearSolver* solver, int worker_count,
    Vector* acceleration, Diagnostic* d, const FiniteMassInput* finite) {
  const int count = space->Count(), free_count = dofs->FreeCount();
  Vector* raw_full = nullptr;
  auto status = Vector::Create(dofs->GlobalCount(), &raw_full, d);
  std::unique_ptr<Vector> full(raw_full);
  if (status != StatusCode::kOk) return status;
  if (free_count == 0) return acceleration->CopyFrom(raw_full, d);
  const int size = free_count + count;
  SparseMatrix* raw_mass = nullptr;
  status = AssembleMass(model, dofs, mode, nullptr, worker_count, size,
                        &raw_mass, d, finite);
  std::unique_ptr<SparseMatrix> mass(raw_mass);
  if (status != StatusCode::kOk) return status;
  if (count != 0) {
    std::unique_ptr<int[]> rows(new (std::nothrow) int[3 * count]);
    std::unique_ptr<int[]> columns(new (std::nothrow) int[3 * count]);
    if (!rows || !columns)
      return model_storage::AllocationFailure("Initial shell constraints", d);
    Vector* raw_values = nullptr;
    status = Vector::Create(3 * count, &raw_values, d);
    std::unique_ptr<Vector> values(raw_values);
    if (status != StatusCode::kOk) return status;
    int entries = 0;
    for (int i = 0; i < count; ++i)
      for (int c = 0; c < 3; ++c) {
        const int row = space->Equation(i, c);
        if (row < 0) continue;
        rows[entries] = row;
        columns[entries] = free_count + i;
        status =
            values->Set(entries++, space->Directions()->Data()[3 * i + c], d);
        if (status != StatusCode::kOk) return status;
      }
    SparseMatrix* raw_constraints = nullptr;
    status = SparseMatrix::Create(size, entries, rows.get(), columns.get(),
                                  raw_values, &raw_constraints, d);
    std::unique_ptr<SparseMatrix> constraints(raw_constraints);
    if (status != StatusCode::kOk) return status;
    SparseMatrix* merged = nullptr;
    status = SparseMatrix::Add(mass.get(), raw_constraints, &merged, d);
    if (status != StatusCode::kOk) return status;
    mass.reset(merged);
  }
  std::unique_ptr<int[]> indices(new (std::nothrow) int[size]);
  if (!indices)
    return model_storage::AllocationFailure("Initial shell indices", d);
  for (int i = 0; i < free_count; ++i)
    indices[i] = dofs->GlobalDofForEquation(i);
  Vector *raw_rhs = nullptr, *raw_solution = nullptr;
  status = Vector::Create(size, &raw_rhs, d);
  std::unique_ptr<Vector> rhs(raw_rhs);
  if (status == StatusCode::kOk)
    status = Vector::Create(size, &raw_solution, d);
  std::unique_ptr<Vector> solution(raw_solution);
  if (status != StatusCode::kOk) return status;
  status = rhs->GatherFrom(net_force, indices.get(), free_count, 0, d);
  for (int i = 0; i < count && status == StatusCode::kOk; ++i)
    status = rhs->Set(free_count + i, null_acceleration->Data()[i], d);
  if (status != StatusCode::kOk) return status;
  EquationSystem* raw_system = nullptr;
  status =
      EquationSystem::Create(mass.get(), raw_rhs,
                             count == 0 ? MatrixType::kSymmetricPositiveDefinite
                                        : MatrixType::kSymmetricIndefinite,
                             &raw_system, d);
  std::unique_ptr<EquationSystem> system(raw_system);
  if (status != StatusCode::kOk) return status;
  mass.release();
  rhs.release();
  status = solver->Solve(raw_system, raw_solution, d);
  if (status != StatusCode::kOk) return status;
  // Discard only auxiliary multipliers; physical acceleration is the first
  // block.
  Vector* raw_free = nullptr;
  status = Vector::Create(free_count, &raw_free, d);
  std::unique_ptr<Vector> free(raw_free);
  if (status != StatusCode::kOk) return status;
  for (int i = 0; i < free_count; ++i) indices[i] = i;
  status = free->GatherFrom(raw_solution, indices.get(), free_count, 0, d);
  for (int i = 0; i < free_count; ++i)
    indices[i] = dofs->GlobalDofForEquation(i);
  if (status == StatusCode::kOk)
    status = full->ScatterAdd(raw_free, indices.get(), free_count, d);
  return status == StatusCode::kOk ? acceleration->CopyFrom(raw_full, d)
                                   : status;
}
}  // namespace
StatusCode ComputeInitialShellAcceleration(
    const AnalysisModel* model, const DofManager* dofs, MassMode mode,
    const Vector* net_force, const Vector* null_acceleration,
    LinearSolver* solver, int worker_count, Vector* acceleration,
    Diagnostic* d) {
  if (!model || !dofs || !solver || worker_count <= 0 || !net_force ||
      !acceleration || net_force == acceleration ||
      null_acceleration == acceleration ||
      net_force->Size() != dofs->GlobalCount() ||
      acceleration->Size() != dofs->GlobalCount() ||
      (mode != MassMode::kConsistent && mode != MassMode::kRowSum))
    return model_storage::Invalid("Initial shell acceleration inputs", d);
  LinearShellNullspace* raw_space = nullptr;
  auto status = LinearShellNullspace::Create(model, dofs, &raw_space, d);
  std::unique_ptr<LinearShellNullspace> space(raw_space);
  if (status != StatusCode::kOk) return status;
  const int count = space->Count(), free_count = dofs->FreeCount();
  if ((count == 0 ? null_acceleration != nullptr
                  : !null_acceleration || null_acceleration->Size() != count) ||
      count > INT_MAX / 3 || free_count > INT_MAX - count)
    return model_storage::Invalid("Initial shell acceleration null constraints",
                                  d);
  for (const auto* input : {net_force, null_acceleration}) {
    if (!input) continue;
    for (int i = 0; i < input->Size(); ++i)
      if (!std::isfinite(input->Data()[i]))
        return model_storage::Invalid(
            "Initial shell acceleration finite inputs", d);
  }
  return SolveInitialShellMass(model, dofs, mode, space.get(), net_force,
                               null_acceleration, solver, worker_count,
                               acceleration, d, nullptr);
}
StatusCode ComputeInitialFiniteShellAcceleration(
    const AnalysisModel* model, const DofManager* dofs,
    const ElementTrialState* pose, const Vector* velocity, MassMode mode,
    const Vector* external, const Vector* external_rate, LinearSolver* solver,
    int worker_count, Vector* acceleration, Diagnostic* d) {
  if (!model || !dofs || !pose || !velocity || !external || !external_rate ||
      !solver || worker_count <= 0 || !acceleration ||
      (mode != MassMode::kConsistent && mode != MassMode::kRowSum) ||
      !pose->translations || !pose->body_spin || !pose->orientations ||
      !pose->committed_orientations || dofs->NodeCount() <= 0 ||
      dofs->NodeCount() > INT_MAX / 6 ||
      dofs->GlobalCount() != 6 * dofs->NodeCount() ||
      pose->translations->Size() != 3 * dofs->NodeCount() ||
      pose->body_spin->Size() != 3 * dofs->NodeCount())
    return model_storage::Invalid("Finite shell initial acceleration inputs",
                                  d);
  for (const auto* q : {pose->orientations, pose->committed_orientations}) {
    if (q->Rows() != dofs->NodeCount() || q->Columns() != 4)
      return model_storage::Invalid("Finite shell initial pose shape", d);
    for (int i = 0; i < 4 * dofs->NodeCount(); ++i)
      if (!std::isfinite(q->Data()[i]))
        return model_storage::Invalid("Finite shell initial pose values", d);
  }
  for (const auto* input : {velocity, external, external_rate})
    if (input == acceleration || input->Size() != dofs->GlobalCount())
      return model_storage::Invalid("Finite shell initial global fields", d);
  for (const auto* input : {velocity, external, external_rate,
                            pose->translations, pose->body_spin}) {
    if (input == acceleration)
      return model_storage::Invalid("Aliased finite shell initial fields", d);
    for (int i = 0; i < input->Size(); ++i)
      if (!std::isfinite(input->Data()[i]))
        return model_storage::Invalid("Finite shell initial finite values", d);
  }
  if (acceleration->Size() != dofs->GlobalCount())
    return model_storage::Invalid("Finite shell initial acceleration shape", d);
  for (int i = 0; i < velocity->Size(); ++i)
    if (dofs->Equation(i) < 0 && velocity->Data()[i] != 0)
      return model_storage::Invalid("Fixed/excluded initial shell velocity", d);
  LinearShellNullspace* raw_space = nullptr;
  auto status = LinearShellNullspace::CreateFinite(
      model, dofs, pose->orientations, &raw_space, d);
  std::unique_ptr<LinearShellNullspace> space(raw_space);
  if (status != StatusCode::kOk) return status;
  const int count = space->Count();
  if (count > INT_MAX / 3 || dofs->FreeCount() > INT_MAX - count)
    return model_storage::Invalid("Finite shell initial augmented size", d);
  // Node-local borrowed data are gathered into call-local MKL-backed vectors.
  std::unique_ptr<Vector> fields[9];
  for (auto& field : fields) {
    Vector* raw = nullptr;
    status = Vector::Create(3, &raw, d);
    field.reset(raw);
    if (status != StatusCode::kOk) return status;
  }
  for (int i = 0; i < count; ++i) {
    const int direction[] = {3 * i, 3 * i + 1, 3 * i + 2};
    const int row = space->GlobalDof(i, 0) / 6;
    const int spin[] = {3 * row, 3 * row + 1, 3 * row + 2};
    const int angular[] = {space->GlobalDof(i, 0), space->GlobalDof(i, 1),
                           space->GlobalDof(i, 2)};
    status = space->ReferenceDirections()->GatherTo(direction, 3,
                                                    fields[0].get(), 0, d);
    if (status == StatusCode::kOk)
      status =
          space->Directions()->GatherTo(direction, 3, fields[1].get(), 0, d);
    if (status == StatusCode::kOk)
      status = fields[2]->GatherFrom(pose->body_spin, spin, 3, 0, d);
    const Vector* inputs[] = {velocity, external, external_rate};
    for (int j = 0; j < 3 && status == StatusCode::kOk; ++j)
      status = fields[3 + j]->GatherFrom(inputs[j], angular, 3, 0, d);
    if (status == StatusCode::kOk)
      status = CheckFiniteShellNodeCompatibility(
          space->NodeId(i), fields[0].get(), fields[1].get(),
          space->DrillingCoefficients()->Data()[i], fields[2].get(),
          fields[3].get(), fields[4].get(), fields[5].get(), d);
    if (status != StatusCode::kOk) return status;
  }
  Vector *raw_net = nullptr, *raw_range = nullptr, *raw_null = nullptr;
  status = Vector::Create(dofs->GlobalCount(), &raw_net, d);
  std::unique_ptr<Vector> net(raw_net);
  if (status == StatusCode::kOk)
    status = Vector::Create(dofs->GlobalCount(), &raw_range, d);
  std::unique_ptr<Vector> range(raw_range);
  if (status == StatusCode::kOk && count > 0)
    status = Vector::Create(count, &raw_null, d);
  std::unique_ptr<Vector> null(raw_null);
  if (status == StatusCode::kOk) status = net->CopyFrom(external, d);
  if (status != StatusCode::kOk) return status;
  FiniteMassInput input;
  input.pose = pose;
  input.velocity = velocity;
  input.net_force = net.get();
  // Z^T eta_r=0 defines the intermediate range solution only. N8 determines
  // the physical null acceleration below; no arbitrary drilling inertia.
  status = SolveInitialShellMass(model, dofs, mode, space.get(), net.get(),
                                 null.get(), solver, worker_count, range.get(),
                                 d, &input);
  if (status != StatusCode::kOk) return status;
  for (int i = 0; i < count; ++i) {
    const int direction[] = {3 * i, 3 * i + 1, 3 * i + 2};
    const int angular[] = {space->GlobalDof(i, 0), space->GlobalDof(i, 1),
                           space->GlobalDof(i, 2)};
    status = space->Directions()->GatherTo(direction, 3, fields[1].get(), 0, d);
    const Vector* inputs[] = {velocity, external, external_rate, range.get()};
    const int destinations[] = {3, 4, 5, 7};
    for (int j = 0; j < 4 && status == StatusCode::kOk; ++j)
      status = fields[destinations[j]]->GatherFrom(inputs[j], angular, 3, 0, d);
    double gamma = 0;
    if (status == StatusCode::kOk)
      status = ComputeFiniteShellNodeNullAcceleration(
          fields[1].get(), space->DrillingCoefficients()->Data()[i],
          fields[3].get(), fields[4].get(), fields[5].get(), fields[6].get(),
          fields[7].get(), &gamma, d);
    if (status == StatusCode::kOk)
      status = fields[8]->CopyFrom(fields[1].get(), d);
    if (status == StatusCode::kOk) status = fields[8]->Scale(gamma, d);
    if (status == StatusCode::kOk)
      status = range->ScatterAdd(fields[8].get(), angular, 3, d);
    if (status != StatusCode::kOk) {
      if (d) d->node_id = space->NodeId(i);
      return status;
    }
  }
  return acceleration->CopyFrom(range.get(), d);
}
}  // namespace fesa
