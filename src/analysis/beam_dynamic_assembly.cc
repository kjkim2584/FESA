#include "analysis/beam_dynamic_assembly.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <memory>
#include <new>

#include "adapters/openmp/thread_runtime.h"
#include "analysis/nonlinear_assembler.h"
#include "analysis/nonlinear_element_pose.h"
#include "analysis/nonlinear_task_failure.h"
#include "fesa/analysis/beam_internal_state.h"
#include "fesa/analysis/rotation_state.h"
#include "fesa/elements/beam3d.h"
#include "model/model_storage.h"

namespace fesa {
namespace {
StatusCode Validate(const AnalysisModel* model, const DofManager* dofs,
                    const AnalysisState* state, bool committed, Diagnostic* d) {
  if (!state || InElementExecution() || !state->Rotations() ||
      !state->FieldHistory() || !state->BeamInternals() ||
      (!committed && !state->HasValidCandidate()))
    return model_storage::Invalid("Beam dynamic state", d);
  auto status = ValidateNonlinearBeamAssemblyCounts(model, dofs, d);
  if (status != StatusCode::kOk) return status;
  status = state->ValidateMapping(dofs, d);
  if (status != StatusCode::kOk) return status;
  const auto* internal = state->BeamInternals();
  if (internal->ElementCount() != dofs->InternalDofCount() / 2 ||
      model->Elements()->size() > static_cast<std::size_t>(INT_MAX / 196))
    return model_storage::Invalid("Beam dynamic internal counts", d);
  const Vector* position =
      committed ? internal->CommittedPosition() : internal->CandidatePosition();
  if (!position || position->Size() != dofs->InternalDofCount())
    return model_storage::Invalid("Beam dynamic internal position", d);
  int row = 0;
  for (const auto& entry : *model->Elements()) {
    if (d) d->element_id = entry.first;
    if (internal->ElementId(row) != entry.first ||
        internal->Offset(entry.first) != 2 * row)
      return model_storage::Invalid("Beam dynamic internal IDs", d);
    const auto* beam = static_cast<const Beam3D*>(entry.second);
    double density = 0;
    status = beam->Section()->GetMaterial()->Density(&density, d);
    if (status != StatusCode::kOk) return status;
    if (!(density > 0) || !std::isfinite(density))
      return model_storage::Invalid("Beam dynamic positive density", d);
    if (!beam->SupportsDynamicEvaluation(DynamicKinematics::kFiniteRotation))
      return ReportFailure(StatusCode::kUnsupported,
                           DiagnosticStage::kValidation,
                           "Beam finite dynamic capability", d);
    ++row;
  }
  return StatusCode::kOk;
}
bool FiniteVector(const Vector* v, int size) {
  if (!v || v->Size() != size) return false;
  for (int i = 0; i < size; ++i)
    if (!std::isfinite(v->Data()[i])) return false;
  return true;
}
// Call-local task owner. Element outputs occupy disjoint item prefixes.
struct DynamicTask final {
  DynamicTask() = default;
  DynamicTask(const DynamicTask&) = delete;
  DynamicTask& operator=(const DynamicTask&) = delete;
  DynamicTask(DynamicTask&&) = delete;
  DynamicTask& operator=(DynamicTask&&) = delete;
  ~DynamicTask() {
    if (mass)
      for (int i = 0; i < count; ++i) delete mass[i];
    delete[] mass;
    delete[] elements;
    delete[] labels;
    delete[] diagnostics;
    delete[] globals;
    delete inertia;
    delete damping;
    delete structural;
  }
  const DofManager* dofs = nullptr;
  const AnalysisState* state = nullptr;
  const Vector* velocity = nullptr;
  const Vector* acceleration = nullptr;
  const Vector* internal_velocity = nullptr;
  const Vector* internal_acceleration = nullptr;
  const ElementWorkspace* const* workspaces = nullptr;  // Borrowed until join.
  bool committed = false;
  int count = 0;
  int* globals = nullptr;              // Owned physical12E indices.
  const Element** elements = nullptr;  // Owned metadata; elements borrowed.
  TrialEvaluationFailure* labels = nullptr;
  Diagnostic* diagnostics = nullptr;
  Vector* inertia = nullptr;     // Owned14E.
  Vector* damping = nullptr;     // Owned14E.
  Vector* structural = nullptr;  // Owned14E, initialization only.
  Matrix** mass = nullptr;       // Owned matrices, initialization only.
};
StatusCode PrepareTask(const AnalysisModel* model, bool initial,
                       DynamicTask* task, Diagnostic* d) {
  task->count = static_cast<int>(model->Elements()->size());
  task->elements = new (std::nothrow) const Element*[task->count];
  task->globals = new (std::nothrow) int[12 * task->count];
  task->labels = new (std::nothrow) TrialEvaluationFailure[task->count]();
  task->diagnostics = new (std::nothrow) Diagnostic[ConfiguredThreadCount()];
  if (initial) task->mass = new (std::nothrow) Matrix*[task->count]();
  if (!task->elements || !task->globals || !task->labels ||
      !task->diagnostics || (initial && !task->mass))
    return model_storage::AllocationFailure("Beam dynamic task metadata", d);
  int row = 0;
  for (const auto& entry : *model->Elements())
    task->elements[row++] = entry.second;
  auto status = Vector::Create(14 * task->count, &task->inertia, d);
  if (status == StatusCode::kOk)
    status = Vector::Create(14 * task->count, &task->damping, d);
  if (status == StatusCode::kOk && initial)
    status = Vector::Create(14 * task->count, &task->structural, d);
  return status;
}
StatusCode Evaluate(int item, int worker, void* context, Diagnostic* d) {
  auto* task = static_cast<DynamicTask*>(context);
  const auto* element = task->elements[item];
  d->element_id = element->Id();
  NonlinearElementPose* raw = nullptr;
  auto status = NonlinearElementPose::Create(element, task->dofs, task->state,
                                             task->committed, &raw, d);
  std::unique_ptr<NonlinearElementPose> pose(raw);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<Vector> fields[7];
  const int sizes[] = {12, 12, 2, 2, 14, 14, 14};
  for (int i = 0; i < (task->mass ? 7 : 6); ++i) {
    Vector* vector = nullptr;
    status = Vector::Create(sizes[i], &vector, d);
    fields[i].reset(vector);
    if (status != StatusCode::kOk) return status;
  }
  const int internal[] = {2 * item, 2 * item + 1};
  status = fields[0]->GatherFrom(task->velocity, pose->GlobalDofs(), 12, 0, d);
  if (status == StatusCode::kOk)
    status = fields[2]->GatherFrom(task->internal_velocity, internal, 2, 0, d);
  if (status == StatusCode::kOk && task->acceleration)
    status =
        fields[1]->GatherFrom(task->acceleration, pose->GlobalDofs(), 12, 0, d);
  if (status == StatusCode::kOk && task->internal_acceleration)
    status =
        fields[3]->GatherFrom(task->internal_acceleration, internal, 2, 0, d);
  if (status != StatusCode::kOk) return status;
  ElementMotionInput input;
  input.kinematics = DynamicKinematics::kFiniteRotation;
  input.finite_state = pose->Trial();
  input.velocity = fields[0].get();
  input.acceleration = fields[1].get();
  input.internal_velocity = fields[2].get();
  input.internal_acceleration = fields[3].get();
  ElementDynamicResult result;
  result.inertia_force = fields[4].get();
  result.damping_force = fields[5].get();
  std::unique_ptr<Matrix> partials[5];
  if (task->mass) {
    for (int i = 0; i < 5; ++i) {
      Matrix* matrix = nullptr;
      status = Matrix::Create(14, 14, &matrix, d);
      partials[i].reset(matrix);
      if (status != StatusCode::kOk) return status;
    }
    result.inertia_displacement = partials[0].get();
    result.inertia_velocity = partials[1].get();
    result.inertia_acceleration = partials[2].get();
    result.damping_displacement = partials[3].get();
    result.damping_velocity = partials[4].get();
  }
  status = element->EvaluateDynamic(
      &input, MassMode::kConsistent,
      task->workspaces ? task->workspaces[worker] : nullptr,
      task->mass ? DynamicEvaluationRequest::kForcesAndLinearization
                 : DynamicEvaluationRequest::kForces,
      &result, task->labels + item, d);
  if (status != StatusCode::kOk) return status;
  int identity[14] = {};
  for (int i = 0; i < 14; ++i) identity[i] = i;
  for (int i = 0; i < 12; ++i)
    task->globals[12 * item + i] = pose->GlobalDofs()[i];
  status =
      task->inertia->GatherFrom(fields[4].get(), identity, 14, 14 * item, d);
  if (status == StatusCode::kOk)
    status =
        task->damping->GatherFrom(fields[5].get(), identity, 14, 14 * item, d);
  if (status == StatusCode::kOk && task->mass) {
    status = element->EvaluateTrial(
        pose->Trial(), task->workspaces ? task->workspaces[worker] : nullptr,
        ElementEvaluationRequest::kInternalForce,
        fields[6].get(), nullptr, task->labels + item, d);
    if (status == StatusCode::kOk)
      status = task->structural->GatherFrom(fields[6].get(), identity, 14,
                                            14 * item, d);
    if (status == StatusCode::kOk) task->mass[item] = partials[2].release();
  }
  return status;
}
StatusCode RunTask(DynamicTask* task, NonlinearFailure* failure,
                   Diagnostic* d) {
  Diagnostic aggregate;
  auto* selected = d ? d : &aggregate;
  const int workers = ConfiguredThreadCount();
  auto status = ExecuteElementTasks(task->count, workers, Evaluate, task,
                                    task->diagnostics, nullptr, selected);
  return SelectNonlinearTaskFailure(status, task->elements, task->count,
                                    task->labels, task->diagnostics, workers,
                                    selected, failure);
}
StatusCode AssembleForces(const DynamicTask* task, const Vector* origins,
                          Vector* physical, Vector* internal, Diagnostic* d) {
  Vector *raw_p = nullptr, *raw_i = nullptr;
  auto status = Vector::Create(12, &raw_p, d);
  std::unique_ptr<Vector> p(raw_p);
  if (status == StatusCode::kOk) status = Vector::Create(2, &raw_i, d);
  std::unique_ptr<Vector> i(raw_i);
  if (status != StatusCode::kOk) return status;
  for (int item = 0; item < task->count; ++item) {
    int indices[14] = {};
    for (int j = 0; j < 14; ++j) indices[j] = 14 * item + j;
    status = p->GatherFrom(origins, indices, 12, 0, d);
    if (status == StatusCode::kOk)
      status = physical->ScatterAdd(p.get(), task->globals + 12 * item, 12, d);
    if (status == StatusCode::kOk)
      status = i->GatherFrom(origins, indices + 12, 2, 0, d);
    const int target[] = {2 * item, 2 * item + 1};
    if (status == StatusCode::kOk)
      status = internal->ScatterAdd(i.get(), target, 2, d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}
bool Owned(const AnalysisState* state, const Vector* v) {
  const auto* i = state->BeamInternals();
  const auto* m = state->Dynamics();
  const auto* r = state->Rotations();
  return v &&
         (v == state->CandidateDisplacement() ||
          v == state->CommittedDisplacement() ||
          v == state->CandidateReaction() || v == state->CommittedReaction() ||
          (i && (v == i->CandidatePosition() || v == i->CommittedPosition() ||
                 v == i->CandidateVelocity() || v == i->CommittedVelocity() ||
                 v == i->CandidateAcceleration() ||
                 v == i->CommittedAcceleration())) ||
          (m && (v == m->CandidateVelocity() || v == m->CommittedVelocity() ||
                 v == m->CandidateAcceleration() ||
                 v == m->CommittedAcceleration())) ||
          (r && (v == r->CandidateBodySpin() || v == r->CommittedBodySpin())));
}
}  // namespace
BeamDynamicForces::~BeamDynamicForces() {
  delete inertia_;
  delete damping_;
  delete internal_inertia_;
  delete internal_damping_;
  delete inertia_origins_;
  delete damping_origins_;
}
StatusCode BeamDynamicForces::Create(const AnalysisModel* model,
                                     const DofManager* dofs,
                                     const AnalysisState* state, MassMode mode,
                                     bool committed, BeamDynamicForces** out,
                                     NonlinearFailure* failure, Diagnostic* d,
                                     const ElementWorkspace* const* workspaces,
                                     int workspace_count) {
  if (out) *out = nullptr;
  if (failure) *failure = NonlinearFailure::kNone;
  if (!out || !failure)
    return model_storage::Invalid("Beam dynamic force outputs", d);
  if (mode != MassMode::kConsistent)
    return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                         "Beam dynamic consistent mass only", d);
  auto status = Validate(model, dofs, state, committed, d);
  if (status != StatusCode::kOk) return status;
  if (workspace_count < 0 || (!workspaces && workspace_count) ||
      (workspaces && workspace_count < ConfiguredThreadCount()))
    return model_storage::Invalid("Beam dynamic worker contexts", d);
  if (workspaces)
    for (int i = 0; i < ConfiguredThreadCount(); ++i)
      if (!workspaces[i]) return model_storage::Invalid("Null Beam worker", d);
  const auto* motion = state->Dynamics();
  const auto* internal = state->BeamInternals();
  DynamicTask task;
  task.workspaces = workspaces;
  task.dofs = dofs;
  task.state = state;
  task.committed = committed;
  if (motion) {
    task.velocity =
        committed ? motion->CommittedVelocity() : motion->CandidateVelocity();
    task.acceleration = committed ? motion->CommittedAcceleration()
                                  : motion->CandidateAcceleration();
  }
  task.internal_velocity =
      committed ? internal->CommittedVelocity() : internal->CandidateVelocity();
  task.internal_acceleration = committed ? internal->CommittedAcceleration()
                                         : internal->CandidateAcceleration();
  if (!FiniteVector(task.velocity, dofs->GlobalCount()) ||
      !FiniteVector(task.acceleration, dofs->GlobalCount()) ||
      !FiniteVector(task.internal_velocity, dofs->InternalDofCount()) ||
      !FiniteVector(task.internal_acceleration, dofs->InternalDofCount()))
    return model_storage::Invalid("Beam dynamic complete source motion", d);
  status = PrepareTask(model, false, &task, d);
  if (status == StatusCode::kOk) status = RunTask(&task, failure, d);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<BeamDynamicForces> packet(new (std::nothrow)
                                                BeamDynamicForces);
  if (!packet)
    return model_storage::AllocationFailure("Beam dynamic forces", d);
  for (auto** target : {&packet->inertia_, &packet->damping_}) {
    status = Vector::Create(dofs->GlobalCount(), target, d);
    if (status != StatusCode::kOk) return status;
  }
  for (auto** target :
       {&packet->internal_inertia_, &packet->internal_damping_}) {
    status = Vector::Create(dofs->InternalDofCount(), target, d);
    if (status != StatusCode::kOk) return status;
  }
  status = AssembleForces(&task, task.inertia, packet->inertia_,
                          packet->internal_inertia_, d);
  if (status == StatusCode::kOk)
    status = AssembleForces(&task, task.damping, packet->damping_,
                            packet->internal_damping_, d);
  if (status != StatusCode::kOk) return status;
  packet->inertia_origins_ = task.inertia;
  packet->damping_origins_ = task.damping;
  task.inertia = nullptr;
  task.damping = nullptr;
  *out = packet.release();
  return StatusCode::kOk;
}
const Vector* BeamDynamicForces::Inertia() const { return inertia_; }
const Vector* BeamDynamicForces::Damping() const { return damping_; }
const Vector* BeamDynamicForces::InternalInertia() const {
  return internal_inertia_;
}
const Vector* BeamDynamicForces::InternalDamping() const {
  return internal_damping_;
}
const Vector* BeamDynamicForces::InertiaOrigins() const {
  return inertia_origins_;
}
const Vector* BeamDynamicForces::DampingOrigins() const {
  return damping_origins_;
}
StatusCode ComputeInitialBeamAcceleration(
    const AnalysisModel* model, const DofManager* dofs,
    const AnalysisState* state, const Vector* physical_velocity,
    const Vector* internal_velocity, const Vector* external,
    LinearSolver* solver, Vector* physical_acceleration,
    Vector* internal_acceleration, NonlinearFailure* failure, Diagnostic* d) {
  if (failure) *failure = NonlinearFailure::kNone;
  if (!failure || !solver)
    return model_storage::Invalid("Beam initial acceleration outputs", d);
  auto status = Validate(model, dofs, state, false, d);
  if (status != StatusCode::kOk) return status;
  if (!FiniteVector(physical_velocity, dofs->GlobalCount()) ||
      !FiniteVector(internal_velocity, dofs->InternalDofCount()) ||
      !FiniteVector(external, dofs->GlobalCount()) || !physical_acceleration ||
      !internal_acceleration ||
      physical_acceleration->Size() != dofs->GlobalCount() ||
      internal_acceleration->Size() != dofs->InternalDofCount() ||
      physical_acceleration == internal_acceleration ||
      Owned(state, physical_acceleration) ||
      Owned(state, internal_acceleration))
    return model_storage::Invalid("Beam initial acceleration shapes/ownership",
                                  d);
  for (auto* output : {physical_acceleration, internal_acceleration})
    if (output == physical_velocity || output == internal_velocity ||
        output == external)
      return model_storage::Invalid("Beam initial acceleration aliases", d);
  for (int i = 0; i < dofs->GlobalCount(); ++i)
    if (dofs->Equation(i) < 0 && physical_velocity->Data()[i] != 0)
      return model_storage::Invalid("Beam initial fixed velocity", d);
  DynamicTask task;
  task.dofs = dofs;
  task.state = state;
  task.velocity = physical_velocity;
  task.internal_velocity = internal_velocity;
  status = PrepareTask(model, true, &task, d);
  if (status == StatusCode::kOk) status = RunTask(&task, failure, d);
  if (status != StatusCode::kOk) return status;
  // Build the local physical residual through library algebra before scatter.
  status = task.structural->AddScaled(1, task.inertia, d);
  if (status == StatusCode::kOk)
    status = task.structural->AddScaled(1, task.damping, d);
  if (status == StatusCode::kOk) status = task.structural->Scale(-1, d);
  if (status != StatusCode::kOk) return status;
  const int count = dofs->EquationCount();
  std::unique_ptr<Vector> rhs, values, solution, staged_p, staged_i;
  const int sizes[] = {count, 105 * task.count, count, dofs->GlobalCount(),
                       dofs->InternalDofCount()};
  std::unique_ptr<Vector>* owners[] = {&rhs, &values, &solution, &staged_p,
                                       &staged_i};
  for (int i = 0; i < 5; ++i) {
    Vector* raw = nullptr;
    status = Vector::Create(sizes[i], &raw, d);
    owners[i]->reset(raw);
    if (status != StatusCode::kOk) return status;
  }
  std::unique_ptr<int[]> rows(new (std::nothrow) int[105 * task.count]);
  std::unique_ptr<int[]> columns(new (std::nothrow) int[105 * task.count]);
  std::unique_ptr<int[]> free_indices(new (std::nothrow) int[count]);
  if (!rows || !columns || !free_indices)
    return model_storage::AllocationFailure("Beam initial COO metadata", d);
  for (int i = 0; i < dofs->FreeCount(); ++i)
    free_indices[i] = dofs->GlobalDofForEquation(i);
  status =
      rhs->GatherFrom(external, free_indices.get(), dofs->FreeCount(), 0, d);
  if (status != StatusCode::kOk) return status;
  int nnz = 0;
  Vector* raw_local = nullptr;
  status = Vector::Create(14, &raw_local, d);
  std::unique_ptr<Vector> local(raw_local);
  if (status != StatusCode::kOk) return status;
  for (int item = 0; item < task.count; ++item) {
    int equations[14] = {}, source[14] = {}, retained[14] = {};
    int n = 0;
    for (int j = 0; j < 14; ++j) {
      equations[j] =
          j < 12 ? dofs->Equation(task.globals[12 * item + j])
                 : dofs->InternalEquation(task.elements[item]->Id(), j - 12);
      if (equations[j] >= 0) {
        source[n] = 14 * item + j;
        retained[n++] = equations[j];
      }
    }
    status = local->GatherFrom(task.structural, source, n, 0, d);
    if (status == StatusCode::kOk)
      status = rhs->ScatterAdd(local.get(), retained, n, d);
    if (status != StatusCode::kOk) return status;
    // Select the computed LOCAL upper entry, then order GLOBAL COO indices.
    // No mirroring, averaging, condensation, or artificial rotational mass.
    for (int j = 0; j < 14; ++j) {
      if (equations[j] < 0) continue;
      for (int k = j; k < 14; ++k) {
        if (equations[k] < 0) continue;
        rows[nnz] = std::min(equations[j], equations[k]);
        columns[nnz] = std::max(equations[j], equations[k]);
        status = values->Set(nnz++, task.mass[item]->Data()[14 * j + k], d);
        if (status != StatusCode::kOk) return status;
      }
    }
  }
  SparseMatrix* raw_mass = nullptr;
  status = SparseMatrix::Create(count, nnz, rows.get(), columns.get(),
                                values.get(), &raw_mass, d);
  std::unique_ptr<SparseMatrix> mass(raw_mass);
  if (status != StatusCode::kOk) return status;
  EquationSystem* raw_system = nullptr;
  status = EquationSystem::Create(mass.get(), rhs.get(),
                                  MatrixType::kSymmetricPositiveDefinite,
                                  &raw_system, d);
  std::unique_ptr<EquationSystem> system(raw_system);
  if (status != StatusCode::kOk) return status;
  mass.release();
  rhs.release();
  status = solver->Solve(system.get(), solution.get(), d);
  if (status != StatusCode::kOk) return status;
  Vector* raw_free = nullptr;
  if (dofs->FreeCount() > 0) {
    status = Vector::Create(dofs->FreeCount(), &raw_free, d);
    std::unique_ptr<Vector> free(raw_free);
    if (status != StatusCode::kOk) return status;
    for (int i = 0; i < dofs->FreeCount(); ++i) free_indices[i] = i;
    status = free->GatherFrom(solution.get(), free_indices.get(),
                              dofs->FreeCount(), 0, d);
    for (int i = 0; i < dofs->FreeCount(); ++i)
      free_indices[i] = dofs->GlobalDofForEquation(i);
    if (status == StatusCode::kOk)
      status = staged_p->ScatterAdd(free.get(), free_indices.get(),
                                    dofs->FreeCount(), d);
    if (status != StatusCode::kOk) return status;
  }
  for (int i = 0; i < dofs->InternalDofCount(); ++i)
    free_indices[i] = dofs->FreeCount() + i;
  status = staged_i->GatherFrom(solution.get(), free_indices.get(),
                                dofs->InternalDofCount(), 0, d);
  if (status == StatusCode::kOk)
    status = physical_acceleration->CopyFrom(staged_p.get(), d);
  if (status == StatusCode::kOk)
    status = internal_acceleration->CopyFrom(staged_i.get(), d);
  return status;
}
}  // namespace fesa
