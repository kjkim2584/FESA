#include "analysis/finite_shell_dynamic_forces.h"

#include <algorithm>
#include <climits>
#include <memory>
#include <new>

#include "adapters/openmp/thread_runtime.h"
#include "analysis/nonlinear_element_pose.h"
#include "analysis/nonlinear_task_failure.h"
#include "analysis/numerical_support.h"
#include "fesa/analysis/rotation_state.h"
#include "fesa/elements/mitc4_shell.h"
#include "model/model_storage.h"

namespace fesa {
namespace {
using numerical_support::Buffers;
struct ForceTask {
  const DofManager* dofs = nullptr;  // Borrowed through join.
  const AnalysisState* state = nullptr;
  const Vector* velocity = nullptr;
  const Vector* acceleration = nullptr;
  Buffers* buffers = nullptr;  // Exclusive bounded output/worker stores.
  Vector* damping = nullptr;   // Exclusive per-item prefix until join.
  TrialEvaluationFailure* labels = nullptr;
  MassMode mode = MassMode::kConsistent;
  bool committed = false;
};
StatusCode ComputeForce(int item, int worker, void* context, Diagnostic* d) {
  auto* task = static_cast<ForceTask*>(context);
  auto* b = task->buffers;
  const auto* element = b->elements[item];
  d->element_id = element->Id();
  NonlinearElementPose* raw_pose = nullptr;
  auto status = NonlinearElementPose::Create(element, task->dofs, task->state,
                                             task->committed, &raw_pose, d);
  std::unique_ptr<NonlinearElementPose> pose(raw_pose);
  if (status != StatusCode::kOk) return status;
  auto* local = b->workers[worker];
  status = local->u->GatherFrom(task->velocity, pose->GlobalDofs(), 24, 0, d);
  if (status == StatusCode::kOk)
    status = local->force->GatherFrom(task->acceleration, pose->GlobalDofs(),
                                      24, 0, d);
  if (status != StatusCode::kOk) return status;
  Vector *raw_i = nullptr, *raw_d = nullptr;
  status = Vector::Create(24, &raw_i, d);
  std::unique_ptr<Vector> inertia(raw_i);
  if (status == StatusCode::kOk) status = Vector::Create(24, &raw_d, d);
  std::unique_ptr<Vector> damping(raw_d);
  if (status != StatusCode::kOk) return status;
  ElementMotionInput input;
  input.kinematics = DynamicKinematics::kFiniteRotation;
  input.finite_state = pose->Trial();
  input.velocity = local->u;
  input.acceleration = local->force;
  ElementDynamicResult output;
  output.inertia_force = raw_i;
  output.damping_force = raw_d;
  status = element->EvaluateDynamic(&input, task->mode, nullptr,
                                    DynamicEvaluationRequest::kForces, &output,
                                    task->labels + item, d);
  if (status != StatusCode::kOk) return status;
  int identity[24] = {};
  for (int i = 0; i < 24; ++i) {
    identity[i] = i;
    b->vector_indices[24 * item + i] = pose->GlobalDofs()[i];
  }
  status = b->vector_values->GatherFrom(raw_i, identity, 24, 24 * item, d);
  if (status == StatusCode::kOk)
    status = task->damping->GatherFrom(raw_d, identity, 24, 24 * item, d);
  return status;
}
}  // namespace
FiniteShellDynamicForces::~FiniteShellDynamicForces() {
  delete inertia_;
  delete damping_;
  delete inertia_origins_;
  delete damping_origins_;
}
StatusCode FiniteShellDynamicForces::Create(
    const AnalysisModel* model, const DofManager* dofs,
    const AnalysisState* state, MassMode mode, bool committed,
    FiniteShellDynamicForces** out, NonlinearFailure* failure, Diagnostic* d) {
  if (out) *out = nullptr;
  if (failure) *failure = NonlinearFailure::kNone;
  if (!out || !failure)
    return model_storage::Invalid("Finite shell force outputs", d);
  std::unique_ptr<FiniteShellDynamicForces> result(
      new (std::nothrow) FiniteShellDynamicForces);
  if (!result)
    return model_storage::AllocationFailure("Finite shell forces", d);
  const auto status =
      result->Initialize(model, dofs, state, mode, committed, failure, d);
  if (status != StatusCode::kOk) return status;
  *out = result.release();
  return StatusCode::kOk;
}
StatusCode FiniteShellDynamicForces::Initialize(const AnalysisModel* model,
                                                const DofManager* dofs,
                                                const AnalysisState* state,
                                                MassMode mode, bool committed,
                                                NonlinearFailure* failure,
                                                Diagnostic* d) {
  if (!model || !dofs || !state || InElementExecution() ||
      !state->FieldHistory() || !state->Rotations() || !state->Dynamics() ||
      dofs->NodeCount() <= 0 || dofs->NodeCount() > INT_MAX / 6 ||
      dofs->GlobalCount() != 6 * dofs->NodeCount() ||
      model->Nodes()->size() != static_cast<std::size_t>(dofs->NodeCount()) ||
      model->Elements()->empty() ||
      model->Elements()->size() > static_cast<std::size_t>(INT_MAX / 24) ||
      (mode != MassMode::kConsistent && mode != MassMode::kRowSum) ||
      (!committed && !state->HasValidCandidate()))
    return model_storage::Invalid("Finite shell force inputs", d);
  auto status = state->ValidateMapping(dofs, d);
  if (status != StatusCode::kOk) return status;
  const auto* dynamics = state->Dynamics();
  const Vector* velocity =
      committed ? dynamics->CommittedVelocity() : dynamics->CandidateVelocity();
  const Vector* acceleration = committed ? dynamics->CommittedAcceleration()
                                         : dynamics->CandidateAcceleration();
  if (!velocity || !acceleration)
    return model_storage::Invalid("Finite shell force source motion", d);
  for (const auto& entry : *model->Nodes())
    if (dofs->NodeDofCount(entry.first) != 6)
      return model_storage::Invalid("Finite shell force node layout", d);
  for (const auto& entry : *model->Elements()) {
    if (d) d->element_id = entry.first;
    const auto* shell = dynamic_cast<const Mitc4Shell*>(entry.second);
    if (!shell)
      return ReportFailure(StatusCode::kUnsupported,
                           DiagnosticStage::kValidation,
                           "Finite shell forces require MITC4", d);
    double density = 0;
    status = shell->material()->Density(&density, d);
    if (status != StatusCode::kOk) return status;
    if (!shell->SupportsDynamicEvaluation(DynamicKinematics::kFiniteRotation) ||
        !shell->SupportsMass(mode))
      return ReportFailure(StatusCode::kUnsupported,
                           DiagnosticStage::kValidation,
                           "Finite shell E0 force capability", d);
  }
  const int entries = 24 * static_cast<int>(model->Elements()->size());
  for (auto** target : {&inertia_, &damping_}) {
    status = Vector::Create(dofs->GlobalCount(), target, d);
    if (status != StatusCode::kOk) return status;
  }
  for (auto** target : {&inertia_origins_, &damping_origins_}) {
    status = Vector::Create(entries, target, d);
    if (status != StatusCode::kOk) return status;
  }
  const int capacity =
      std::min(static_cast<int>(model->Elements()->size()), Buffers::kElements);
  const int workers = ConfiguredThreadCount();
  Buffers* raw_buffers = nullptr;
  status = Buffers::Create(false, capacity, &raw_buffers, d, nullptr, workers,
                           24, false);
  std::unique_ptr<Buffers> buffers(raw_buffers);
  if (status != StatusCode::kOk) return status;
  Vector* raw_damping = nullptr;
  status = Vector::Create(24 * capacity, &raw_damping, d);
  std::unique_ptr<Vector> damping_chunk(raw_damping);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<int[]> identity(new (std::nothrow) int[24 * capacity]);
  std::unique_ptr<TrialEvaluationFailure[]> labels(
      new (std::nothrow) TrialEvaluationFailure[capacity]);
  if (!identity || !labels)
    return model_storage::AllocationFailure("Finite force chunk metadata", d);
  for (int i = 0; i < 24 * capacity; ++i) identity[i] = i;
  ForceTask task;
  task.dofs = dofs;
  task.state = state;
  task.velocity = velocity;
  task.acceleration = acceleration;
  task.buffers = raw_buffers;
  task.damping = raw_damping;
  task.labels = labels.get();
  task.mode = mode;
  task.committed = committed;
  auto next = model->Elements()->begin();
  int base = 0;
  while (next != model->Elements()->end()) {
    int count = 0;
    while (count < capacity && next != model->Elements()->end()) {
      raw_buffers->elements[count] = next->second;
      labels[count++] = TrialEvaluationFailure::kNone;
      ++next;
    }
    Diagnostic aggregate;
    Diagnostic* selected = d ? d : &aggregate;
    status = ExecuteElementTasks(count, workers, ComputeForce, &task,
                                 raw_buffers->diagnostics, nullptr, selected);
    status = SelectNonlinearTaskFailure(status, raw_buffers->elements, count,
                                        labels.get(), raw_buffers->diagnostics,
                                        workers, selected, failure);
    if (status != StatusCode::kOk) return status;
    status = inertia_origins_->GatherFrom(
        raw_buffers->vector_values, identity.get(), 24 * count, 24 * base, d);
    if (status == StatusCode::kOk)
      status = damping_origins_->GatherFrom(raw_damping, identity.get(),
                                            24 * count, 24 * base, d);
    if (status == StatusCode::kOk)
      status = inertia_->ScatterAdd(raw_buffers->vector_values,
                                    raw_buffers->vector_indices, 24 * count, d);
    if (status == StatusCode::kOk)
      status = damping_->ScatterAdd(raw_damping, raw_buffers->vector_indices,
                                    24 * count, d);
    if (status != StatusCode::kOk) return status;
    base += count;
  }
  return StatusCode::kOk;
}
const Vector* FiniteShellDynamicForces::Inertia() const { return inertia_; }
const Vector* FiniteShellDynamicForces::Damping() const { return damping_; }
const Vector* FiniteShellDynamicForces::InertiaOrigins() const {
  return inertia_origins_;
}
const Vector* FiniteShellDynamicForces::DampingOrigins() const {
  return damping_origins_;
}
}  // namespace fesa
