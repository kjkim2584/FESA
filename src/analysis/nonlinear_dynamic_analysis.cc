#include "fesa/analysis/nonlinear_dynamic_analysis.h"

#include <climits>
#include <cmath>
#include <memory>
#include <new>
#include <set>
#include <tuple>
#include <vector>

#include "adapters/mkl/pardiso_solver.h"
#include "adapters/openmp/thread_runtime.h"
#include "analysis/beam_dynamic_assembly.h"
#include "analysis/dynamic_initialization.h"
#include "analysis/dynamic_load.h"
#include "analysis/finite_dynamic_kinematics.h"
#include "analysis/nonlinear_assembler.h"
#include "analysis/nonlinear_scaling.h"
#include "analysis/nonlinear_step_conditions.h"
#include "base/diagnostic_copy.h"
#include "fesa/analysis/beam_internal_state.h"
#include "fesa/analysis/nonlinear_field_history.h"
#include "fesa/analysis/rotation_state.h"
#include "fesa/analysis/static_workspace.h"
#include "fesa/elements/beam3d.h"
#include "fesa/elements/mitc4_shell.h"
#include "fesa/model/amplitude.h"
#include "model/model_storage.h"

namespace fesa {
namespace {
// Caller-side physical CLOAD typing. Zero/cancelling/fixed origins retain
// their identities; no weighted or assembled resultant is used for R11.
StatusCode ReadLoads(const DynamicStepContext* context, double time,
                     double total, Vector** external, Vector** force,
                     Vector** moment, Diagnostic* d) {
  *external = *force = *moment = nullptr;
  Vector *raw = nullptr, *raw_origins = nullptr;
  auto status = EvaluateDynamicClLoads(context->conditions, context->dofs, time,
                                       total, &raw, &raw_origins, d);
  std::unique_ptr<Vector> global(raw), origins(raw_origins), fields[2];
  if (status != StatusCode::kOk) return status;
  std::vector<int> indices[2];
  int offset = 0;
  for (const auto* load : *context->conditions->Loads())
    indices[(load->definition->Dof() - 1) / 3].push_back(offset++);
  for (int field = 0; field < 2; ++field) {
    if (indices[field].empty()) continue;
    status = Vector::Create(static_cast<int>(indices[field].size()), &raw, d);
    fields[field].reset(raw);
    if (status == StatusCode::kOk)
      status = raw->GatherFrom(origins.get(), indices[field].data(),
                               raw->Size(), 0, d);
    if (status != StatusCode::kOk) return status;
  }
  *external = global.release();
  *force = fields[0].release();
  *moment = fields[1].release();
  return StatusCode::kOk;
}
// Only scalar input provenance is compared here, not a mechanical force
// cache. Constants may be redefined equivalently. TABULAR contexts must keep
// their origin coefficients/amplitude and STEP clock or TOTAL clock identity.
bool SameLoadFunctions(const EffectiveConditions* previous,
                       const EffectiveConditions* current) {
  using Signature = std::tuple<int, int, int, int, double>;
  std::multiset<Signature> signatures[2];
  const EffectiveConditions* sources[] = {previous, current};
  for (int side = 0; side < 2; ++side)
    for (const auto* load : *sources[side]->Loads()) {
      const auto* definition = load->definition;
      const auto* amplitude = definition->DynamicAmplitude();
      if (!load->active || !amplitude || amplitude->PointCount() == 1) continue;
      const int clock = amplitude->TimeBasis() == AmplitudeTimeBasis::kStepTime
                            ? definition->Step()
                            : 0;
      signatures[side].insert(std::make_tuple(load->node_id, definition->Dof(),
                                              amplitude->Id(), clock,
                                              load->value));
    }
  return signatures[0] == signatures[1];
}
StatusCode ClassifyContinuation(const DynamicStepContext* context,
                                const EffectiveConditions* previous,
                                const Vector* current_external,
                                bool* continuous, Diagnostic* d) {
  *continuous = false;
  if (!SameLoadFunctions(previous, context->conditions)) return StatusCode::kOk;
  const auto* old_step =
      previous->SourceDomain()->Steps()->at(previous->StepNumber());
  Vector *raw = nullptr, *raw_origins = nullptr;
  auto status = EvaluateDynamicClLoads(
      previous, context->dofs, old_step->Period(),
      context->total_time_before_step, &raw, &raw_origins, d);
  std::unique_ptr<Vector> old(raw), origins(raw_origins);
  double change = 0;
  if (status == StatusCode::kOk)
    status = old->AddScaled(-1, current_external, d);
  // Exact zero predicates only; this is not a force/moment convergence norm.
  if (status == StatusCode::kOk) status = old->Norm2(&change, d);
  if (status != StatusCode::kOk || change != 0) return status;
  Vector *raw_left = nullptr, *left_origins = nullptr;
  status = EvaluateDynamicClLoadRates(
      previous, context->dofs, std::nextafter(old_step->Period(), 0.),
      std::nextafter(context->total_time_before_step, 0.), &raw_left,
      &left_origins, d);
  std::unique_ptr<Vector> left(raw_left), left_owner(left_origins);
  Vector *raw_right = nullptr, *right_origins = nullptr;
  if (status == StatusCode::kOk)
    status = EvaluateDynamicClLoadRates(context->conditions, context->dofs, 0,
                                        context->total_time_before_step,
                                        &raw_right, &right_origins, d);
  std::unique_ptr<Vector> right(raw_right), right_owner(right_origins);
  if (status == StatusCode::kOk) status = left->AddScaled(-1, right.get(), d);
  if (status == StatusCode::kOk) status = left->Norm2(&change, d);
  if (status == StatusCode::kOk) *continuous = change == 0;
  return status;
}
StatusCode InitialAcceleration(const DynamicStepContext* context,
                               const Vector* velocity, MassMode mode,
                               const Vector* external, LinearSolver* solver,
                               Vector* acceleration, Diagnostic* d) {
  const auto* state = context->state;
  const int nodes = context->dofs->NodeCount();
  Vector *raw_u = nullptr, *raw_node = nullptr;
  auto status = Vector::Create(3 * nodes, &raw_u, d);
  std::unique_ptr<Vector> u(raw_u);
  if (status == StatusCode::kOk) status = Vector::Create(3, &raw_node, d);
  std::unique_ptr<Vector> node(raw_node);
  if (status != StatusCode::kOk) return status;
  for (int row = 0; row < nodes; ++row) {
    const int indices[] = {3 * row, 3 * row + 1, 3 * row + 2};
    status = state->GatherCandidateTranslation(row, node.get(), d);
    if (status == StatusCode::kOk)
      status = u->ScatterAdd(node.get(), indices, 3, d);
    if (status != StatusCode::kOk) return status;
  }
  Vector *raw_rate = nullptr, *raw_origins = nullptr;
  status = EvaluateDynamicClLoadRates(context->conditions, context->dofs, 0,
                                      context->total_time_before_step,
                                      &raw_rate, &raw_origins, d);
  std::unique_ptr<Vector> rate(raw_rate), origins(raw_origins);
  if (status != StatusCode::kOk) return status;
  ElementTrialState pose;
  pose.translations = u.get();
  pose.orientations = state->Rotations()->CandidateOrientations();
  pose.committed_orientations = state->Rotations()->CommittedOrientations();
  pose.body_spin = state->Rotations()->CandidateBodySpin();
  return ComputeInitialFiniteShellAcceleration(
      context->model, context->dofs, &pose, velocity, mode, external,
      rate.get(), solver, ConfiguredThreadCount(), acceleration, d);
}
}  // namespace
StatusCode NonlinearDynamicAnalysis::Create(const DynamicStepSettings* settings,
                                            MassMode mode,
                                            const EffectiveConditions* previous,
                                            const Vector* initial_velocity,
                                            NonlinearDynamicAnalysis** out,
                                            Diagnostic* d) {
  if (!out) return model_storage::Invalid("NonlinearDynamicAnalysis output", d);
  *out = nullptr;
  if (mode != MassMode::kConsistent && mode != MassMode::kRowSum)
    return model_storage::Invalid("Finite dynamic mass mode", d);
  try {
    std::unique_ptr<NonlinearDynamicAnalysis> result(
        new NonlinearDynamicAnalysis);
    auto status = result->InitializeIntegration(settings, d);
    if (status == StatusCode::kOk)
      status = NewtonRaphson::Create(&result->newton_, d);
    if (status == StatusCode::kOk)
      status = NewtonReport::Create(&result->newton_report_, d);
    if (status != StatusCode::kOk) return status;
    result->solver_ = new PardisoSolver;
    result->mode_ = mode;
    result->previous_ = previous;
    result->initial_velocity_ = initial_velocity;
    *out = result.release();
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return model_storage::AllocationFailure("NonlinearDynamicAnalysis", d);
  }
}
NonlinearDynamicAnalysis::~NonlinearDynamicAnalysis() {
  ReleaseStepResources();
  delete solver_;
  delete newton_report_;
  delete newton_;
}
StatusCode NonlinearDynamicAnalysis::ValidateStep(
    const AnalysisModel* model, const DofManager* dofs, MassMode mode,
    const EffectiveConditions* previous, Diagnostic* d) {
  if (!model || !dofs || !model->Conditions() || model->Elements()->empty())
    return model_storage::Invalid("Finite dynamic model preflight", d);
  const bool beam_role = dofs->InternalDofCount() > 0;
  for (const auto& entry : *model->Elements())
    if (entry.second->family() !=
        (beam_role ? ElementFamily::kBeam : ElementFamily::kShell)) {
      if (d) d->element_id = entry.first;
      return ReportFailure(StatusCode::kUnsupported,
                           DiagnosticStage::kValidation,
                           "Finite dynamics element family/mapping", d);
    }
  if (dofs->NodeCount() <= 0 || dofs->NodeCount() > INT_MAX / 9 ||
      dofs->GlobalCount() != 6 * dofs->NodeCount() ||
      (mode != MassMode::kConsistent && mode != MassMode::kRowSum))
    return model_storage::Invalid("Finite dynamic model preflight", d);
  const StepDefinition* step = nullptr;
  auto status = model->Conditions()->SourceDomain()->FindStepDefinition(
      model->Conditions()->StepNumber(), &step, d);
  if (status != StatusCode::kOk) return status;
  if (!step->IsDynamic() || !step->HasFiniteGeometry())
    return model_storage::Invalid("Finite dynamic step identity", d);
  if (beam_role) {
    status = ValidateNonlinearBeamAssemblyCounts(model, dofs, d);
    if (status != StatusCode::kOk) return status;
    if (mode != MassMode::kConsistent)
      return ReportFailure(StatusCode::kUnsupported,
                           DiagnosticStage::kValidation,
                           "Finite Beam requires consistent mass", d);
  }
  try {
    status =
        ValidateDynamicConstraintTransition(previous, model->Conditions(), d);
    if (status != StatusCode::kOk) return status;
    for (int row = 0; row < dofs->NodeCount(); ++row) {
      const int id = dofs->NodeId(row);
      if (dofs->NodeDofCount(id) != 6 ||
          model->Nodes()->find(id) == model->Nodes()->end())
        return model_storage::Invalid("Finite dynamic node layout", d);
      int fixed = 0;
      for (int c = 4; c <= 6; ++c)
        if (dofs->Equation(dofs->GlobalDof(id, c)) == -1) ++fixed;
      if (fixed != 0 && fixed != 3)
        return ReportFailure(
            StatusCode::kUnsupported, DiagnosticStage::kValidation,
            "Finite dynamic requires full rotation constraints", d);
    }
    for (const auto& entry : *model->Elements()) {
      if (d) d->element_id = entry.first;
      const auto* shell = dynamic_cast<const Mitc4Shell*>(entry.second);
      const auto* beam = dynamic_cast<const Beam3D*>(entry.second);
      if ((beam_role ? !beam : !shell || !shell->SupportsMass(mode)) ||
          !entry.second->SupportsDynamicEvaluation(
              DynamicKinematics::kFiniteRotation))
        return ReportFailure(StatusCode::kUnsupported,
                             DiagnosticStage::kValidation,
                             "Finite dynamic element capability", d);
      double density = 0;
      status = (beam_role ? beam->Section()->GetMaterial() : shell->material())
                   ->Density(&density, d);
      if (status != StatusCode::kOk) return status;
    }
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return model_storage::AllocationFailure("Finite dynamic transition", d);
  }
}
StatusCode NonlinearDynamicAnalysis::BeginStep(
    const DynamicStepContext* context, Diagnostic* d) {
  ReleaseStepResources();
  newton_report_->Reset();
  element_started_ = ElementDispatchSeconds();
  if (!context->controls->newton ||
      context->model->Conditions() != context->conditions ||
      !context->state->Rotations() || !context->state->FieldHistory())
    return model_storage::Invalid("Finite dynamic execution context", d);
  auto status = context->state->ValidateMapping(context->dofs, d);
  if (status == StatusCode::kOk)
    status = ValidateStep(context->model, context->dofs, mode_, previous_, d);
  if (status != StatusCode::kOk || !context->dofs->InternalDofCount()) return status;
  const int count = ConfiguredThreadCount();
  status = PrepareWorkerWorkspaces(count, 12, &workers_, d);
  if (status == StatusCode::kOk) worker_count_ = count;
  return status;
}
StatusCode NonlinearDynamicAnalysis::InitializeCandidate(
    const DynamicStepContext* context, Diagnostic* d) {
  auto* state = context->state;
  const bool has_motion =
      state->Dynamics() && state->Dynamics()->CommittedVelocity();
  if ((state->StepNumber() != 0 && initial_velocity_) ||
      (has_motion &&
       context->total_time_before_step != state->Dynamics()->TotalTime()))
    return model_storage::Invalid("Finite dynamic initial history", d);
  Vector *raw_load = nullptr, *raw_force = nullptr, *raw_moment = nullptr;
  auto status = ReadLoads(context, 0, context->total_time_before_step,
                          &raw_load, &raw_force, &raw_moment, d);
  std::unique_ptr<Vector> external(raw_load), force_origins(raw_force),
      moment_origins(raw_moment);
  if (status != StatusCode::kOk) return status;
  bool continuous = has_motion;
  if (has_motion) {
    status = ClassifyContinuation(context, previous_, external.get(),
                                  &continuous, d);
    if (status != StatusCode::kOk) return status;
  }
  status = state->EnableFiniteDynamics(d);
  if (status == StatusCode::kOk && !has_motion)
    status = state->EnableDynamicConvergenceHistory(d);
  if (status == StatusCode::kOk) status = state->PrepareCandidate(d);
  if (status != StatusCode::kOk) return status;
  NonlinearFailure failure = NonlinearFailure::kNone;
  if (state->StepNumber() == 0) {
    NonlinearStepConditionInput path;
    path.conditions = context->conditions;
    path.dofs = context->dofs;
    path.start_displacement = state->CommittedDisplacement();
    path.start_reaction = state->CommittedReaction();
    path.start_orientations = state->Rotations()->CommittedOrientations();
    path.lambda = 1;
    status = ApplyNonlinearStepPrescriptions(&path, state, &failure, d);
    if (status != StatusCode::kOk) return status;
  }
  Vector *raw_v = nullptr, *raw_a = nullptr;
  status = Vector::Create(state->GlobalCount(), &raw_v, d);
  std::unique_ptr<Vector> velocity(raw_v);
  if (status == StatusCode::kOk)
    status = Vector::Create(state->GlobalCount(), &raw_a, d);
  std::unique_ptr<Vector> acceleration(raw_a);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<Vector> internal_velocity, internal_acceleration;
  const int internal_count = context->dofs->InternalDofCount();
  if (internal_count > 0) {
    Vector* raw = nullptr;
    status = Vector::Create(internal_count, &raw, d);
    internal_velocity.reset(raw);
    if (status == StatusCode::kOk) {
      status = Vector::Create(internal_count, &raw, d);
      internal_acceleration.reset(raw);
    }
    if (status != StatusCode::kOk) return status;
    const auto* internal = state->BeamInternals();
    if (!internal || (has_motion && (!internal->CommittedVelocity() ||
                                     !internal->CommittedAcceleration())))
      return model_storage::Invalid("Finite Beam committed joint motion", d);
    if (has_motion)
      status = internal_velocity->CopyFrom(internal->CommittedVelocity(), d);
    if (status != StatusCode::kOk) return status;
  }
  if (has_motion)
    status = velocity->CopyFrom(state->Dynamics()->CommittedVelocity(), d);
  else if (initial_velocity_)
    status = velocity->CopyFrom(initial_velocity_, d);
  if (status != StatusCode::kOk) return status;
  for (int i = 0; i < context->dofs->GlobalCount(); ++i)
    if (context->dofs->Equation(i) < 0 &&
        (velocity->Data()[i] != 0 ||
         (has_motion &&
          state->Dynamics()->CommittedAcceleration()->Data()[i] != 0))) {
      if (d) {
        d->node_id = context->dofs->NodeId(i / 6);
        d->condition_dof = i % 6 + 1;
      }
      return ReportFailure(initial_velocity_ ? StatusCode::kInputError
                                             : StatusCode::kUnsupported,
                           DiagnosticStage::kValidation,
                           "Nonzero fixed/excluded finite initial motion", d);
    }
  const auto calls = PardisoSolveCallCount();
  if (continuous) {
    status =
        acceleration->CopyFrom(state->Dynamics()->CommittedAcceleration(), d);
    if (status == StatusCode::kOk && internal_count > 0)
      status = internal_acceleration->CopyFrom(
          state->BeamInternals()->CommittedAcceleration(), d);
  } else if (internal_count > 0) {
    status = ComputeInitialBeamAcceleration(
        context->model, context->dofs, state, velocity.get(),
        internal_velocity.get(), external.get(), solver_, acceleration.get(),
        internal_acceleration.get(), &failure, d);
  } else {
    status = InitialAcceleration(context, velocity.get(), mode_, external.get(),
                                 solver_, acceleration.get(), d);
  }
  context->report->linear_solver_calls +=
      static_cast<int>(PardisoSolveCallCount() - calls);
  if (status == StatusCode::kOk)
    status =
        state->StageDynamicKinematics(velocity.get(), acceleration.get(), 0,
                                      context->total_time_before_step, d);
  if (status == StatusCode::kOk && internal_count > 0)
    status = state->StageBeamInternalKinematics(internal_velocity.get(),
                                                internal_acceleration.get(), d);
  if (status != StatusCode::kOk) return status;
  Vector* raw_r = nullptr;
  status = Vector::Create(state->GlobalCount(), &raw_r, d);
  std::unique_ptr<Vector> residual(raw_r);
  std::unique_ptr<Vector> internal_residual;
  if (status == StatusCode::kOk && internal_count > 0) {
    Vector* raw = nullptr;
    status = Vector::Create(internal_count, &raw, d);
    internal_residual.reset(raw);
  }
  if (status != StatusCode::kOk) return status;
  NonlinearResidualInput input;
  input.workspaces = workers_;
  input.workspace_count = worker_count_;
  input.model = context->model;
  input.dofs = context->dofs;
  input.state = state;
  input.controls = context->controls;
  input.external = external.get();
  input.force_origins = force_origins.get();
  input.moment_origins = moment_origins.get();
  NonlinearResidualOutput output;
  FieldResidual force, moment;
  output.residual = residual.get();
  output.reaction = state->CandidateReaction();
  output.force = &force;
  output.moment = &moment;
  output.internal_residual = internal_residual.get();
  status = EvaluateFiniteHhtResidual(&input, external.get(), 0, mode_, &output,
                                     &failure, d);
  if (status != StatusCode::kOk) return status;
  context->report->characteristic_length = output.characteristic_length;
  return state->PrepareDynamicConvergenceHistory(
      context->conditions->StepNumber(), context->controls, true,
      force.maximum_flux, moment.maximum_flux, d);
}
StatusCode NonlinearDynamicAnalysis::AdvanceCandidate(
    const DynamicStepContext* context, double size, double time, double total,
    Diagnostic* d) {
  ReleaseAttemptResources();
  context_ = context;
  numerical_diagnostic_ = d;
  size_ = size;
  time_ = time;
  total_ = total;
  newton_report_->Reset();
  auto status = context->state->PrepareCandidate(d);
  bool branch = false;
  if (status == StatusCode::kOk)
    status = StageFiniteHhtKinematics(context->dofs, Integrator(), size, time,
                                      total, true, context->state, &branch, d);
  if (status == StatusCode::kOk)
    status = ReadLoads(context, time, total, &external_, &force_origins_,
                       &moment_origins_, d);
  if (status == StatusCode::kOk) {
    Vector* old_origins = nullptr;
    status = EvaluateDynamicClLoads(context->conditions, context->dofs,
                                    context->state->Dynamics()->StepTime(),
                                    context->state->Dynamics()->TotalTime(),
                                    &old_external_, &old_origins, d);
    delete old_origins;
  }
  bool ran_newton = false;
  if (status == StatusCode::kOk) {
    ran_newton = true;
    status = newton_->Solve(this, solver_, context->controls,
                            context->controls->newton, newton_report_, this, d,
                            Settings()->IsAutomatic()
                                ? NewtonExecution::kAutomatic
                                : NewtonExecution::kDirect);
    UpdateReport();
  }
  if (ran_newton) {
    // Scalar metadata only; no physical state/buffer copy or fake iteration.
    *context->newton = *newton_report_->Info();
  } else {
    context->newton->status = context->newton->numerical_status = status;
    context->newton->outcome =
        status == StatusCode::kEquilibriumFailure && branch
            ? NewtonOutcome::kRotationBranch : NewtonOutcome::kFatal;
  }
  ReleaseAttemptResources();
  context_ = nullptr;
  numerical_diagnostic_ = nullptr;
  return status;
}
StatusCode NonlinearDynamicAnalysis::EvaluateResidual(
    int completed, FieldResidual* force, FieldResidual* moment, double* length,
    NonlinearFailure* failure, Diagnostic* d) {
  *failure = NonlinearFailure::kNone;
  bool branch = false;
  auto* state = context_->state;
  auto status =
      StageFiniteHhtKinematics(context_->dofs, Integrator(), size_, time_,
                               total_, false, state, &branch, d);
  if (status != StatusCode::kOk) {
    if (status == StatusCode::kEquilibriumFailure && branch)
      *failure = NonlinearFailure::kRotationBranch;
    return status;
  }
  status = state->PrepareDynamicConvergenceHistory(
      context_->conditions->StepNumber(), context_->controls, false, 0, 0, d);
  delete residual_;
  residual_ = nullptr;
  delete internal_residual_;
  internal_residual_ = nullptr;
  if (status == StatusCode::kOk)
    status = Vector::Create(state->GlobalCount(), &residual_, d);
  if (status == StatusCode::kOk && context_->dofs->InternalDofCount() > 0)
    status = Vector::Create(context_->dofs->InternalDofCount(),
                            &internal_residual_, d);
  if (status != StatusCode::kOk) return status;
  NonlinearResidualInput input;
  input.workspaces = workers_;
  input.workspace_count = worker_count_;
  input.model = context_->model;
  input.dofs = context_->dofs;
  input.state = state;
  input.controls = context_->controls;
  input.external = external_;
  input.force_origins = force_origins_;
  input.moment_origins = moment_origins_;
  input.actual_correction = completed ? correction_ : nullptr;
  input.internal_actual_correction = completed ? internal_correction_ : nullptr;
  input.previous_force_max_mean = state->FieldHistory()->PreviousMaximumMean(0);
  input.previous_moment_max_mean =
      state->FieldHistory()->PreviousMaximumMean(1);
  NonlinearResidualOutput output;
  output.residual = residual_;
  output.reaction = state->CandidateReaction();
  output.force = force;
  output.moment = moment;
  output.internal_residual = internal_residual_;
  status = EvaluateFiniteHhtResidual(&input, old_external_,
                                     Settings()->IntegratorParameters()->alpha,
                                     mode_, &output, failure, d);
  if (status != StatusCode::kOk) return status;
  *length = length_ = output.characteristic_length;
  status = state->PreviewDynamicFluxHistory(
      force->spatial_mean, moment->spatial_mean, force->maximum_flux,
      moment->maximum_flux, &force->raw_mean, &moment->raw_mean, d);
  return status == StatusCode::kOk
             ? ConvertNonlinearFieldMeans(context_->controls, length_, force,
                                          moment, d)
             : status;
}
StatusCode NonlinearDynamicAnalysis::AssembleTangent(
    const EquationSystem** system, Vector** correction,
    NonlinearFailure* failure, Diagnostic* d) {
  *system = nullptr;
  *correction = nullptr;
  Matrix *raw_h = nullptr, *raw_b = nullptr;
  auto status = Matrix::Create(context_->dofs->NodeCount(), 9, &raw_h, d);
  std::unique_ptr<Matrix> h(raw_h);
  if (status == StatusCode::kOk)
    status = Matrix::Create(context_->dofs->NodeCount(), 9, &raw_b, d);
  std::unique_ptr<Matrix> b(raw_b);
  if (status != StatusCode::kOk) return status;
  RotationalLinearization rotation;
  rotation.velocity_per_rotation = h.get();
  rotation.acceleration_per_rotation = b.get();
  TranslationalLinearization translation;
  bool branch = false;
  status = EvaluateFiniteHhtKinematicLinearization(
      context_->dofs, Integrator(), size_, context_->state, &translation,
      &rotation, &branch, d);
  if (status != StatusCode::kOk) {
    if (status == StatusCode::kEquilibriumFailure && branch)
      *failure = NonlinearFailure::kRotationBranch;
    return status;
  }
  FiniteHhtLinearization linearization;
  linearization.alpha = Settings()->IntegratorParameters()->alpha;
  linearization.mode = mode_;
  linearization.translation = &translation;
  linearization.rotation = &rotation;
  status = AssembleFiniteHhtTangent(
      context_->model, context_->dofs, context_->state, residual_,
      &linearization, &equations_, failure, d, internal_residual_,
      workers_, worker_count_);
  delete residual_;
  residual_ = nullptr;
  delete internal_residual_;
  internal_residual_ = nullptr;
  if (status == StatusCode::kOk && context_->dofs->EquationCount() > 0)
    status =
        Vector::Create(context_->dofs->EquationCount(), &free_correction_, d);
  if (status != StatusCode::kOk) {
    ReleaseEquationSystem();
    return status;
  }
  *system = equations_;
  *correction = free_correction_;
  return StatusCode::kOk;
}
StatusCode NonlinearDynamicAnalysis::ApplyCorrection(const Vector* correction,
                                                     NonlinearFailure* failure,
                                                     Diagnostic* d) {
  *failure = NonlinearFailure::kNone;
  ReleaseCorrectionAssessment();
  auto status = Vector::Create(context_->dofs->GlobalCount(), &correction_, d);
  if (status != StatusCode::kOk) return status;
  std::vector<int> indices(context_->dofs->FreeCount());
  for (int i = 0; i < context_->dofs->FreeCount(); ++i)
    indices[i] = context_->dofs->GlobalDofForEquation(i);
  if (context_->dofs->FreeCount() > 0) {
    Vector* raw = nullptr;
    status = Vector::Create(context_->dofs->FreeCount(), &raw, d);
    std::unique_ptr<Vector> physical(raw);
    std::vector<int> equations(context_->dofs->FreeCount());
    for (int i = 0; i < context_->dofs->FreeCount(); ++i) equations[i] = i;
    if (status == StatusCode::kOk)
      status = physical->GatherFrom(correction, equations.data(),
                                    physical->Size(), 0, d);
    if (status == StatusCode::kOk)
      status = correction_->ScatterAdd(physical.get(), indices.data(),
                                       physical->Size(), d);
  }
  if (status == StatusCode::kOk && context_->dofs->InternalDofCount() > 0) {
    const int count = context_->dofs->InternalDofCount();
    status = Vector::Create(count, &internal_correction_, d);
    std::vector<int> equations(count);
    for (int i = 0; i < count; ++i)
      equations[i] = context_->dofs->FreeCount() + i;
    if (status == StatusCode::kOk)
      status = internal_correction_->GatherFrom(correction, equations.data(),
                                                count, 0, d);
  }
  if (status != StatusCode::kOk) return status;
  bool branch = false;
  status = context_->state->ApplyNonlinearCorrection(correction_, &branch, d);
  if (status == StatusCode::kOk && internal_correction_)
    status =
        context_->state->ApplyBeamInternalCorrection(internal_correction_, d);
  if (status == StatusCode::kEquilibriumFailure && branch)
    *failure = NonlinearFailure::kRotationBranch;
  return status;
}
StatusCode NonlinearDynamicAnalysis::OnIteration(const NewtonReport* report,
                                                 Diagnostic* d) {
  context_->attempt->characteristic_length = length_;
  const auto* info = report->Info();
  return context_->progress->OnIncrementProgress(
      IncrementEvent::kIteration, context_->attempt, info,
      info->residual_valid ? report->Displacement() : nullptr,
      info->residual_valid ? report->Rotation() : nullptr,
      info->numerical_status != StatusCode::kOk ? numerical_diagnostic_
                                                : nullptr,
      d);
}
void NonlinearDynamicAnalysis::UpdateReport() noexcept {
  const auto* info = newton_report_->Info();
  auto* report = context_->report;
  report->iteration_count = info->iteration;
  report->linear_solver_calls += info->linear_solver_calls;
  report->assembly_seconds += info->tangent_seconds;
  report->solve_seconds += info->solve_seconds;
  report->equilibrium_recovery_seconds +=
      info->residual_seconds + info->criteria_seconds;
  report->residual_valid = info->residual_valid;
  report->characteristic_length = length_;
  report->displacement_field =
      info->residual_valid ? newton_report_->Displacement() : nullptr;
  report->rotation_field =
      info->residual_valid ? newton_report_->Rotation() : nullptr;
  report->logging_diagnostic = newton_report_->LoggingDiagnostic();
}
StatusCode NonlinearDynamicAnalysis::RecoverCandidate(
    const DynamicStepContext* context, Diagnostic* d) {
  auto* state = context->state;
  auto status = state->PrepareNonlinearOutput(d);
  NonlinearFailure failure = NonlinearFailure::kNone;
  if (status == StatusCode::kOk) {
    state->MarkResultsUnconfirmed(context->conditions->StepNumber(),
                                  context->report->increment_number);
    status = RecoverNonlinearResults(context->model, context->dofs, state,
                                     workers_, worker_count_, &failure, d);
  }
  context->report->element_dispatch_seconds =
      ElementDispatchSeconds() - element_started_;
  context->report->observed_element_team_max = ObservedElementTeamMaximum();
  return status;
}
void NonlinearDynamicAnalysis::ReleaseEquationSystem() noexcept {
  delete equations_;
  equations_ = nullptr;
}
void NonlinearDynamicAnalysis::ReleaseFreeCorrection() noexcept {
  delete free_correction_;
  free_correction_ = nullptr;
}
void NonlinearDynamicAnalysis::ReleaseCorrectionAssessment() noexcept {
  delete correction_;
  correction_ = nullptr;
  delete internal_correction_;
  internal_correction_ = nullptr;
}
void NonlinearDynamicAnalysis::ReleaseAttemptResources() noexcept {
  ReleaseEquationSystem();
  ReleaseFreeCorrection();
  ReleaseCorrectionAssessment();
  delete residual_;
  residual_ = nullptr;
  delete internal_residual_;
  internal_residual_ = nullptr;
  delete external_;
  external_ = nullptr;
  delete old_external_;
  old_external_ = nullptr;
  delete force_origins_;
  force_origins_ = nullptr;
  delete moment_origins_;
  moment_origins_ = nullptr;
}
void NonlinearDynamicAnalysis::ReleaseStepResources() noexcept {
  for (int i = 0; i < worker_count_; ++i) delete workers_[i];
  delete[] workers_;
  workers_ = nullptr;
  worker_count_ = 0;
  ReleaseAttemptResources();
  context_ = nullptr;
  numerical_diagnostic_ = nullptr;
}
}  // namespace fesa
