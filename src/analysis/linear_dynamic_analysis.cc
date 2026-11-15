#include "fesa/analysis/linear_dynamic_analysis.h"

#include <algorithm>
#include <chrono>
#include <climits>
#include <memory>
#include <new>
#include <vector>

#include "adapters/mkl/pardiso_solver.h"
#include "adapters/openmp/thread_runtime.h"
#include "analysis/dynamic_initialization.h"
#include "analysis/dynamic_load.h"
#include "analysis/equilibrium_check.h"
#include "analysis/field_measurement.h"
#include "analysis/linear_dynamic_element.h"
#include "analysis/linear_result_recovery.h"
#include "analysis/solid_dynamic_reaction.h"
#include "analysis/static_assembler.h"
#include "fesa/analysis/dynamic_state.h"
#include "fesa/analysis/nonlinear_field_history.h"
#include "fesa/elements/solid_element.h"
#include "model/model_storage.h"

namespace fesa {
namespace {
StatusCode NewVector(int size, std::unique_ptr<Vector>* owner, Diagnostic* d) {
  Vector* raw = nullptr;
  const auto status = Vector::Create(size, &raw, d);
  owner->reset(raw);
  return status;
}
double Seconds() {
  return std::chrono::duration<double>(
             std::chrono::steady_clock::now().time_since_epoch())
      .count();
}
}  // namespace
StatusCode LinearDynamicAnalysis::Create(const DynamicStepSettings* settings,
                                         MassMode mode,
                                         const EffectiveConditions* previous,
                                         const Vector* initial_velocity,
                                         LinearDynamicAnalysis** out,
                                         Diagnostic* d) {
  if (!out) return model_storage::Invalid("LinearDynamicAnalysis output", d);
  *out = nullptr;
  if (mode != MassMode::kConsistent && mode != MassMode::kRowSum &&
      mode != MassMode::kDefault)
    return model_storage::Invalid("Dynamic mass mode", d);
  try {
    std::unique_ptr<LinearDynamicAnalysis> result(new LinearDynamicAnalysis);
    auto status = result->InitializeIntegration(settings, d);
    if (status != StatusCode::kOk) return status;
    result->solver_ = new PardisoSolver;
    status = LinearAlgorithm::Create(&result->algorithm_, d);
    if (status != StatusCode::kOk) return status;
    result->mode_ = mode;
    result->previous_ = previous;
    result->initial_velocity_ = initial_velocity;
    *out = result.release();
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return model_storage::AllocationFailure("LinearDynamicAnalysis", d);
  }
}
LinearDynamicAnalysis::~LinearDynamicAnalysis() {
  ReleaseStepResources();
  delete force_;
  delete moment_;
  delete algorithm_;
  delete solver_;
}
void LinearDynamicAnalysis::ReleaseEquationSystem() noexcept {
  delete equations_;
  equations_ = nullptr;
}
void LinearDynamicAnalysis::ReleaseSolution() noexcept {
  if (workspace_) workspace_->ReleaseReducedBuffers();
}
void LinearDynamicAnalysis::ReleaseAttemptResources() noexcept {
  ReleaseEquationSystem();
  delete workspace_;
  workspace_ = nullptr;
  delete predictor_;
  predictor_ = nullptr;
  delete velocity_;
  velocity_ = nullptr;
  delete acceleration_;
  acceleration_ = nullptr;
  delete external_;
  external_ = nullptr;
  delete external_origins_;
  external_origins_ = nullptr;
  delete weighted_;
  weighted_ = nullptr;
}
void LinearDynamicAnalysis::ReleaseStepResources() noexcept {
  ReleaseAttemptResources();
  for (int i = 0; i < worker_count_; ++i) delete workers_[i];
  delete[] workers_;
  workers_ = nullptr;
  worker_count_ = 0;
  for (int i = 0; i < reference_count_; ++i) delete references_[i];
  delete[] references_;
  references_ = nullptr;
  reference_count_ = 0;
  context_ = nullptr;
}
StatusCode LinearDynamicAnalysis::BeginStep(const DynamicStepContext* context,
                                            Diagnostic* d) {
  ReleaseStepResources();
  delete force_;
  force_ = nullptr;
  delete moment_;
  moment_ = nullptr;
  if (context->model->Conditions() != context->conditions ||
      context->model->Elements()->empty() ||
      context->model->Elements()->size() > INT_MAX)
    return model_storage::Invalid("Dynamic model/condition binding", d);
  auto status = context->state->ValidateMapping(context->dofs, d);
  if (status == StatusCode::kOk)
    status = ValidateStep(context->model, mode_, previous_, d);
  if (status != StatusCode::kOk) return status;
  shell_ = context->model->Elements()->begin()->second->family() ==
           ElementFamily::kShell;
  if (shell_ && mode_ == MassMode::kDefault) mode_ = MassMode::kRowSum;
  if (!shell_) {
    std::vector<SolidFormulation> forms;
    for (const auto& entry : *context->model->Elements()) {
      const auto form =
          static_cast<const SolidElement*>(entry.second)->formulation();
      if (std::find(forms.begin(), forms.end(), form) == forms.end())
        forms.push_back(form);
    }
    const bool reduced =
        std::find(forms.begin(), forms.end(), SolidFormulation::kHexa8Reduced) !=
        forms.end();
    const int table_count =
        2 * static_cast<int>(forms.size()) + (reduced ? 1 : 0);
    references_ = new SolidReferenceData* [table_count] {};
    reference_count_ = table_count;
    for (int i = 0; i < static_cast<int>(forms.size()); ++i) {
      status = SolidReferenceData::Create(forms[i],
                                          IntegrationPurpose::kStiffness,
                                          &references_[2 * i], d);
      if (status == StatusCode::kOk)
        status = SolidReferenceData::Create(forms[i], IntegrationPurpose::kMass,
                                            &references_[2 * i + 1], d);
      if (status != StatusCode::kOk) return status;
    }
    if (reduced) {
      status = SolidReferenceData::Create(
          SolidFormulation::kHexa8Reduced, IntegrationPurpose::kGeometry,
          &references_[reference_count_ - 1], d);
      if (status != StatusCode::kOk) return status;
    }
  }
  const int count = ConfiguredThreadCount();
  workers_ = new ElementWorkspace* [count] {};
  worker_count_ = count;
  for (int i = 0; !shell_ && i < count; ++i) {
    status = ElementWorkspace::CreateForSolidCatalog(
        references_, reference_count_, &workers_[i], d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}
StatusCode LinearDynamicAnalysis::ValidateStep(
    const AnalysisModel* model, MassMode mode,
    const EffectiveConditions* previous, Diagnostic* d) {
  if (!model || !model->Conditions() || model->Elements()->empty())
    return model_storage::Invalid("Dynamic preflight model", d);
  auto status = ValidateDynamicConstraintTransition(previous, model->Conditions(), d);
  if (status != StatusCode::kOk) return status;
  bool shell = false;
  return ValidateLinearDynamicModel(model, mode, &shell, d);
}
StatusCode LinearDynamicAnalysis::InitializeCandidate(
    const DynamicStepContext* context, Diagnostic* d) {
  Vector *raw = nullptr, *raw_origins = nullptr;
  auto status = EvaluateDynamicClLoads(context->conditions, context->dofs, 0,
                                       context->total_time_before_step, &raw,
                                       &raw_origins, d);
  std::unique_ptr<Vector> external(raw), origins(raw_origins);
  if (status != StatusCode::kOk) return status;
  auto* state = context->state;
  const bool continuous =
      state->Dynamics() && state->Dynamics()->CommittedVelocity();
  auto start = continuous                 ? SolidDynamicStart::kContinuous
               : state->StepNumber() == 0 ? SolidDynamicStart::kDirect
                                          : SolidDynamicStart::kStaticPreload;
  if (continuous) {
    const auto* old_step =
        previous_->SourceDomain()->Steps()->at(previous_->StepNumber());
    Vector *old_raw = nullptr, *old_origins = nullptr;
    status = EvaluateDynamicClLoads(
        previous_, context->dofs, old_step->Period(),
        context->total_time_before_step, &old_raw, &old_origins, d);
    std::unique_ptr<Vector> old(old_raw), old_origin_owner(old_origins);
    if (status != StatusCode::kOk) return status;
    status = old->AddScaled(-1, external.get(), d);
    double difference = 0;
    if (status == StatusCode::kOk) status = old->Norm2(&difference, d);
    if (status != StatusCode::kOk) return status;
    if (difference != 0) start = SolidDynamicStart::kLoadJump;
  }
  status = shell_ ? state->EnableLinearDynamics(d)
                  : state->EnableTranslationalDynamics(d);
  if (status == StatusCode::kOk && !continuous)
    status = state->EnableDynamicConvergenceHistory(d);
  const auto calls = PardisoSolveCallCount();
  if (status == StatusCode::kOk && shell_) {
    Vector *rate_raw = nullptr, *origin_rate_raw = nullptr;
    if (start != SolidDynamicStart::kContinuous)
      status = EvaluateDynamicClLoadRates(context->conditions, context->dofs, 0,
                                          context->total_time_before_step,
                                          &rate_raw, &origin_rate_raw, d);
    std::unique_ptr<Vector> rate(rate_raw), origin_rate(origin_rate_raw);
    if (status == StatusCode::kOk)
      status = PrepareInitialShellDynamicState(
          start, context->model, context->dofs, context->conditions, mode_,
          external.get(), rate_raw, initial_velocity_, solver_, worker_count_,
          context->total_time_before_step, state, d);
  } else {
    if (status == StatusCode::kOk)
      status = PrepareInitialSolidDynamicState(
          start, context->model, context->dofs, context->conditions, mode_,
          external.get(), initial_velocity_, solver_, workers_, worker_count_,
          context->total_time_before_step, state, d);
  }
  context->report->linear_solver_calls =
      static_cast<int>(PardisoSolveCallCount() - calls);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<Vector> residual;
  status = NewVector(context->dofs->GlobalCount(), &residual, d);
  if (status != StatusCode::kOk) return status;
  HhtSolidForceInput input;
  input.alpha = 0;
  input.previous_displacement = state->CandidateDisplacement();
  input.external_origins = origins.get();
  Vector* all_raw = nullptr;
  status = ComputeSolidHhtResidual(
      context->model, context->dofs, mode_, state->CandidateDisplacement(),
      state->Dynamics()->CandidateVelocity(),
      state->Dynamics()->CandidateAcceleration(), external.get(), &input,
      workers_, worker_count_, residual.get(), &all_raw, d);
  std::unique_ptr<Vector> all(all_raw);
  Vector *force_raw = nullptr, *moment_raw = nullptr;
  if (status == StatusCode::kOk)
    status =
        SplitLinearDynamicFluxOrigins(context->model, context->conditions,
                                      all.get(), &force_raw, &moment_raw, d);
  std::unique_ptr<Vector> force_values(force_raw), moment_values(moment_raw);
  double maxima[2] = {};
  Vector* values[] = {force_raw, moment_raw};
  for (int field = 0; field < 2 && status == StatusCode::kOk; ++field) {
    if (!values[field]) continue;
    int index = 0;
    status = values[field]->MaximumAbs(&maxima[field], &index, d);
  }
  return status == StatusCode::kOk
             ? state->PrepareDynamicConvergenceHistory(
                   context->conditions->StepNumber(), context->controls, true,
                   maxima[0], maxima[1], d)
             : status;
}
StatusCode LinearDynamicAnalysis::AdvanceCandidate(
    const DynamicStepContext* context, double size, double time, double total,
    Diagnostic* d) {
  ReleaseAttemptResources();
  context_ = context;
  size_ = size;
  time_ = time;
  total_ = total;
  const int count = context->dofs->GlobalCount();
  auto* state = context->state;
  auto status = state->PrepareCandidate(d);
  if (status == StatusCode::kOk)
    status = StaticWorkspace::Create(context->dofs, &workspace_, d);
  if (status == StatusCode::kOk) status = Vector::Create(count, &predictor_, d);
  if (status == StatusCode::kOk) status = Vector::Create(count, &velocity_, d);
  if (status == StatusCode::kOk)
    status = Vector::Create(count, &acceleration_, d);
  if (status == StatusCode::kOk) status = Vector::Create(count, &weighted_, d);
  if (status == StatusCode::kOk) {
    TranslationalStateView committed;
    committed.displacement = state->CommittedDisplacement();
    committed.velocity = state->Dynamics()->CommittedVelocity();
    committed.acceleration = state->Dynamics()->CommittedAcceleration();
    TranslationalCandidate candidate;
    candidate.displacement = state->CandidateDisplacement();
    candidate.velocity = velocity_;
    candidate.acceleration = acceleration_;
    status = Integrator()->PrepareKinematics(&committed, size, &candidate, d);
    if (status == StatusCode::kOk)
      status = predictor_->CopyFrom(candidate.displacement, d);
  }
  if (status == StatusCode::kOk)
    status = EvaluateDynamicClLoads(context->conditions, context->dofs, time,
                                    total, &external_, &external_origins_, d);
  if (status == StatusCode::kOk) {
    Vector *old_raw = nullptr, *old_origins = nullptr;
    status = EvaluateDynamicClLoads(
        context->conditions, context->dofs, state->Dynamics()->StepTime(),
        state->Dynamics()->TotalTime(), &old_raw, &old_origins, d);
    std::unique_ptr<Vector> old(old_raw), old_origin_owner(old_origins);
    const double alpha = Settings()->IntegratorParameters()->alpha;
    if (status == StatusCode::kOk) status = weighted_->CopyFrom(external_, d);
    if (status == StatusCode::kOk) status = weighted_->Scale(1 + alpha, d);
    if (status == StatusCode::kOk)
      status = weighted_->AddScaled(-alpha, old.get(), d);
  }
  const auto calls = PardisoSolveCallCount();
  if (status == StatusCode::kOk) {
    status = algorithm_->Solve(this, solver_, d);
    context->report->solve_seconds += algorithm_->SolveSeconds();
  }
  context->report->linear_solver_calls +=
      static_cast<int>(PardisoSolveCallCount() - calls);
  ReleaseAttemptResources();
  context_ = nullptr;
  return status;
}
StatusCode LinearDynamicAnalysis::AssembleSystem(const EquationSystem** system,
                                                 Vector** solution,
                                                 bool* solve_required,
                                                 Diagnostic* d) {
  *system = nullptr;
  *solution = nullptr;
  *solve_required = false;
  TranslationalLinearization linearization;
  auto status =
      Integrator()->EvaluateKinematicLinearization(size_, &linearization, d);
  if (status != StatusCode::kOk) return status;
  HhtSolidAssemblyInput input;
  input.committed = context_->state->CommittedDisplacement();
  input.predictor = predictor_;
  input.external = weighted_;
  input.alpha = Settings()->IntegratorParameters()->alpha;
  input.c0 = linearization.acceleration_per_displacement;
  input.c1 = linearization.velocity_per_displacement;
  input.previous_velocity = context_->state->Dynamics()->CommittedVelocity();
  input.predicted_velocity = velocity_;
  input.mode = mode_;
  const double started = Seconds();
  status = AssembleReducedSystem(context_->model, context_->dofs,
                                 context_->conditions, workspace_, &equations_,
                                 d, workers_, worker_count_, &input);
  context_->report->assembly_seconds += Seconds() - started;
  if (status != StatusCode::kOk) return status;
  *system = equations_;
  *solution = workspace_->ReducedSolution();
  *solve_required = context_->dofs->FreeCount() > 0;
  return StatusCode::kOk;
}
StatusCode LinearDynamicAnalysis::ApplySolution(const Vector* solution,
                                                Diagnostic* d) {
  const auto* dofs = context_->dofs;
  std::vector<int> indices(dofs->FreeCount());
  for (int i = 0; i < dofs->FreeCount(); ++i)
    indices[i] = dofs->GlobalDofForEquation(i);
  std::unique_ptr<Vector> old;
  auto status = NewVector(dofs->FreeCount(), &old, d);
  auto* displacement = context_->state->CandidateDisplacement();
  if (status == StatusCode::kOk)
    status =
        old->GatherFrom(displacement, indices.data(), dofs->FreeCount(), 0, d);
  if (status == StatusCode::kOk) status = old->Scale(-1, d);
  if (status == StatusCode::kOk)
    status = displacement->ScatterAdd(old.get(), indices.data(),
                                      dofs->FreeCount(), d);
  return status == StatusCode::kOk
             ? displacement->ScatterAdd(solution, indices.data(),
                                        dofs->FreeCount(), d)
             : status;
}
StatusCode LinearDynamicAnalysis::EvaluateResidual(Diagnostic* d) {
  auto* state = context_->state;
  TranslationalStateView committed;
  committed.displacement = state->CommittedDisplacement();
  committed.velocity = state->Dynamics()->CommittedVelocity();
  committed.acceleration = state->Dynamics()->CommittedAcceleration();
  TranslationalCandidate candidate;
  candidate.displacement = state->CandidateDisplacement();
  candidate.velocity = velocity_;
  candidate.acceleration = acceleration_;
  auto status =
      Integrator()->UpdateKinematics(&committed, size_, &candidate, d);
  if (status == StatusCode::kOk)
    status = state->StageDynamicKinematics(velocity_, acceleration_, time_,
                                           total_, d);
  if (status == StatusCode::kOk)
    status = state->PrepareDynamicConvergenceHistory(
        context_->conditions->StepNumber(), context_->controls, false, 0, 0, d);
  HhtSolidForceInput input;
  input.previous_displacement = committed.displacement;
  input.previous_velocity = committed.velocity;
  input.external_origins = external_origins_;
  input.alpha = Settings()->IntegratorParameters()->alpha;
  Vector* raw = nullptr;
  if (status == StatusCode::kOk)
    status = ComputeSolidHhtResidual(
        context_->model, context_->dofs, mode_, candidate.displacement,
        velocity_, acceleration_, weighted_, &input, workers_, worker_count_,
        workspace_->InternalForce(), &raw, d);
  std::unique_ptr<Vector> origins(raw);
  if (status != StatusCode::kOk) return status;
  delete force_;
  force_ = nullptr;
  delete moment_;
  moment_ = nullptr;
  force_ = new FieldResidual;
  moment_ = new FieldResidual;
  Vector *force_raw = nullptr, *moment_raw = nullptr;
  status =
      SplitLinearDynamicFluxOrigins(context_->model, context_->conditions,
                                    origins.get(), &force_raw, &moment_raw, d);
  std::unique_ptr<Vector> force_values(force_raw), moment_values(moment_raw);
  Vector* values[] = {force_raw, moment_raw};
  FieldResidual* fields[] = {force_, moment_};
  const FieldControls* controls[] = {context_->controls->displacement,
                                     context_->controls->rotation};
  for (int field = 0; field < 2 && status == StatusCode::kOk; ++field) {
    status = ComputeSpatialFluxMean(
        values[field], controls[field],
        state->Dynamics()->ConvergenceHistory()->PreviousMaximumMean(field),
        field, fields[field], d);
    if (status == StatusCode::kOk)
      status = MeasureLinearDynamicField(
          context_->dofs, field, workspace_->InternalForce(),
          committed.displacement, candidate.displacement, predictor_,
          fields[field], d);
  }
  if (status == StatusCode::kOk)
    status = state->PreviewDynamicFluxHistory(
        force_->spatial_mean, moment_->spatial_mean, force_->maximum_flux,
        moment_->maximum_flux, &force_->raw_mean, &moment_->raw_mean, d);
  std::unique_ptr<Vector> lengths;
  const int count = static_cast<int>(context_->model->Elements()->size());
  if (status == StatusCode::kOk) status = NewVector(count, &lengths, d);
  int index = 0;
  for (const auto& entry : *context_->model->Elements()) {
    if (status != StatusCode::kOk) break;
    double length = 0;
    status = entry.second->CharacteristicLength(&length, d);
    if (status == StatusCode::kOk) status = lengths->Set(index++, length, d);
  }
  double length = 0;
  if (status == StatusCode::kOk) status = lengths->AbsoluteSum(&length, d);
  if (status != StatusCode::kOk) return status;
  length /= count;
  force_->used_mean = force_->raw_mean;
  moment_->used_mean = moment_->raw_mean;
  auto* report = context_->report;
  report->characteristic_length = length;
  report->displacement_field = force_;
  report->rotation_field = moment_;
  report->residual_valid = true;
  return EvaluateJointFieldCriteria(context_->controls, length, force_, moment_,
                                    d);
}
StatusCode LinearDynamicAnalysis::RecoverCandidate(
    const DynamicStepContext* context, Diagnostic* d) {
  const double started = Seconds();
  auto* state = context->state;
  Vector *raw = nullptr, *raw_origins = nullptr;
  auto status = EvaluateDynamicClLoads(context->conditions, context->dofs,
                                       state->Dynamics()->CandidateStepTime(),
                                       state->Dynamics()->CandidateTotalTime(),
                                       &raw, &raw_origins, d);
  std::unique_ptr<Vector> external(raw), origins(raw_origins);
  if (status == StatusCode::kOk)
    status = ComputeSolidDynamicReaction(
        context->model, context->dofs, mode_, state->CandidateDisplacement(),
        state->Dynamics()->CandidateVelocity(),
        state->Dynamics()->CandidateAcceleration(), external.get(), workers_,
        worker_count_, state->CandidateReaction(), d);
  if (status == StatusCode::kOk)
    status = RecoverLinearResults(
        context->model, context->dofs, state, context->conditions->StepNumber(),
        context->report->increment_number, workers_, worker_count_, d);
  context->report->equilibrium_recovery_seconds += Seconds() - started;
  return status;
}
}  // namespace fesa
