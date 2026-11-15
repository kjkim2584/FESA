#include "fesa/analysis/nonlinear_static_analysis.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <memory>
#include <new>
#include <utility>
#include <vector>

#include "adapters/mkl/pardiso_solver.h"
#include "adapters/openmp/thread_runtime.h"
#include "analysis/increment_control.h"
#include "analysis/nonlinear_assembler.h"
#include "analysis/nonlinear_scaling.h"
#include "analysis/nonlinear_step_conditions.h"
#include "base/diagnostic_copy.h"
#include "elements/solid_element_read_view.h"
#include "fesa/analysis/beam_internal_state.h"
#include "fesa/analysis/incremental_solid_history.h"
#include "fesa/analysis/nonlinear_field_history.h"
#include "fesa/analysis/rotation_state.h"
#include "fesa/analysis/static_workspace.h"
#include "fesa/elements/beam3d.h"
#include "fesa/math/rotation.h"
#include "model/model_storage.h"

namespace fesa {
namespace {
double Seconds() {
  return std::chrono::duration<double>(
             std::chrono::steady_clock::now().time_since_epoch())
      .count();
}

StatusCode Notify(AnalysisProgress* progress, AnalysisEvent event,
                  const AnalysisReport* report, Diagnostic* d) {
  try {
    return progress->OnProgress(event, report, d);
  } catch (const std::bad_alloc&) {
    return model_storage::AllocationFailure("Nonlinear progress allocation", d);
  }
}
StatusCode PhysicalIncrement(const DofManager* dofs, const AnalysisState* state,
                             Vector* out, NonlinearFailure* failure,
                             Diagnostic* d) {
  auto status = out->CopyFrom(state->CandidateDisplacement(), d);
  if (status == StatusCode::kOk)
    status = out->AddScaled(-1, state->CommittedDisplacement(), d);
  if (status != StatusCode::kOk) return status;
  if (!state->Rotations()) return StatusCode::kOk;
  Vector *raw_ref = nullptr, *raw_q = nullptr, *raw_spin = nullptr;
  status = Vector::Create(4, &raw_ref, d);
  std::unique_ptr<Vector> ref(raw_ref);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(4, &raw_q, d);
  std::unique_ptr<Vector> q(raw_q);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(3, &raw_spin, d);
  std::unique_ptr<Vector> spin(raw_spin);
  if (status != StatusCode::kOk) return status;
  for (int node = 0; node < dofs->NodeCount(); ++node) {
    const int indices[4] = {4 * node, 4 * node + 1, 4 * node + 2, 4 * node + 3};
    status = state->Rotations()->CommittedOrientations()->GatherTo(
        indices, 4, ref.get(), 0, d);
    if (status == StatusCode::kOk)
      status = state->Rotations()->CandidateOrientations()->GatherTo(
          indices, 4, q.get(), 0, d);
    if (status != StatusCode::kOk) return status;
    bool branch = false;
    status =
        SpatialRotationIncrement(ref.get(), q.get(), spin.get(), &branch, d);
    if (status != StatusCode::kOk) {
      if (branch && status == StatusCode::kEquilibriumFailure)
        *failure = NonlinearFailure::kRotationBranch;
      return status;
    }
    const int axes[3] = {0, 1, 2};
    status = out->GatherFrom(spin.get(), axes, 3, 6 * node + 3, d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}
const NonlinearStaticParameters* Parameters(
    const EffectiveConditions* conditions) {
  const StepDefinition* step = nullptr;
  if (conditions->SourceDomain()->FindStepDefinition(
          conditions->StepNumber(), &step, nullptr) != StatusCode::kOk)
    return nullptr;
  return step->NonlinearParameters();
}
}  // namespace

StatusCode NonlinearStaticAnalysis::Create(NonlinearStaticAnalysis** out,
                                           Diagnostic* d,
                                           bool b33_material_first) {
  if (!out) return model_storage::Invalid("Nonlinear Analysis output", d);
  *out = nullptr;
  try {
    std::unique_ptr<NonlinearStaticAnalysis> result(
        new NonlinearStaticAnalysis);
    result->b33_material_first_ = b33_material_first;
    auto status = NewtonRaphson::Create(&result->newton_, d);
    if (status != StatusCode::kOk) return status;
    status = NewtonReport::Create(&result->newton_report_, d);
    if (status != StatusCode::kOk) return status;
    result->solver_ = new PardisoSolver;
    *out = result.release();
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return model_storage::AllocationFailure("Nonlinear Analysis creation", d);
  }
}
NonlinearStaticAnalysis::~NonlinearStaticAnalysis() {
  ClearExecution();
  delete solver_;
  delete newton_report_;
  delete newton_;
}
void NonlinearStaticAnalysis::ReleaseEquationSystem() noexcept {
  delete equations_;
  equations_ = nullptr;
}
void NonlinearStaticAnalysis::ReleaseFreeCorrection() noexcept {
  delete free_correction_;
  free_correction_ = nullptr;
}
void NonlinearStaticAnalysis::ReleaseCorrectionAssessment() noexcept {
  delete correction_;
  correction_ = nullptr;
  delete internal_correction_;
  internal_correction_ = nullptr;
}
void NonlinearStaticAnalysis::ReleaseAttemptResources() noexcept {
  ReleaseEquationSystem();
  ReleaseFreeCorrection();
  ReleaseCorrectionAssessment();
  delete residual_;
  residual_ = nullptr;
  delete internal_residual_;
  internal_residual_ = nullptr;
  delete external_;
  external_ = nullptr;
  delete force_origins_;
  force_origins_ = nullptr;
  delete moment_origins_;
  moment_origins_ = nullptr;
}
void NonlinearStaticAnalysis::ClearExecution() noexcept {
  if (attempt_) attempt_->geometry_diagnostic = nullptr;
  ReleaseAttemptResources();
  for (int i = 0; i < worker_count_; ++i) delete workers_[i];
  delete[] workers_;
  workers_ = nullptr;
  worker_count_ = 0;
  for (int i = 0; i < solid_reference_count_; ++i) delete solid_references_[i];
  delete[] solid_references_;
  solid_references_ = nullptr;
  solid_reference_count_ = 0;
  delete start_u_;
  start_u_ = nullptr;
  delete start_reaction_;
  start_reaction_ = nullptr;
  delete start_q_;
  start_q_ = nullptr;
  delete predictor_;
  predictor_ = nullptr;
  delete next_predictor_;
  next_predictor_ = nullptr;
  delete internal_predictor_;
  internal_predictor_ = nullptr;
  delete next_internal_predictor_;
  next_internal_predictor_ = nullptr;
  delete clock_;
  clock_ = nullptr;
  delete next_clock_;
  next_clock_ = nullptr;
  delete target_;
  target_ = nullptr;
  delete attempt_;
  attempt_ = nullptr;
  delete path_;
  path_ = nullptr;
  model_ = nullptr;
  dofs_ = nullptr;
  conditions_ = nullptr;
  controls_ = nullptr;
  state_ = nullptr;
  report_ = nullptr;
  progress_ = nullptr;
  attempt_diagnostic_ = nullptr;
  length_ = 0;
  element_started_ = 0;
}
StatusCode NonlinearStaticAnalysis::PrepareSolidWorkspaces(Diagnostic* d) {
  if (b33_material_first_) {
    if (model_->Elements()->empty())
      return model_storage::Invalid("B33 material-first empty model", d);
    for (const auto& entry : *model_->Elements())
      if (!dynamic_cast<const Beam3D*>(entry.second))
        return ReportFailure(StatusCode::kUnsupported,
                             DiagnosticStage::kValidation,
                             "B33 material-first requires pure B33 model", d);
  }
  if (state_->Rotations()) {
    if (!dofs_->InternalDofCount()) return StatusCode::kOk;
    const int count = ConfiguredThreadCount();
    const auto status = PrepareWorkerWorkspaces(count, 12, &workers_, d);
    if (status == StatusCode::kOk) worker_count_ = count;
    return status;
  }
  std::vector<SolidFormulation> formulations;
  for (const auto& entry : *model_->Elements()) {
    if (dynamic_cast<const SolidElement*>(entry.second)) {
      if (d) d->element_id = entry.first;
      return ReportFailure(StatusCode::kUnsupported,
          DiagnosticStage::kValidation,
          "StVK nonlinear Solid is deprecated; use incremental_elastic", d);
    }
    SolidElementReadView solid_view;
    const auto view_status = SolidElementReadView::Bind(
        entry.second, &solid_view, d);
    if (view_status != StatusCode::kOk) return view_status;
    const auto* solid = &solid_view;
    if (!solid->element()->SupportsTrialEvaluation())
      return ReportFailure(StatusCode::kUnsupported,
                           DiagnosticStage::kValidation,
                           "Unsupported nonlinear Solid execution", d);
    if (std::find(formulations.begin(), formulations.end(),
                  solid->formulation()) == formulations.end())
      formulations.push_back(solid->formulation());
  }
  const int reference_count = 2 * static_cast<int>(formulations.size());
  solid_references_ = new SolidReferenceData* [reference_count] {};
  solid_reference_count_ = reference_count;
  for (int i = 0; i < solid_reference_count_; ++i) {
    const auto status = SolidReferenceData::Create(
        formulations[i / 2],
        i % 2 == 0 ? IntegrationPurpose::kStiffness
                   : IntegrationPurpose::kTrialGeometry,
        &solid_references_[i], d);
    if (status != StatusCode::kOk) return status;
  }
  const int count = ConfiguredThreadCount();
  workers_ = new ElementWorkspace* [count] {};
  worker_count_ = count;
  for (int i = 0; i < worker_count_; ++i) {
    const auto status = state_->SolidHistory()
        ? ElementWorkspace::CreateForIncrementalSolidCatalog(
              solid_references_, solid_reference_count_, &workers_[i], d)
        : ElementWorkspace::CreateForSolidCatalog(
              solid_references_, solid_reference_count_, &workers_[i], d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}
StatusCode NonlinearStaticAnalysis::PrepareStep(Diagnostic* d) {
  const int n = dofs_->GlobalCount();
  auto status = PrepareSolidWorkspaces(d);
  if (status == StatusCode::kOk) status = Vector::Create(n, &start_u_, d);
  if (status == StatusCode::kOk)
    status = start_u_->CopyFrom(state_->CommittedDisplacement(), d);
  if (status == StatusCode::kOk)
    status = Vector::Create(n, &start_reaction_, d);
  if (status == StatusCode::kOk)
    status = start_reaction_->CopyFrom(state_->CommittedReaction(), d);
  if (status == StatusCode::kOk && state_->Rotations()) {
    status = Matrix::Create(dofs_->NodeCount(), 4, &start_q_, d);
    if (status == StatusCode::kOk)
      status =
          start_q_->CopyFrom(state_->Rotations()->CommittedOrientations(), d);
  }
  if (status == StatusCode::kOk) status = Vector::Create(n, &predictor_, d);
  if (status == StatusCode::kOk)
    status = Vector::Create(n, &next_predictor_, d);
  if (status == StatusCode::kOk && dofs_->InternalDofCount() > 0)
    status = Vector::Create(dofs_->InternalDofCount(), &internal_predictor_, d);
  if (status == StatusCode::kOk && dofs_->InternalDofCount() > 0)
    status = Vector::Create(dofs_->InternalDofCount(), &next_internal_predictor_, d);
  if (status != StatusCode::kOk) return status;
  clock_ = new IncrementClock;
  next_clock_ = new IncrementClock;
  target_ = new IncrementTarget;
  attempt_ = new IncrementAttemptInfo;
  path_ = new NonlinearStepConditionInput;
  path_->conditions = conditions_;
  path_->dofs = dofs_;
  path_->start_displacement = start_u_;
  path_->start_orientations = start_q_;
  path_->start_reaction = start_reaction_;
  report_->equation_count = dofs_->EquationCount();
  return InitializeIncrementClock(Parameters(conditions_), clock_, d);
}
StatusCode NonlinearStaticAnalysis::NotifyIncrement(
    IncrementEvent event, const NewtonIterationInfo* iteration,
    const Diagnostic* primary, Diagnostic* d) {
  Diagnostic error;
  CopyDiagnosticMetadata(d, &error);
  StatusCode status;
  try {
    status = progress_->OnIncrementProgress(
        event, attempt_, iteration,
        iteration && iteration->residual_valid ? newton_report_->Displacement()
                                               : nullptr,
        iteration && iteration->residual_valid ? newton_report_->Rotation()
                                               : nullptr,
        primary, &error);
  } catch (const std::bad_alloc&) {
    status = model_storage::AllocationFailure(
        "Nonlinear increment progress allocation", &error);
  }
  if (status != StatusCode::kOk) {
    newton_report_->RecordLoggingFailure(status, &error);
    report_->logging_diagnostic = newton_report_->LoggingDiagnostic();
    CopyDiagnosticMetadata(&error, d);
  }
  return status;
}
StatusCode NonlinearStaticAnalysis::NotifyProgress(AnalysisEvent event,
                                                   Diagnostic* d) {
  Diagnostic error;
  CopyDiagnosticMetadata(d, &error);
  const auto status = Notify(progress_, event, report_, &error);
  if (status != StatusCode::kOk) {
    newton_report_->RecordLoggingFailure(status, &error);
    report_->logging_diagnostic = newton_report_->LoggingDiagnostic();
    CopyDiagnosticMetadata(&error, d);
  }
  return status;
}
StatusCode NonlinearStaticAnalysis::OnIteration(const NewtonReport* report,
                                                Diagnostic* d) {
  if (!model_ || report != newton_report_)
    return model_storage::Invalid("Inactive nonlinear observation", d);
  attempt_->characteristic_length = length_;
  return NotifyIncrement(
      IncrementEvent::kIteration, report->Info(),
      attempt_diagnostic_ && attempt_diagnostic_->code != StatusCode::kOk
          ? attempt_diagnostic_
          : nullptr,
      d);
}
StatusCode NonlinearStaticAnalysis::PrepareAttempt(NonlinearFailure* failure,
                                                   Diagnostic* d) {
  *failure = NonlinearFailure::kNone;
  auto status = state_->PrepareNonlinearCandidate(conditions_->StepNumber(),
                                                  controls_, d);
  if (status != StatusCode::kOk) return status;
  path_->lambda = target_->time / Parameters(conditions_)->period;
  path_->committed_lambda = clock_->time / Parameters(conditions_)->period;
  delete external_;
  external_ = nullptr;
  delete force_origins_;
  force_origins_ = nullptr;
  delete moment_origins_;
  moment_origins_ = nullptr;
  status = EvaluateNonlinearStepLoads(path_, &external_, &force_origins_,
                                      &moment_origins_, d);
  if (status != StatusCode::kOk) return status;
  if (target_->use_extrapolation) {
    Vector *raw_free = nullptr, *raw_global = nullptr;
    if (dofs_->FreeCount() > 0) {
      status = Vector::Create(dofs_->FreeCount(), &raw_free, d);
      std::unique_ptr<Vector> free(raw_free);
      if (status != StatusCode::kOk) return status;
      status = Vector::Create(dofs_->GlobalCount(), &raw_global, d);
      std::unique_ptr<Vector> global(raw_global);
      if (status != StatusCode::kOk) return status;
      std::vector<int> indices(dofs_->FreeCount());
      for (int i = 0; i < dofs_->FreeCount(); ++i)
        indices[i] = dofs_->GlobalDofForEquation(i);
      status = free->GatherFrom(predictor_, indices.data(), free->Size(), 0, d);
      if (status == StatusCode::kOk)
        status = free->Scale(target_->extrapolation_ratio, d);
      if (status == StatusCode::kOk)
        status =
            global->ScatterAdd(free.get(), indices.data(), free->Size(), d);
      if (status != StatusCode::kOk) return status;
      bool branch = false;
      status = state_->ApplyNonlinearCorrection(global.get(), &branch, d);
      if (status != StatusCode::kOk) {
        if (branch && status == StatusCode::kEquilibriumFailure)
          *failure = NonlinearFailure::kRotationBranch;
        return status;
      }
      attempt_->predictor_used = true;
      attempt_->predictor_ratio = target_->extrapolation_ratio;
    }
    if (dofs_->InternalDofCount() > 0) {
      Vector* raw = nullptr;
      status = Vector::Create(dofs_->InternalDofCount(), &raw, d);
      std::unique_ptr<Vector> internal(raw);
      if (status == StatusCode::kOk)
        status = internal->CopyFrom(internal_predictor_, d);
      if (status == StatusCode::kOk)
        status = internal->Scale(target_->extrapolation_ratio, d);
      if (status == StatusCode::kOk)
        status = state_->ApplyBeamInternalCorrection(internal.get(), d);
      if (status != StatusCode::kOk) return status;
      attempt_->predictor_used = true;
      attempt_->predictor_ratio = target_->extrapolation_ratio;
    }
  }
  return ApplyNonlinearStepPrescriptions(path_, state_, failure, d);
}
StatusCode NonlinearStaticAnalysis::EvaluateResidual(
    int completed, FieldResidual* force, FieldResidual* moment, double* length,
    NonlinearFailure* failure, Diagnostic* d) {
  if (!model_) return model_storage::Invalid("Inactive nonlinear residual", d);
  delete residual_;
  residual_ = nullptr;
  auto status = Vector::Create(dofs_->GlobalCount(), &residual_, d);
  delete internal_residual_;
  internal_residual_ = nullptr;
  if (status == StatusCode::kOk && dofs_->InternalDofCount() > 0)
    status = Vector::Create(dofs_->InternalDofCount(), &internal_residual_, d);
  if (status != StatusCode::kOk) return status;
  NonlinearResidualInput input;
  input.model = model_;
  input.dofs = dofs_;
  input.state = state_;
  input.controls = controls_;
  input.workspaces = workers_;
  input.workspace_count = worker_count_;
  input.external = external_;
  input.force_origins = force_origins_;
  input.moment_origins = moment_origins_;
  input.actual_correction = completed ? correction_ : nullptr;
  input.internal_actual_correction = completed ? internal_correction_ : nullptr;
  input.previous_force_max_mean =
      state_->FieldHistory()->PreviousMaximumMean(0);
  input.previous_moment_max_mean =
      state_->FieldHistory()->PreviousMaximumMean(1);
  NonlinearResidualOutput output;
  output.residual = residual_;
  output.reaction = state_->CandidateReaction();
  output.force = force;
  output.moment = moment;
  output.internal_residual = internal_residual_;
  status = EvaluateNonlinearResidual(&input, &output, failure, d);
  if (status != StatusCode::kOk) return status;
  length_ = output.characteristic_length;
  *length = length_;
  status = state_->PreviewNonlinearFluxHistory(
      force->spatial_mean, moment->spatial_mean, force->maximum_flux,
      moment->maximum_flux, &force->raw_mean, &moment->raw_mean, d);
  if (status != StatusCode::kOk) return status;
  if (!state_->Rotations()) {
    // A translation-only model has no moment origin or cross-field conversion.
    force->used_mean = force->raw_mean;
    return StatusCode::kOk;
  }
  return ConvertNonlinearFieldMeans(controls_, length_, force, moment, d);
}
StatusCode NonlinearStaticAnalysis::AssembleTangent(
    const EquationSystem** system, Vector** correction,
    NonlinearFailure* failure, Diagnostic* d) {
  if (system) *system = nullptr;
  if (correction) *correction = nullptr;
  if (!model_ || !system || !correction)
    return model_storage::Invalid("Inactive nonlinear tangent", d);
  delete equations_;
  equations_ = nullptr;
  delete free_correction_;
  free_correction_ = nullptr;
  // IE5-R: lag only geometric stress on the first extrapolated correction.
  // Keep this a call argument so retries/workers cannot retain a prior choice.
  const bool use_committed_geometric_stress =
      state_->SolidHistory() && attempt_->predictor_used &&
      newton_report_->Info()->iteration == 0;
  auto status = AssembleNonlinearTangent(
      model_, dofs_, state_, residual_, workers_, worker_count_,
      state_->Rotations() ? MatrixType::kGeneral
                          : MatrixType::kSymmetricIndefinite,
      &equations_, failure, d, use_committed_geometric_stress,
      internal_residual_,
      b33_material_first_ && newton_report_->Info()->iteration == 0);
  delete residual_;
  residual_ = nullptr;
  delete internal_residual_;
  internal_residual_ = nullptr;
  if (status == StatusCode::kOk && dofs_->EquationCount() > 0)
    status = Vector::Create(dofs_->EquationCount(), &free_correction_, d);
  if (status != StatusCode::kOk) {
    delete equations_;
    equations_ = nullptr;
    return status;
  }
  *system = equations_;
  *correction = free_correction_;
  return StatusCode::kOk;
}
StatusCode NonlinearStaticAnalysis::ApplyCorrection(const Vector* correction,
                                                    NonlinearFailure* failure,
                                                    Diagnostic* d) {
  if (failure) *failure = NonlinearFailure::kNone;
  if (!model_ || !failure || correction != free_correction_)
    return model_storage::Invalid("Inactive nonlinear correction", d);
  delete correction_;
  correction_ = nullptr;
  auto status = Vector::Create(dofs_->GlobalCount(), &correction_, d);
  if (status != StatusCode::kOk) return status;
  std::vector<int> indices(dofs_->FreeCount());
  for (int i = 0; i < dofs_->FreeCount(); ++i)
    indices[i] = dofs_->GlobalDofForEquation(i);
  if (dofs_->FreeCount() > 0) {
    Vector* raw = nullptr;
    status = Vector::Create(dofs_->FreeCount(), &raw, d);
    std::unique_ptr<Vector> physical(raw);
    std::vector<int> equations(dofs_->FreeCount());
    for (int i = 0; i < dofs_->FreeCount(); ++i) equations[i] = i;
    if (status == StatusCode::kOk)
      status = physical->GatherFrom(correction, equations.data(),
                                    physical->Size(), 0, d);
    if (status == StatusCode::kOk)
      status = correction_->ScatterAdd(physical.get(), indices.data(),
                                       physical->Size(), d);
  }
  if (status != StatusCode::kOk) return status;
  delete internal_correction_;
  internal_correction_ = nullptr;
  if (dofs_->InternalDofCount() > 0) {
    status =
        Vector::Create(dofs_->InternalDofCount(), &internal_correction_, d);
    std::vector<int> equations(dofs_->InternalDofCount());
    for (int i = 0; i < dofs_->InternalDofCount(); ++i)
      equations[i] = dofs_->FreeCount() + i;
    if (status == StatusCode::kOk)
      status = internal_correction_->GatherFrom(
          correction, equations.data(), internal_correction_->Size(), 0, d);
    if (status != StatusCode::kOk) return status;
  }
  bool branch = false;
  status = state_->ApplyNonlinearCorrection(correction_, &branch, d);
  if (status == StatusCode::kOk && internal_correction_)
    status = state_->ApplyBeamInternalCorrection(internal_correction_, d);
  if (status == StatusCode::kEquilibriumFailure && branch)
    *failure = NonlinearFailure::kRotationBranch;
  return status;
}
void NonlinearStaticAnalysis::UpdateReport() noexcept {
  const auto* info = newton_report_->Info();
  report_->linear_solver_calls += info->linear_solver_calls;
  report_->iteration_count += info->iteration;
  report_->assembly_seconds += info->tangent_seconds;
  report_->solve_seconds += info->solve_seconds;
  report_->equilibrium_recovery_seconds +=
      info->residual_seconds + info->criteria_seconds;
  report_->residual_valid = info->residual_valid;
  report_->displacement_field =
      info->residual_valid ? newton_report_->Displacement() : nullptr;
  report_->rotation_field =
      info->residual_valid ? newton_report_->Rotation() : nullptr;
  report_->characteristic_length = length_;
  report_->calculation_status = info->numerical_status;
  report_->logging_diagnostic = newton_report_->LoggingDiagnostic();
}
StatusCode NonlinearStaticAnalysis::FinishAttempt(
    ResultWriter* writer, std::int64_t frame, double total_time,
    bool had_cutback, NonlinearFailure* failure, Diagnostic* d) {
  const double began = Seconds();
  auto status = AcceptIncrement(Parameters(conditions_), controls_->increments,
                                clock_, target_, newton_report_->Info(),
                                had_cutback, next_clock_, d);
  if (status == StatusCode::kOk) status = state_->PrepareNonlinearOutput(d);
  if (status == StatusCode::kOk)
    status = PhysicalIncrement(dofs_, state_, next_predictor_, failure, d);
  if (status == StatusCode::kOk && next_internal_predictor_) {
    status = next_internal_predictor_->CopyFrom(
        state_->BeamInternals()->CandidatePosition(), d);
    if (status == StatusCode::kOk)
      status = next_internal_predictor_->AddScaled(
          -1, state_->BeamInternals()->CommittedPosition(), d);
  }
  if (status != StatusCode::kOk) return status;
  state_->MarkResultsUnconfirmed(conditions_->StepNumber(),
                                 clock_->successful + 1);
  const double recovery = Seconds();
  status = RecoverNonlinearResults(model_, dofs_, state_, workers_,
                                   worker_count_, failure, d);
  if (status == StatusCode::kOk && state_->SolidHistory())
    status = state_->SolidHistory()->ValidateCandidate(d);
  report_->equilibrium_recovery_seconds += Seconds() - recovery;
  report_->calculation_seconds += Seconds() - began;
  report_->calculation_finished = true;
  report_->calculation_status = status;
  if (status != StatusCode::kOk) return status;
  report_->element_dispatch_seconds =
      ElementDispatchSeconds() - element_started_;
  report_->observed_element_team_max = ObservedElementTeamMaximum();
  status = NotifyProgress(AnalysisEvent::kCalculationFinished, d);
  if (status != StatusCode::kOk) return status;
  const StepDefinition* step = nullptr;
  status = conditions_->SourceDomain()->FindStepDefinition(
      conditions_->StepNumber(), &step, d);
  if (status != StatusCode::kOk) return status;
  FrameInfo output;
  output.geometry = FrameGeometry::kNonlinear;
  output.frame_number = frame;
  output.step_number = conditions_->StepNumber();
  output.increment_number = clock_->successful + 1;
  output.step_name = step->Name();
  output.step_time = target_->time;
  output.total_time = total_time + target_->time;
  output.completes_step = target_->completes_step;
  const double started = Seconds();
  status = writer->WriteFrame(&output, conditions_->SourceDomain(), state_, d);
  if (status == StatusCode::kOk) {
    state_->Commit(conditions_->StepNumber(), clock_->successful + 1);
    std::swap(clock_, next_clock_);
    std::swap(predictor_, next_predictor_);
    std::swap(internal_predictor_, next_internal_predictor_);
    report_->step_completed = target_->completes_step;
    report_->last_committed_step = conditions_->StepNumber();
    report_->last_committed_increment = clock_->successful;
    report_->last_committed_frame = frame;
    report_->last_committed_total_time = output.total_time;
  }
  report_->output_seconds += Seconds() - started;
  if (status != StatusCode::kOk) return status;
  status = NotifyProgress(AnalysisEvent::kIncrementCommitted, d);
  if (status != StatusCode::kOk) return status;
  attempt_->proposed_available = true;
  attempt_->proposed_increment = clock_->nominal;
  return NotifyIncrement(IncrementEvent::kCommitted, newton_report_->Info(),
                         nullptr, d);
}
StatusCode NonlinearStaticAnalysis::InspectFailedCandidate(
    NewtonIterationInfo* failed, Diagnostic* d) {
  if (failed->logging_failed || failed->outcome == NewtonOutcome::kFatal ||
      failed->status != StatusCode::kEquilibriumFailure)
    return StatusCode::kOk;
  attempt_->original_outcome = failed->outcome;
  if (!state_->HasValidCandidate()) {
    attempt_->geometry_inspection = GeometryInspection::kUnavailableCandidate;
    return StatusCode::kOk;
  }
  bool supported = false;
  for (const auto& entry : *model_->Elements())
    supported = supported || entry.second->SupportsTrialGeometryCheck();
  if (!supported) {
    attempt_->geometry_inspection = GeometryInspection::kUnsupported;
    return StatusCode::kOk;
  }
  NonlinearFailure failure = NonlinearFailure::kNone;
  const auto status = CheckNonlinearTrialGeometry(
      model_, dofs_, state_, workers_, worker_count_, &failure, d);
  if (status == StatusCode::kOk) {
    attempt_->geometry_inspection = GeometryInspection::kPassed;
    return StatusCode::kOk;
  }
  attempt_->geometry_diagnostic = d;
  if (status == StatusCode::kEquilibriumFailure &&
      failure == NonlinearFailure::kGeometry) {
    attempt_->geometry_inspection = GeometryInspection::kGeometry;
    failed->outcome = NewtonOutcome::kGeometry;
    return StatusCode::kOk;
  }
  attempt_->geometry_inspection = GeometryInspection::kFailed;
  failed->outcome = NewtonOutcome::kFatal;
  failed->status = status;
  return status;
}

StatusCode NonlinearStaticAnalysis::RunIncrements(ResultWriter* writer,
                                                  std::int64_t first_frame,
                                                  double total_time,
                                                  Diagnostic* d) {
  const auto* parameters = Parameters(conditions_);
  while (clock_->time < parameters->period) {
    double nominal = clock_->nominal;
    bool had_cutback = false;
    for (int attempt_number = 1;; ++attempt_number) {
      Diagnostic numeric, inspection;
      attempt_->geometry_inspection = GeometryInspection::kNotRequested;
      attempt_->geometry_diagnostic = nullptr;
      attempt_->original_outcome = NewtonOutcome::kNotStarted;
      numeric.step = conditions_->StepNumber();
      numeric.source_file = conditions_->SourceDomain()->SourceFile();
      attempt_diagnostic_ = &numeric;
      auto status = ProposeIncrement(parameters, controls_->increments, clock_,
                                     nominal, target_, &numeric);
      if (status != StatusCode::kOk) {
        CopyDiagnosticMetadata(&numeric, d);
        attempt_diagnostic_ = nullptr;
        return status;
      }
      attempt_->step_number = conditions_->StepNumber();
      attempt_->increment_number = clock_->successful + 1;
      attempt_->attempt_number = attempt_number;
      attempt_->step_time = target_->time;
      attempt_->total_time = total_time + target_->time;
      attempt_->increment_size = target_->size;
      attempt_->nominal_increment = target_->nominal;
      attempt_->completes_step = target_->completes_step;
      attempt_->predictor_used = false;
      attempt_->predictor_ratio = 0;
      attempt_->proposed_available = false;
      report_->increment_number = clock_->successful + 1;
      ++report_->attempt_count;
      report_->calculation_finished = false;
      report_->residual_valid = false;
      report_->displacement_field = report_->rotation_field = nullptr;
      newton_report_->Reset();
      status = NotifyIncrement(IncrementEvent::kAttemptStarted, nullptr,
                               nullptr, &numeric);
      if (status != StatusCode::kOk) {
        CopyDiagnosticMetadata(&numeric, d);
        attempt_diagnostic_ = nullptr;
        return status;
      }
      NonlinearFailure failure = NonlinearFailure::kNone;
      const double started = Seconds();
      status = PrepareAttempt(&failure, &numeric);
      bool ran_newton = false;
      if (status == StatusCode::kOk) {
        ran_newton = true;
        status = newton_->Solve(this, solver_, controls_, controls_->newton,
                                newton_report_, this, &numeric);
        UpdateReport();
      }
      report_->calculation_seconds += Seconds() - started;
      report_->calculation_finished = true;
      report_->calculation_status = status;
      const bool converged = status == StatusCode::kOk;
      if (converged)
        status = FinishAttempt(writer, first_frame + clock_->successful,
                               total_time, had_cutback, &failure, &numeric);
      if (status == StatusCode::kOk) {
        attempt_diagnostic_ = nullptr;
        break;
      }
      NewtonIterationInfo failed;
      failed.outcome = NewtonOutcome::kFatal;
      if (ran_newton) {
        // Local scalar metadata copy, not a physical vector/state copy.
        failed = *newton_report_->Info();
      }
      failed.status = failed.numerical_status = status;
      // P06: after convergence, recovery/output failures terminate this step.
      // A typed geometry diagnostic describes the failed operation, not a
      // request to retry an already converged increment with DH.
      if (converged) {
        failed.outcome = NewtonOutcome::kFatal;
      } else if (status == StatusCode::kEquilibriumFailure) {
        if (failure == NonlinearFailure::kGeometry)
          failed.outcome = NewtonOutcome::kGeometry;
        else if (failure == NonlinearFailure::kRotationBranch)
          failed.outcome = NewtonOutcome::kRotationBranch;
      }
      inspection.step = numeric.step;
      inspection.source_file = numeric.source_file;
      const auto inspection_status =
          InspectFailedCandidate(&failed, &inspection);
      state_->DiscardCandidate();
      Diagnostic log_error;
      const auto logged = NotifyIncrement(IncrementEvent::kAttemptFailed,
                                          &failed, &numeric, &log_error);
      if (logged != StatusCode::kOk) failed.logging_failed = true;
      if (inspection_status != StatusCode::kOk) {
        CopyDiagnosticMetadata(&inspection, d);
        attempt_diagnostic_ = nullptr;
        return inspection_status;
      }
      // A postcommit observer failure never retries an already accepted frame.
      if (clock_->successful != target_->base_successful) {
        CopyDiagnosticMetadata(&numeric, d);
        attempt_diagnostic_ = nullptr;
        return status;
      }
      double retry = 0;
      const auto cutback =
          CutbackIncrement(parameters, controls_->increments, clock_, target_,
                           &failed, attempt_number, &retry, &numeric);
      if (cutback != StatusCode::kOk) {
        CopyDiagnosticMetadata(&numeric, d);
        attempt_diagnostic_ = nullptr;
        return cutback;
      }
      attempt_->proposed_available = true;
      attempt_->proposed_increment = retry;
      status = NotifyIncrement(IncrementEvent::kCutback, &failed, nullptr,
                               &log_error);
      if (status != StatusCode::kOk) {
        CopyDiagnosticMetadata(&log_error, d);
        attempt_diagnostic_ = nullptr;
        return status;
      }
      nominal = retry;
      had_cutback = true;
      attempt_diagnostic_ = nullptr;
    }
  }
  return StatusCode::kOk;
}
StatusCode NonlinearStaticAnalysis::SolveStep(
    const AnalysisModel* model, const DofManager* dofs,
    const EffectiveConditions* conditions, const ControlSettings* controls,
    AnalysisState* state, ResultWriter* writer, std::int64_t first_frame,
    double total_time, AnalysisReport* report, AnalysisProgress* progress,
    Diagnostic* d) {
  ClearExecution();
  newton_report_->Reset();
  if (report) report->Reset();
  const double element_started = ElementDispatchSeconds();
  StatusCode status = StatusCode::kOk;
  try {
    if (!model || !dofs || !conditions || !controls || !controls->newton ||
        !controls->increments || !controls->displacement ||
        !controls->rotation || !state || !state->FieldHistory() || !writer ||
        !report || !progress || model->Conditions() != conditions ||
        !Parameters(conditions) ||
        conditions->StepNumber() !=
            static_cast<std::int64_t>(state->StepNumber()) + 1 ||
        dofs->GlobalCount() != state->GlobalCount() || first_frame <= 0 ||
        !std::isfinite(total_time) || total_time < 0 ||
        !std::isfinite(total_time + Parameters(conditions)->period) ||
        first_frame > std::numeric_limits<std::int64_t>::max() -
                          Parameters(conditions)->increment_limit)
      status = model_storage::Invalid("Nonlinear SolveStep context", d);
    else {
      model_ = model;
      dofs_ = dofs;
      conditions_ = conditions;
      controls_ = controls;
      state_ = state;
      report_ = report;
      element_started_ = element_started;
      progress_ = progress;
      report->step_number = conditions->StepNumber();
      report->last_committed_step = state->StepNumber();
      report->last_committed_increment = state->IncrementNumber();
      report->last_committed_frame = first_frame - 1;
      report->last_committed_total_time = total_time;
      status = PrepareStep(d);
      if (status == StatusCode::kOk)
        status = RunIncrements(writer, first_frame, total_time, d);
    }
  } catch (const std::bad_alloc&) {
    status = model_storage::AllocationFailure(
        "Nonlinear Analysis execution storage", d);
  }
  if (status != StatusCode::kOk && state) state->DiscardCandidate();
  if (report) {
    report->status = status;
    report->element_dispatch_seconds =
        ElementDispatchSeconds() - element_started;
    report->observed_element_team_max = ObservedElementTeamMaximum();
    if (status != StatusCode::kOk &&
        status != StatusCode::kEquilibriumFailure) {
      report->residual_valid = false;
      report->displacement_field = report->rotation_field = nullptr;
    }
  }
  ClearExecution();
  return status;
}
}  // namespace fesa
