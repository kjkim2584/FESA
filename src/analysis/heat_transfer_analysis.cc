#include "fesa/analysis/heat_transfer_analysis.h"

#include <chrono>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <new>
#include <set>
#include <vector>

#include "adapters/mkl/pardiso_solver.h"
#include "adapters/openmp/thread_runtime.h"
#include "analysis/thermal_assembler.h"
#include "base/diagnostic_copy.h"
#include "fesa/analysis/thermal_increment_controller.h"
#include "fesa/model/film_property_table.h"
#include "model/model_storage.h"

namespace fesa {
namespace {
double WallTime() {
  return std::chrono::duration<double>(
             std::chrono::steady_clock::now().time_since_epoch())
      .count();
}
void RefreshDispatchReport(AnalysisReport* report, double baseline) {
  if (!report) return;
  report->element_dispatch_seconds = ElementDispatchSeconds() - baseline;
  report->observed_element_team_max = ObservedElementTeamMaximum();
}
StatusCode Error(const char* operation, Diagnostic* d,
                 StatusCode code = StatusCode::kInvalidArgument) {
  return ReportFailure(code,
                       (code == StatusCode::kEquilibriumFailure ||
                        code == StatusCode::kLibraryFailure)
                           ? DiagnosticStage::kComputation
                           : DiagnosticStage::kValidation,
                       operation, d);
}
int Root(const std::map<int, int>* parents, int id) {
  while (parents->at(id) != id) id = parents->at(id);
  return id;
}
// A call-local problem owns numerical storage; no binding survives SolveStep.
class ThermalStepProblem final : public ThermalNonlinearProblem {
 public:
  ThermalStepProblem(const ThermalModel* model, const ThermalDofMap* dofs,
                     const ThermalConditions* conditions,
                     const FieldControls* field, ThermalState* state,
                     const FirstOrderTimeIntegrator* integrator, double period,
                     AnalysisReport* report)
      : model_(model),
        dofs_(dofs),
        conditions_(conditions),
        field_(field),
        state_(state),
        integrator_(integrator),
        period_(period),
        report_(report) {}
  ~ThermalStepProblem() override {
    ReleaseAttemptResources();
    delete assembly_;
    delete free_globals_;
    delete increment_;
    delete free_increment_;
  }
  ThermalStepProblem(const ThermalStepProblem&) = delete;
  ThermalStepProblem& operator=(const ThermalStepProblem&) = delete;
  ThermalStepProblem(ThermalStepProblem&&) = delete;
  ThermalStepProblem& operator=(ThermalStepProblem&&) = delete;
  StatusCode Initialize(const Vector* step_start, Diagnostic* d) {
    auto status = ThermalAssembly::CreateForStep(model_, dofs_, conditions_,
                                                 step_start, &assembly_, d);
    if (status != StatusCode::kOk) return status;
    free_globals_ = new std::vector<int>(0);
    for (int e = 0; e < dofs_->FreeCount(); ++e)
      free_globals_->push_back(dofs_->GlobalDofForEquation(e));
    status = Vector::Create(dofs_->GlobalCount(), &increment_, d);
    if (status == StatusCode::kOk && dofs_->FreeCount() > 0)
      status = Vector::Create(dofs_->FreeCount(), &free_increment_, d);
    return status;
  }
  void BeginAttempt(double time, double dt, bool extrapolated) {
    extrapolated_start_ = extrapolated;
    extrapolated_evaluation_ = false;
    time_ = time;
    dt_ = dt;
    correction_max_ = 0;
  }
  StatusCode EvaluateResidual(int completed, HeatFluxResidual* out,
                              Diagnostic* d) override {
    extrapolated_evaluation_ = completed == 0 && extrapolated_start_;
    double start = WallTime();
    auto status = ThermalAssembler::Assemble(
        model_, dofs_, conditions_, state_, time_, period_, integrator_,
        integrator_ ? dt_ : 0, assembly_, d, extrapolated_evaluation_);
    report_->assembly_seconds += WallTime() - start;
    if (status != StatusCode::kOk) return status;
    start = WallTime();
    status = MeasureHeatFlux(assembly_, state_->FluxHistory(), field_, out, d);
    report_->equilibrium_recovery_seconds += WallTime() - start;
    if (status == StatusCode::kOk)
      out->correction_max = completed > 0 ? correction_max_ : 0;
    return status;
  }
  StatusCode AssembleTangent(const EquationSystem** out, Vector** correction,
                             Diagnostic* d) override {
    *out = nullptr;
    *correction = nullptr;
    if (!assembly_->IsValid() || dofs_->FreeCount() <= 0)
      return Error("Thermal tangent without valid free equations", d);
    const double start = WallTime();
    SparseMatrix* raw_matrix = nullptr;
    auto status = SparseMatrix::Create(
        dofs_->FreeCount(), static_cast<int>(assembly_->Rows()->size()),
        assembly_->Rows()->data(), assembly_->Columns()->data(),
        assembly_->Values(), &raw_matrix, d);
    std::unique_ptr<SparseMatrix> matrix(raw_matrix);
    if (status != StatusCode::kOk) return status;
    Vector* raw_rhs = nullptr;
    status = Vector::Create(dofs_->FreeCount(), &raw_rhs, d);
    std::unique_ptr<Vector> rhs(raw_rhs);
    if (status == StatusCode::kOk)
      status = rhs->CopyFrom(assembly_->FreeResidual(), d);
    if (status != StatusCode::kOk) return status;
    MatrixType matrix_type = MatrixType::kSymmetricPositiveDefinite;
    for (const auto* row : *conditions_->ElementLoads())
      if ((row->origin->specification()->kind == ThermalConditionKind::kFilm &&
           row->origin->specification()->film_property) ||
          (extrapolated_evaluation_ && row->origin->specification()->kind ==
                                           ThermalConditionKind::kRadiation))
        matrix_type = MatrixType::kSymmetricIndefinite;
    status = EquationSystem::Create(matrix.get(), rhs.get(), matrix_type,
                                    &system_, d);
    if (status != StatusCode::kOk) return status;
    matrix.release();
    rhs.release();
    status = Vector::Create(dofs_->FreeCount(), &correction_, d);
    report_->assembly_seconds += WallTime() - start;
    if (status != StatusCode::kOk) return status;
    *out = system_;
    *correction = correction_;
    return StatusCode::kOk;
  }
  StatusCode ApplyCorrection(const Vector* correction, Diagnostic* d) override {
    int index = 0;
    auto status = correction->MaximumAbs(&correction_max_, &index, d);
    if (status == StatusCode::kOk)
      status = state_->CandidateTemperature()->ScatterAdd(
          correction, free_globals_->data(), dofs_->FreeCount(), d);
    return status;
  }
  StatusCode FreeChange(double* out, Diagnostic* d) {
    *out = 0;
    if (!free_increment_) return StatusCode::kOk;
    auto status = increment_->CopyFrom(state_->CandidateTemperature(), d);
    if (status == StatusCode::kOk)
      status = increment_->AddScaled(-1, state_->CommittedTemperature(), d);
    if (status == StatusCode::kOk)
      status = free_increment_->GatherFrom(increment_, free_globals_->data(),
                                           dofs_->FreeCount(), 0, d);
    int index = 0;
    if (status == StatusCode::kOk)
      status = free_increment_->MaximumAbs(out, &index, d);
    return status;
  }
  const Vector* Reaction() const { return assembly_->Reaction(); }
  void ReleaseEquationSystem() noexcept override {
    delete system_;
    system_ = nullptr;
  }
  void ReleaseFreeCorrection() noexcept override {
    delete correction_;
    correction_ = nullptr;
  }
  void ReleaseCorrectionAssessment() noexcept override { correction_max_ = 0; }
  void ReleaseAttemptResources() noexcept override {
    ReleaseEquationSystem();
    ReleaseFreeCorrection();
    ReleaseCorrectionAssessment();
  }

 private:
  const ThermalModel* model_;            // Borrowed for SolveStep only.
  const ThermalDofMap* dofs_;            // Borrowed from model during call.
  const ThermalConditions* conditions_;  // Borrowed for call.
  const FieldControls* field_;           // Borrowed for call.
  ThermalState* state_;                  // Borrowed for call.
  const FirstOrderTimeIntegrator* integrator_;  // Borrowed for call or null.
  double period_, time_ = 0, dt_ = 0, correction_max_ = 0;
  bool extrapolated_start_ = false;       ///< Current attempt has a predictor.
  bool extrapolated_evaluation_ = false;  ///< Last evaluation was that start.
  AnalysisReport* report_;               // Borrowed for call.
  ThermalAssembly* assembly_ = nullptr;  // Owned, preserves converged reaction.
  EquationSystem* system_ = nullptr;     // Owned temporary tangent.
  Vector* correction_ = nullptr;         // Owned temporary correction.
  std::vector<int>* free_globals_ = nullptr;  // Owned indices.
  Vector* increment_ = nullptr;       // Owned nominal difference scratch.
  Vector* free_increment_ = nullptr;  // Owned optional free scratch.
};
class ThermalTimedSolver final : public LinearSolver {
 public:
  ThermalTimedSolver(LinearSolver* solver, AnalysisReport* report)
      : solver_(solver), report_(report) {}
  StatusCode Solve(const EquationSystem* system, Vector* out,
                   Diagnostic* d) override {
    const double start = WallTime();
    const auto status = solver_->Solve(system, out, d);
    report_->solve_seconds += WallTime() - start;
    return status;
  }

 private:
  LinearSolver* solver_;    // Borrowed for enclosing SolveStep.
  AnalysisReport* report_;  // Borrowed for enclosing SolveStep.
};
class ThermalStepObserver final : public ThermalNewtonProgress {
 public:
  ThermalStepObserver(AnalysisProgress* progress,
                      const IncrementAttemptInfo* attempt,
                      const Diagnostic* primary, AnalysisReport* report,
                      double baseline)
      : progress_(progress),
        attempt_(attempt),
        primary_(primary),
        report_(report),
        baseline_(baseline) {}
  StatusCode OnIteration(const ThermalIterationReport* report,
                         Diagnostic* d) override {
    if (!progress_) return StatusCode::kOk;
    RefreshDispatchReport(report_, baseline_);
    return progress_->OnThermalIncrementProgress(
        IncrementEvent::kIteration, attempt_, report->Info(),
        report->Info()->residual_valid ? report->Residual() : nullptr,
        report->Info()->numerical_status == StatusCode::kOk ? nullptr
                                                            : primary_,
        d);
  }

 private:
  AnalysisProgress* progress_;           // Borrowed for this SolveStep.
  const IncrementAttemptInfo* attempt_;  // Borrowed current attempt metadata.
  AnalysisReport* report_;               ///< Borrowed for this SolveStep.
  double baseline_;
  const Diagnostic* primary_;  // Borrowed current numerical output.
};
StatusCode WriteNodes(const ThermalDofMap* dofs, ThermalState* state,
                      Diagnostic* d) {
  for (int row = 0; row < dofs->NodeCount(); ++row) {
    const int id = dofs->NodeId(row);
    if (!dofs->Participates(id)) continue;
    const int count = dofs->NodeDofCount(id);
    std::vector<int> globals(static_cast<std::size_t>(count));
    for (int j = 0; j < count; ++j) globals[j] = dofs->GlobalDof(id, 11 + j);
    Vector* raw = nullptr;
    auto status = Vector::Create(count, &raw, d);
    if (status != StatusCode::kOk) return status;
    std::unique_ptr<Vector> nt(raw);
    status = Vector::Create(count, &raw, d);
    if (status != StatusCode::kOk) return status;
    std::unique_ptr<Vector> rfl(raw);
    status = nt->GatherFrom(state->CandidateTemperature(), globals.data(),
                            count, 0, d);
    if (status == StatusCode::kOk)
      status = rfl->GatherFrom(state->CandidateReactionFlux(), globals.data(),
                               count, 0, d);
    if (status == StatusCode::kOk)
      status = dofs->SourceDomain()->Nodes()->at(id)->UpdateTemperatureResults(
          nt.get(), rfl.get(), d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}
// D6 applies only to loaded faces, even for explicit zero magnitude. This
// preflight runs for every step before the solver creates/truncates output.
StatusCode ValidateLoadedFaces(const ThermalModel* model,
                               const ThermalConditions* conditions,
                               Diagnostic* d) {
  static constexpr const char* errors[] = {
      "Loaded DC3D8 face S1 requires planar geometry",
      "Loaded DC3D8 face S2 requires planar geometry",
      "Loaded DC3D8 face S3 requires planar geometry",
      "Loaded DC3D8 face S4 requires planar geometry",
      "Loaded DC3D8 face S5 requires planar geometry",
      "Loaded DC3D8 face S6 requires planar geometry"};
  const StepDefinition* step = nullptr;
  auto lookup = model->SourceDomain()->FindStepDefinition(
      conditions->StepNumber(), &step, d);
  if (lookup != StatusCode::kOk) return lookup;
  const bool ramp =
      step->HeatTransfer()->DefaultAmplitude() == ThermalStepAmplitude::kRamp;
  for (const auto* row : *conditions->ElementLoads()) {
    const auto* spec = row->origin->specification();
    if (spec->kind == ThermalConditionKind::kBodyFlux) continue;
    const ThermalElement* element = nullptr;
    auto status =
        model->SourceDomain()->FindThermalElement(row->element_id, &element, d);
    if (status != StatusCode::kOk) return status;
    if (d) {
      d->source_file = spec->source_file;
      d->source_line = spec->source_line;
      d->step = conditions->StepNumber();
      d->element_id = row->element_id;
      d->node_id = 0;
      d->condition_dof = 0;
    }
    if (spec->kind == ThermalConditionKind::kRadiation) {
      ThermalLoadContext values;
      status =
          conditions->PopulateExchangeContext(row, step->Period(), &values, d);
      if (status != StatusCode::kOk) return status;
      if (ramp) {
        for (const auto* other : *conditions->ElementLoads()) {
          if (other->element_id == row->element_id &&
              other->origin->specification()->kind ==
                  ThermalConditionKind::kFilm)
            return ReportFailure(StatusCode::kUnsupported,
                                 DiagnosticStage::kValidation,
                                 "D9 same-element film/radiation RAMP", d);
        }
      }
    }
    Vector* raw = nullptr;
    status = Vector::Create(element->dof_count(), &raw, d);
    std::unique_ptr<Vector> weights(raw);
    if (status != StatusCode::kOk) return status;
    for (int i = 0; i < spec->expanded_count; ++i) {
      if (spec->expanded_element_ids[i] != row->element_id) continue;
      const int face = spec->expanded_face_codes[i];
      if (spec->kind == ThermalConditionKind::kFilm &&
          spec->family == ThermalLoadFamily::kSurfaceFilm && ramp &&
          element->type() == ThermalElementType::kDs4 &&
          (!spec->film_amplitude || spec->film_property)) {
        if (d)
          d->condition_dof = face == -1 ? 11 : 10 + element->dofs_per_node();
        return ReportFailure(
            StatusCode::kUnsupported, DiagnosticStage::kValidation,
            face == -1
                ? "D9 DS4 SFILM SNEG RAMP without effective FILM AMPLITUDE"
                : "D9 DS4 SFILM SPOS RAMP without effective FILM AMPLITUDE",
            d);
      }
      status = element->ComputeSurfaceWeights(NormalizedThermalFaceLabel(face),
                                              raw, d);
      if (status == StatusCode::kUnsupported &&
          element->type() == ThermalElementType::kDc3d8 && face >= 1 &&
          face <= 6)
        return ReportFailure(status, DiagnosticStage::kValidation,
                             errors[face - 1], d);
      if (status != StatusCode::kOk) return status;
    }
  }
  return StatusCode::kOk;
}
}  // namespace
StatusCode HeatTransferAnalysis::Create(const HeatTransferStepParameters* step,
                                        const ThermalConditions*,
                                        HeatTransferAnalysis** out,
                                        Diagnostic* d) {
  if (!out) return Error("HeatTransferAnalysis output", d);
  *out = nullptr;
  try {
    std::unique_ptr<HeatTransferAnalysis> a(new HeatTransferAnalysis);
    auto status = CloneHeatTransferStepParameters(step, &a->settings_, d);
    if (status == StatusCode::kOk)
      status = ThermalIterationReport::Create(&a->iterations_, d);
    if (status == StatusCode::kOk)
      status = ThermalNewton::Create(&a->newton_, d);
    if (status != StatusCode::kOk) return status;
    a->logging_error_ = new Diagnostic;
    *out = a.release();
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return model_storage::AllocationFailure("HeatTransferAnalysis allocation",
                                            d);
  }
}
HeatTransferAnalysis::~HeatTransferAnalysis() {
  delete settings_;
  delete iterations_;
  delete newton_;
  delete logging_error_;
}
StatusCode HeatTransferAnalysis::ValidateStep(
    const ThermalModel* model, const ThermalDofMap* dofs,
    const ThermalConditions* conditions, const HeatTransferStepParameters* step,
    Diagnostic* d) {
  if (!model || !dofs || !conditions || !step || model->DofMap() != dofs ||
      model->Conditions() != conditions ||
      model->SourceDomain() != conditions->SourceDomain() ||
      dofs->StepNumber() != conditions->StepNumber() ||
      model->Elements()->empty())
    return Error("Thermal step bindings", d);
  const StepDefinition* definition = nullptr;
  if (model->SourceDomain()->FindStepDefinition(
          conditions->StepNumber(), &definition, d) != StatusCode::kOk ||
      !definition->IsHeatTransfer() || definition->Period() != step->period ||
      definition->HeatTransfer()->steady != step->steady)
    return Error("Thermal step procedure/period mismatch", d);
  auto status = ValidateHeatTransferStepParameters(step, d);
  if (status != StatusCode::kOk) return status;
  try {
    std::map<int, int> parents;
    for (const auto* element : *model->Elements()) {
      status = element->ValidateGeometry(d);
      if (status != StatusCode::kOk) return status;
      if (!step->steady && !element->material()->HasCapacity())
        return Error(
            "Transient heat transfer requires density and specific heat", d,
            StatusCode::kInputError);
      for (int j = 0; j < element->node_count(); ++j) {
        const int id = element->GetNode(j)->Id();
        parents.emplace(id, id);
      }
      for (int j = 1; j < element->node_count(); ++j)
        parents[Root(&parents, element->GetNode(j)->Id())] =
            Root(&parents, element->GetNode(0)->Id());
    }
    status = ValidateLoadedFaces(model, conditions, d);
    if (status != StatusCode::kOk) return status;
    if (step->steady) {
      std::set<int> anchored;
      for (const auto* row : *conditions->Temperatures())
        anchored.insert(Root(&parents, row->node_id));
      for (const auto* row : *conditions->ElementLoads()) {
        if (row->origin->specification()->kind != ThermalConditionKind::kFilm &&
            row->origin->specification()->kind !=
                ThermalConditionKind::kRadiation)
          continue;
        ThermalLoadContext values;
        status =
            conditions->PopulateExchangeContext(row, step->period, &values, d);
        if (status != StatusCode::kOk) return status;
        const bool positive = values.kind == ThermalConditionKind::kRadiation
                                  ? values.emissivity > 0
                              : values.film_property
                                  ? values.film_property->HasHeatExchange()
                                  : values.coefficient > 0;
        if (!positive) continue;
        const ThermalElement* element = nullptr;
        status = model->SourceDomain()->FindThermalElement(row->element_id,
                                                           &element, d);
        if (status != StatusCode::kOk) return status;
        anchored.insert(Root(&parents, element->GetNode(0)->Id()));
      }
      for (const auto& item : parents)
        if (anchored.find(Root(&parents, item.first)) == anchored.end())
          return Error(
              "Steady thermal component lacks prescribed temperature or "
              "positive exchange",
              d, StatusCode::kSingular);
    }
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return model_storage::AllocationFailure("Thermal step preflight indices",
                                            d);
  }
}
StatusCode HeatTransferAnalysis::SolveStep(
    const ThermalModel* model, const ThermalDofMap* dofs,
    const ThermalConditions* conditions, const ControlSettings* controls,
    ThermalState* state, ResultWriter* writer, std::int64_t first_frame,
    double total_before, AnalysisReport* report, AnalysisProgress* progress,
    Diagnostic* diagnostic) {
  Diagnostic local;
  Diagnostic* d = diagnostic ? diagnostic : &local;
  if (report) report->Reset();
  const double dispatch_baseline = ElementDispatchSeconds();
  Diagnostic empty;
  CopyDiagnosticMetadata(&empty, logging_error_);
  auto finish = [&](StatusCode status) {
    RefreshDispatchReport(report, dispatch_baseline);
    if (status != StatusCode::kOk && state) state->DiscardCandidate();
    if (report) {
      report->status = status;
      if (logging_error_->code != StatusCode::kOk)
        report->logging_diagnostic = logging_error_;
    }
    return status;
  };
  if (!report || !state || !writer || !controls || !controls->temperature ||
      !controls->newton || first_frame < 1 || !std::isfinite(total_before) ||
      total_before < 0)
    return finish(Error("HeatTransferAnalysis SolveStep inputs", d));
  auto status = ValidateStep(model, dofs, conditions, settings_, d);
  if (status != StatusCode::kOk) return finish(status);
  status = state->ValidateMapping(dofs, d);
  if (status != StatusCode::kOk) return finish(status);
  report->step_number = conditions->StepNumber();
  report->equation_count = dofs->FreeCount();
  report->last_committed_frame = first_frame - 1;
  report->last_committed_step = state->ResultStepNumber();
  report->last_committed_increment = state->ResultIncrementNumber();
  report->last_committed_total_time = total_before;
  status = state->BeginStep(report->step_number, controls->temperature, d);
  if (status != StatusCode::kOk) return finish(status);
  try {
    Vector* raw = nullptr;
    status = Vector::Create(state->GlobalCount(), &raw, d);
    std::unique_ptr<Vector> step_start(raw);
    if (status == StatusCode::kOk)
      status = step_start->CopyFrom(state->CommittedTemperature(), d);
    if (status != StatusCode::kOk) return finish(status);
    Vector* raw_delta = nullptr;
    status = Vector::Create(state->GlobalCount(), &raw_delta, d);
    std::unique_ptr<Vector> accepted_delta(raw_delta);
    if (status != StatusCode::kOk) return finish(status);
    status = Vector::Create(state->GlobalCount(), &raw_delta, d);
    std::unique_ptr<Vector> staged_delta(raw_delta);
    if (status != StatusCode::kOk) return finish(status);
    double previous_dt = 0;
    BackwardEulerIntegrator* raw_integrator = nullptr;
    if (!settings_->steady)
      status = BackwardEulerIntegrator::Create(&raw_integrator, d);
    std::unique_ptr<BackwardEulerIntegrator> integrator(raw_integrator);
    if (status != StatusCode::kOk) return finish(status);
    ThermalIncrementController* raw_controller = nullptr;
    status = ThermalIncrementController::Create(settings_, &raw_controller, d);
    std::unique_ptr<ThermalIncrementController> controller(raw_controller);
    if (status != StatusCode::kOk) return finish(status);
    ThermalStepProblem problem(model, dofs, conditions, controls->temperature,
                               state, integrator.get(), settings_->period,
                               report);
    status = problem.Initialize(step_start.get(), d);
    if (status != StatusCode::kOk) return finish(status);
    PardisoSolver solver;
    ThermalTimedSolver timed_solver(&solver, report);
    IncrementAttemptInfo attempt;
    attempt.step_number = report->step_number;
    auto notify = [&](IncrementEvent event, StatusCode primary) {
      if (!progress) return primary;
      RefreshDispatchReport(report, dispatch_baseline);
      Diagnostic secondary;
      CopyDiagnosticMetadata(d, &secondary);
      StatusCode log_status = StatusCode::kOk;
      try {
        log_status = progress->OnThermalIncrementProgress(
            event, &attempt, nullptr, nullptr,
            primary == StatusCode::kOk ? nullptr : d, &secondary);
      } catch (const std::bad_alloc&) {
        log_status = model_storage::AllocationFailure(
            "Thermal increment progress allocation", &secondary);
      }
      if (log_status != StatusCode::kOk) {
        if (logging_error_->code == StatusCode::kOk) {
          CopyDiagnosticMetadata(&secondary, logging_error_);
          logging_error_->code = log_status;
        }
        if (primary == StatusCode::kOk) CopyDiagnosticMetadata(&secondary, d);
      }
      return primary != StatusCode::kOk ? primary : log_status;
    };
    auto common_progress = [&](AnalysisEvent event) {
      if (!progress) return StatusCode::kOk;
      RefreshDispatchReport(report, dispatch_baseline);
      Diagnostic secondary;
      StatusCode logged = StatusCode::kOk;
      try {
        logged = progress->OnProgress(event, report, &secondary);
      } catch (const std::bad_alloc&) {
        logged = model_storage::AllocationFailure(
            "Thermal common progress allocation", &secondary);
      }
      if (logged != StatusCode::kOk) {
        if (logging_error_->code == StatusCode::kOk) {
          CopyDiagnosticMetadata(&secondary, logging_error_);
          logging_error_->code = logged;
        }
        CopyDiagnosticMetadata(&secondary, d);
      }
      return logged;
    };
    while (controller->StepTime() < settings_->period) {
      if (controller->AcceptedCount() >= settings_->max_increments)
        return finish(Error("Heat transfer maximum increments reached", d,
                            StatusCode::kEquilibriumFailure));
      if (report->last_committed_frame ==
              std::numeric_limits<std::int64_t>::max() ||
          report->attempt_count == std::numeric_limits<int>::max())
        return finish(Error("Thermal frame/attempt counter overflow", d));
      const double time = controller->StepTime();
      double dt = 0, end = 0;
      status = controller->PlanIncrement(time, &dt, &end, d);
      if (status != StatusCode::kOk) return finish(status);
      if (!std::isfinite(end) || !(end > time) ||
          !std::isfinite(total_before + end))
        return finish(Error("Thermal physical time overflow", d,
                            StatusCode::kLibraryFailure));
      Diagnostic context;
      context.source_file = model->SourceDomain()->SourceFile();
      context.step = report->step_number;
      CopyDiagnosticMetadata(&context, d);
      ++report->attempt_count;
      report->increment_number = controller->AcceptedCount() + 1;
      attempt.increment_number = report->increment_number;
      attempt.attempt_number = controller->AttemptNumber();
      attempt.step_time = end;
      attempt.total_time = total_before + end;
      attempt.increment_size = dt;
      attempt.completes_step = end == settings_->period;
      attempt.proposed_available = false;
      attempt.thermal_change_available = false;
      attempt.thermal_free_change = 0;
      status = notify(IncrementEvent::kAttemptStarted, StatusCode::kOk);
      if (status != StatusCode::kOk) return finish(status);
      status = state->PrepareCandidate(d);
      if (status != StatusCode::kOk) return finish(status);
      double extrapolation = 0;
      status = controller->ExtrapolationFactor(previous_dt, &extrapolation, d);
      if (status == StatusCode::kOk && extrapolation != 0)
        status = state->CandidateTemperature()->AddScaled(
            extrapolation, accepted_delta.get(), d);
      if (status != StatusCode::kOk) return finish(status);
      for (const auto* row : *conditions->Temperatures()) {
        const int global = dofs->GlobalDof(row->node_id, row->label);
        double theta = 0;
        status = conditions->EvaluateTemperature(
            row, end, step_start->Data()[global], &theta, d);
        if (status == StatusCode::kOk)
          status = state->CandidateTemperature()->Set(global, theta, d);
        if (status != StatusCode::kOk) return finish(status);
      }
      problem.BeginAttempt(end, dt, extrapolation != 0);
      ThermalStepObserver observer(progress, &attempt, d, report,
                                   dispatch_baseline);
      const double numerical_start = WallTime();
      status = newton_->Solve(&problem, &timed_solver, controls->temperature,
                              controls->newton, iterations_, d,
                              progress ? &observer : nullptr);
      report->calculation_seconds += WallTime() - numerical_start;
      report->calculation_finished = true;
      report->calculation_status = iterations_->Info()->numerical_status;
      report->iteration_count = iterations_->Info()->iteration;
      if (iterations_->Info()->linear_solver_calls >
          std::numeric_limits<int>::max() - report->linear_solver_calls)
        return finish(Error("Thermal solve counter overflow", d));
      report->linear_solver_calls += iterations_->Info()->linear_solver_calls;
      report->residual_valid = iterations_->Info()->residual_valid;
      report->heat_flux_field =
          report->residual_valid ? iterations_->Residual() : nullptr;
      if (iterations_->Info()->logging_failed) {
        CopyDiagnosticMetadata(iterations_->LoggingDiagnostic(),
                               logging_error_);
        return finish(status);
      }
      double change = 0;
      if (status == StatusCode::kOk) {
        status = problem.FreeChange(&change, d);
        if (status == StatusCode::kOk) {
          attempt.thermal_change_available = true;
          attempt.thermal_free_change = change;
        }
      }
      const bool accuracy = status == StatusCode::kOk &&
                            settings_->has_deltmx && change > settings_->deltmx;
      if (accuracy)
        status = Error("Heat transfer DELTMX exceeded", d,
                       StatusCode::kEquilibriumFailure);
      if (status != StatusCode::kOk) {
        report->calculation_status = status;
        const auto notified = notify(IncrementEvent::kAttemptFailed, status);
        if (logging_error_->code != StatusCode::kOk) return finish(notified);
        if (status != StatusCode::kEquilibriumFailure) return finish(status);
        IncrementRejection cause = IncrementRejection::kTimeAccuracy;
        if (!accuracy) {
          const auto outcome = iterations_->Info()->outcome;
          cause = outcome == ThermalNewtonOutcome::kDivergence
                      ? IncrementRejection::kDivergence
                  : outcome == ThermalNewtonOutcome::kSlowConvergence
                      ? IncrementRejection::kSlowConvergence
                      : IncrementRejection::kIterationCap;
        }
        state->DiscardCandidate();
        IncrementRetry retry;
        status = controller->OnRejected(cause, dt, change, &retry, d);
        if (status != StatusCode::kOk) return finish(status);
        attempt.proposed_available = true;
        attempt.proposed_increment = retry.dt;
        status = notify(IncrementEvent::kCutback, StatusCode::kOk);
        if (status != StatusCode::kOk) return finish(status);
        continue;
      }
      status = common_progress(AnalysisEvent::kCalculationFinished);
      if (status != StatusCode::kOk) return finish(status);
      status = state->CandidateReactionFlux()->CopyFrom(problem.Reaction(), d);
      if (status == StatusCode::kOk)
        status = state->StageFirstAcceptedTemperature(d);
      if (status == StatusCode::kOk && !settings_->steady)
        status = state->IncrementRates()->Preview(change, dt, d);
      if (status != StatusCode::kOk) return finish(status);
      // Prepare the next predictor before any result writes. Rejected trials
      // never promote this buffer; successful flush/commit uses pointer swap.
      status = staged_delta->CopyFrom(state->CandidateTemperature(), d);
      if (status == StatusCode::kOk)
        status = staged_delta->AddScaled(-1, state->CommittedTemperature(), d);
      if (status != StatusCode::kOk) return finish(status);
      state->MarkResultsUnconfirmed(report->step_number,
                                    report->increment_number);
      status = WriteNodes(dofs, state, d);
      if (status != StatusCode::kOk) return finish(status);
      FrameInfo frame;
      frame.frame_number = report->last_committed_frame + 1;
      frame.step_number = report->step_number;
      frame.increment_number = report->increment_number;
      frame.step_name =
          model->SourceDomain()->Steps()->at(report->step_number)->Name();
      frame.step_time = end;
      frame.total_time = total_before + end;
      frame.completes_step = attempt.completes_step;
      frame.procedure = settings_->steady
                            ? FrameProcedure::kHeatTransferSteady
                            : FrameProcedure::kHeatTransferTransient;
      const double output_start = WallTime();
      status =
          writer->WriteThermalFrame(&frame, model->SourceDomain(), state, d);
      if (status != StatusCode::kOk) return finish(status);
      state->Commit(report->step_number, report->increment_number);
      accepted_delta.swap(staged_delta);
      previous_dt = dt;
      report->output_seconds += WallTime() - output_start;
      report->last_committed_step = report->step_number;
      report->last_committed_increment = report->increment_number;
      report->last_committed_frame = frame.frame_number;
      report->last_committed_total_time = frame.total_time;
      report->step_completed = frame.completes_step;
      status = common_progress(AnalysisEvent::kIncrementCommitted);
      if (status != StatusCode::kOk) return finish(status);
      status = controller->OnAccepted(
          dt, change, iterations_->Info()->iteration,
          iterations_->Residual()->branch, state->IncrementRates(), d);
      if (status != StatusCode::kOk) return finish(status);
      status = notify(IncrementEvent::kCommitted, StatusCode::kOk);
      if (status != StatusCode::kOk) return finish(status);
    }
    return finish(StatusCode::kOk);
  } catch (const std::bad_alloc&) {
    return finish(
        model_storage::AllocationFailure("Thermal SolveStep allocation", d));
  }
}
}  // namespace fesa
