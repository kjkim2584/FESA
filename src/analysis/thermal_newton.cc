#include "fesa/analysis/thermal_newton.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <new>

#include "analysis/field_decision.h"
#include "analysis/nonlinear_convergence.h"
#include "base/diagnostic_copy.h"
#include "fesa/io/control_settings.h"
#include "model/model_storage.h"
namespace fesa {
namespace {
class ThermalAttemptCleanup final {
 public:
  explicit ThermalAttemptCleanup(ThermalNonlinearProblem* problem)
      : problem_(problem) {}
  ~ThermalAttemptCleanup() {
    if (problem_) problem_->ReleaseAttemptResources();
  }
  ThermalAttemptCleanup(const ThermalAttemptCleanup&) = delete;
  ThermalAttemptCleanup& operator=(const ThermalAttemptCleanup&) = delete;
  ThermalAttemptCleanup(ThermalAttemptCleanup&&) = delete;
  ThermalAttemptCleanup& operator=(ThermalAttemptCleanup&&) = delete;

 private:
  ThermalNonlinearProblem* problem_;  ///< Borrowed until enclosing Solve ends.
};
StatusCode Emit(ThermalNewtonProgress* progress, ThermalIterationReport* report,
                ThermalIterationInfo* info, Diagnostic* d) {
  if (!progress) return info->status;
  Diagnostic secondary;
  const bool numeric_failure = info->numerical_status != StatusCode::kOk;
  if (!numeric_failure) CopyDiagnosticMetadata(d, &secondary);
  StatusCode status = StatusCode::kOk;
  try {
    status = progress->OnIteration(report, &secondary);
  } catch (const std::bad_alloc&) {
    status = model_storage::AllocationFailure(
        "ThermalNewton progress allocation", &secondary);
  }
  if (status != StatusCode::kOk) {
    report->RecordLoggingFailure(status, &secondary);
    if (!numeric_failure) {
      CopyDiagnosticMetadata(&secondary, d);
      info->status = status;
      info->outcome = ThermalNewtonOutcome::kFatal;
    }
  }
  return info->status;
}
bool Nonnegative(double value) { return std::isfinite(value) && value >= 0; }
bool ValidControls(const FieldControls* c) {
  if (!c) return false;
  const double ratios[] = {c->rn, c->rp, c->cn,        c->epsilon,  c->cepsilon,
                           c->rl, c->cf, c->epsilon_l, c->epsilon_d};
  for (double value : ratios)
    if (!Nonnegative(value)) return false;
  return c->qu_explicit ? std::isfinite(c->qu) && c->qu > 0
                        : std::isfinite(c->q0) && c->q0 > 0;
}
bool ValidMeasurements(const HeatFluxResidual* f) {
  return f->dof_count >= 0 && Nonnegative(f->maximum) &&
         Nonnegative(f->spatial_mean) && Nonnegative(f->increment_max) &&
         Nonnegative(f->correction_max) &&
         (f->dof_count == 0 ||
          (std::isfinite(f->time_mean_used) && f->time_mean_used > 0));
}
// formulation7: correction estimate a_i/min(a_(i-1),a_(i-2))*c_i.
// At i1 only a0 is available. Exponent-separated scalar arithmetic avoids
// spurious intermediate overflow; a positive underflow is unavailable, not0.
bool Estimate(double current, double previous, double older, double correction,
              double* out) {
  const double denominator = std::min(previous, older);
  if (!(denominator > 0)) return false;
  if (current == 0 || correction == 0) {
    *out = 0;
    return true;
  }
  int ec = 0, er = 0, ed = 0;
  const double mc = std::frexp(current, &ec), mr = std::frexp(correction, &er),
               md = std::frexp(denominator, &ed);
  const double value = std::scalbn(mc * mr / md, ec + er - ed);
  if (!std::isfinite(value) || !(value > 0)) return false;
  *out = value;
  return true;
}
double Predict(int i, double current, double previous, double limit,
               double used_mean) {
  if (!(current > 0) || !(previous > 0) || current >= previous || limit == 0)
    return std::numeric_limits<double>::infinity();
  const double rate = current > previous / 2
                          ? std::log1p((current - previous) / previous)
                          : std::log(current) - std::log(previous);
  if (!(rate < 0)) return std::numeric_limits<double>::infinity();
  const double target =
      std::log(limit) + std::log(used_mean) - std::log(current);
  return static_cast<double>(i) + std::ceil(std::max(0.0, target / rate));
}
ThermalResidualBranch ThermalBranch(ScalarFieldBranch branch) {
  switch (branch) {
    case ScalarFieldBranch::kNoDofs:
      return ThermalResidualBranch::kNoDofs;
    case ScalarFieldBranch::kZeroFlux:
      return ThermalResidualBranch::kZeroFlux;
    case ScalarFieldBranch::kLinear:
      return ThermalResidualBranch::kLinear;
    case ScalarFieldBranch::kGeneral:
      return ThermalResidualBranch::kGeneral;
  }
  return ThermalResidualBranch::kNoDofs;
}
}  // namespace
StatusCode ThermalIterationReport::Create(ThermalIterationReport** out,
                                          Diagnostic* d) {
  if (!out) return model_storage::Invalid("Thermal report output", d);
  *out = nullptr;
  std::unique_ptr<ThermalIterationReport> owner(new (std::nothrow)
                                                    ThermalIterationReport);
  if (!owner) return model_storage::AllocationFailure("Thermal report", d);
  owner->info_ = new (std::nothrow) ThermalIterationInfo;
  owner->residual_ = new (std::nothrow) HeatFluxResidual;
  owner->logging_diagnostic_ = new (std::nothrow) Diagnostic;
  if (!owner->info_ || !owner->residual_ || !owner->logging_diagnostic_)
    return model_storage::AllocationFailure("Thermal report records", d);
  *out = owner.release();
  return StatusCode::kOk;
}
ThermalIterationReport::~ThermalIterationReport() {
  delete info_;
  delete residual_;
  delete logging_diagnostic_;
}
void ThermalIterationReport::Reset() noexcept {
  Diagnostic empty;
  CopyDiagnosticMetadata(&empty, logging_diagnostic_);
  info_->numerical_status = info_->logging_status = StatusCode::kOk;
  info_->logging_failed = false;
  info_->iteration = info_->linear_solver_calls = info_->residual_evaluations =
      0;
  info_->phi = info_->correction_max = info_->increment_max = 0;
  info_->branch = ThermalResidualBranch::kNoDofs;
  info_->residual_valid = info_->estimate_available =
      info_->prediction_available = false;
  info_->estimate_max = info_->predicted_total = 0;
  info_->status = StatusCode::kOk;
  info_->outcome = ThermalNewtonOutcome::kNotStarted;
  residual_->dof_count = residual_->contribution_count =
      residual_->active_contribution_count = 0;
  residual_->maximum_node = residual_->maximum_label = 0;
  residual_->maximum = residual_->spatial_mean = residual_->maximum_flux = 0;
  residual_->phi = residual_->increment_max = residual_->correction_max = 0;
  residual_->time_mean_used = .01;
  residual_->branch = ThermalResidualBranch::kNoDofs;
  residual_->passed = false;
}
const ThermalIterationInfo* ThermalIterationReport::Info() const {
  return info_;
}
const HeatFluxResidual* ThermalIterationReport::Residual() const {
  return residual_;
}
const Diagnostic* ThermalIterationReport::LoggingDiagnostic() const {
  return info_->logging_failed ? logging_diagnostic_ : nullptr;
}
void ThermalIterationReport::RecordLoggingFailure(
    StatusCode status, const Diagnostic* d) noexcept {
  if (status == StatusCode::kOk || info_->logging_failed) return;
  CopyDiagnosticMetadata(d, logging_diagnostic_);
  logging_diagnostic_->code = status;
  info_->logging_status = status;
  info_->logging_failed = true;
}
StatusCode ThermalNewton::Create(ThermalNewton** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("ThermalNewton output", d);
  *out = new (std::nothrow) ThermalNewton;
  return *out ? StatusCode::kOk
              : model_storage::AllocationFailure("ThermalNewton", d);
}
StatusCode ThermalNewton::Solve(ThermalNonlinearProblem* problem,
                                LinearSolver* solver,
                                const FieldControls* controls,
                                const NewtonControls* counts,
                                ThermalIterationReport* report, Diagnostic* d,
                                ThermalNewtonProgress* progress) {
  ThermalAttemptCleanup cleanup(problem);
  if (report) report->Reset();
  auto fail = [&](StatusCode status) {
    if (report) {
      report->info_->status = report->info_->numerical_status = status;
      report->info_->outcome = ThermalNewtonOutcome::kFatal;
    }
    return report ? Emit(progress, report, report->info_, d) : status;
  };
  if (!problem || !solver || !ValidControls(controls) || !counts || !report ||
      counts->divergence_start < 3 || counts->slow_start <= 0 ||
      counts->nonquadratic_start <= 0 || counts->prediction_limit <= 0 ||
      counts->iteration_limit <= 0)
    return fail(model_storage::Invalid("ThermalNewton bindings/controls", d));
  auto* info = report->info_;
  auto* flux = report->residual_;
  info->outcome = ThermalNewtonOutcome::kInProgress;
  double previous = 0, older = 0;
  try {
    for (int i = 0;; ++i) {
      info->iteration = i;
      info->residual_valid = info->estimate_available =
          info->prediction_available = false;
      info->estimate_max = info->predicted_total = 0;
      ++info->residual_evaluations;
      auto status = problem->EvaluateResidual(i, flux, d);
      problem->ReleaseCorrectionAssessment();
      if (status != StatusCode::kOk) return fail(status);
      if (!ValidMeasurements(flux))
        return fail(model_storage::Invalid("ThermalNewton measurements", d));
      FieldConvergenceOptions options;
      options.correction_available = i > 0;
      options.residual_ratio =
          i >= counts->nonquadratic_start ? controls->rp : controls->rn;
      if (i > 0 && flux->dof_count > 0) {
        options.estimate_available =
            Estimate(flux->maximum, previous, i > 1 ? older : previous,
                     flux->correction_max, &options.estimate_max);
      }
      ScalarFieldDecisionInput input;
      input.dof_count = flux->dof_count;
      input.maximum = flux->maximum;
      input.spatial_mean = flux->spatial_mean;
      input.used_mean = flux->time_mean_used;
      input.increment_max = flux->increment_max;
      input.correction_max = flux->correction_max;
      ScalarFieldDecision decision;
      status = EvaluateScalarFieldDecision(controls, &options, false, 0, &input,
                                           &decision, d);
      if (status != StatusCode::kOk &&
          status != StatusCode::kEquilibriumFailure)
        return fail(status);
      flux->phi = decision.phi;
      flux->branch = ThermalBranch(decision.branch);
      flux->passed = decision.passed;
      info->phi = flux->phi;
      info->correction_max = flux->correction_max;
      info->increment_max = flux->increment_max;
      info->branch = flux->branch;
      info->residual_valid = true;
      info->estimate_available = options.estimate_available;
      info->estimate_max = options.estimate_max;
      if (decision.passed && (i > 0 || flux->dof_count == 0)) {
        info->outcome = ThermalNewtonOutcome::kConverged;
        return Emit(progress, report, info, d);
      }
      const bool diverging = i >= counts->divergence_start &&
                             std::min(flux->maximum, previous) > older;
      bool slow = false;
      if (i >= counts->slow_start && flux->phi > decision.residual_limit) {
        info->prediction_available = true;
        info->predicted_total =
            Predict(i, flux->maximum, previous, decision.residual_limit,
                    flux->time_mean_used);
        slow = info->predicted_total > counts->prediction_limit;
      }
      if (diverging || slow || i >= counts->iteration_limit) {
        info->outcome = diverging ? ThermalNewtonOutcome::kDivergence
                        : slow    ? ThermalNewtonOutcome::kSlowConvergence
                                  : ThermalNewtonOutcome::kIterationLimit;
        info->status = info->numerical_status = ReportFailure(
            StatusCode::kEquilibriumFailure, DiagnosticStage::kComputation,
            diverging ? "ThermalNewton diverging residual"
            : slow    ? "ThermalNewton slow residual convergence"
                      : "ThermalNewton iteration limit",
            d);
        return Emit(progress, report, info, d);
      }
      const auto notified = Emit(progress, report, info, d);
      if (notified != StatusCode::kOk) return notified;
      older = previous;
      previous = flux->maximum;
      const EquationSystem* system = nullptr;
      Vector* correction = nullptr;
      status = problem->AssembleTangent(&system, &correction, d);
      if (status != StatusCode::kOk) return fail(status);
      if (!system || !correction || correction == system->Rhs() ||
          correction->Size() != system->Rhs()->Size())
        return fail(
            model_storage::Invalid("ThermalNewton correction system", d));
      ++info->linear_solver_calls;
      status = solver->Solve(system, correction, d);
      problem->ReleaseEquationSystem();
      if (status != StatusCode::kOk) return fail(status);
      status = problem->ApplyCorrection(correction, d);
      problem->ReleaseFreeCorrection();
      if (status != StatusCode::kOk) return fail(status);
    }
  } catch (const std::bad_alloc&) {
    return fail(
        model_storage::AllocationFailure("ThermalNewton allocation", d));
  }
}
}  // namespace fesa
