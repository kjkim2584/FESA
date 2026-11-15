#include "fesa/analysis/thermal_increment_controller.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <new>

#include "analysis/thermal_increment_policy.h"
#include "model/model_storage.h"
namespace fesa {
namespace {
bool Positive(double value) { return std::isfinite(value) && value > 0; }
StatusCode Invalid(Diagnostic* d) {
  return model_storage::Invalid("Thermal increment metadata", d);
}
StatusCode Arithmetic(Diagnostic* d) {
  return ReportFailure(StatusCode::kLibraryFailure,
                       DiagnosticStage::kComputation,
                       "Thermal increment arithmetic", d);
}
}  // namespace
StatusCode PlanThermalIncrement(double period, double time, double proposal,
                                int accepted, double* dt, double* end,
                                Diagnostic* d) {
  if (!dt || !end || dt == end) return Invalid(d);
  if (!Positive(period) || !Positive(proposal) || !std::isfinite(time) ||
      time < 0 || time > period || accepted < 0)
    return Arithmetic(d);
  const double remaining = period - time;
  if (remaining == 0) {
    *dt = 0;
    *end = period;
    return StatusCode::kOk;
  }
  // Higham recursive-sum bound, with input representation terms; scalar
  // coefficient arithmetic only. int accepted keeps n*u strictly below1.
  const double u = std::numeric_limits<double>::epsilon() / 2;
  const double nu = (static_cast<double>(accepted) + 1) * u;
  if (remaining <= (nu / (1 - nu)) * period) return Arithmetic(d);
  double actual = std::min(proposal, remaining);
  double next = actual == remaining ? period : time + actual;
  if (!std::isfinite(next) || !(next > time) || next > period)
    return Arithmetic(d);
  const double next_nu = (static_cast<double>(accepted) + 2) * u;
  if (next < period && period - next <= (next_nu / (1 - next_nu)) * period) {
    actual = remaining;
    next = period;
  }
  *dt = actual;
  *end = next;
  return StatusCode::kOk;
}
StatusCode CreateThermalIncrementControllerWithPolicy(
    const HeatTransferStepParameters* p, const ThermalIncrementPolicy* policy,
    ThermalIncrementController** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("Thermal controller output", d);
  *out = nullptr;
  if (!policy || !std::isfinite(policy->dl) || policy->dl <= 0 ||
      policy->dl > 1)
    return model_storage::Invalid("Thermal increment verification policy", d);
  std::unique_ptr<ThermalIncrementController> owner(
      new (std::nothrow) ThermalIncrementController);
  if (!owner) return model_storage::AllocationFailure("Thermal controller", d);
  const auto status = CloneHeatTransferStepParameters(p, &owner->settings_, d);
  if (status != StatusCode::kOk) return status;
  owner->nominal_ = p->initial_increment;
  owner->dl_ = policy->dl;
  *out = owner.release();
  return StatusCode::kOk;
}
StatusCode ThermalIncrementController::Create(
    const HeatTransferStepParameters* p, ThermalIncrementController** out,
    Diagnostic* d) {
  ThermalIncrementPolicy policy;
  return CreateThermalIncrementControllerWithPolicy(p, &policy, out, d);
}
ThermalIncrementController::~ThermalIncrementController() { delete settings_; }
bool ThermalIncrementController::IsAutomatic() const {
  return settings_->steady
             ? !(settings_->has_minimum && settings_->has_maximum &&
                 settings_->initial_increment == settings_->minimum &&
                 settings_->initial_increment == settings_->maximum)
             : settings_->has_deltmx;
}
double ThermalIncrementController::NextIncrement(double step_time) const {
  double dt = 0, end = 0;
  return PlanIncrement(step_time, &dt, &end, nullptr) == StatusCode::kOk
             ? dt
             : std::numeric_limits<double>::quiet_NaN();
}
StatusCode ThermalIncrementController::PlanIncrement(double step_time,
                                                     double* dt, double* end,
                                                     Diagnostic* d) const {
  if (!std::isfinite(step_time) || step_time != time_) return Invalid(d);
  return PlanThermalIncrement(settings_->period, time_, nominal_, accepted_, dt,
                              end, d);
}
double ThermalIncrementController::StepTime() const { return time_; }
int ThermalIncrementController::AcceptedCount() const { return accepted_; }
int ThermalIncrementController::AttemptNumber() const { return attempt_; }
StatusCode ThermalIncrementController::ExtrapolationFactor(
    double previous_dt, double* factor, Diagnostic* d) const {
  const double actual = NextIncrement(time_);
  if (!factor || !std::isfinite(previous_dt) || previous_dt < 0 ||
      !Positive(actual))
    return Invalid(d);
  double value = 0;
  if (previous_dt > 0 && nominal_ / previous_dt >= .1)
    value = actual / previous_dt;
  if (!std::isfinite(value)) return Arithmetic(d);
  *factor = value;
  return StatusCode::kOk;
}
StatusCode ThermalIncrementController::OnAccepted(double dt, double change,
                                                  int iterations,
                                                  ThermalResidualBranch branch,
                                                  const RateHistory* committed,
                                                  Diagnostic* d) {
  if (!Positive(dt) || dt != NextIncrement(time_) || !std::isfinite(change) ||
      change < 0 || iterations < 0 ||
      accepted_ == std::numeric_limits<int>::max() ||
      (branch != ThermalResidualBranch::kNoDofs &&
       branch != ThermalResidualBranch::kZeroFlux &&
       branch != ThermalResidualBranch::kLinear &&
       branch != ThermalResidualBranch::kGeneral))
    return Invalid(d);
  const bool automatic_transient = IsAutomatic() && !settings_->steady;
  if (automatic_transient &&
      (change > settings_->deltmx || !committed || committed->CanCommit() ||
       committed->AcceptedCount() != accepted_ + 1))
    return Invalid(d);
  double actual = 0, next_time = 0;
  const auto plan_status = PlanIncrement(time_, &actual, &next_time, d);
  if (plan_status != StatusCode::kOk) return plan_status;
  double proposal = nominal_;
  bool next_floor_attempt = false;
  const int accuracy_accepted =
      accuracy_rejected_ ? 1 : std::min(3, accuracy_accepted_ + 1);
  if (next_time < settings_->period && IsAutomatic()) {
    if (iterations > 10) {
      // KEY-CT IL=10 / DB=.75: an expensive successful increment reduces
      // the next proposal, independently of growth eligibility.
      const double minimum = settings_->EffectiveMinimum();
      proposal = std::min(.75 * dt, settings_->period - next_time);
      if (settings_->has_maximum)
        proposal = std::min(proposal, settings_->maximum);
      next_floor_attempt = proposal <= minimum;
      // Clamp the stored proposal last. PlanIncrement caps actual dt, so DE
      // still sees the nominal proposal and R<nominal preserves cap retries.
      proposal = std::max(proposal, minimum);
    } else if (settings_->steady) {
      // Two actual accepted iterations, including a legitimate no-DOF0.
      // Task16: IG bounds each count, not their sum; pairs slide.
      if (accepted_ > 0 && previous_iterations_ <= 4 && iterations <= 4) {
        proposal = std::min(1.5 * dt, settings_->period - next_time);
        if (settings_->has_maximum)
          proposal = std::min(proposal, settings_->maximum);
      }
    } else {
      bool eligible = false;
      auto status =
          committed->AllRecentBelow(dt, .75 * settings_->deltmx, &eligible, d);
      if (status != StatusCode::kOk) return status;
      if (eligible && accepted_ > 0 && previous_iterations_ <= 4 &&
          iterations <= 4 && accuracy_accepted >= 3) {
        // G-NL-FACTOR: accepted attempt n>=2 uses min(DD,DMdiff)=1.5;
        // n<=1 (including no-DOF0) retains the diffusion bound2.
        // An infinite accuracy ratio or S=0 has the same bounded limit.
        const double maximum_factor = iterations >= 2 ? 1.5 : 2.;
        const double factor =
            change > 0
                ? std::min(.8 * (settings_->deltmx / change), maximum_factor)
                : maximum_factor;
        const double unbounded = factor * dt;
        double bounded = std::min(unbounded, settings_->period - next_time);
        if (settings_->has_maximum)
          bounded = std::min(bounded, settings_->maximum);
        if (iterations >= 2 || branch == ThermalResidualBranch::kGeneral ||
            bounded >= dl_ * 2 * dt ||
            (settings_->has_maximum && unbounded >= settings_->maximum))
          proposal = bounded;
      }
    }
  }
  if (!Positive(proposal)) return Arithmetic(d);
  time_ = next_time;
  nominal_ = proposal;
  ++accepted_;
  attempt_ = 1;
  floor_attempt_ = next_floor_attempt;
  previous_iterations_ = iterations;
  accuracy_accepted_ = accuracy_accepted;
  accuracy_rejected_ = false;
  return StatusCode::kOk;
}
StatusCode ThermalIncrementController::OnRejected(IncrementRejection cause,
                                                  double dt, double measure,
                                                  IncrementRetry* out,
                                                  Diagnostic* d) {
  if (!out || !Positive(dt) || dt != NextIncrement(time_) ||
      !std::isfinite(measure) || measure < 0)
    return Invalid(d);
  double factor = 0;
  switch (cause) {
    case IncrementRejection::kTimeAccuracy:
      if (settings_->steady || !settings_->has_deltmx ||
          measure <= settings_->deltmx)
        return Invalid(d);
      factor = .85 * (settings_->deltmx / measure);
      break;
    case IncrementRejection::kDivergence:
      factor = .25;
      break;
    case IncrementRejection::kSlowConvergence:
    case IncrementRejection::kIterationCap:
      factor = .5;
      break;
    default:
      return Invalid(d);
  }
  auto terminal = [&](IncrementFailure reason, const char* operation) {
    out->dt = 0;
    out->failure = reason;
    return ReportFailure(StatusCode::kEquilibriumFailure,
                         DiagnosticStage::kComputation, operation, d);
  };
  if (!IsAutomatic())
    return terminal(IncrementFailure::kFixed,
                    "Fixed thermal increment rejected");
  const double minimum = settings_->EffectiveMinimum();
  double proposal = factor * dt;
  const bool below = proposal < minimum;
  const bool capped = settings_->period - time_ < nominal_;
  if (below && floor_attempt_ && !capped)
    return terminal(IncrementFailure::kMinimum,
                    "Thermal minimum increment exhausted");
  if (attempt_ >= 5)
    return terminal(IncrementFailure::kAttemptLimit,
                    "Thermal increment attempt limit");
  if (below) proposal = minimum;
  if (!Positive(proposal)) return Arithmetic(d);
  double retry = 0, end = 0;
  const auto status = PlanThermalIncrement(settings_->period, time_, proposal,
                                           accepted_, &retry, &end, d);
  if (status != StatusCode::kOk) return status;
  out->dt = retry;
  out->failure = IncrementFailure::kNone;
  nominal_ = proposal;
  floor_attempt_ = below;
  if (cause == IncrementRejection::kTimeAccuracy) accuracy_rejected_ = true;
  ++attempt_;
  return StatusCode::kOk;
}
}  // namespace fesa
