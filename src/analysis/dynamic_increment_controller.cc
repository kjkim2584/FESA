#include "fesa/analysis/dynamic_increment_controller.h"

#include <algorithm>
#include <cmath>

#include "analysis/direct_increment.h"

namespace fesa {
namespace {
bool Positive(double value) { return std::isfinite(value) && value > 0; }
bool Reduction(double value) { return Positive(value) && value < 1; }
StatusCode Invalid(Diagnostic* d) {
  return ReportFailure(StatusCode::kInvalidArgument,
                       DiagnosticStage::kValidation,
                       "Automatic increment metadata", d);
}
StatusCode Exhausted(const char* operation, Diagnostic* d) {
  return ReportFailure(StatusCode::kEquilibriumFailure,
                       DiagnosticStage::kComputation, operation, d);
}
StatusCode Primary(const NewtonIterationInfo* attempt, Diagnostic* d) {
  if (attempt->numerical_status != StatusCode::kOk)
    return attempt->numerical_status;
  if (attempt->status != StatusCode::kOk) return attempt->status;
  return Invalid(d);
}
void CopyClock(const IncrementClock* from, IncrementClock* to) {
  to->time = from->time;
  to->compensation = from->compensation;
  to->nominal = from->nominal;
  to->last_increment = from->last_increment;
  to->successful = from->successful;
  to->quick_successes = from->quick_successes;
}
}  // namespace

DynamicIncrementController::DynamicIncrementController(
    const DynamicStepSettings* settings) : settings_(settings) {}
StatusCode DynamicIncrementController::Initialize(IncrementClock* out,
                                                   Diagnostic* d) const {
  if (!settings_ || !settings_->IsAutomatic() || !out) return Invalid(d);
  IncrementClock initial;
  initial.nominal = settings_->NominalIncrement();
  CopyClock(&initial, out);
  return StatusCode::kOk;
}
StatusCode DynamicIncrementController::Validate(
    const IncrementControls* controls, const IncrementClock* clock,
    Diagnostic* d) const {
  if (!settings_ || !settings_->IsAutomatic() || !controls || !clock ||
      controls->attempt_limit <= 0 || controls->slow_success_iterations < 0 ||
      controls->quick_iterations < 0 ||
      !Reduction(controls->divergence_cutback) ||
      !Reduction(controls->slow_cutback) ||
      !Reduction(controls->geometry_cutback) ||
      !Positive(controls->success_scale) || controls->success_scale > 1 ||
      !std::isfinite(controls->growth) || controls->growth < 1 ||
      !std::isfinite(controls->maximum_growth) ||
      controls->maximum_growth < 1 || !std::isfinite(clock->time) ||
      clock->time < 0 || clock->time >= settings_->Parameters()->period ||
      !std::isfinite(clock->compensation) || !Positive(clock->nominal) ||
      !std::isfinite(clock->last_increment) || clock->last_increment < 0 ||
      clock->successful < 0 || clock->quick_successes < 0 ||
      clock->quick_successes > 2)
    return Invalid(d);
  return StatusCode::kOk;
}
StatusCode DynamicIncrementController::Propose(
    const IncrementControls* controls, const IncrementClock* clock,
    double nominal, IncrementTarget* out, Diagnostic* d) const {
  auto status = Validate(controls, clock, d);
  if (status != StatusCode::kOk) return status;
  if (!out || !Positive(nominal)) return Invalid(d);
  nominal = std::min(nominal, settings_->MaximumIncrement());
  if (nominal < settings_->MinimumIncrement())
    return Exhausted("Automatic minimum increment exhausted", d);
  IncrementClock bounded;
  CopyClock(clock, &bounded);
  // min may exceed the period; endpoint clipping remains legal in that case.
  bounded.nominal = std::min(nominal, settings_->Parameters()->period);
  status = ProposeDirectIncrement(settings_->Parameters()->period,
                                  settings_->Parameters()->increment_limit,
                                  &bounded, out, d);
  if (status == StatusCode::kOk) out->nominal = nominal;
  return status;
}
StatusCode DynamicIncrementController::ValidateTarget(
    const IncrementControls* controls, const IncrementClock* clock,
    const IncrementTarget* target, Diagnostic* d) const {
  if (!target) return Invalid(d);
  IncrementTarget expected;
  auto status = Propose(controls, clock, target->nominal, &expected, d);
  if (status != StatusCode::kOk) return status;
  if (target->base_time != expected.base_time ||
      target->base_compensation != expected.base_compensation ||
      target->base_successful != expected.base_successful ||
      target->nominal != expected.nominal || target->size != expected.size ||
      target->time != expected.time ||
      target->compensation != expected.compensation ||
      target->completes_step != expected.completes_step ||
      target->use_extrapolation || target->extrapolation_ratio != 0)
    return Invalid(d);
  return StatusCode::kOk;
}
StatusCode DynamicIncrementController::Cutback(
    const IncrementControls* controls, const IncrementClock* clock,
    const IncrementTarget* target, const NewtonIterationInfo* attempt,
    int attempt_number, double* retry_nominal, Diagnostic* d) const {
  auto status = ValidateTarget(controls, clock, target, d);
  if (status != StatusCode::kOk) return status;
  if (!attempt || !retry_nominal || attempt_number <= 0) return Invalid(d);
  if (attempt->logging_failed || attempt->outcome == NewtonOutcome::kFatal ||
      attempt->status != StatusCode::kEquilibriumFailure ||
      attempt->numerical_status != StatusCode::kEquilibriumFailure)
    return Primary(attempt, d);
  double factor = 0;
  switch (attempt->outcome) {
    case NewtonOutcome::kDivergence:
      factor = controls->divergence_cutback;
      break;
    case NewtonOutcome::kSlowConvergence:
    case NewtonOutcome::kIterationLimit:
      factor = controls->slow_cutback;
      break;
    case NewtonOutcome::kGeometry:
    case NewtonOutcome::kRotationBranch:
      factor = controls->geometry_cutback;
      break;
    default:
      return Invalid(d);
  }
  if (attempt_number >= controls->attempt_limit)
    return Exhausted("Automatic attempt limit exhausted", d);
  const double retry = target->size * factor;
  if (!Positive(retry) || retry < settings_->MinimumIncrement())
    return Exhausted("Automatic cutback below minimum increment", d);
  *retry_nominal = retry;
  return StatusCode::kOk;
}
StatusCode DynamicIncrementController::Accept(
    const IncrementControls* controls, const IncrementClock* clock,
    const IncrementTarget* target, const NewtonIterationInfo* attempt,
    IncrementClock* out, Diagnostic* d) const {
  auto status = ValidateTarget(controls, clock, target, d);
  if (status != StatusCode::kOk) return status;
  if (!out || out == clock || (attempt && attempt->iteration < 0))
    return Invalid(d);
  if (attempt &&
      (attempt->logging_failed || attempt->status != StatusCode::kOk ||
       attempt->numerical_status != StatusCode::kOk ||
       attempt->outcome != NewtonOutcome::kConverged))
    return Primary(attempt, d);
  const int quick_limit = controls->quick_iterations == 0
                              ? 4 : controls->quick_iterations;
  const int slow_limit = controls->slow_success_iterations == 0
                             ? 10 : controls->slow_success_iterations;
  const bool quick = !attempt || attempt->iteration <= quick_limit;
  const bool slow = attempt && attempt->iteration > slow_limit;
  IncrementClock next;
  next.time = target->time;
  next.compensation = target->compensation;
  next.last_increment = target->size;
  next.successful = clock->successful + 1;
  next.quick_successes = quick ? std::min(2, clock->quick_successes + 1) : 0;
  next.nominal = target->size;
  if (slow) {
    next.nominal *= controls->success_scale;
  } else if (next.quick_successes == 2) {
    const double growth = std::min(controls->growth, controls->maximum_growth);
    const double cap = settings_->MaximumIncrement();
    next.nominal = next.nominal >= cap / growth ? cap : next.nominal * growth;
  }
  if (!Positive(next.nominal))
    return ReportFailure(StatusCode::kLibraryFailure,
                         DiagnosticStage::kComputation,
                         "Unrepresentable next automatic increment", d);
  CopyClock(&next, out);
  return StatusCode::kOk;
}
}  // namespace fesa
