#include "analysis/increment_control.h"

#include <algorithm>
#include <cmath>

#include "model/model_storage.h"

namespace fesa {
namespace {
bool Positive(double value) { return std::isfinite(value) && value > 0; }
bool Reduction(double value) { return Positive(value) && value < 1; }
bool Growth(double value) { return std::isfinite(value) && value >= 1; }
StatusCode Exhausted(const char* operation, Diagnostic* d) {
  return ReportFailure(StatusCode::kEquilibriumFailure,
                       DiagnosticStage::kComputation, operation, d);
}
StatusCode Arithmetic(const char* operation, Diagnostic* d) {
  return ReportFailure(StatusCode::kLibraryFailure,
                       DiagnosticStage::kComputation, operation, d);
}
StatusCode Validate(const NonlinearStaticParameters* p,
                    const IncrementControls* c, const IncrementClock* clock,
                    Diagnostic* d) {
  auto status = ValidateNonlinearStaticParameters(p, d);
  if (status != StatusCode::kOk) return status;
  if (!c || !clock || c->attempt_limit <= 0 || c->slow_success_iterations < 0 ||
      c->quick_iterations < 0 || !Reduction(c->divergence_cutback) ||
      !Reduction(c->slow_cutback) || !Reduction(c->geometry_cutback) ||
      !Positive(c->success_scale) || c->success_scale > 1 ||
      !Growth(c->growth) || !Growth(c->maximum_growth) ||
      !Growth(c->endpoint_stretch) || !std::isfinite(c->extrapolation_cutoff) ||
      c->extrapolation_cutoff < 0 || !std::isfinite(clock->time) ||
      clock->time < 0 || clock->time > p->period ||
      !std::isfinite(clock->compensation) || !Positive(clock->nominal) ||
      !std::isfinite(clock->last_increment) || clock->last_increment < 0 ||
      clock->successful < 0 || clock->quick_successes < 0 ||
      clock->quick_successes > 2) {
    return model_storage::Invalid("Increment controls/clock", d);
  }
  return StatusCode::kOk;
}
bool Bound(const IncrementClock* clock, const IncrementTarget* target) {
  return target && target->base_time == clock->time &&
         target->base_compensation == clock->compensation &&
         target->base_successful == clock->successful &&
         Positive(target->nominal) && Positive(target->size) &&
         std::isfinite(target->time) && target->time > clock->time &&
         std::isfinite(target->compensation);
}
StatusCode Primary(const NewtonIterationInfo* attempt, Diagnostic* d) {
  if (attempt->numerical_status != StatusCode::kOk)
    return attempt->numerical_status;
  if (attempt->status != StatusCode::kOk) return attempt->status;
  return model_storage::Invalid("Inconsistent Newton attempt outcome", d);
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

StatusCode InitializeIncrementClock(const NonlinearStaticParameters* p,
                                    IncrementClock* out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("Initialize increment clock", d);
  auto status = ValidateNonlinearStaticParameters(p, d);
  if (status != StatusCode::kOk) return status;
  IncrementClock initial;
  initial.nominal = p->initial_increment;
  CopyClock(&initial, out);
  return StatusCode::kOk;
}
StatusCode ProposeIncrement(const NonlinearStaticParameters* p,
                            const IncrementControls* c,
                            const IncrementClock* clock, double nominal,
                            IncrementTarget* out, Diagnostic* d) {
  auto status = Validate(p, c, clock, d);
  if (status != StatusCode::kOk) return status;
  if (!out || !Positive(nominal) || clock->time == p->period)
    return model_storage::Invalid("Propose increment", d);
  if (clock->successful >= p->increment_limit)
    return Exhausted("STATIC INC exhausted", d);
  if (p->has_maximum) nominal = std::min(nominal, p->maximum_increment);
  if (nominal < p->minimum_increment)
    return Exhausted("STATIC minimum increment exhausted", d);

  IncrementTarget target;
  target.base_time = clock->time;
  target.base_compensation = clock->compensation;
  target.base_successful = clock->successful;
  target.nominal = nominal;
  target.size = nominal;
  const double remaining = p->period - clock->time;
  // Division avoids overflow of DT*h. An overflowing ratio cannot pass a
  // finite DT, so a very small h still follows the ordinary compensated sum.
  const bool stretch = remaining / nominal <= c->endpoint_stretch &&
                       (!p->has_maximum || remaining <= p->maximum_increment);
  const double adjusted = nominal - clock->compensation;
  const double proposed = clock->time + adjusted;
  if (remaining <= nominal || stretch || proposed >= p->period) {
    target.time = p->period;
    target.size = remaining;
    target.completes_step = true;
  } else {
    target.time = proposed;
    target.compensation = (proposed - clock->time) - adjusted;
  }
  if (!std::isfinite(target.time) || target.time <= clock->time ||
      !std::isfinite(target.compensation))
    return Arithmetic("STATIC time cannot advance representably", d);
  if (p->extrapolate && clock->successful > 0) {
    if (!Positive(clock->last_increment))
      return model_storage::Invalid("Missing accepted predictor interval", d);
    const double ratio = target.size / clock->last_increment;
    if (!Positive(ratio))
      return Arithmetic("Unrepresentable STATIC predictor ratio", d);
    target.use_extrapolation = ratio > c->extrapolation_cutoff;
    if (target.use_extrapolation) target.extrapolation_ratio = ratio;
  }
  out->base_time = target.base_time;
  out->base_compensation = target.base_compensation;
  out->base_successful = target.base_successful;
  out->nominal = target.nominal;
  out->size = target.size;
  out->time = target.time;
  out->compensation = target.compensation;
  out->completes_step = target.completes_step;
  out->use_extrapolation = target.use_extrapolation;
  out->extrapolation_ratio = target.extrapolation_ratio;
  return StatusCode::kOk;
}
StatusCode CutbackIncrement(const NonlinearStaticParameters* p,
                            const IncrementControls* c,
                            const IncrementClock* clock,
                            const IncrementTarget* target,
                            const NewtonIterationInfo* attempt,
                            int attempt_number, double* retry_nominal,
                            Diagnostic* d) {
  auto status = Validate(p, c, clock, d);
  if (status != StatusCode::kOk) return status;
  if (!retry_nominal || !attempt || !Bound(clock, target) ||
      attempt_number <= 0)
    return model_storage::Invalid("Cutback increment", d);
  if (attempt->logging_failed || attempt->outcome == NewtonOutcome::kFatal)
    return Primary(attempt, d);
  if (attempt->numerical_status != StatusCode::kEquilibriumFailure ||
      attempt->status != StatusCode::kEquilibriumFailure)
    return Primary(attempt, d);
  double factor = 0;
  switch (attempt->outcome) {
    case NewtonOutcome::kDivergence:
      factor = c->divergence_cutback;
      break;
    case NewtonOutcome::kSlowConvergence:
    case NewtonOutcome::kIterationLimit:
      factor = c->slow_cutback;
      break;
    case NewtonOutcome::kGeometry:
    case NewtonOutcome::kRotationBranch:
      factor = c->geometry_cutback;
      break;
    default:
      return model_storage::Invalid("Non-retryable Newton outcome", d);
  }
  if (attempt_number >= c->attempt_limit)
    return Exhausted("STATIC attempt limit exhausted", d);
  const double retry = target->size * factor;
  if (retry < p->minimum_increment || !Positive(retry))
    return Exhausted("STATIC cutback below minimum increment", d);
  *retry_nominal = retry;
  return StatusCode::kOk;
}
StatusCode AcceptIncrement(const NonlinearStaticParameters* p,
                           const IncrementControls* c,
                           const IncrementClock* clock,
                           const IncrementTarget* target,
                           const NewtonIterationInfo* attempt, bool had_cutback,
                           IncrementClock* out, Diagnostic* d) {
  auto status = Validate(p, c, clock, d);
  if (status != StatusCode::kOk) return status;
  if (!out || out == clock || !attempt || !Bound(clock, target) ||
      target->time > p->period || clock->successful >= p->increment_limit ||
      attempt->iteration < 0)
    return model_storage::Invalid("Accept increment", d);
  if (attempt->logging_failed || attempt->status != StatusCode::kOk ||
      attempt->numerical_status != StatusCode::kOk ||
      attempt->outcome != NewtonOutcome::kConverged)
    return Primary(attempt, d);
  IncrementClock next;
  next.time = target->time;
  next.compensation = target->compensation;
  next.nominal = target->nominal;
  next.last_increment = target->size;
  next.successful = clock->successful + 1;
  if (attempt->iteration > c->slow_success_iterations) {
    next.nominal *= c->success_scale;
  } else if (!had_cutback && attempt->iteration <= c->quick_iterations) {
    next.quick_successes = std::min(2, clock->quick_successes + 1);
    if (next.quick_successes == 2) {
      const double cap = p->has_maximum ? p->maximum_increment : p->period;
      const double growth = std::min(c->growth, c->maximum_growth);
      next.nominal = next.nominal >= cap / growth ? cap : next.nominal * growth;
    }
  }
  // A subminimum next nominal is not silently clamped. The completed output
  // remains legal; Propose reports exhaustion if another increment is needed.
  if (!Positive(next.nominal))
    return Arithmetic("Unrepresentable next STATIC increment", d);
  CopyClock(&next, out);
  return StatusCode::kOk;
}
}  // namespace fesa
