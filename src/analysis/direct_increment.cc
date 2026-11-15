#include "analysis/direct_increment.h"

#include <algorithm>
#include <cmath>

namespace fesa {
StatusCode InitializeDirectIncrementClock(double requested, double period,
                                          double minimum, double maximum,
                                          IncrementClock* out, Diagnostic* d) {
  if (!out)
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation, "DIRECT clock output",
                         d);
  double nominal = 0;
  const auto status =
      ResolveDirectIncrement(requested, period, minimum, maximum, &nominal, d);
  if (status != StatusCode::kOk) return status;
  out->time = out->compensation = out->last_increment = 0;
  out->nominal = nominal;
  out->successful = out->quick_successes = 0;
  return StatusCode::kOk;
}
StatusCode ProposeDirectIncrement(double period, int increment_limit,
                                  const IncrementClock* clock,
                                  IncrementTarget* out, Diagnostic* d) {
  if (!clock || !out || !std::isfinite(period) || period <= 0 ||
      increment_limit <= 0 || !std::isfinite(clock->time) || clock->time < 0 ||
      clock->time >= period || !std::isfinite(clock->nominal) ||
      clock->nominal <= 0 || clock->nominal > period ||
      !std::isfinite(clock->compensation) || clock->successful < 0)
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation, "DIRECT proposal inputs",
                         d);
  if (clock->successful >= increment_limit)
    return ReportFailure(StatusCode::kEquilibriumFailure,
                         DiagnosticStage::kComputation,
                         "DIRECT increment limit exhausted", d);
  const double remaining = period - clock->time;
  const bool clipped = remaining <= clock->nominal;
  // A fixed h=T/N can round-trip exactly while fl(N*h) is one ULP below T.
  // Recognize only the matching uniform clock, not an arbitrary success count.
  // The actual last size remains T-t, using the existing terminal contract.
  // @see docs/HHT/direct_endpoint_roundtrip.md.
  const double next_count = static_cast<double>(clock->successful) + 1;
  const bool closes_uniform_grid =
      !clipped && period / clock->nominal == next_count &&
      period / next_count == clock->nominal &&
      clock->time == (next_count - 1) * clock->nominal;
  const double adjusted =
      clipped ? remaining : clock->nominal - clock->compensation;
  const double proposed = clipped ? period : clock->time + adjusted;
  if (!std::isfinite(adjusted) || !std::isfinite(proposed))
    return ReportFailure(StatusCode::kLibraryFailure,
                         DiagnosticStage::kComputation,
                         "DIRECT compensated time overflow", d);
  const bool terminal = clipped || closes_uniform_grid || proposed >= period;
  const double time = terminal ? period : proposed;
  const double compensation =
      terminal ? 0 : (proposed - clock->time) - adjusted;
  if (time <= clock->time || !std::isfinite(compensation))
    return ReportFailure(StatusCode::kLibraryFailure,
                         DiagnosticStage::kComputation,
                         "DIRECT time cannot advance representably", d);
  out->base_time = clock->time;
  out->base_compensation = clock->compensation;
  out->base_successful = clock->successful;
  out->nominal = clock->nominal;
  out->size = terminal ? remaining : clock->nominal;
  out->time = time;
  out->compensation = compensation;
  out->completes_step = terminal;
  out->use_extrapolation = false;
  out->extrapolation_ratio = 0;
  return StatusCode::kOk;
}
StatusCode AcceptDirectIncrement(double period, int increment_limit,
                                 const IncrementClock* clock,
                                 const IncrementTarget* target,
                                 StatusCode numerical_status,
                                 IncrementClock* out, Diagnostic* d) {
  if (numerical_status != StatusCode::kOk) return numerical_status;
  if (!out || !target || out == clock)
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation,
                         "DIRECT acceptance output/target", d);
  IncrementTarget expected;
  auto status =
      ProposeDirectIncrement(period, increment_limit, clock, &expected, d);
  if (status != StatusCode::kOk) return status;
  if (target->base_time != expected.base_time ||
      target->base_compensation != expected.base_compensation ||
      target->base_successful != expected.base_successful ||
      target->nominal != expected.nominal || target->size != expected.size ||
      target->time != expected.time ||
      target->compensation != expected.compensation ||
      target->completes_step != expected.completes_step ||
      target->use_extrapolation || target->extrapolation_ratio != 0)
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation,
                         "DIRECT stale or modified proposal", d);
  out->time = target->time;
  out->compensation = target->compensation;
  out->nominal = clock->nominal;
  out->last_increment = target->size;
  out->successful = clock->successful + 1;
  out->quick_successes = 0;
  return StatusCode::kOk;
}
StatusCode ResolveDirectIncrement(double requested, double period,
                                  double minimum, double maximum,
                                  double* nominal, Diagnostic* d) {
  if (!nominal)
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation, "DIRECT nominal output",
                         d);
  if (!std::isfinite(requested) || requested <= 0 || !std::isfinite(period) ||
      period <= 0 || !std::isfinite(minimum) || minimum < 0 ||
      !std::isfinite(maximum) || maximum < 0)
    return ReportFailure(StatusCode::kInputError, DiagnosticStage::kValidation,
                         "DIRECT finite time settings", d);
  // Check original inputs before terminal clipping. An initial request above
  // T can legitimately finish the whole step with a duration below minimum.
  if (minimum > requested || (maximum > 0 && minimum > maximum))
    return ReportFailure(StatusCode::kInputError, DiagnosticStage::kValidation,
                         "DIRECT minimum exceeds initial or maximum", d);
  *nominal = std::min(requested, period);
  if (maximum > 0) *nominal = std::min(*nominal, maximum);
  return StatusCode::kOk;
}
}  // namespace fesa
