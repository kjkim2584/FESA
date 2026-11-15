/// @file
/// @brief Scalar increment proposals/precommit clock updates; Analysis owns
/// loop.
#ifndef FESA_ANALYSIS_INCREMENT_CONTROL_H_
#define FESA_ANALYSIS_INCREMENT_CONTROL_H_
#include "fesa/analysis/increment_clock.h"
#include "fesa/analysis/increment_controls.h"
#include "fesa/analysis/newton_report.h"
#include "fesa/model/nonlinear_static_parameters.h"
namespace fesa {
/// Sets initial clock after parameter validation. All objects borrowed,
/// output creation-only metadata; preflight preserves it. No allocation.
StatusCode InitializeIncrementClock(const NonlinearStaticParameters* parameters,
                                    IncrementClock* out, Diagnostic* d);
/// Proposes using provided nominal (clock.nominal first, local cutback retry).
/// Kahan sum compensation is only promoted on accepted output/physical commit.
/// Exact T used at endpoint; no epsilon skips unfinished time. No-progress,
/// unrepresentable predictor arithmetic fatal; INC/min exhaustion EqFailure.
/// Outputs distinct, preflight/any failure preserves out; no input mutation.
StatusCode ProposeIncrement(const NonlinearStaticParameters* parameters,
                            const IncrementControls* controls,
                            const IncrementClock* clock, double nominal,
                            IncrementTarget* out, Diagnostic* d);
/// Failure-count IA includes first attempt. Only Df/DC/DH outcomes retry;
/// logging_failed/fatal stop with existing primary status, no new nominal.
/// EqFailure on exhausted attempts/minimum; caller aborts on any non-OK.
/// Output scalar preserved on failure; no committed metadata changes.
StatusCode CutbackIncrement(const NonlinearStaticParameters* parameters,
                            const IncrementControls* controls,
                            const IncrementClock* clock,
                            const IncrementTarget* target,
                            const NewtonIterationInfo* attempt,
                            int attempt_number, double* retry_nominal,
                            Diagnostic* d);
/// Precomputes candidate clock only for converged status/no logging failure.
/// No commit here: Analysis swaps clock pointers in the same no-fail gate as
/// physical commit AFTER writer flush. Out distinct from input clock; any
/// failure preserves out. Cutback breaks quick streak; shrink wins overlap.
StatusCode AcceptIncrement(const NonlinearStaticParameters* parameters,
                           const IncrementControls* controls,
                           const IncrementClock* clock,
                           const IncrementTarget* target,
                           const NewtonIterationInfo* attempt, bool had_cutback,
                           IncrementClock* out, Diagnostic* d);
}  // namespace fesa
#endif  // FESA_ANALYSIS_INCREMENT_CONTROL_H_
