/// @file
/// @brief Convergence-history AGGRESSIVE proposals, without physical state.
#ifndef FESA_ANALYSIS_DYNAMIC_INCREMENT_CONTROLLER_H_
#define FESA_ANALYSIS_DYNAMIC_INCREMENT_CONTROLLER_H_

#include "fesa/analysis/dynamic_step_settings.h"
#include "fesa/analysis/increment_clock.h"
#include "fesa/analysis/increment_controls.h"
#include "fesa/analysis/newton_report.h"

namespace fesa {
/// Borrows immutable automatic settings for its lifetime. No allocations,
/// model/state/writer ownership, physical mutation or implicit promotion.
/// All method inputs are call-borrowed; failure preserves outputs and inputs.
/// Outputs must use storage distinct from all borrowed input metadata.
/// @see docs/HHT/automatic_increment_control.md for evidence and equations.
class DynamicIncrementController final {
 public:
  /// settings must outlive this controller; invalid/nonautomatic is rejected
  /// by methods, so construction itself has no fallible side effects.
  explicit DynamicIncrementController(const DynamicStepSettings* settings);
  DynamicIncrementController(const DynamicIncrementController&) = delete;
  DynamicIncrementController& operator=(const DynamicIncrementController&) =
      delete;
  DynamicIncrementController(DynamicIncrementController&&) = delete;
  DynamicIncrementController& operator=(DynamicIncrementController&&) = delete;
  /// Creates time0 clock/q0; initial frame consumes no increment count.
  StatusCode Initialize(IncrementClock* out, Diagnostic* d) const;
  /// Bounded nominal then compensated DIRECT time proposal; endpoint-only
  /// remainder may be subminimum. No static extrapolation/stability stretch.
  StatusCode Propose(const IncrementControls* controls,
                     const IncrementClock* clock, double nominal,
                     IncrementTarget* out, Diagnostic* d) const;
  /// Typed failure shrink of actual attempted h; IA includes first attempt.
  /// Logging/fatal status remains primary, without retry. No quick reset.
  StatusCode Cutback(const IncrementControls* controls,
                     const IncrementClock* clock,
                     const IncrementTarget* target,
                     const NewtonIterationInfo* attempt, int attempt_number,
                     double* retry_nominal, Diagnostic* d) const;
  /// Stages next clock after numerical success; caller promotes only together
  /// with physical/history commit after flush. Null attempt means successful
  /// direct linear/no-equation path, quick-class without invented iterations.
  /// Non-null must describe actual converged Newton; output distinct from
  /// clock. Slow wins overlap but still contributes to the quick streak.
  StatusCode Accept(const IncrementControls* controls,
                    const IncrementClock* clock,
                    const IncrementTarget* target,
                    const NewtonIterationInfo* attempt, IncrementClock* out,
                    Diagnostic* d) const;

 private:
  /// Validates borrowed settings/controls and analysis-owned scalar clock.
  StatusCode Validate(const IncrementControls* controls,
                      const IncrementClock* clock, Diagnostic* d) const;
  /// Recomputes proposal to reject stale/modified candidate metadata.
  StatusCode ValidateTarget(const IncrementControls* controls,
                            const IncrementClock* clock,
                            const IncrementTarget* target,
                            Diagnostic* d) const;
  const DynamicStepSettings* settings_ = nullptr;  ///< Borrowed immutable.
};
}  // namespace fesa
#endif  // FESA_ANALYSIS_DYNAMIC_INCREMENT_CONTROLLER_H_
