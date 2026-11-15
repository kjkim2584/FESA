/// @file
/// @brief Thermal fixed/automatic increment proposals without physical state.
#ifndef FESA_ANALYSIS_THERMAL_INCREMENT_CONTROLLER_H_
#define FESA_ANALYSIS_THERMAL_INCREMENT_CONTROLLER_H_
#include "fesa/analysis/heat_flux_residual.h"
#include "fesa/analysis/thermal_history.h"
#include "fesa/model/heat_transfer_step.h"
namespace fesa {
struct ThermalIncrementPolicy;
/// Typed cause; never infer equilibrium failure from a temperature change.
enum class IncrementRejection {
  kTimeAccuracy,
  kDivergence,
  kSlowConvergence,
  kIterationCap
};
/// Terminal retry reason, distinct from fatal input/arithmetic errors.
enum class IncrementFailure { kNone, kFixed, kMinimum, kAttemptLimit };
/// Caller-owned scalar output; dt is actual next attempted size, time units.
struct IncrementRetry {
  double dt = 0;
  IncrementFailure failure = IncrementFailure::kNone;
};
/// Owns immutable settings copy and scalar clock/policy bookkeeping. Never
/// owns or retains model/state/rate pointers. Formulation6 defines the rules.
class ThermalIncrementController final {
 public:
  /// Clones borrowed settings; caller owns success. Failure nulls *out and
  /// frees partial storage. M1 coefficients are fixed, not input-configurable.
  static StatusCode Create(const HeatTransferStepParameters* parameters,
                           ThermalIncrementController** out, Diagnostic* d);
  ~ThermalIncrementController();
  ThermalIncrementController(const ThermalIncrementController&) = delete;
  ThermalIncrementController& operator=(const ThermalIncrementController&) =
      delete;
  ThermalIncrementController(ThermalIncrementController&&) = delete;
  ThermalIncrementController& operator=(ThermalIncrementController&&) = delete;
  /// Requires the exact current accepted step_time. Returns remaining-capped
  /// proposal;0 only at completion, NaN for stale/nonfinite time. No mutation.
  double NextIncrement(double step_time) const;
  /// Computes the actual dt and end together, without mutation. Formulation6.5
  /// closes only the accumulated gamma_(accepted+2)*period roundoff remainder.
  /// Completion writes (0,period). Stale input is invalid; a broken clock
  /// invariant is kLibraryFailure. Any failure preserves both scalar outputs.
  StatusCode PlanIncrement(double step_time, double* dt, double* end,
                           Diagnostic* d) const;
  double StepTime() const;  ///< Accepted time, initially0, no physical state.
  int AcceptedCount() const;
  int AttemptNumber() const;  ///< IA includes the first attempt, initially1.
  /// Computes the default LINEAR predictor multiplier for the next attempt:
  /// actual dt / previous accepted dt, provided nominal/previous >= DE=.1.
  /// previous_dt=0 means no same-step history and writes0 (carry). Borrows the
  /// scalar output only; no mutation. Invalid/completed/numeric failure
  /// preserves *factor and controller state.
  StatusCode ExtrapolationFactor(double previous_dt, double* factor,
                                 Diagnostic* d) const;
  /// Call after output/State commit. Automatic transient requires committed
  /// rate ring including this increment (not staged). Other modes allow null.
  /// Borrows history for this call, never Preview/Commit/Discard. Computes in
  /// locals then promotes clock/proposal; every failure preserves controller.
  /// A failure here is fatal after physical commit, never a cutback request.
  /// This overload retains M1 linear-driver semantics; physical analysis uses
  /// the explicit-branch overload below.
  StatusCode OnAccepted(double dt, double max_free_change, int iterations,
                        const RateHistory* committed, Diagnostic* d) {
    return OnAccepted(dt, max_free_change, iterations,
                      ThermalResidualBranch::kLinear, committed, d);
  }
  /// Same lifetime/failure contract; branch is the current accepted residual
  /// branch, not the IT window aggregate. General convergence or iterations>=2
  /// bypasses DL. Transient growth uses factor2 for iterations<=1, otherwise
  /// 1.5; both modes require each of the last two accepted counts <= IG.
  /// A successful accuracy retry counts first of IT3 without clearing rates.
  StatusCode OnAccepted(double dt, double max_free_change, int iterations,
                        ThermalResidualBranch branch,
                        const RateHistory* committed, Diagnostic* d);
  /// Actual attempted dt must match NextIncrement(StepTime()). Accuracy cause
  /// requires S>DELTMX; others use public .25/.5 factors, with minimum-floor
  /// behaviour from formulation6.4/8.7. Retry kOk writes dt and records any
  /// accuracy rejection for this increment's subsequent accepted counter;
  /// terminal kEquilibriumFailure writes reason. Invalid/numeric failure
  /// preserves out; all failure preserves controller and physical/rate data.
  StatusCode OnRejected(IncrementRejection cause, double attempted_dt,
                        double accuracy_measure, IncrementRetry* out,
                        Diagnostic* d);

 private:
  friend StatusCode CreateThermalIncrementControllerWithPolicy(
      const HeatTransferStepParameters*, const ThermalIncrementPolicy*,
      ThermalIncrementController**, Diagnostic*);
  ThermalIncrementController() = default;
  bool IsAutomatic() const;
  HeatTransferStepParameters* settings_ = nullptr;  ///< Owned immutable clone.
  double nominal_ = 0, time_ = 0, dl_ = .95;
  int accepted_ = 0, attempt_ = 1, previous_iterations_ = 0;
  int accuracy_accepted_ =
      3;  ///< Saturated IT3 count; eligible before cutback.
  bool accuracy_rejected_ = false;  ///< Pending since last accepted increment.
  bool floor_attempt_ = false;
};
}  // namespace fesa
#endif  // FESA_ANALYSIS_THERMAL_INCREMENT_CONTROLLER_H_
