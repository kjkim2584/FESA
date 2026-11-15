/// @file
/// @brief Immutable relative TABULAR load amplitude in model time units.
#ifndef FESA_MODEL_AMPLITUDE_H_
#define FESA_MODEL_AMPLITUDE_H_

#include "fesa/math/vector.h"

namespace fesa {
/// Caller selects the matching physical clock; no time shift is implicit.
enum class AmplitudeTimeBasis { kStepTime, kTotalTime };

/// Borrowed scalar context for an active dynamic load's relative amplitude.
/// The load resolver supplies the last definition step (not the table's scope)
/// and its period. It retains them across unchanged steps and replaces them
/// upon redefinition. Removed loads must not be evaluated. Times use model T.
struct AmplitudeEvaluationContext {
  int definition_step = 0;  ///< One-based step of the active load definition.
  int current_step = 0;     ///< One-based evaluation step, >= definition_step.
  double definition_period = 0;  ///< Finite positive defining step duration.
  double step_time = 0;     ///< Finite nonnegative time in the current step.
  double total_time = 0;    ///< Includes prior static periods; >= step_time.
};

/// Owns its name and table. No mutable state or borrowed inputs are retained.
/// The creating caller deletes it; future Domain registration transfers
/// ownership.
/// @see docs/HHT/io_contract.md for equations and supported semantics.
class Amplitude final {
 public:
  /// Explicitly copies borrowed times (T) and relative values (dimensionless).
  /// Requires positive ID, nonempty name, valid basis, equal positive sizes,
  /// finite entries and strictly increasing times. One point is constant.
  /// Success transfers ownership to caller. Failure nulls *out and releases
  /// partial resources. out is required and must not hold an existing object.
  static StatusCode Create(int id, const char* name, AmplitudeTimeBasis basis,
                           const Vector* times, const Vector* values,
                           Amplitude** out, Diagnostic* d);
  ~Amplitude();
  Amplitude(const Amplitude&) = delete;
  Amplitude& operator=(const Amplitude&) = delete;
  Amplitude(Amplitude&&) = delete;
  Amplitude& operator=(Amplitude&&) = delete;
  /// Original positive model ID.
  int Id() const;
  /// Borrowed original spelling, valid until destruction.
  const char* Name() const;
  /// Clock used by the caller when selecting evaluation time.
  AmplitudeTimeBasis TimeBasis() const;
  /// Number of immutable tabulated points.
  int PointCount() const;
  /// Linearly interpolates a finite time in the selected clock, holding the
  /// first/last value outside the table. No smoothing or derivatives.
  /// Mutates only *value on success. Invalid input preserves *value and this.
  /// No arguments are retained; value required, diagnostic optional.
  StatusCode Evaluate(double time, double* value, Diagnostic* d) const;
  /// Exact minimum/maximum of the held, piecewise-linear table on the closed
  /// interval [start,end] in its own clock. Extrema occur at endpoints or
  /// enclosed knots. Requires finite ordered times and distinct nonnull
  /// outputs. Failure preserves outputs; no allocation or retained pointers.
  /// Used to validate a physical ambient over a whole step before output.
  StatusCode EvaluateRange(double start, double end, double* minimum,
                           double* maximum, Diagnostic* d) const;
  /// Returns the relative coefficient for the supplied active load context.
  /// STEP TIME uses step_time in the definition step (<= definition_period),
  /// then holds its defining end value; TOTAL TIME always uses total_time.
  /// Inputs are borrowed for this call; no allocation, retention or mutation.
  /// Success overwrites only *value. Invalid context/null output preserves it;
  /// output aliasing any context scalar is rejected before writing.
  /// diagnostic source fields are preserved. Does not resolve OP/targets,
  /// multiply the reference load, or apply static ramp semantics.
  /// @see docs/HHT/io_contract.md for continuation and clock contracts.
  StatusCode EvaluateInStep(const AmplitudeEvaluationContext* context,
                            double* value, Diagnostic* d) const;
  /// Right-hand TABULAR slope in1/T. Uses the following segment at a knot;
  /// before the first point/at or after the last, the held-value slope is0.
  /// No smoothing or distributional derivative. Mutates only *rate on success.
  /// Invalid inputs preserve it; overflow/nonzero underflow is kLibraryFailure.
  /// No inputs retained or allocations; rate required, diagnostic optional.
  StatusCode EvaluateRightDerivative(double time, double* rate,
                                     Diagnostic* d) const;
  /// Same context/alias contract as EvaluateInStep, returning its right slope.
  /// STEP TIME returns0 at the defining step end and in propagated steps;
  /// TOTAL TIME continues evaluating its clock. Does not alter motion/loads.
  StatusCode EvaluateRightDerivativeInStep(
      const AmplitudeEvaluationContext* context, double* rate,
      Diagnostic* d) const;

 private:
  Amplitude() = default;
  int id_ = 0;
  char* name_ = nullptr;  ///< Owned original name.
  AmplitudeTimeBasis basis_ = AmplitudeTimeBasis::kStepTime;
  Vector* times_ = nullptr;   ///< Owned strictly increasing model times.
  Vector* values_ = nullptr;  ///< Owned relative scalar coefficients.
};
}  // namespace fesa
#endif  // FESA_MODEL_AMPLITUDE_H_
