/// @file
/// @brief Scalar heat-transfer step settings and Abaqus default resolution.
#ifndef FESA_MODEL_HEAT_TRANSFER_STEP_H_
#define FESA_MODEL_HEAT_TRANSFER_STEP_H_

#include "fesa/base/diagnostic.h"

namespace fesa {
/// Default step amplitude, distinct from an explicit tabular Amplitude.
enum class ThermalStepAmplitude { kStep, kRamp };
/// Caller-owned scalar definition; StepDefinition clones and owns its settings.
/// Durations use consistent model time, deltmx the input temperature unit.
struct HeatTransferStepParameters {
  HeatTransferStepParameters() = default;
  HeatTransferStepParameters(const HeatTransferStepParameters&) = delete;
  HeatTransferStepParameters& operator=(const HeatTransferStepParameters&) =
      delete;
  HeatTransferStepParameters(HeatTransferStepParameters&&) = delete;
  HeatTransferStepParameters& operator=(HeatTransferStepParameters&&) = delete;
  bool steady = false;
  double initial_increment = 1;
  double period = 1;
  bool has_minimum = false;
  double minimum = 0;  ///< Ignored when absent; resolved by EffectiveMinimum.
  bool has_maximum = false;
  double maximum = 0;  ///< No period substitution when absent.
  bool has_deltmx = false;
  double deltmx = 0;  ///< Positive automatic transient temperature limit.
  bool has_default_amplitude = false;  ///< Preserves an explicit STEP override.
  ThermalStepAmplitude default_amplitude = ThermalStepAmplitude::kStep;
  int max_increments = 100;
  /// Formulation 6.2: min(.8*dt0, explicit min or 1e-5*period).
  double EffectiveMinimum() const;
  /// Whether a maximum was explicitly provided; values are validated at
  /// creation.
  bool HasMaximum() const;
  /// Explicit override, else steady RAMP or transient STEP; no mutation.
  ThermalStepAmplitude DefaultAmplitude() const;
};
/// Validates finite positive time/limits and procedure-compatible DELTMX.
/// Explicit maximum below EffectiveMinimum is kInputError in every procedure
/// (Task11 T2/T8/T14/T15), even when fixed integration ignores its value.
/// Borrows settings; failure changes only optional diagnostic.
StatusCode ValidateHeatTransferStepParameters(
    const HeatTransferStepParameters* parameters, Diagnostic* d);
/// Explicit scalar clone; caller owns success, failure nulls *out.
/// Borrowed input and any previously created objects remain unchanged.
StatusCode CloneHeatTransferStepParameters(
    const HeatTransferStepParameters* parameters,
    HeatTransferStepParameters** out, Diagnostic* d);
}  // namespace fesa

#endif  // FESA_MODEL_HEAT_TRANSFER_STEP_H_
