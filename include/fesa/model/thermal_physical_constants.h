/// @file
/// @brief Optional immutable heat-transfer physical constants.
#ifndef FESA_MODEL_THERMAL_PHYSICAL_CONSTANTS_H_
#define FESA_MODEL_THERMAL_PHYSICAL_CONSTANTS_H_
#include "fesa/base/diagnostic.h"
namespace fesa {
/// Scalar input definition cloned by Domain. Temperature origin is unchanged;
/// sigma uses model power/(area*temperature^4). No M1 numerical contribution.
struct ThermalPhysicalConstants {
  ThermalPhysicalConstants() = default;
  ThermalPhysicalConstants(const ThermalPhysicalConstants&) = delete;
  ThermalPhysicalConstants& operator=(const ThermalPhysicalConstants&) = delete;
  ThermalPhysicalConstants(ThermalPhysicalConstants&&) = delete;
  ThermalPhysicalConstants& operator=(ThermalPhysicalConstants&&) = delete;
  bool has_absolute_zero = false;
  double absolute_zero = 0;
  bool has_stefan_boltzmann = false;
  double stefan_boltzmann = 0;
  int source_line = 0;  ///< Positive input line, source is Domain::SourceFile.
};
/// Validates present finite constants (sigma positive), leaving input intact.
StatusCode ValidateThermalPhysicalConstants(
    const ThermalPhysicalConstants* values, Diagnostic* d);
/// Explicit scalar clone; caller owns success, failure nulls *out.
StatusCode CloneThermalPhysicalConstants(const ThermalPhysicalConstants* values,
                                         ThermalPhysicalConstants** out,
                                         Diagnostic* d);
}  // namespace fesa
#endif  // FESA_MODEL_THERMAL_PHYSICAL_CONSTANTS_H_
