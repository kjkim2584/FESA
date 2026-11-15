/// @file
/// @brief Immutable scalar thermal-property interpolation and integration.
#ifndef FESA_MODEL_TEMPERATURE_TABLE_H_
#define FESA_MODEL_TEMPERATURE_TABLE_H_
#include "fesa/math/vector.h"
namespace fesa {
/// Owns explicit copies of values and strictly increasing model temperatures.
/// One row is constant. Interpolation is linear; extrapolation is constant.
class TemperatureTable final {
 public:
  /// Borrows inputs during call; caller owns success, failure nulls *out.
  static StatusCode Create(const Vector* values, const Vector* temperatures,
                           TemperatureTable** out, Diagnostic* d);
  ~TemperatureTable();
  TemperatureTable(const TemperatureTable&) = delete;
  TemperatureTable& operator=(const TemperatureTable&) = delete;
  TemperatureTable(TemperatureTable&&) = delete;
  TemperatureTable& operator=(TemperatureTable&&) = delete;
  /// Independent caller-owned deep copy; same factory failure contract.
  StatusCode Clone(TemperatureTable** out, Diagnostic* d) const;
  /// Finite theta; outputs required/distinct, preserved on failure.
  StatusCode Evaluate(double theta, double* value, double* slope,
                      Diagnostic* d) const;
  /// Exact signed integral in property*temperature units; preserves *out on
  /// failure. Borrows d and never changes the table. Bounds must be finite.
  StatusCode Integrate(double from, double to, double* out,
                       Diagnostic* d) const;
  /// True for one row; no assumptions about equality of multiple row values.
  bool IsConstant() const;
  /// Tests all tabulated values for positivity without changing the table.
  bool HasPositiveValues() const;
  /// True when every tabulated value is nonnegative, including zero.
  bool HasNonnegativeValues() const;
  /// True only when every tabulated value equals zero (no exchange anchor).
  bool IsZero() const;

 private:
  TemperatureTable(Vector* values, Vector* temperatures);
  Vector* values_;        ///< Owned fixed-size property data.
  Vector* temperatures_;  ///< Owned model-temperature coordinates.
};
}  // namespace fesa
#endif  // FESA_MODEL_TEMPERATURE_TABLE_H_
