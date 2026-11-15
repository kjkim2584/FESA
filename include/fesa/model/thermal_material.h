/// @file
/// @brief Owned thermal property tables with optional transient capacity.
#ifndef FESA_MODEL_THERMAL_MATERIAL_H_
#define FESA_MODEL_THERMAL_MATERIAL_H_
#include "fesa/model/material.h"
#include "fesa/model/temperature_table.h"
namespace fesa {
/// Domain-owned immutable thermal material. Clones the borrowed input tables.
/// Units follow formulation.md: k [E/T/L/theta], rho [M/L^3], c [E/M/theta].
class ThermalMaterial final : public Material {
 public:
  /// Caller owns *out; failure nulls it and releases partial resources.
  /// Requires positive table values, valid ID/name, and positive finite rho
  /// when present. Absent density is ignored. Missing rho/c means steady only.
  static StatusCode Create(int id, const char* name,
                           const TemperatureTable* conductivity,
                           bool has_density, double density,
                           const TemperatureTable* specific_heat,
                           ThermalMaterial** out, Diagnostic* d);
  ~ThermalMaterial() override;
  ThermalMaterial(const ThermalMaterial&) = delete;
  ThermalMaterial& operator=(const ThermalMaterial&) = delete;
  ThermalMaterial(ThermalMaterial&&) = delete;
  ThermalMaterial& operator=(ThermalMaterial&&) = delete;
  /// Borrowed this pointer, valid until material destruction.
  const ThermalMaterial* AsThermal() const override;
  /// Evaluates k; required output is preserved on failure. Table unchanged.
  StatusCode Conductivity(double theta, double* k, Diagnostic* d) const;
  /// Evaluates rho*c; missing rho/c is kInputError. Output preserved on
  /// failure.
  StatusCode Capacity(double theta, double* rho_c, Diagnostic* d) const;
  /// Returns rho*integral(c dtheta), signed; same failure/ownership contract.
  StatusCode EnergyDifference(double from, double to, double* rho_du,
                              Diagnostic* d) const;
  /// Capacity availability, independent of the current analysis procedure.
  bool HasCapacity() const;
  /// True if any material property table has more than one row.
  bool IsTemperatureDependent() const;
  /// True if k has a temperature table, independent of the specific-heat table.
  /// Elements use this to preserve their constant-k quadrature/evaluation path.
  bool IsConductivityTemperatureDependent() const;

 private:
  ThermalMaterial(int id, char* name, TemperatureTable* conductivity,
                  bool has_density, double density, TemperatureTable* heat);
  TemperatureTable* conductivity_;   ///< Owned clone.
  TemperatureTable* specific_heat_;  ///< Owned optional clone.
  bool has_density_;
  double density_;
};
}  // namespace fesa
#endif  // FESA_MODEL_THERMAL_MATERIAL_H_
