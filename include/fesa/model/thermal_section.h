/// @file
/// @brief Thermal solid section with borrowed immutable material.
#ifndef FESA_MODEL_THERMAL_SECTION_H_
#define FESA_MODEL_THERMAL_SECTION_H_

#include "fesa/model/thermal_material.h"

namespace fesa {
/// Domain-owned solid section. Its material must outlive this section.
class ThermalSolidSection final {
 public:
  /// Positive ID and non-null material required. Caller owns success;
  /// failure nulls *out. Domain registration checks material identity.
  static StatusCode Create(int id, const ThermalMaterial* material,
                           ThermalSolidSection** out, Diagnostic* d);
  ~ThermalSolidSection() = default;
  ThermalSolidSection(const ThermalSolidSection&) = delete;
  ThermalSolidSection& operator=(const ThermalSolidSection&) = delete;
  ThermalSolidSection(ThermalSolidSection&&) = delete;
  ThermalSolidSection& operator=(ThermalSolidSection&&) = delete;
  /// Positive original section identifier.
  int Id() const;
  /// Borrowed material, valid while the owning Domain lives.
  const ThermalMaterial* material() const;

 private:
  ThermalSolidSection(int id, const ThermalMaterial* material);
  int id_;
  const ThermalMaterial* material_;  ///< Borrowed; Domain owns the material.
};
/// Domain-owned planar, single-layer Simpson shell section, formulation 3.3.
/// Material is borrowed and must outlive this immutable definition.
class ThermalShellSection final {
 public:
  /// Positive finite thickness [L], odd ns in 1..19, positive ID and material
  /// required. Caller owns success; failure nulls *out. No input mutation.
  static StatusCode Create(int id, double thickness, int section_points,
                           const ThermalMaterial* material,
                           ThermalShellSection** out, Diagnostic* d);
  ~ThermalShellSection() = default;
  ThermalShellSection(const ThermalShellSection&) = delete;
  ThermalShellSection& operator=(const ThermalShellSection&) = delete;
  ThermalShellSection(ThermalShellSection&&) = delete;
  ThermalShellSection& operator=(ThermalShellSection&&) = delete;
  /// Original positive section ID.
  int Id() const;
  /// Uniform thickness [L].
  double thickness() const;
  /// Number of temperatures at each node, labels 11..10+ns.
  int section_points() const;
  /// Borrowed immutable material until Domain destruction.
  const ThermalMaterial* material() const;

 private:
  ThermalShellSection(int id, double thickness, int section_points,
                      const ThermalMaterial* material);
  int id_;
  double thickness_;
  int section_points_;
  const ThermalMaterial* material_;  ///< Borrowed from Domain.
};
}  // namespace fesa

#endif  // FESA_MODEL_THERMAL_SECTION_H_
