/// @file
/// @brief Immutable purpose-selected natural interpolation and quadrature data.
#ifndef FESA_ELEMENTS_SOLID_REFERENCE_DATA_H_
#define FESA_ELEMENTS_SOLID_REFERENCE_DATA_H_
#include "fesa/math/matrix.h"
namespace fesa {
/// Immutable formulation chosen by the input element TYPE.
enum class SolidFormulation {
  kTetra4Full,
  kTetra10Full,
  kTetra10Reduced,  ///< Retired identifier; factories return kUnsupported.
  kHexa8Full,
  kHexa8Reduced,
  kWedge6,
  kPyramid5Full,
  kPyramid5Reduced
};
/// Selects only the reference data needed by a requested operation.
enum class IntegrationPurpose {
  kStiffness, kMass, kGeometry,
  kTrialGeometry  ///< Q23/Q53 ordered current-shape samples; no integration.
};
/// Owns read-only natural coordinates g by 3, weights g, shape values g by n
/// and natural derivatives 3 by n per point. Analysis/caller owns this object.
class SolidReferenceData {
 public:
  /// Caller owns success; failure nulls output and releases partial
  /// allocations. Unsupported rules return kUnsupported without a substitute
  /// rule.
  static StatusCode Create(SolidFormulation formulation,
                           IntegrationPurpose purpose, SolidReferenceData** out,
                           Diagnostic* d);
  ~SolidReferenceData();
  SolidReferenceData(const SolidReferenceData&) = delete;
  SolidReferenceData& operator=(const SolidReferenceData&) = delete;
  SolidReferenceData(SolidReferenceData&&) = delete;
  SolidReferenceData& operator=(SolidReferenceData&&) = delete;
  /// Immutable input formulation identity.
  SolidFormulation formulation() const;
  /// Immutable requested table purpose, not a computed result.
  IntegrationPurpose purpose() const;
  /// Immutable interpolation node count.
  int node_count() const;
  /// Number of physical reference integration points.
  int point_count() const;
  /// Borrowed stable rule identifier, valid until destruction.
  const char* rule_id() const;
  /// Borrowed natural coordinates g by 3, dimensionless.
  const Matrix* natural_coordinates() const;
  /// Borrowed weights including reference-domain measure. TrialGeometry
  /// weights are zero: these samples must never be used for integration.
  const Vector* weights() const;
  /// Borrowed interpolation values g by n, dimensionless.
  const Matrix* shape_values() const;
  /// Borrowed natural gradients 3 by n; invalid point returns null.
  const Matrix* natural_gradients(int point) const;

 private:
  SolidReferenceData() = default;
  /// Builds Q23/Q53 samples after immutable formulation/node count selection.
  /// Owns partial allocations on failure; factory deletes the entire object.
  StatusCode InitializeTrialGeometry(Diagnostic* d);
  SolidFormulation formulation_ = SolidFormulation::kTetra4Full;
  IntegrationPurpose purpose_ = IntegrationPurpose::kStiffness;
  int node_count_ = 0;
  int point_count_ = 0;
  const char* rule_id_ = nullptr;  ///< Borrowed static identifier.
  Matrix* coordinates_ = nullptr;  ///< Owned g by 3 natural coordinates.
  Vector* weights_ = nullptr;      ///< Owned reference weights.
  Matrix* shapes_ = nullptr;       ///< Owned g by n interpolation coefficients.
  Matrix** gradients_ = nullptr;   ///< Owned array and gradient matrices.
};
}  // namespace fesa
#endif  // FESA_ELEMENTS_SOLID_REFERENCE_DATA_H_
