/// @file
/// @brief Linear-static model ownership and immutable input definitions.
#ifndef FESA_MODEL_BEAM_SECTION_H_
#define FESA_MODEL_BEAM_SECTION_H_

#include "fesa/math/vector.h"
#include "fesa/model/elastic_material.h"

namespace fesa {

/// Immutable beam coefficients; borrows material and owns direction copy.
/// Units: area L^2; I11/I12/I22/J L^4; approximate n1 dimensionless.
/// Direction validity/normalization and stiffness belong to S1/S4.
class BeamSection {
 public:
  /// Caller owns *out; failure sets null and preserves all borrowed inputs.
  /// Requires positive ID/A/I11/I22/J, I12=0 and a 3-component direction.
  /// Material must outlive this section.
  static StatusCode Create(int id, const ElasticMaterial* material, double area,
                           double i11, double i12, double i22, double torsion_j,
                           const Vector* approximate_n1, BeamSection** out,
                           Diagnostic* d);
  ~BeamSection();
  BeamSection(const BeamSection&) = delete;
  BeamSection& operator=(const BeamSection&) = delete;
  BeamSection(BeamSection&&) = delete;
  BeamSection& operator=(BeamSection&&) = delete;
  /// Section ID.
  int Id() const;
  /// Borrowed material, valid while the owning Domain lives.
  const ElasticMaterial* GetMaterial() const;
  /// Borrowed explicitly copied direction, valid until section destruction.
  const Vector* ApproximateN1() const;
  /// Scalar section coefficients, in the units specified above.
  double Area() const;
  /// Second moment about n1.
  double I11() const;
  /// Product of inertia (currently zero only).
  double I12() const;
  /// Second moment about n2.
  double I22() const;
  /// Torsion constant.
  double TorsionJ() const;

 private:
  BeamSection(int id, const ElasticMaterial* material, double area, double i11,
              double i12, double i22, double torsion_j, Vector* direction);
  int id_;
  const ElasticMaterial* material_;  ///< Borrowed; Domain owns material.
  double area_;
  double i11_;
  double i12_;
  double i22_;
  double torsion_j_;
  Vector* approximate_n1_;  ///< Owned explicit copy.
};

}  // namespace fesa

#endif  // FESA_MODEL_BEAM_SECTION_H_
