/// @file
/// @brief Call-owned native-compatible C3D8R hourglass operator (F14).
#ifndef FESA_ELEMENTS_HEXA_HOURGLASS_H_
#define FESA_ELEMENTS_HEXA_HOURGLASS_H_
#include "fesa/math/matrix.h"
namespace fesa {
class ElasticMaterial;
/// Owns H3/C and purpose-specific temporary products, never retained by
/// Element. All global vectors are node-major XYZ. Coordinates are
/// (X-Xanchor)/h; gradient is the physical volume average3x8. The empirical
/// .2357 multiplier is user-adopted and documented in
/// SolidElement/formulation.md, not derived from the public .005G default
/// factor.
class HexaHourglass {
 public:
  /// Borrows inputs only during creation. Caller owns success; failure nulls
  /// output and releases partial resources. Builds no full24x24 HG matrix.
  /// stiffness selects CH3 storage; otherwise only q/Q/resisting-force storage.
  static StatusCode Create(const Matrix* coordinates, double scale,
                           const Matrix* gradient, double volume,
                           const ElasticMaterial* material, bool stiffness,
                           HexaHourglass** out, Diagnostic* d);
  ~HexaHourglass();
  HexaHourglass(const HexaHourglass&) = delete;
  HexaHourglass& operator=(const HexaHourglass&) = delete;
  HexaHourglass(HexaHourglass&&) = delete;
  HexaHourglass& operator=(HexaHourglass&&) = delete;
  /// Adds H3^T C H3 to borrowed global24x24 K. No physical K/D evaluation.
  /// Invalid arguments preserve out; later backend failure invalidates out.
  StatusCode AddStiffness(Matrix* out, Diagnostic* d);
  /// Adds H3^T C H3*u through q/Q; no K is constructed. u/out24 disjoint.
  /// Invalid arguments preserve out; later backend failure invalidates out.
  StatusCode AddForce(const Vector* u, Vector* out, Diagnostic* d);

 private:
  HexaHourglass() = default;
  StatusCode Initialize(const Matrix* coordinates, double scale,
                        const Matrix* gradient, double volume,
                        const ElasticMaterial* material, bool stiffness,
                        Diagnostic* d);
  Matrix *h_ = nullptr, *c_ = nullptr,
         *weighted_h_ = nullptr;  ///< Owned H3/C/CH3.
  Vector *q_ = nullptr, *resistance_ = nullptr,
         *force_ = nullptr;  ///< Owned temporary.
};
}  // namespace fesa
#endif  // FESA_ELEMENTS_HEXA_HOURGLASS_H_
