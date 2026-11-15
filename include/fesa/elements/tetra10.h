/// @file
/// @brief Straight quadratic tetrahedron with full stiffness integration.
#ifndef FESA_ELEMENTS_TETRA10_H_
#define FESA_ELEMENTS_TETRA10_H_
#include "fesa/elements/solid_element.h"
namespace fesa {
/// Ten translations-only nodes; midsides follow Abaqus edges 12,23,31,14,24,34.
/// Borrows Domain definitions; owns official results through SolidElement.
class Tetra10 final : public SolidElement {
 public:
  /// Requires ten distinct 3-DOF nodes, a section and kTetra10Full.
  /// Retired kTetra10Reduced returns kUnsupported; use C3D10.
  /// Straight midsides satisfy geometry_contract.md; no coordinate adjustment.
  /// Caller owns success; failure nulls output and cleans partial storage.
  static StatusCode Create(int id, const Node* const* nodes, int node_count,
                           const SolidSection* section,
                           SolidFormulation formulation, Tetra10** out,
                           Diagnostic* d);
  ~Tetra10() override = default;
  Tetra10(const Tetra10&) = delete;
  Tetra10& operator=(const Tetra10&) = delete;
  Tetra10(Tetra10&&) = delete;
  Tetra10& operator=(Tetra10&&) = delete;

 private:
  Tetra10(int id, const SolidSection* section, SolidFormulation formulation);
};
}  // namespace fesa
#endif  // FESA_ELEMENTS_TETRA10_H_
