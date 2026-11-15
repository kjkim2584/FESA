/// @file
/// @brief Affine four-node linear-elastic tetrahedron, F1/F2/F6-F10.
#ifndef FESA_ELEMENTS_TETRA4_H_
#define FESA_ELEMENTS_TETRA4_H_
#include "fesa/elements/solid_element.h"
namespace fesa {
/// C3D4, TET1 stiffness and TET4 consistent mass; no hourglass terms.
class Tetra4 final : public SolidElement {
 public:
  /// Requires four distinct three-component nodes, nonnull section and
  /// kTetra4Full. Borrows inputs through destruction, caller owns success.
  /// Geometry must pass geometry_contract.md; failure nulls and cleans output.
  static StatusCode Create(int id, const Node* const* nodes, int node_count,
                           const SolidSection* section,
                           SolidFormulation formulation, Tetra4** out,
                           Diagnostic* d);
  ~Tetra4() override = default;
  Tetra4(const Tetra4&) = delete;
  Tetra4& operator=(const Tetra4&) = delete;
  Tetra4(Tetra4&&) = delete;
  Tetra4& operator=(Tetra4&&) = delete;

 private:
  Tetra4(int id, const SolidSection* section);
};
}  // namespace fesa
#endif  // FESA_ELEMENTS_TETRA4_H_
