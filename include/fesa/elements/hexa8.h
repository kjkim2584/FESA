/// @file
/// @brief Trilinear eight-node hexahedron with full or stabilized one-point K.
#ifndef FESA_ELEMENTS_HEXA8_H_
#define FESA_ELEMENTS_HEXA8_H_
#include "fesa/elements/solid_element.h"
namespace fesa {
/// Eight translations-only nodes in formulation F3/Abaqus sign order.
/// Borrows Domain definitions; owns official results through SolidElement.
class Hexa8 final : public SolidElement {
 public:
  /// Requires eight distinct 3-DOF nodes, a section and Hexa8 full/reduced.
  /// Geometry guards follow geometry_contract.md; reduced HG follows F14.
  /// Caller owns success; failure nulls output and cleans partial storage.
  static StatusCode Create(int id, const Node* const* nodes, int node_count,
                           const SolidSection* section,
                           SolidFormulation formulation, Hexa8** out,
                           Diagnostic* d);
  ~Hexa8() override = default;
  Hexa8(const Hexa8&) = delete;
  Hexa8& operator=(const Hexa8&) = delete;
  Hexa8(Hexa8&&) = delete;
  Hexa8& operator=(Hexa8&&) = delete;

 private:
  Hexa8(int id, const SolidSection* section, SolidFormulation formulation);
};
}  // namespace fesa
#endif  // FESA_ELEMENTS_HEXA8_H_
