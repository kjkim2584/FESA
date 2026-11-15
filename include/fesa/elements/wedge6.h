/// @file
/// @brief Six-node wedge with WEDGE2 stiffness and WEDGE18 consistent mass.
#ifndef FESA_ELEMENTS_WEDGE6_H_
#define FESA_ELEMENTS_WEDGE6_H_
#include "fesa/elements/solid_element.h"
namespace fesa {
/// Borrows six Domain nodes and section/material; owns official result stores.
/// Node order is the lower triangle followed by the corresponding upper one.
class Wedge6 final : public SolidElement {
 public:
  /// Requires six distinct 3-DOF nodes and kWedge6; no stabilized fallback.
  /// Geometry guards follow geometry_contract.md. Caller owns success;
  /// failure nulls output and releases partial storage without definition
  /// changes.
  static StatusCode Create(int id, const Node* const* nodes, int node_count,
                           const SolidSection* section,
                           SolidFormulation formulation, Wedge6** out,
                           Diagnostic* d);
  ~Wedge6() override = default;
  Wedge6(const Wedge6&) = delete;
  Wedge6& operator=(const Wedge6&) = delete;
  Wedge6(Wedge6&&) = delete;
  Wedge6& operator=(Wedge6&&) = delete;

 private:
  Wedge6(int id, const SolidSection* section);
};
}  // namespace fesa
#endif  // FESA_ELEMENTS_WEDGE6_H_
