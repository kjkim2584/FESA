/// @file
/// @brief Isolated C3D6 factory for incremental elastic static response.
#ifndef FESA_ELEMENTS_INCREMENTAL_WEDGE6_H_
#define FESA_ELEMENTS_INCREMENTAL_WEDGE6_H_
#include "fesa/elements/incremental_solid_element.h"
namespace fesa {
/// Six-node WEDGE2 solid. Existing Wedge6/StVK classes are not base classes.
class IncrementalWedge6 final : public IncrementalSolidElement {
 public:
  /// Borrows nodes/section through destruction. Caller owns success; Domain
  /// registration transfers ownership only on success. Failure cleans partial
  /// resources and nulls out. Requires distinct3DOF nodes and kWedge6.
  static StatusCode Create(int id, const Node* const* nodes, int node_count,
      const SolidSection* section, SolidFormulation formulation,
      IncrementalWedge6** out, Diagnostic* d);
  ~IncrementalWedge6() override = default;
  IncrementalWedge6(const IncrementalWedge6&) = delete;
  IncrementalWedge6& operator=(const IncrementalWedge6&) = delete;
  IncrementalWedge6(IncrementalWedge6&&) = delete;
  IncrementalWedge6& operator=(IncrementalWedge6&&) = delete;
 private:
  /// Factory-only construction; borrows validated section.
  IncrementalWedge6(int id, const SolidSection* section);
};
}  // namespace fesa
#endif  // FESA_ELEMENTS_INCREMENTAL_WEDGE6_H_
