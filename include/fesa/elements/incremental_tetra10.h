/// @file
/// @brief Isolated C3D10 factory for incremental elastic static response.
#ifndef FESA_ELEMENTS_INCREMENTAL_TETRA10_H_
#define FESA_ELEMENTS_INCREMENTAL_TETRA10_H_
#include "fesa/elements/incremental_solid_element.h"
namespace fesa {
/// Ten-node TET4 solid. Existing Tetra10/StVK classes are not base classes.
class IncrementalTetra10 final : public IncrementalSolidElement {
 public:
  /// Borrows nodes/section through destruction. Caller owns success; Domain
  /// registration transfers ownership only on success. Failure cleans partial
  /// resources and nulls out. Requires ten distinct3DOF nodes and kTetra10Full.
  static StatusCode Create(int id, const Node* const* nodes, int node_count,
      const SolidSection* section, SolidFormulation formulation,
      IncrementalTetra10** out, Diagnostic* d);
  ~IncrementalTetra10() override = default;
  IncrementalTetra10(const IncrementalTetra10&) = delete;
  IncrementalTetra10& operator=(const IncrementalTetra10&) = delete;
  IncrementalTetra10(IncrementalTetra10&&) = delete;
  IncrementalTetra10& operator=(IncrementalTetra10&&) = delete;
 private:
  /// Factory-only construction; borrows validated section.
  IncrementalTetra10(int id, const SolidSection* section);
};
}  // namespace fesa
#endif  // FESA_ELEMENTS_INCREMENTAL_TETRA10_H_
