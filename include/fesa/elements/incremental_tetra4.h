/// @file
/// @brief Isolated C3D4 factory for incremental elastic static response.
#ifndef FESA_ELEMENTS_INCREMENTAL_TETRA4_H_
#define FESA_ELEMENTS_INCREMENTAL_TETRA4_H_
#include "fesa/elements/incremental_solid_element.h"
namespace fesa {
/// Four-node TET1 solid. Existing Tetra4/StVK classes are not base classes.
class IncrementalTetra4 final : public IncrementalSolidElement {
 public:
  /// Borrows nodes/section through destruction. Caller owns success; Domain
  /// registration transfers ownership only on success. Failure cleans partial
  /// resources and nulls out. Requires distinct3DOF nodes and kTetra4Full.
  static StatusCode Create(int id, const Node* const* nodes, int node_count,
      const SolidSection* section, SolidFormulation formulation,
      IncrementalTetra4** out, Diagnostic* d);
  ~IncrementalTetra4() override = default;
  IncrementalTetra4(const IncrementalTetra4&) = delete;
  IncrementalTetra4& operator=(const IncrementalTetra4&) = delete;
  IncrementalTetra4(IncrementalTetra4&&) = delete;
  IncrementalTetra4& operator=(IncrementalTetra4&&) = delete;
 private:
  /// Factory-only construction; borrows validated section.
  IncrementalTetra4(int id, const SolidSection* section);
};
}  // namespace fesa
#endif  // FESA_ELEMENTS_INCREMENTAL_TETRA4_H_
