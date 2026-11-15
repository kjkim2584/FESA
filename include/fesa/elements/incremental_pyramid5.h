/// @file
/// @brief Isolated C3D5 factory for incremental elastic static response.
#ifndef FESA_ELEMENTS_INCREMENTAL_PYRAMID5_H_
#define FESA_ELEMENTS_INCREMENTAL_PYRAMID5_H_
#include "fesa/elements/incremental_solid_element.h"
namespace fesa {
/// Five-node PYR5_FELIPPA solid. Existing Pyramid5/StVK classes are not base classes.
class IncrementalPyramid5 final : public IncrementalSolidElement {
 public:
  /// Borrows nodes/section through destruction. Caller owns success; Domain
  /// registration transfers ownership only on success. Failure cleans partial
  /// resources and nulls out. Requires distinct3DOF nodes and kPyramid5Full.
  static StatusCode Create(int id, const Node* const* nodes, int node_count,
      const SolidSection* section, SolidFormulation formulation,
      IncrementalPyramid5** out, Diagnostic* d);
  ~IncrementalPyramid5() override = default;
  IncrementalPyramid5(const IncrementalPyramid5&) = delete;
  IncrementalPyramid5& operator=(const IncrementalPyramid5&) = delete;
  IncrementalPyramid5(IncrementalPyramid5&&) = delete;
  IncrementalPyramid5& operator=(IncrementalPyramid5&&) = delete;
 private:
  /// Factory-only construction; borrows validated section.
  IncrementalPyramid5(int id, const SolidSection* section);
};
}  // namespace fesa
#endif  // FESA_ELEMENTS_INCREMENTAL_PYRAMID5_H_
