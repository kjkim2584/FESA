/// @file
/// @brief Isolated C3D8 factory for incremental elastic static response.
#ifndef FESA_ELEMENTS_INCREMENTAL_HEXA8_H_
#define FESA_ELEMENTS_INCREMENTAL_HEXA8_H_
#include "fesa/elements/incremental_solid_element.h"
namespace fesa {
/// Eight-node coupled HEX8 solid. Existing Hexa8/StVK classes are not base
/// classes.
class IncrementalHexa8 final : public IncrementalSolidElement {
 public:
  /// Borrows nodes/section through destruction. Caller owns success; Domain
  /// registration transfers ownership only on success. Failure cleans partial
  /// resources and nulls out. Requires distinct3DOF nodes and kHexa8Full.
  static StatusCode Create(int id, const Node* const* nodes, int node_count,
                           const SolidSection* section,
                           SolidFormulation formulation, IncrementalHexa8** out,
                           Diagnostic* d);
  ~IncrementalHexa8() override = default;
  IncrementalHexa8(const IncrementalHexa8&) = delete;
  IncrementalHexa8& operator=(const IncrementalHexa8&) = delete;
  IncrementalHexa8(IncrementalHexa8&&) = delete;
  IncrementalHexa8& operator=(IncrementalHexa8&&) = delete;

 private:
  /// Factory-only construction; borrows validated section.
  IncrementalHexa8(int id, const SolidSection* section);
};
}  // namespace fesa
#endif  // FESA_ELEMENTS_INCREMENTAL_HEXA8_H_
