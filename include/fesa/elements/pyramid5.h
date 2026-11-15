/// @file
/// @brief Rational five-node pyramid with full or raw reduced stiffness.
#ifndef FESA_ELEMENTS_PYRAMID5_H_
#define FESA_ELEMENTS_PYRAMID5_H_
#include "fesa/elements/solid_element.h"
namespace fesa {
/// Borrows four base nodes and apex plus Domain section/material definitions.
/// Full stiffness/results use Felippa five points; mass uses eight points.
/// SolidElement owns official results. Raw reduced has no hourglass term.
class Pyramid5 final : public SolidElement {
 public:
  /// Requires five distinct 3-DOF nodes and Pyramid5 full/reduced formulation.
  /// Geometry guards exclude the natural apex gradient, never physical
  /// collapse. Caller owns success; failure nulls output and releases partial
  /// storage.
  static StatusCode Create(int id, const Node* const* nodes, int node_count,
                           const SolidSection* section,
                           SolidFormulation formulation, Pyramid5** out,
                           Diagnostic* d);
  ~Pyramid5() override = default;
  Pyramid5(const Pyramid5&) = delete;
  Pyramid5& operator=(const Pyramid5&) = delete;
  Pyramid5(Pyramid5&&) = delete;
  Pyramid5& operator=(Pyramid5&&) = delete;

 private:
  Pyramid5(int id, const SolidSection* section, SolidFormulation formulation);
};
}  // namespace fesa
#endif  // FESA_ELEMENTS_PYRAMID5_H_
