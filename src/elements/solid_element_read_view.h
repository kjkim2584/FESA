/// @file
/// @brief Call-lifetime read-only adapter for both isolated Solid families.
#ifndef FESA_ELEMENTS_SOLID_ELEMENT_READ_VIEW_H_
#define FESA_ELEMENTS_SOLID_ELEMENT_READ_VIEW_H_
#include "fesa/elements/incremental_solid_element.h"
#include "fesa/elements/solid_element.h"
namespace fesa {
/// Borrows one element; owns no definition/state/result. Bind before getters.
/// Element must outlive the view. State, not getters, determines result validity.
class SolidElementReadView {
 public:
  /// Supports existing and incremental Solid only. Failure preserves out.
  static StatusCode Bind(const Element* element, SolidElementReadView* out,
                         Diagnostic* d);
  /// Borrowed polymorphic source, valid until its destruction.
  const Element* element() const;
  /// Immutable formulation, connectivity and physical IP layout.
  SolidFormulation formulation() const;
  int node_count() const;
  int integration_point_count() const;
  /// Borrowed node; invalid index returns null.
  const Node* GetNode(int local) const;
  /// Borrowed immutable Domain definitions.
  const SolidSection* section() const;
  const ElasticMaterial* material() const;
  /// Borrowed official3n global resisting force, F.
  const Vector* internal_force() const;
  /// Borrowed ng×6 global S/LE, tensor order [11,22,33,12,13,23].
  const Matrix* ip_stress() const;
  const Matrix* ip_strain() const;
  /// Borrowed n×6 element-node tensors, never averaged across elements.
  const Matrix* nodal_stress() const;
  const Matrix* nodal_strain() const;
 private:
  const SolidElement* legacy_ = nullptr;  ///< Borrowed, mutually exclusive.
  const IncrementalSolidElement* incremental_ = nullptr;  ///< Borrowed.
};
}  // namespace fesa
#endif  // FESA_ELEMENTS_SOLID_ELEMENT_READ_VIEW_H_
