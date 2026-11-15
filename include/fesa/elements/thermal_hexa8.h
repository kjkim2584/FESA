/// @file
/// @brief DC3D8 conduction and nearest-Gauss capacity-station weights.
#ifndef FESA_ELEMENTS_THERMAL_HEXA8_H_
#define FESA_ELEMENTS_THERMAL_HEXA8_H_

#include "fesa/model/thermal_element.h"
#include "fesa/model/thermal_section.h"

namespace fesa {
/// Immutable eight-node solid, using formulation §3.1 and Abaqus P12 evidence.
class ThermalHexa8 final : public ThermalElement {
 public:
  /// Caller owns success; failure nulls *out and frees partial resources.
  /// Requires eight distinct one-temperature nodes and a valid section.
  /// Section is borrowed during creation only; nodes/material must outlive
  /// the element. Invalid geometry is kInputError, with element ID in d.
  static StatusCode Create(int id, const Node* const* nodes, int count,
                           const ThermalSolidSection* section,
                           ThermalHexa8** out, Diagnostic* d);
  ~ThermalHexa8() override = default;
  /// DC3D8 input identity.
  ThermalElementType type() const override;
  const char* type_name() const override;
  /// Eight nodes, one temperature per node.
  int node_count() const override;
  int dofs_per_node() const override;
  /// Borrowed connectivity, or null for an invalid local index.
  const Node* GetNode(int local) const override;
  /// Borrowed material until Domain destruction.
  const ThermalMaterial* material() const override;
  /// Positive detJ at every 2x2x2 Gauss station; no model mutation.
  StatusCode ValidateGeometry(Diagnostic* d) const override;
  /// Overwrites 8x8 K using constant k and 2x2x2 integration. M2 material
  /// evaluation remains unsupported until the corresponding gate passes.
  StatusCode ComputeConduction(const Vector* element_theta, Matrix* out,
                               Diagnostic* d) const override;
  /// Overwrites eight detJ values at the Gauss station nearest each node.
  StatusCode ComputeCapacityWeights(Vector* out, Diagnostic* d) const override;

  /// Formulation8.4 face-area weights; full DOF order, zero off this face.
  /// Borrows normalized face label/output for this call. No model mutation;
  /// unsupported label preserves out, numerical failure may invalidate out.
  StatusCode ComputeSurfaceWeights(const char* face_label, Vector* out,
                                   Diagnostic* d) const override;

 private:
  ThermalHexa8(int id, const ThermalMaterial* material);
  const Node* nodes_[8] = {};        ///< Borrowed, input connectivity order.
  const ThermalMaterial* material_;  ///< Borrowed; Domain owns material.
};
}  // namespace fesa

#endif  // FESA_ELEMENTS_THERMAL_HEXA8_H_
