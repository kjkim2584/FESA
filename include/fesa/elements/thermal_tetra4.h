/// @file
/// @brief Linear DC3D4 conduction and equal V/4 capacity-station weights.
#ifndef FESA_ELEMENTS_THERMAL_TETRA4_H_
#define FESA_ELEMENTS_THERMAL_TETRA4_H_

#include "fesa/model/thermal_element.h"
#include "fesa/model/thermal_section.h"

namespace fesa {
/// Immutable linear tetrahedron, formulation §3.2 and Abaqus P8 evidence.
class ThermalTetra4 final : public ThermalElement {
 public:
  /// Caller owns success; failure nulls *out and cleans partial resources.
  /// Requires four distinct one-temperature nodes and positive volume.
  /// Section is borrowed for creation; nodes/material must outlive this.
  static StatusCode Create(int id, const Node* const* nodes, int count,
                           const ThermalSolidSection* section,
                           ThermalTetra4** out, Diagnostic* d);
  ~ThermalTetra4() override = default;
  /// DC3D4 keyword identity.
  ThermalElementType type() const override;
  const char* type_name() const override;
  /// Four nodes with one temperature component each.
  int node_count() const override;
  int dofs_per_node() const override;
  /// Borrowed node, null for an invalid local index.
  const Node* GetNode(int local) const override;
  /// Borrowed material, valid until Domain destruction.
  const ThermalMaterial* material() const override;
  /// Positive volume; kInputError for degenerate/inverted geometry.
  StatusCode ValidateGeometry(Diagnostic* d) const override;
  /// Overwrites 4x4 constant-property K [E/T/theta]; failure may invalidate
  /// out. Optional theta has four finite components. M2 remains gated.
  StatusCode ComputeConduction(const Vector* element_theta, Matrix* out,
                               Diagnostic* d) const override;
  /// Overwrites four equal V/4 [L^3] weights; failure may invalidate out.
  StatusCode ComputeCapacityWeights(Vector* out, Diagnostic* d) const override;

  /// Formulation8.4 face-area weights; full DOF order, zero off this face.
  /// Borrows normalized face label/output for this call. No model mutation;
  /// unsupported label preserves out, numerical failure may invalidate out.
  StatusCode ComputeSurfaceWeights(const char* face_label, Vector* out,
                                   Diagnostic* d) const override;

 private:
  ThermalTetra4(int id, const ThermalMaterial* material);
  /// Computes volume and optional 3x4 physical gradients via LU. Failure
  /// preserves scalar volume, may invalidate gradient; no model mutation.
  StatusCode Geometry(Matrix* gradient, double* volume, Diagnostic* d) const;
  const Node* nodes_[4] = {};        ///< Borrowed connectivity.
  const ThermalMaterial* material_;  ///< Borrowed; Domain owns material.
};
}  // namespace fesa

#endif  // FESA_ELEMENTS_THERMAL_TETRA4_H_
