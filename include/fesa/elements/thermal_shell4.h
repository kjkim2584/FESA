/// @file
/// @brief Planar DS4 conduction and Simpson station volume weights.
#ifndef FESA_ELEMENTS_THERMAL_SHELL4_H_
#define FESA_ELEMENTS_THERMAL_SHELL4_H_

#include "fesa/model/thermal_element.h"
#include "fesa/model/thermal_section.h"

namespace fesa {
/// Immutable four-node shell with node-major thickness temperatures.
/// Uses formulation 3.3 and the ROT/MESH gates; no retained trial state.
class ThermalShell4 final : public ThermalElement {
 public:
  /// Four distinct nodes with section_points thermal slots required.
  /// Caller owns success; failure nulls *out and releases partial resources.
  /// Nodes/section/material are borrowed and must outlive this element.
  /// Non-planar geometry is kUnsupported; degenerate geometry kInputError.
  static StatusCode Create(int id, const Node* const* nodes, int count,
                           const ThermalShellSection* section,
                           ThermalShell4** out, Diagnostic* d);
  ~ThermalShell4() override = default;
  /// DS4 identity, four nodes and ns slots per node.
  ThermalElementType type() const override;
  const char* type_name() const override;
  int node_count() const override;
  int dofs_per_node() const override;
  /// Borrowed connectivity in input order; null for invalid local index.
  const Node* GetNode(int local) const override;
  /// Borrowed immutable material and scalar section thickness [L].
  const ThermalMaterial* material() const override;
  StatusCode GetShellThickness(double* out, Diagnostic* d) const override;
  /// Checks normalized planarity and positive 2x2 Gauss Jacobians.
  /// No state mutation; diagnostics identify this element on failure.
  StatusCode ValidateGeometry(Diagnostic* d) const override;
  /// Overwrites 4ns-square symmetric K [E/T/theta] with the weighted Gram.
  /// Null theta selects M1; non-null finite theta must have 4ns entries.
  /// M2 table evaluation remains gated. Failure may invalidate out only.
  /// Nonfinite K or a non-normal positive diagonal is kLibraryFailure;
  /// split-root coefficients include sqrt(k), formulation 3.3 range contract.
  StatusCode ComputeConduction(const Vector* theta, Matrix* out,
                               Diagnostic* d) const override;
  /// Overwrites 4ns volume weights A_a W_i [L^3], without rho*c.
  /// Caller owns output; failure may invalidate it; model never changes.
  /// Non-normal Gauss volume contributions report kLibraryFailure explicitly.
  StatusCode ComputeCapacityWeights(Vector* out, Diagnostic* d) const override;

  /// Formulation8.4 face-area weights; full DOF order, zero off this face.
  /// Borrows normalized face label/output for this call. No model mutation;
  /// unsupported label preserves out, numerical failure may invalidate out.
  StatusCode ComputeSurfaceWeights(const char* face_label, Vector* out,
                                   Diagnostic* d) const override;

 private:
  ThermalShell4(int id, const ThermalShellSection* section);
  const Node* nodes_[4] = {};           ///< Borrowed until element destruction.
  const ThermalShellSection* section_;  ///< Borrowed from Domain.
};
}  // namespace fesa
#endif  // FESA_ELEMENTS_THERMAL_SHELL4_H_
