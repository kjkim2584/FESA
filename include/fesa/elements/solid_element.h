/// @file
/// @brief Three-translation solid definitions and typed official result
/// storage.
#ifndef FESA_ELEMENTS_SOLID_ELEMENT_H_
#define FESA_ELEMENTS_SOLID_ELEMENT_H_
#include "fesa/elements/solid_reference_data.h"
#include "fesa/model/element.h"
#include "fesa/model/solid_section.h"
namespace fesa {
/// Borrows nodes, SolidSection and identical ElasticMaterial from Domain.
/// Derived factories validate geometry. Owns one official result set; common
/// AnalysisState controls validity. No persistent geometry, D or K cache.
class SolidElement : public Element {
 public:
  ~SolidElement() override;
  SolidElement(const SolidElement&) = delete;
  SolidElement& operator=(const SolidElement&) = delete;
  SolidElement(SolidElement&&) = delete;
  SolidElement& operator=(SolidElement&&) = delete;
  /// Immutable solid family and TYPE identity.
  ElementFamily family() const override;
  /// Borrowed TYPE name, valid until destruction.
  const char* type_name() const override;
  /// Immutable TYPE-selected mathematical formulation.
  SolidFormulation formulation() const;
  /// Immutable node count.
  int node_count() const override;
  /// Three physical translations, no rotations.
  int dofs_per_node() const override;
  /// Node-major [U1,U2,U3] local size 3n.
  int dof_count() const override;
  /// Borrowed node; invalid local index returns null.
  const Node* GetNode(int local) const override;
  /// User-adopted Tetra/Hexa/Wedge/Pyramid compatible convergence-length
  /// estimate, L. Recomputed each call; Hexa/Wedge need center J, no LU/D/B/K/M
  /// or results. Pyramid uses four bilinear base samples and apex distances, no
  /// apex J. Tetra uses ordered corner edge norms without J. Success writes
  /// *out; validation/backend failure preserves it. See geometry_contract.md.
  StatusCode CharacteristicLength(double* out, Diagnostic* d) const override;
  /// Overwrites global 3n square K, F/L. Requires matching stiffness context.
  /// Preflight preserves output; backend failure invalidates output.
  StatusCode ComputeStiffness(ElementWorkspace* context, Matrix* out,
                              Diagnostic* d) const override;
  /// Direct quadrature resisting force in global axes, F. No K/M or official
  /// recovery is performed. Input/output disjoint; same failure contract.
  StatusCode ComputeInternalForce(const Vector* u, ElementWorkspace* context,
                                  Vector* out, Diagnostic* d) const override;
  /// Nonlinear trial is implemented for all five full Solid types.
  /// Reduced and retired formulations remain unsupported.
  bool SupportsTrialEvaluation() const override;
  /// M1-M8 total-Lagrangian StVK, with initial-volume F-bar for C3D8.
  /// Borrows state/context for this call only. Requires translations3n in
  /// global length units and null rotation/history fields, matching immutable
  /// stiffness reference data, and nonnull failure. Requested force3n has force
  /// units; tangent3nx3n is the full derivative in force/length units.
  /// Unrequested outputs must be null; input/output buffers must be disjoint.
  /// Preflight preserves outputs; computational failure invalidates requested
  /// outputs. Nonpositive det(F) sets kGeometry; resource/backend failures do
  /// not. No official results/state/cache change. Force-only skips tangent,
  /// mass, logarithmic strain and recovery. See nonlinear_formulation M1-M8.
  StatusCode EvaluateTrial(const ElementTrialState* state,
                           const ElementWorkspace* workspace,
                           ElementEvaluationRequest request, Vector* force,
                           Matrix* tangent, TrialEvaluationFailure* failure,
                           Diagnostic* d) const override;
  /// Q23/Q53 sample checks are available for all five full solid types.
  bool SupportsTrialGeometryCheck() const override;
  /// Borrows translations-only state and exact TrialGeometry catalog. Checks
  /// original X+u with current h; Tetra10 midsides are never straightened.
  /// No material/force/results/state changes. Reports the first failed sample,
  /// typed geometry vs fatal backend failure, and measured determinant only
  /// when available. A successful sample check is not global noninversion proof.
  StatusCode CheckTrialGeometry(const ElementTrialState* state,
                                const ElementWorkspace* workspace,
                                TrialEvaluationFailure* failure,
                                Diagnostic* d) const override;
  /// Optional consistent/row-sum modes. Tetra10 is consistent-only.
  bool SupportsMass(MassMode mode) const override;
  /// Overwrites global mass, F T^2/L. Requires density and matching mass data;
  /// no D/B/HG or official recovery. Preflight preserves output.
  StatusCode ComputeMass(MassMode mode, ElementWorkspace* context, Matrix* out,
                         Diagnostic* d) const override;
  /// Linear dynamics for all seven current Solid formulations; finite
  /// kinematics and the retired reduced Tetra10 identifier are unsupported.
  bool SupportsDynamicEvaluation(DynamicKinematics kinematics) const override;
  /// E0 constant reference mass: fI=M*a, MI=M; fD=(alpha_R*M+beta_R*K_el)*v.
  /// K_el uses the existing physical linear elastic quadrature, excluding
  /// C3D8R hourglass stabilization. CD is this material damping matrix; other
  /// partials are zero. Official stress/force results remain unchanged.
  /// Borrows global node-major 3n motion and reference context for this call.
  /// Overwrites requested forces (F) and derivatives per element_dynamic.h.
  /// No recovery, official state change or HHT weighting. Preflight
  /// preserves outputs; computational failure invalidates requested outputs.
  /// Requires failure; mass/backend errors leave it kNone. See
  /// docs/HHT/linear_solid_dynamics.md for units and supported mass modes.
  StatusCode EvaluateDynamic(const ElementMotionInput* input,
                             MassMode mass_mode,
                             const ElementWorkspace* workspace,
                             DynamicEvaluationRequest request,
                             ElementDynamicResult* out,
                             TrialEvaluationFailure* failure,
                             Diagnostic* d) const override;
  /// Row sums into a separate diagonal output without new geometry/integration.
  /// Unsupported/shape/alias failures preserve both matrices.
  StatusCode LumpFromConsistentMass(const Matrix* consistent, Matrix* lumped,
                                    Diagnostic* d) const override;
  /// Completes internal force, IP and element-node stress/strain. Caller marks
  /// common state unconfirmed first. Failure may leave partial unconfirmed
  /// results; no state commit or cross-call factor reuse occurs here.
  StatusCode UpdateResults(const Vector* u, ElementWorkspace* context,
                           Diagnostic* d) override;
  /// Overwrites official force and global Cauchy S/logarithmic LE at IPs and
  /// element nodes for all five full types. Uses the same state/context as EvaluateTrial;
  /// translations cannot alias owned force. Caller marks results unconfirmed
  /// first. Preflight preserves payload; calculation failure may leave partial
  /// unconfirmed results. Never commits state or writes files. Recovery failure
  /// is fatal to the caller phase even if failure reports geometry. No tangent
  /// or mass is computed. Tensor ordering/shear follows the result getters.
  StatusCode UpdateTrialResults(const ElementTrialState* state,
                                const ElementWorkspace* workspace,
                                TrialEvaluationFailure* failure,
                                Diagnostic* d) override;
  /// Physical IP count, distinct from geometry quadrature count.
  int integration_point_count() const;
  /// Borrowed global node-major internal resisting force, F, size 3n.
  const Vector* internal_force() const;
  /// Borrowed g by 6 global stress [11,22,33,12,13,23], F/L^2.
  /// Linear update uses small-strain stress; trial update uses Cauchy S.
  const Matrix* ip_stress() const;
  /// Borrowed g by 6 global tensor strain, same order; dimensionless.
  /// Linear update uses infinitesimal strain; trial update uses logarithmic LE.
  const Matrix* ip_strain() const;
  /// Borrowed n by 6 element-node stress, never shared-node averaged.
  const Matrix* nodal_stress() const;
  /// Borrowed n by 6 element-node tensor strain, never averaged.
  const Matrix* nodal_strain() const;
  /// Borrowed immutable section and direct identical material reference.
  const SolidSection* section() const;
  /// Borrowed ElasticMaterial, which must outlive this element.
  const ElasticMaterial* material() const;

 protected:
  /// Derived factory owns construction and validates immutable inputs.
  SolidElement(int id, int node_count, int point_count,
               const SolidSection* section, SolidFormulation formulation);
  /// Copies borrowed connectivity and allocates official result storage. Failed
  /// partial storage remains owned for destruction by the derived factory.
  StatusCode Initialize(const Node* const* nodes, Diagnostic* d);
  /// Checks purpose/boundary geometry without D/K/M; factory discards failures.
  StatusCode ValidateGeometry(Diagnostic* d) const;

 private:
  /// Call-local stiffness evaluation borrowing immutable reference context.
  /// Same output/failure contract as ComputeStiffness; no context mutation.
  /// include_hourglass=false omits HG construction/contribution for material
  /// beta damping only; the default preserves the structural operator.
  StatusCode ComputeStiffnessReadOnly(const ElementWorkspace* context,
                                      Matrix* out, Diagnostic* d,
                                      bool include_hourglass = true) const;
  /// Shared mass implementation; borrows immutable reference data only.
  /// Same output and failure contract as ComputeMass.
  StatusCode ComputeMassReadOnly(MassMode mode, const ElementWorkspace* context,
                                 Matrix* out, Diagnostic* d) const;
  int node_count_;
  int point_count_;
  SolidFormulation formulation_;
  const Node** nodes_ = nullptr;      ///< Owned pointer array, borrowed nodes.
  const SolidSection* section_;       ///< Borrowed from Domain.
  const ElasticMaterial* material_;   ///< Direct borrowed identical material.
  Vector* internal_force_ = nullptr;  ///< Owned global official force.
  Matrix* ip_stress_ = nullptr;  ///< Owned global tensor-order official stress.
  Matrix* ip_strain_ = nullptr;  ///< Owned global tensor-order official strain.
  Matrix* nodal_stress_ = nullptr;  ///< Owned element-node official stress.
  Matrix* nodal_strain_ =
      nullptr;  ///< Owned element-node official tensor strain.
};
}  // namespace fesa
#endif  // FESA_ELEMENTS_SOLID_ELEMENT_H_
