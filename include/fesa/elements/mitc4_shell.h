/// @file
/// @brief Four-node MITC4 definition, linear and finite-state operators.
#ifndef FESA_ELEMENTS_MITC4_SHELL_H_
#define FESA_ELEMENTS_MITC4_SHELL_H_
#include "fesa/elements/mitc4_results.h"
#include "fesa/model/element.h"
#include "fesa/model/shell_section.h"
namespace fesa {
/// Borrows four nodes, section and identical material through destruction.
/// No reference/workspace/cache storage. Initial kd is an immutable scalar
/// model coefficient. Parser/writer are separate stages.
/// Linear operators use global node-major [U1,U2,U3,UR1,UR2,UR3] input and
/// [F1,F2,F3,M1,M2,M3] output. U has length units, UR radians, F force,
/// M force*length. See docs/MITC4/initial_kernel_contract.md for formulation.
class Mitc4Shell final : public Element {
 public:
  /// Positive ID, four distinct6DOF nodes and nonnull section required.
  /// Initial flat geometry must pass the documented guards. Computes fixed kd
  /// from initial physical K, discards temporary K. Caller owns successful*out;
  /// failure nulls output and releases partial resources, preserving inputs.
  static StatusCode Create(int id, const Node* const* nodes, int count,
                           const ShellSection* section, Mitc4Shell** out,
                           Diagnostic* d);
  ~Mitc4Shell() override;
  Mitc4Shell(const Mitc4Shell&) = delete;
  Mitc4Shell& operator=(const Mitc4Shell&) = delete;
  Mitc4Shell(Mitc4Shell&&) = delete;
  Mitc4Shell& operator=(Mitc4Shell&&) = delete;
  /// Initial shell family and supported input TYPE identity.
  ElementFamily family() const override;
  /// Static borrowed S4 string; not full Abaqus formulation compatibility.
  const char* type_name() const override;
  /// Four connected nodes; six physical components per node,24 node-major.
  int node_count() const override;
  int dofs_per_node() const override;
  int dof_count() const override;
  /// Borrowed original node, invalid local index returns null.
  const Node* GetNode(int local) const override;
  /// Initial mean edge length Lc; failure preserves scalar. No current state.
  StatusCode CharacteristicLength(double* out, Diagnostic* d) const override;
  /// Full initial global24x24 K; workspace unused and may be null.
  /// Preflight/geometry errors preserve out, backend failure invalidates it.
  StatusCode ComputeStiffness(ElementWorkspace* workspace, Matrix* out,
                              Diagnostic* d) const override;
  /// Direct linear force from borrowed global24 U/UR and fixed initial kd.
  /// Separate out24 required; no K/M/results. Same failure contract as K.
  StatusCode ComputeInternalForce(const Vector* u, ElementWorkspace* workspace,
                                  Vector* out, Diagnostic* d) const override;
  /// Supports finite translations/unit quaternion/body-spin trial evaluation.
  bool SupportsTrialEvaluation() const override;
  /// Recomputes requested force24/tangent24x24 with frozen initial kd, no
  /// density/mass/official recovery. Input3n translation/bodyspin and n×4 q/qc
  /// follow element_trial.h. Tangent is the unsymmetrized spatial Jacobian.
  /// Same output/failure contract as Element; initial flat definition is fixed,
  /// current nonplanar geometry is validated at Gauss/tying points.
  /// Workspace is unused, may be null, and is never retained.
  StatusCode EvaluateTrial(const ElementTrialState* state,
                           const ElementWorkspace* workspace,
                           ElementEvaluationRequest request, Vector* force,
                           Matrix* tangent, TrialEvaluationFailure* failure,
                           Diagnostic* d) const override;
  /// Consistent and local-component row-sum mass modes.
  bool SupportsMass(MassMode mode) const override;
  /// Linear and finite-rotation inertia/Rayleigh evaluation.
  /// density and unit-quaternion validity are checked on evaluation.
  bool SupportsDynamicEvaluation(DynamicKinematics kinematics) const override;
  /// Linear E0: fI=M0*A, MI=M0, fD=C*v, CD=C; other derivatives zero.
  /// C=alpha_R*M0+beta_R*K_el, excluding transverse shear/drilling (RD12).
  /// Finite E0 uses current directors, gyro moment and independent GI/CI/MI
  /// (I2/I3). RD1--RD10 provide current-plane fD and spatial GD/CD.
  /// Global24
  /// node-major [translation,small rotation] or [translation,spatial rotation]
  /// with F/FL force rows. No physical
  /// drilling inertia is added. Uses the existing consistent/local-row-sum M0.
  /// Borrows all inputs/outputs during this call; workspace is unused and not
  /// retained. No structural force, HHT weighting, results/state or I/O change.
  /// Common dynamic shape/alias/request checks apply. Preflight preserves out;
  /// calculation failure invalidates all requested outputs. failure is required
  /// and set to kNone first. Reference-area inertia/mass drag do not evaluate
  /// current volume/Log chart. Nonzero beta_R checks current volume at Gauss
  /// points and reports typed geometry; structural tying checks stay separate.
  /// Quaternion/resource/backend errors retain their actual fatal status.
  /// @see docs/HHT/implementation_contract.md, rotation_inertia_damping.md,
  /// mitc4_rayleigh.md.
  StatusCode EvaluateDynamic(const ElementMotionInput* input, MassMode mode,
                             const ElementWorkspace* workspace,
                             DynamicEvaluationRequest request,
                             ElementDynamicResult* out,
                             TrialEvaluationFailure* failure,
                             Diagnostic* d) const override;
  /// Initial global24x24 mass; missing density fails only this request.
  /// Workspace unused, preflight preserves out, backend invalidates it.
  StatusCode ComputeMass(MassMode mode, ElementWorkspace* workspace,
                         Matrix* out, Diagnostic* d) const override;
  /// Same local-component lumping contract as initial kernels; preserves C.
  StatusCode LumpFromConsistentMass(const Matrix* consistent, Matrix* lumped,
                                    Diagnostic* d) const override;
  /// Recalculates official linear payload in one owned set. Caller marks
  /// common state unconfirmed first; backend/recovery failure may leave partial
  /// unconfirmed results. No global commit or file I/O here.
  StatusCode UpdateResults(const Vector* u, ElementWorkspace* workspace,
                           Diagnostic* d) override;
  /// Recalculates finite-trial force and current S/LE payload with fixed kd.
  /// No K/M, state/history/commit or I/O. Caller marks results unconfirmed.
  /// Borrowed state; required typed failure. Same validity/failure lifetime as
  /// Element::UpdateTrialResults; no principal-UR reconstruction or cache.
  /// Workspace is unused, may be null, and is never retained.
  StatusCode UpdateTrialResults(const ElementTrialState* state,
                                const ElementWorkspace* workspace,
                                TrialEvaluationFailure* failure,
                                Diagnostic* d) override;
  /// Immutable initial output-model axes3x3; borrowed through destruction.
  /// Published model metadata, not used as a numerical operator cache.
  const Matrix* initial_output_axes() const;
  /// Recomputes the initial global-to-local3x3 row basis [a,b,n] from immutable
  /// nodes, for call-local numerical consumers. No stiffness, mass, material,
  /// state/results or retained scratch. Caller owns required out3x3; invalid
  /// shape preserves it, backend failure invalidates it. Diagnostic optional.
  StatusCode ComputeInitialFrame(Matrix* out, Diagnostic* d) const;
  /// Borrowed official output payload; validity comes from common state.
  const Mitc4Results* results() const;
  /// Borrowed immutable definition references, valid through destruction.
  const ShellSection* section() const;
  const ElasticMaterial* material() const;
  /// Fixed initial drilling coefficient in moment/rad.
  double drilling_stiffness() const;

 private:
  Mitc4Shell(int id, const ShellSection* section);
  const Node* nodes_[4] = {};        ///< Borrowed; no nodal-normal averaging.
  const ShellSection* section_;      ///< Borrowed.
  const ElasticMaterial* material_;  ///< Same borrowed pointer as section.
  double drilling_ = 0.0;
  Matrix* initial_axes_ = nullptr;  ///< Owned immutable model-output metadata.
  Mitc4Results* results_ =
      nullptr;  ///< Owned official result set, not scratch.
};
}  // namespace fesa
#endif  // FESA_ELEMENTS_MITC4_SHELL_H_
