/// @file
/// @brief Linear12 and finite selected-D14 Green/EB beam responses.
#ifndef FESA_ELEMENTS_BEAM3D_H_
#define FESA_ELEMENTS_BEAM3D_H_
#include "fesa/model/beam_section.h"
#include "fesa/model/element.h"
namespace fesa {
/// B1-B5 in docs/Linear static/beam3d_linear.md: small displacement/rotation,
/// homogeneous isotropic elasticity and constant section. Loaded evaluation
/// accepts caller-prepared equivalent vectors under the L1/L2 line-load
/// contract. Borrows nodes, section and its material from Domain, which must
/// outlive this. Owns its initial frame and result set; owns no evaluation
/// workspace.
class Beam3D final : public Element {
 public:
  /// Validates B1-B3 in docs/Linear static/beam3d_linear.md using MKL
  /// operations. Requires distinct endpoint IDs, finite positive length, and
  /// direction sine > sqrt(eps64). No axis substitution or geometry correction
  /// occurs. Caller owns *out; failure nulls it without changing any borrowed
  /// input.
  static StatusCode Create(int id, const Node* first, const Node* second,
                           const ElasticMaterial* material,
                           const BeamSection* section, Beam3D** out,
                           Diagnostic* d);
  ~Beam3D() override;
  Beam3D(const Beam3D&) = delete;
  Beam3D& operator=(const Beam3D&) = delete;
  Beam3D(Beam3D&&) = delete;
  Beam3D& operator=(Beam3D&&) = delete;

  /// Immutable beam layout: two nodes, six components each, type B33.
  ElementFamily family() const override;
  /// Borrowed original type name.
  const char* type_name() const override;
  /// Two connected endpoints.
  int node_count() const override;
  /// Six physical components [U1,U2,U3,UR1,UR2,UR3] per endpoint.
  int dofs_per_node() const override;
  /// Twelve local components, endpoint order then physical component order.
  int dof_count() const override;
  /// Two eta components for 14-component finite local evaluation.
  /// Physical layout stays 12; finite dynamic evaluation appends these tails.
  int internal_dof_count() const override;
  /// Selected-D Green/EB structural evaluation, with caller-owned eta2.
  bool SupportsTrialEvaluation() const override;
  /// Exact physical14 inertia and spatial Rayleigh evaluation, finite only.
  /// Legacy SupportsMass/ComputeMass remain unsupported in every mode.
  bool SupportsDynamicEvaluation(DynamicKinematics kinematics) const override;
  /// Borrows physical V/A12 and internal eta/V/A2 for this call only. Outputs
  /// are global spatial force14 and independent q/V/A partials14x14, without
  /// HHT weights. Consistent mass only; default/row-sum/linear unsupported.
  /// Optional worker-exclusive workspace lazily owns overwritten kinematic/
  /// time scratch. Null uses call-local ownership. No model pointer is retained.
  /// Required failure is cleared first. Supported finite preflight preserves
  /// outputs and ignores their old values; subsequent failure invalidates all
  /// requested blocks. Geometry/rotation and original backend diagnostics
  /// propagate. No state/history/official force or recovery storage changes.
  /// @see docs/Linear static/b33_dynamic_metric_contract.md.
  StatusCode EvaluateDynamic(
      const ElementMotionInput* input, MassMode mass_mode,
      const ElementWorkspace* workspace, DynamicEvaluationRequest request,
      ElementDynamicResult* out, TrialEvaluationFailure* failure,
      Diagnostic* d) const override;
  /// Borrows translations6/unit orientations2x4/eta2 for this call only.
  /// Requested global spatial force14 and tangent14x14 use [ua,wa,ub,wb,eta].
  /// Optional worker-exclusive workspace owns overwritten kinematic scratch;
  /// null uses call-local ownership. No
  /// state/history/result mutation or dynamic evaluation. Preflight preserves
  /// outputs; later failure invalidates them. Required failure distinguishes
  /// geometry/rotation branch from fatal resource/input/backend errors.
  StatusCode EvaluateTrial(const ElementTrialState* state,
                           const ElementWorkspace* workspace,
                           ElementEvaluationRequest request, Vector* force,
                           Matrix* tangent, TrialEvaluationFailure* failure,
                           Diagnostic* d) const override;
  /// Geometry-only K1 Gauss4 domain check, independent of material/density.
  bool SupportsTrialGeometryCheck() const override;
  /// Same borrowed pose, scratch ownership and typed failure as EvaluateTrial.
  /// Uses value-only K1 and endpoint transport; no official result/state change.
  StatusCode CheckTrialGeometry(const ElementTrialState* state,
                                const ElementWorkspace* workspace,
                                TrialEvaluationFailure* failure,
                                Diagnostic* d) const override;
  /// Recovers structural force12 in each CURRENT endpoint material frame and
  /// EndResults16=[-fa6,fb6,E*Green_a,E*Green_b,Green_a,Green_b], where
  /// Green=((1+eta/L)^2-1)/2. Initial LocalAxes remain immutable. All numerical
  /// values are computed in call-local scratch before copying official
  /// payloads. Failure leaves results unconfirmed; no AnalysisState commit is
  /// performed.
  StatusCode UpdateTrialResults(const ElementTrialState* state,
                                const ElementWorkspace* workspace,
                                TrialEvaluationFailure* failure,
                                Diagnostic* d) override;
  /// Borrowed endpoint; invalid index returns null.
  const Node* GetNode(int end) const override;
  /// Immutable initial length in model units L. Null output preserves state.
  StatusCode CharacteristicLength(double* out, Diagnostic* d) const override;
  /// Borrowed immutable section; Domain owns it and its material.
  const BeamSection* Section() const;
  /// Borrowed stored LOCAL forces [Fx,Fy,Fz,Mx,My,Mz] per endpoint.
  /// Linear results use initial axes; finite trials use CURRENT material axes.
  /// Valid only with confirmed AnalysisState, units F and F L.
  const Vector* InternalForce() const;
  /// Borrowed B5 [sa6,sb6,sigma_a,sigma_b,epsilon_a,epsilon_b].
  /// sa=-ga, sb=gb, order [N,V1,V2,Ttor,M1,M2]. Linear axial S/strain use
  /// signed N/A and N/(EA); finite results use E*Green and Green from eta.
  /// Valid only with confirmed AnalysisState; units F/F L, F/L^2, unitless.
  const Vector* EndResults() const;
  /// B3-B4 evaluation; follows Element size/coordinate/ownership contracts.
  StatusCode ComputeStiffness(ElementWorkspace* workspace, Matrix* global_k,
                              Diagnostic* d) const override;
  /// Computes B4 from the candidate only; stored results remain unchanged.
  StatusCode ComputeInternalForce(const Vector* global_u,
                                  ElementWorkspace* workspace,
                                  Vector* global_force,
                                  Diagnostic* d) const override;
  /// Stores B4 local generalized force. AnalysisState must be unconfirmed;
  /// includes B5 recovery before success. Failure does not commit.
  StatusCode UpdateResults(const Vector* global_u, ElementWorkspace* workspace,
                           Diagnostic* d) override;
  /// L2 actual global end force K*u-element_load. null load preserves B4.
  /// Load is borrowed Vector(12), distinct from output. Validation preserves
  /// output; later failure invalidates it. Stored results remain unchanged.
  StatusCode ComputeLoadedInternalForce(
      const Vector* global_u, const Vector* element_load,
      ElementWorkspace* workspace, Vector* global_force,
      Diagnostic* d) const override;
  /// Stores local Kl*ul-local_element_load and completes B5 recovery. Load must
  /// be finite Vector(12), distinct from owned InternalForce. Requires common
  /// unconfirmed state. Validation preserves storage; later failure can leave
  /// partial unconfirmed results. No state commit or history mutation.
  StatusCode UpdateLoadedResults(const Vector* global_u,
                                const Vector* local_element_load,
                                ElementWorkspace* workspace,
                                Diagnostic* d) override;
  /// Recovers B5 [sa6,sb6,sigma_a,b,epsilon_a,b] from stored local force.
  /// Requires successful current UpdateResults and unconfirmed AnalysisState.
  /// Failure leaves the whole result set unconfirmed, with no restoration.
  /// Unrepresentable B5 map coefficients (1/A or 1/(EA)) fail explicitly;
  /// this is the current recovery operator's range limit, not invalid physics
  /// or input. It can fail even if final stress/strain would be representable;
  /// coefficient underflow is not accepted as zero stress/strain.
  StatusCode RecoverEndForces(ElementWorkspace* workspace, Diagnostic* d);
  /// Immutable B1 frame, rows t,n1,n2; borrowed until destruction.
  const Matrix* LocalAxes() const;

 private:
  Beam3D(int id, const Node* first, const Node* second,
         const BeamSection* section);
  /// Fills the fixed B3/B4 coefficient patterns in borrowed worker storage.
  /// No allocation; failure may overwrite scratch only.
  StatusCode Prepare(ElementWorkspace* workspace, Diagnostic* d) const;
  /// Computes local_u and local_force into scratch; no stored-result mutation.
  StatusCode EvaluateLocal(const Vector* global_u, ElementWorkspace* workspace,
                           Diagnostic* d) const;
  /// Allocates owned result storage during creation; failure discards this.
  StatusCode CreateResults(Diagnostic* d);
  /// Stores candidate local force without recovery; used by both update paths.
  StatusCode StoreLocalForce(const Vector* global_u,
                             ElementWorkspace* workspace, Diagnostic* d);
  /// Mutable owned force/result storage, confined to beam implementation.
  Vector* MutableInternalForce();
  /// Mutable owned B5 storage.
  Vector* MutableEndResults();
  const Node* first_;                 ///< Borrowed from Domain.
  const Node* second_;                ///< Borrowed from Domain.
  const BeamSection* section_;        ///< Borrowed from Domain.
  Vector* internal_force_ = nullptr;  ///< Owned local force, size 12.
  Vector* end_results_ = nullptr;     ///< Owned B5 results, size 16.
  double length_ = 0;  ///< Immutable initial centerline length, model units L.
  Matrix* local_axes_ = nullptr;  ///< Owned immutable B1 frame, nine doubles.
};
}  // namespace fesa
#endif  // FESA_ELEMENTS_BEAM3D_H_
