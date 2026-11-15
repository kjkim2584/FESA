/// @file
/// @brief Family-independent element identity, connectivity and optional math.
#ifndef FESA_MODEL_ELEMENT_H_
#define FESA_MODEL_ELEMENT_H_

#include "fesa/elements/element_workspace.h"
#include "fesa/model/element_dynamic.h"
#include "fesa/model/element_trial.h"
#include "fesa/model/node.h"

namespace fesa {

/// Physical family; result payloads retain family-specific contracts.
enum class ElementFamily { kBeam, kShell, kSolid };
/// Optional element mass construction, independent of static stiffness.
/// kDefault is an Analysis request only: resolve each element's policy before
/// calling physical Element mass/dynamic operations. Those retain explicit
/// consistent/row-sum capability contracts.
enum class MassMode { kConsistent, kRowSum, kDefault };

/// Domain-owned polymorphic element. Derived objects own their connectivity
/// references and official results. AnalysisState controls result validity.
class Element {
 public:
  /// Enables Domain deletion through the common interface.
  virtual ~Element() = default;
  Element(const Element&) = delete;
  Element& operator=(const Element&) = delete;
  Element(Element&&) = delete;
  Element& operator=(Element&&) = delete;
  /// Immutable original element ID.
  int Id() const;
  /// Immutable physical family.
  virtual ElementFamily family() const = 0;
  /// Borrowed immutable input type name, valid until destruction.
  virtual const char* type_name() const = 0;
  /// Number of connected nodes.
  virtual int node_count() const = 0;
  /// Physical components per node.
  virtual int dofs_per_node() const = 0;
  /// Immutable physical local components, connectivity then component order.
  virtual int dof_count() const = 0;
  /// Optional future finite local components beyond immutable physical DOFs.
  /// Default zero; metadata alone claims no trial/dynamic/mass capability.
  virtual int internal_dof_count() const;
  /// Borrowed connected node; invalid local index returns null.
  virtual const Node* GetNode(int local_node) const = 0;
  /// Writes all original node IDs in connectivity order. Insufficient capacity,
  /// null output or missing connectivity preserves the entire output.
  virtual StatusCode GetDofNodes(int* node_ids, int capacity,
                                 Diagnostic* d) const;
  /// Initial convergence length in model units L; default unsupported preserves
  /// scalar output. No generic length formula is assumed.
  virtual StatusCode CharacteristicLength(double* out, Diagnostic* d) const;
  /// Optional mass modes. Unsupported operations preserve output matrices.
  virtual bool SupportsMass(MassMode mode) const;
  /// Overwrites mass only for supported modes; default returns kUnsupported.
  virtual StatusCode ComputeMass(MassMode mode, ElementWorkspace* workspace,
                                 Matrix* mass, Diagnostic* d) const;
  /// Row-sums existing consistent mass into a separate output without new
  /// integration. Default unsupported preserves both matrices.
  virtual StatusCode LumpFromConsistentMass(const Matrix* consistent,
                                            Matrix* lumped,
                                            Diagnostic* d) const;
  /// Overwrites a dof_count square global matrix; derived formulation defines
  /// units. Validation preserves output, backend failure invalidates output.
  /// No stored results change. Workspace is borrowed for this call only.
  virtual StatusCode ComputeStiffness(ElementWorkspace* workspace,
                                      Matrix* global_k, Diagnostic* d) const;
  /// Evaluates global resisting force from the current candidate without
  /// mutating official results. Derived formulation defines units/components.
  /// Input/output must be separate; validation preserves output, later failure
  /// invalidates it. Default unsupported never changes output.
  virtual StatusCode ComputeInternalForce(const Vector* global_u,
                                          ElementWorkspace* workspace,
                                          Vector* global_force,
                                          Diagnostic* d) const;
  /// Optional finite-state evaluation; default false, linear families
  /// unchanged.
  virtual bool SupportsTrialEvaluation() const;
  /// Borrowed finite translations/q/a input; requested global resisting force
  /// and its spatial-correction Jacobian. Unrequested outputs null, outputs
  /// disjoint from all inputs. No mass, official recovery, state mutation or
  /// retained numerical scratch. Workspace is immutable, borrowed for the call
  /// only; Shell permits null. Tangent may be nonsymmetric. Derived contract
  /// defines dimensions; preflight preserves outputs, later failure discards
  /// them. Required failure flag distinguishes geometry/chart cutback from
  /// fatal input/resource/library status. Default clears flag if provided,
  /// returns unsupported and preserves every numerical output.
  /// kMaterialTangent is a B33-only approximate static Jacobian request:
  /// force must be null, tangent is integral B^T C B, not the consistent
  /// derivative. Other families reject it with non-OK status before modifying
  /// numerical outputs (including caller-supplied candidate stress). Derived
  /// implementations must explicitly validate recognized requests; never infer
  /// force from request!=kTangent or silently reinterpret an unknown value.
  virtual StatusCode EvaluateTrial(const ElementTrialState* state,
                                   const ElementWorkspace* workspace,
                                   ElementEvaluationRequest request,
                                   Vector* force, Matrix* tangent,
                                   TrialEvaluationFailure* failure,
                                   Diagnostic* d) const;
  /// Optional geometry-only trial validation; default false. Independent of
  /// force evaluation and result recovery, with no stored state mutation.
  virtual bool SupportsTrialGeometryCheck() const;
  /// Borrows state and immutable workspace for this call only. Derived
  /// formulation defines geometry points/criteria. Required failure separates
  /// geometry from fatal failures. Default clears it if provided and returns
  /// unsupported without modifying input, numerical output or official results.
  virtual StatusCode CheckTrialGeometry(
      const ElementTrialState* state, const ElementWorkspace* workspace,
      TrialEvaluationFailure* failure, Diagnostic* d) const;
  /// Optional physical dynamic evaluation for the given pose representation;
  /// default false. A true capability does not validate mass mode or inputs.
  virtual bool SupportsDynamicEvaluation(DynamicKinematics kinematics) const;
  /// Borrows physical motion and immutable workspace for this call only;
  /// neither is retained. Writes physical forces/partial derivatives defined
  /// in element_dynamic.h, without HHT weighting, commit or official recovery.
  /// Supported implementations require failure and validate mass mode/context.
  /// Preflight preserves output; computational failure invalidates requested
  /// blocks. Default clears failure if provided and returns unsupported,
  /// preserving all input/output buffers and official results, even for null
  /// input/output. Geometry/rotation flags do not override caller phase policy.
  virtual StatusCode EvaluateDynamic(
      const ElementMotionInput* input, MassMode mass_mode,
      const ElementWorkspace* workspace,
      DynamicEvaluationRequest request, ElementDynamicResult* out,
      TrialEvaluationFailure* failure, Diagnostic* d) const;
  /// Optional official finite-trial recovery, separate from numerical
  /// EvaluateTrial. Analysis marks common validity unconfirmed first. Borrows
  /// state/immutable workspace for this call only (Shell allows null);
  /// no ownership/commit/I/O. Default unsupported preserves results and
  /// clears failure if provided. Supported implementations require failure;
  /// preflight preserves payload, later failure may leave unconfirmed partials.
  virtual StatusCode UpdateTrialResults(const ElementTrialState* state,
                                        const ElementWorkspace* workspace,
                                        TrialEvaluationFailure* failure,
                                        Diagnostic* d);
  /// Completes this element's official result evaluation AND recovery.
  /// Caller marks AnalysisState unconfirmed first. Failure can partially change
  /// derived storage, which stays unconfirmed; no global commit occurs here.
  virtual StatusCode UpdateResults(const Vector* global_u,
                                   ElementWorkspace* workspace, Diagnostic* d);
  /// Optional element-load force evaluation. Null load delegates to the
  /// original virtual method; nonnull load is unsupported by default.
  /// Inputs are borrowed, validation preserves output; later failure discards
  /// it.
  virtual StatusCode ComputeLoadedInternalForce(const Vector* global_u,
                                                const Vector* element_load,
                                                ElementWorkspace* workspace,
                                                Vector* global_force,
                                                Diagnostic* d) const;
  /// Optional load-aware official update including recovery. Same validity and
  /// failure contract as UpdateResults. Null load delegates to UpdateResults.
  virtual StatusCode UpdateLoadedResults(const Vector* global_u,
                                         const Vector* local_element_load,
                                         ElementWorkspace* workspace,
                                         Diagnostic* d);

 protected:
  /// Derived factory validates the original ID before construction.
  explicit Element(int id);

 private:
  int id_;
};

}  // namespace fesa

#endif  // FESA_MODEL_ELEMENT_H_
