/// @file
/// @brief Physical nonlinear assembly with explicit Beam internal equations.
#ifndef FESA_ANALYSIS_NONLINEAR_ASSEMBLER_H_
#define FESA_ANALYSIS_NONLINEAR_ASSEMBLER_H_
#include "fesa/analysis/analysis_model.h"
#include "fesa/analysis/analysis_state.h"
#include "fesa/analysis/dof_manager.h"
#include "fesa/analysis/nonlinear_problem.h"
#include "fesa/analysis/time_integrator.h"
#include "fesa/io/control_settings.h"
#include "fesa/model/element_trial.h"
namespace fesa {
/// Test-only read-only observation of the actual B33 joined tangent request.
/// Registers a nonthrowing callback/context on this coordinator thread until
/// reset with null. Called once after preflight, before dispatch; no worker
/// invocation or numerical buffer is exposed. The caller owns context and must
/// clear the registration before destroying it. May not mutate solver state.
void SetBeamTangentObserverForTest(void (*callback)(ElementEvaluationRequest,
                                                    void*),
                                   void* context);
/// Read-only preflight for the explicit nonlinear Beam mapping. Checks actual
/// Beam12+2 connectivity, matching internal equation offsets, and full-general
/// CSR/COO bounds before output creation. No trial state, numerical allocation
/// or element evaluation. Model/mapping borrowed and unchanged on all outcomes.
/// Ordinary mappings are invalid for this role; use ValidateAssemblyCounts.
StatusCode ValidateNonlinearBeamAssemblyCounts(const AnalysisModel* model,
                                              const DofManager* dofs,
                                              Diagnostic* d);
/// All objects borrowed for one residual evaluation. Definition/view/mapping
/// remain fixed; the finite state is read only. Units are consistent F/L/rad.
struct NonlinearResidualInput {
  const AnalysisModel* model = nullptr;
  const DofManager* dofs = nullptr;
  const AnalysisState* state = nullptr;
  const ControlSettings* controls = nullptr;
  /// Borrowed worker-indexed contexts, live through join; null/0 for Shell.
  const ElementWorkspace* const* workspaces = nullptr;
  int workspace_count = 0;
  const Vector* external =
      nullptr;  ///< Physical global force/moment at current target.
  const Vector* force_origins =
      nullptr;  ///< Individual current F; static requires nonzero, or null.
  const Vector* moment_origins =
      nullptr;  ///< Individual current F*L; static requires nonzero, or null.
  /// Physical global[du,dw], null at i0.
  const Vector* actual_correction = nullptr;
  double previous_force_max_mean =
      0;  ///< THIS step prior-success qmax average.
  double previous_moment_max_mean = 0;
  /// Internal eta correction [L], null at the initial iteration.
  const Vector* internal_actual_correction = nullptr;
};
/// Caller-owned global outputs and borrowed Newton field summaries. No output
/// aliases other outputs or numerical inputs. reaction may be the state's
/// candidate reaction buffer; physical pose and committed state never change.
struct NonlinearResidualOutput {
  /// Global residual including prescribed entries; evaluator defines weights.
  Vector* residual = nullptr;
  Vector* reaction = nullptr;  ///< Current physical reaction on constrained
                               ///< DOFs, zero elsewhere.
  FieldResidual* force = nullptr;
  FieldResidual* moment = nullptr;
  double characteristic_length = 0;  ///< Mean initial element Lc.
  /// Required 2E internal residual for the explicit nonlinear Beam role.
  Vector* internal_residual = nullptr;
};
/// Common physical structural-force evaluation, with a prepared candidate
/// and matching bound model/state/mapping. Borrows all inputs through joined
/// dispatch; Shell workspaces null/0, Solid workspaces as for residual below.
/// committed selects the copied committed pose/history rather than candidate.
/// Creates global f_int [F,FL] and unassembled packed element force/moment
/// [original element ID, connectivity, actual component] order, including
/// stabilization; caller deletes both. Empty elements give null packed/zero
/// global/Lc0. Initial characteristic_length is the element Lc arithmetic mean.
/// Distinct creation output locations required; failure nulls both and
/// preserves the scalar length. Required failure is cleared/typed as below.
/// A committed-source failure identifies that fixed source; retry eligibility
/// remains with Analysis and is not implied by the typed element cause.
/// No HHT weighting, inertia, damping, field/time statistics, official results
/// or physical state mutation; current/committed input state remains unchanged.
/// Release copied evaluation data after its last consumer; never cache across
/// iterations/increments.
/// The nonlinear Beam role requires global_internal_force and returns separate
/// internal2E [F], physical6N [F,FL], and packed14E outputs. No internal global
/// node index is invented. All creation output locations must differ; failure
/// nulls all outputs. Other roles may supply the optional output and get null.
StatusCode EvaluateNonlinearStructuralForce(
    const AnalysisModel* model, const DofManager* dofs,
    const AnalysisState* state, bool committed,
    const ElementWorkspace* const* workspaces, int workspace_count,
    Vector** global_force, Vector** element_force,
    double* characteristic_length, NonlinearFailure* failure, Diagnostic* d,
    Vector** global_internal_force = nullptr);
/// Static r=p-fint and reaction=fint-p on constrained DOFs (zero elsewhere).
/// Beam additionally requires internal_residual2E=-fint_eta, with no eta load.
/// Internal residual/motion join the force field; moments remain physical.
/// Evaluates element force only, unassembled field flux with origin
/// multiplicity, free residual norms/ID maxima and all-participating physical
/// increment maxima. No K/M/official recovery/convergence/time-history commit
/// or field conversion. force/moment spatial_mean and maximum_flux are filled;
/// raw/used means and branch decisions remain unset for Analysis's time preview
/// then Newton. Nodes have three or six physical components; elements support
/// trial evaluation,
/// mapping/state sizes consistent and finite candidate prepared. Origin entries
/// must be nonzero. Element tasks use joined OpenMP with MKLlocal1 and disjoint
/// per-call output; no cross-call numerical cache/workspace. Deterministic
/// original element order assembly. failure required/cleared; typed
/// geometry/chart maps to recoverable status, resources/library/unsupported
/// fatal. Preflight preserves outputs; later failure discards all outputs.
/// Caller discards state candidate on error.
/// @see docs/MITC4/nonlinear_assembly_contract.md.
StatusCode EvaluateNonlinearResidual(const NonlinearResidualInput* input,
                                     NonlinearResidualOutput* output,
                                     NonlinearFailure* failure, Diagnostic* d);
/// Finite MITC4/Beam HHT residual with Rayleigh, alpha in [-1/3,0]. Borrows
/// prepared mechanical inputs plus previous physical global external load.
/// Requires finite dynamics with staged current motion; nonzero alpha also
/// requires committed motion. alpha0 permits physical initial-frame evaluation
/// before the first dynamic commit, without replacing physical initialization.
/// Computes r=(1+alpha)(p-fint-fD)-alpha(p_old-fint_old-fD_old)-FI_current,
/// including current gyro. Reaction estimate is current fint+FI+fD-p on
/// constrained DOFs only, using algorithmic inertia per P18. R11 origins use
/// unweighted current element [fint,FI] plus an independent fD block only for
/// nonzero material Rayleigh, then typed current external origins. Retains
/// zero external origins as in linear Shell
/// dynamics. No old inertia, extra gyro origin, history update or convergence
/// decision. Same output/alias/failure ownership contract as static evaluator.
/// Also forbids alias with previous external and all dynamic motion stores.
/// Evaluates current structure first; alpha0 does not evaluate old structure.
/// Starting damping is reevaluated at the committed pose/V. Does not mutate
/// physical state/history/official results.
/// @see docs/HHT/implementation_contract.md.
StatusCode EvaluateFiniteHhtResidual(const NonlinearResidualInput* input,
                                     const Vector* previous_external,
                                     double alpha, MassMode mode,
                                     NonlinearResidualOutput* output,
                                     NonlinearFailure* failure, Diagnostic* d);
/// Call-scoped H6 coefficients for the same current trial and dt. Borrows
/// positive Cartesian c1/c0 and physical-node-order n-by-9 angular H/B blocks.
/// No ownership, state mutation or retained numerical data.
struct FiniteHhtLinearization {
  double alpha = 0;  ///< HHT alpha in [-1/3,0].
  MassMode mode = MassMode::kConsistent;  ///< Approved physical mass mode.
  const TranslationalLinearization* translation = nullptr;
  const RotationalLinearization* rotation = nullptr;
};
/// Uses common sparse assembly with finite Kt/physical partials/R2 to
/// create full-general H6 Kff*delta=r_free. Borrows matching prepared finite
/// pose/staged motion, physical residual and global-node linearization through
/// joined execution. Beam requires matching internal_residual2E and14 blocks.
/// No old force/physical state/history/
/// recovery/solve. Fixed corrections are0. Caller deletes successful *out;
/// EquationCount0 gives null success. All failures null output, clean scratch and
/// preserve inputs. Typed failure follows the shared joined selector.
/// @see docs/HHT/implementation_contract.md.
/// Optional worker contexts are borrowed until join and overwritten as scratch.
StatusCode AssembleFiniteHhtTangent(
    const AnalysisModel* model, const DofManager* dofs,
    const AnalysisState* state, const Vector* current_residual,
    const FiniteHhtLinearization* linearization, EquationSystem** out,
    NonlinearFailure* failure, Diagnostic* d,
    const Vector* internal_residual = nullptr,
    const ElementWorkspace* const* workspaces = nullptr,
    int workspace_count = 0);
/// Recomputes tangent only at current pose and creates caller-typed CSR reduced
/// Kff*delta=r_free. Prescribed corrections are0: no Kfc*delta_uc term, since
/// trial already meets target. Zero free count succeeds with null system.
/// Analysis owns
/// successful*out; failure nulls it/cleans partial allocations, inputs
/// preserved. No density/mass, final element force, official recovery or
/// physical mutation. Borrowed current_residual is physical global from the
/// same current trial.
/// Optional call-local IE5-R selector changes incremental-solid geometric
/// stress coefficients only. Copied through joined worker tasks; default false
/// uses candidate stress. No selector/history state survives the call.
/// Nonlinear Beam requires internal_residual2E from the same trial and general
/// CSR. Its system order is EquationCount=FreeCount+2E, including when all
/// physical DOFs are prescribed. No internal static condensation is performed.
/// Optional material-only selector requires a pure B33 model. Uses B^T C B
/// for this joined call only, with an unchanged residual. HHT never selects it.
StatusCode AssembleNonlinearTangent(
    const AnalysisModel* model, const DofManager* dofs,
    const AnalysisState* state, const Vector* current_residual,
    const ElementWorkspace* const* workspaces, int workspace_count,
    MatrixType type, EquationSystem** out, NonlinearFailure* failure,
    Diagnostic* d, bool use_committed_geometric_stress = false,
    const Vector* internal_residual = nullptr,
    bool beam_material_tangent = false);
/// Q21/Q37 read-only sampling of a valid failed candidate before discard.
/// Skips elements without the optional geometry capability; a claimed
/// capability returning unsupported is fatal. Uses joined workers and the same
/// physical trial gather as residual/tangent. No force, tangent, recovery,
/// output projection, or state mutation. Caller owns the separate inspection
/// diagnostic and retains the original Newton failure. Sample success is not
/// convergence. Solid workspaces must contain immutable kTrialGeometry
/// reference data and remain alive through the synchronous join. Failure
/// clears/sets typed output; candidate and committed state remain unchanged
/// even on failure.
StatusCode CheckNonlinearTrialGeometry(
    const AnalysisModel* model, const DofManager* dofs,
    const AnalysisState* state, const ElementWorkspace* const* workspaces,
    int workspace_count, NonlinearFailure* failure, Diagnostic* d);
/// Recovers official Node and optional element trial results from the same
/// finite candidate. Requires projected outputUR and unconfirmed result tags.
/// Uses joined element dispatch/MKLlocal1 and the same trial gathering as f/K;
/// does not assemble K/M/residual, write/commit or change physical state.
/// Mutates only Node/Element result payloads. Any failure invalidates the whole
/// unconfirmed payload set; no rollback. Caller discards candidate on error.
StatusCode RecoverNonlinearResults(const AnalysisModel* model,
                                   const DofManager* dofs,
                                   const AnalysisState* state,
                                   const ElementWorkspace* const* workspaces,
                                   int workspace_count,
                                   NonlinearFailure* failure, Diagnostic* d);
}  // namespace fesa
#endif  // FESA_ANALYSIS_NONLINEAR_ASSEMBLER_H_
