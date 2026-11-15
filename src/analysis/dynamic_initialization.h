/// @file
/// @brief Physical initialization for supported Solid/Shell dynamics.
#ifndef FESA_ANALYSIS_DYNAMIC_INITIALIZATION_H_
#define FESA_ANALYSIS_DYNAMIC_INITIALIZATION_H_

#include "fesa/analysis/analysis_model.h"
#include "fesa/analysis/analysis_state.h"
#include "fesa/analysis/dof_manager.h"
#include "fesa/analysis/linear_solver.h"
#include "fesa/elements/element_workspace.h"

namespace fesa {
/// Read-only same-Domain sequential-step constraint identity/INPUT-value
/// check, shared by linear/finite dynamics. Null previous requires step1.
/// No comparison with output principal UR, physical mutation or retained data.
/// Requires unchanged node/DOF keys and command values. Failure reports the
/// existing unsupported transition; inputs remain unchanged.
StatusCode ValidateDynamicConstraintTransition(
    const EffectiveConditions* previous, const EffectiveConditions* current,
    Diagnostic* d);
/// N4 acceleration solve for constrained linear positive-density S4.
/// Caller first verifies P14 initial configuration/velocity compatibility.
/// Borrows global physical net_force [F,FL] and derived null_acceleration
/// [rad/T^2], ordered by LinearShellNullspace. The latter is required exactly
/// when kernel Count()>0; it is not an arbitrary zero-acceleration convention.
/// Solves the sparse augmented mass system described in initialization.md.
/// worker_count>0; S4 uses no element workspace. No state/results change.
/// Writes global [L/T^2,rad/T^2], fixed/excluded components zero. Inputs and
/// output must be distinct. Failure preserves output until final CopyFrom,
/// whose backend failure invalidates output. No mass/solver fallback or retry.
StatusCode ComputeInitialShellAcceleration(
    const AnalysisModel* model, const DofManager* dofs, MassMode mode,
    const Vector* net_force, const Vector* null_acceleration,
    LinearSolver* solver, int worker_count, Vector* acceleration, Diagnostic* d);
/// N6/N7/P17 and N8 physical initialization for finite-rotation S4 with Rayleigh.
/// Borrows node-ID-ordered pose (3N translations/body_spin, N-by-4 unit wxyz
/// orientations/committed_orientations), global6N velocity [v,omega], physical
/// external [F,FL] and right external_rate [F/T,FL/T]. Constant prescribed
/// motion requires zero fixed/excluded velocity; all3 rotations free or fixed.
/// Constant/right TABULAR segments have Mddot=0. Validates P17 before solving
/// current MI with a zero-null range constraint, then derives the physical
/// null acceleration from N8; no artificial mass or input projection.
/// The physical range RHS is p-fint-gyro-fD at this same pose/velocity.
/// Solver is borrowed through the call; exclusive scratch is call-local and
/// released after joined dispatch (worker_count>0).
/// Output is distinct global6N [a,eta], fixed/excluded zero. Preserves output
/// on failure until final CopyFrom (whose backend failure invalidates it).
/// Initial geometry/recovery/frame checks remain the caller responsibility;
/// this helper alone does not activate nonlinear step execution. No state,
/// official results or retained scratch are changed.
StatusCode ComputeInitialFiniteShellAcceleration(
    const AnalysisModel* model, const DofManager* dofs,
    const ElementTrialState* pose, const Vector* velocity, MassMode mode,
    const Vector* external, const Vector* external_rate, LinearSolver* solver,
    int worker_count, Vector* acceleration, Diagnostic* d);
/// Caller selects the approved physical start or unchanged-load continuation.
enum class SolidDynamicStart { kDirect, kStaticPreload, kContinuous, kLoadJump };
/// Linear S4 counterpart of PrepareInitialSolidDynamicState, with Rayleigh.
/// Same start, mapping, prescribed-motion and transactional failure contract;
/// 6N fields contain [translation,rotation] in global axes. external/external_rate
/// are physical [F,FL]/[F/T,FL/T], required for physical starts. Supported CLOAD
/// is constant or right-piecewise-linear TABULAR: its right second derivative
/// is zero, so N4 derives zero null acceleration when Kzr=0. Checks P14 before
/// mass solve without changing supplied/default velocity or preload. Continuous
/// start ignores external/rate, preserves discrete motion and skips P14/solve.
/// No worker contexts, result recovery, output, commit or retained scratch.
StatusCode PrepareInitialShellDynamicState(
    SolidDynamicStart start, const AnalysisModel* model, const DofManager* dofs,
    const EffectiveConditions* conditions, MassMode mode,
    const Vector* external, const Vector* external_rate,
    const Vector* initial_velocity, LinearSolver* solver, int worker_count,
    double total_time, AnalysisState* state, Diagnostic* d);
/// Prepares initial u/v/a (step_time=0), without recovery, output or commit.
/// Borrows model/mapping/conditions/solver/workers through the call. state owns
/// all staged copies. Requires an unprepared linear 3N state and the next step.
/// Direct start copies optional global initial_velocity [L/T], otherwise zero;
/// static preload preserves u and starts v=0; continuous preserves u/v/a and
/// performs no mass solve. Non-direct starts require null initial_velocity and
/// unchanged prescribed displacements. Caller establishes continuous loads.
/// Load jump preserves committed u/v and recomputes a from the right-side N1
/// physical force; it requires committed motion at the same total time.
/// external [F] is evaluated at the exact new start time, required for physical
/// starts and ignored for continuous. Initial structural force is recomputed
/// from the prepared u; N1 uses p-f_int-C*v, without HHT weighting or convergence
/// flux. Fixed v/a are zero; nonzero supplied fixed velocity is an input error.
/// Any failure discards candidates and preserves committed data. Activation
/// may leave allocated but uncommitted dynamic storage. Node/element official
/// results and tags are unchanged. Initial frame still requires recovery.
StatusCode PrepareInitialSolidDynamicState(
    SolidDynamicStart start, const AnalysisModel* model, const DofManager* dofs,
    const EffectiveConditions* conditions, MassMode mode,
    const Vector* external, const Vector* initial_velocity,
    LinearSolver* solver, ElementWorkspace** workers, int worker_count,
    double total_time, AnalysisState* state, Diagnostic* d);
/// Solves N1 in docs/HHT/initialization.md for supported linear 3N elements.
/// Borrows all inputs for the call; workers are exclusive through joined
/// dispatch. net_force is physical p-f_int-C*v at the initial state/time, global
/// node-ID order [F1,F2,F3], without HHT weighting. Fixed prescribed motion
/// must be time-independent: its acceleration is zero. Excluded components
/// also remain zero. Caller handles initial velocity, recovery and output.
/// Writes acceleration [L/T^2] in the same layout, distinct from net_force.
/// Uses bounded sparse mass assembly and the supplied solver (SPD); no
/// regularization, mass substitution or retry. Zero free DOFs skip assembly
/// and solve. Failure preserves output until final CopyFrom (whose backend
/// failure invalidates output), never changes model/results/analysis state.
StatusCode ComputeInitialSolidAcceleration(
    const AnalysisModel* model, const DofManager* dofs, MassMode mode,
    const Vector* net_force, LinearSolver* solver, ElementWorkspace** workers,
    int worker_count, Vector* acceleration, Diagnostic* d);
}  // namespace fesa

#endif  // FESA_ANALYSIS_DYNAMIC_INITIALIZATION_H_
