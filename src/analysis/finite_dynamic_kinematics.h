/// @file
/// @brief Call-local bridge from physical six-DOF state to R1/R2 integration.
#ifndef FESA_ANALYSIS_FINITE_DYNAMIC_KINEMATICS_H_
#define FESA_ANALYSIS_FINITE_DYNAMIC_KINEMATICS_H_
#include "fesa/analysis/analysis_state.h"
#include "fesa/analysis/dof_manager.h"
#include "fesa/analysis/time_integrator.h"
namespace fesa {
/// Borrows Domain-layout6N mapping/integrator and prepared finite state with
/// valid committed motion. dt/step_time/total_time are physical T. predict
/// writes predicted u/q through state setters (including body-spin); false
/// preserves the caller-corrected pose and recomputes v/omega/a/eta. Stages
/// global6N motion only after both R1 roles succeed. For a matching exclusive
/// Beam mapping/state, predicts eta [L] before physical staging and recomputes
/// eta velocity/acceleration [L/T,L/T^2] after it via the same integrator.
/// Requires explicitly committed eta motion; never infers zero motion.
/// predict=false preserves candidate physical pose and eta. No outputUR
/// reconstruction, prescriptions, history, recovery, convergence or commit.
/// Required branch cleared/set for existing chart failure. Every failure
/// discards candidate, preserving committed state; all numerical snapshots die
/// on return.
/// @see docs/HHT/rotational_integrator.md AnalysisState connection.
StatusCode StageFiniteHhtKinematics(const DofManager* dofs,
                                    const TimeIntegrator* integrator, double dt,
                                    double step_time, double total_time,
                                    bool predict, AnalysisState* state,
                                    bool* branch, Diagnostic* d);
/// Read-only R2 bridge using matching prepared/staged current and committed
/// physical motion. Writes caller-owned distinct n-by-9 angular H/B and
/// Cartesian c1/c0 for dt, also applying to Beam eta. Beam role requires
/// matching IDs/offsets and valid joint candidate eta velocity/acceleration.
/// All pointers borrowed for call; no retained scratch.
/// Required branch follows R2. Preflight preserves outputs; later failure
/// invalidates angular blocks but preserves scalar outputs and all state.
StatusCode EvaluateFiniteHhtKinematicLinearization(
    const DofManager* dofs, const TimeIntegrator* integrator, double dt,
    const AnalysisState* state, TranslationalLinearization* translation,
    RotationalLinearization* rotation, bool* branch, Diagnostic* d);
}  // namespace fesa
#endif  // FESA_ANALYSIS_FINITE_DYNAMIC_KINEMATICS_H_
