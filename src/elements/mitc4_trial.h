/// @file
/// @brief Recomputed nonlinear MITC4 resisting force and spatial Jacobian.
#ifndef FESA_ELEMENTS_MITC4_TRIAL_H_
#define FESA_ELEMENTS_MITC4_TRIAL_H_
#include "fesa/model/element_trial.h"
#include "fesa/model/node.h"
#include "fesa/model/shell_section.h"
namespace fesa {
/// Four distinct sixDOF nodes and homogeneous section borrowed. State
/// dimensions are translations/body_spin12 and orientations/committed4×4; q
/// unit wxyz. fixed_drilling>=0 is frozen initial kd (zero allowed for isolated
/// kernel verification). No initial K evaluation, mass or official recovery.
/// Call-local basis scratch and eight raw tying strain/gradient entries are
/// owned only until return/failure, never reused across states or workers.
/// Caller maintains body_spin consistency from accepted history; routine
/// checks relative chart but cannot reconstruct a_c from these inputs alone.
/// Requested force24 is global node-major[Fx,Fy,Fz,Mx,My,Mz]; tangent24×24 is
/// its derivative under spatial Exp(dw) updates and the same body-spin rule.
/// Tangent generally nonsymmetric; never symmetrize. Unrequested pointers null.
/// Outputs disjoint from inputs; readonly inputs may alias each other.
/// Failure flag required, cleared on entry.
/// Preflight preserves outputs. Later failure invalidates requested outputs.
/// Initial invalid geometry is fatal; only valid-initial trial geometry/chart
/// failure returns kEquilibriumFailure plus typed failure. Other errors fatal.
/// No inputs/definition/official results retained or changed; d optional.
/// Geometry failure records its actual Gauss/tying natural location in d.
StatusCode EvaluateMitc4Trial(const Node* const* nodes, int count,
                              const ShellSection* section,
                              const ElementTrialState* state,
                              ElementEvaluationRequest request,
                              double fixed_drilling, Vector* force,
                              Matrix* tangent, TrialEvaluationFailure* failure,
                              Diagnostic* d);
}  // namespace fesa
#endif  // FESA_ELEMENTS_MITC4_TRIAL_H_
