/// @file
/// @brief Q23/Q53 current mapping checks on immutable ordered sample tables.
#ifndef FESA_ELEMENTS_SOLID_TRIAL_GEOMETRY_H_
#define FESA_ELEMENTS_SOLID_TRIAL_GEOMETRY_H_
#include "fesa/elements/solid_element.h"
namespace fesa {
class SolidElementReadView;
/// Borrows element, total translations and exact-key TrialGeometry catalog.
/// Uses original X+u at all nodes, including Tetra10 midsides, and current h.
/// No force/tangent/material/results/state mutation. First failed sample is
/// reported with normalized mapping determinant when actually computed.
/// Singular LU/collapsed h is geometry without an invented measured value;
/// allocation/backend/range errors remain fatal. failure is required and reset.
/// Success means only these samples passed, not global noninversion proof.
StatusCode CheckSolidTrialGeometry(const SolidElement* element,
    const ElementTrialState* state, const ElementWorkspace* workspace,
    TrialEvaluationFailure* failure, Diagnostic* d);
/// Same contract for either Solid family through a borrowed call-lifetime view.
/// Includes every physical midside displacement; never applies initial-only
/// straight-edge constraints to a candidate.
StatusCode CheckSolidTrialGeometry(const SolidElementReadView* element,
    const ElementTrialState* state, const ElementWorkspace* workspace,
    TrialEvaluationFailure* failure, Diagnostic* d);
}  // namespace fesa
#endif  // FESA_ELEMENTS_SOLID_TRIAL_GEOMETRY_H_
