/// @file
/// @brief Internal fixed-gradient total-Lagrangian integration, M1-M8.
#ifndef FESA_ELEMENTS_SOLID_NONLINEAR_OPS_H_
#define FESA_ELEMENTS_SOLID_NONLINEAR_OPS_H_

#include "fesa/elements/solid_reference_data.h"
#include "fesa/model/elastic_material.h"
#include "fesa/model/element_trial.h"

namespace fesa {
/// Borrows all inputs for this call; supports the five full Solid types.
/// C3D4/TET1, C3D10/TET4, C3D6/WEDGE2, C3D5/PYR5 use actual F;
/// C3D8/HEX8 uses F-bar.
/// Each gradient is 3xn, dN/dX in 1/length; positive weights are initial
/// integration measures. All actual det(F)>0 before any F-bar average.
/// translations is node-major [u1,u2,u3] in global length units, size 3n.
/// material is immutable. Caller owns force3n (force units), tangent3nx3n
/// (force/length), optional energy (force*length), required failure, optional
/// diagnostic.
/// Unrequested force/tangent must be null. Inputs/outputs must not overlap.
/// Preflight failure preserves numerical outputs; computation failure makes
/// requested outputs invalid. failure is cleared on entry and set to kGeometry
/// only for singular/nonpositive det(F). No state/cache survives the call.
/// Force-only does not evaluate/allocate either point or element tangent.
/// @see docs/SolidElement/nonlinear_formulation_detail_draft.md, M1-M8.
StatusCode EvaluateSolidNonlinearPoints(SolidFormulation formulation,
    const Matrix* const* gradients, const Vector* weights,
    const Vector* translations, const ElasticMaterial* material,
    ElementEvaluationRequest request, double* energy, Vector* force,
    Matrix* tangent, TrialEvaluationFailure* failure, Diagnostic* d);
/// Computes force and material-point F/Fbar together, then M9/F11 outputs.
/// Borrows the matching stiffness table and the same G/q/u/material contract as
/// EvaluateSolidNonlinearPoints. Caller owns force3n, IP ng-by-6 and node
/// n-by-6
/// tensors in [11,22,33,12,13,23] order, tensor shear, global axes. Stress has
/// force/length^2 units; logarithmic strain is dimensionless. All outputs are
/// required and mutually disjoint. Preflight preserves values; later failure
/// invalidates every output. No tangent, mass, state/commit or stored cache.
/// Ordinary TL types use actual F; C3D8 uses the same-call modified Fbar.
StatusCode RecoverSolidNonlinearPoints(const SolidReferenceData* data,
    const Matrix* const* gradients, const Vector* weights,
    const Vector* translations, const ElasticMaterial* material, Vector* force,
    Matrix* ip_stress, Matrix* ip_strain, Matrix* nodal_stress,
    Matrix* nodal_strain, TrialEvaluationFailure* failure, Diagnostic* d);
}  // namespace fesa
#endif  // FESA_ELEMENTS_SOLID_NONLINEAR_OPS_H_
