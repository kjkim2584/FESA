/// @file
/// @brief Total-Lagrangian StVK material-point response, M2/M3.
#ifndef FESA_ELEMENTS_SOLID_STVK_H_
#define FESA_ELEMENTS_SOLID_STVK_H_
#include "fesa/math/matrix.h"
#include "fesa/model/elastic_material.h"
namespace fesa {
/// Borrows finite dimensionless F3x3 and material for this call; det(F)>0.
/// Optional caller-owned outputs: W in force/length^2, PK2/P3x3 in the same
/// units, A9x9=d(vec9 P)/d(vec9 F) in row-major [11,12,13,21,...,33] order.
/// At least one output required. All matrices must be distinct from F and
/// each other. Preflight failure preserves outputs; computation failure
/// invalidates requested outputs. No state/cache is retained. Singular or
/// nonpositive determinant returns kEquilibriumFailure; backend errors remain
/// fatal. Unrequested tangent is neither allocated nor evaluated.
/// @see docs/SolidElement/nonlinear_formulation_detail_draft.md, M2/M3.
StatusCode EvaluateStvkPoint(const Matrix* deformation_gradient,
    const ElasticMaterial* material, double* energy_density,
    Matrix* second_piola, Matrix* first_piola, Matrix* first_piola_tangent,
    Diagnostic* d);
}  // namespace fesa
#endif  // FESA_ELEMENTS_SOLID_STVK_H_
