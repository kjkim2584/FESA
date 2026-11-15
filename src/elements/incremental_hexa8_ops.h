/// @file
/// @brief Coupled corrected-strain HEX8 force/exact and symmetric rate API.
#ifndef FESA_ELEMENTS_INCREMENTAL_HEXA8_OPS_H_
#define FESA_ELEMENTS_INCREMENTAL_HEXA8_OPS_H_
#include "elements/incremental_solid_ops.h"
namespace fesa {
/// H[8]3x8/q8/Un,U24/Sn8x6 borrowed for call, disjoint from outputs/scratch.
/// Same request/output/null/selector/failure contract as ordinary ops, but
/// only8/8 workspace. All averages rebuilt. Preflight preserves outputs;
/// computational failure invalidates all candidates/results and reports IP.
/// No history/result commit. Geometry uses physical F, LE uses separate Fbar.
/// @see docs/SolidElement/incremental_hexa8_contracts.md.
StatusCode EvaluateIncrementalHexa8Points(
    const Matrix* const* gradients, const Vector* weights,
    const Vector* committed_u, const Vector* candidate_u,
    const Matrix* committed_stress, const ElasticMaterial* material,
    IncrementalElasticWorkspace* workspace, ElementEvaluationRequest request,
    Vector* force, Matrix* exact_tangent, Matrix* symmetric_tangent,
    Matrix* candidate_stress, TrialEvaluationFailure* failure, Diagnostic* d,
    bool use_committed_geometric_stress = false, int* failed_point = nullptr);
}  // namespace fesa
#endif  // FESA_ELEMENTS_INCREMENTAL_HEXA8_OPS_H_
