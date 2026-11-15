/// @file
/// @brief IE1/IE5 ordinary full-solid IP resisting force and rate correction.
#ifndef FESA_ELEMENTS_INCREMENTAL_SOLID_OPS_H_
#define FESA_ELEMENTS_INCREMENTAL_SOLID_OPS_H_
#include "elements/incremental_elastic_point.h"
namespace fesa {
/// Borrows ng initial gradients3×n in1/L, positive weights ng entries in L³,
/// committed/candidate U3n/global L and committed stress ng×6/global tensor
/// shear F/L².
/// n/ng must match the worker:4/1,10/4,6/2 or5/5. Hexa8 requires coupled ops.
/// Candidate stress ng×6 is required and separate from committed storage.
/// Requested force3n/F and rate tangent3n×3n/F/L are caller-owned;
/// unrequested pointers must be null. Optional exact tangent3n×3n is computed
/// only with a tangent request. All outputs and inputs must be disjoint.
/// Preflight preserves outputs; later failure invalidates all requested outputs.
/// Worker scratch owns temporaries; no official result/history is changed.
/// Internal point geometry failures are labeled for the Element adapter;
/// committed/backend faults stay fatal. @see IE1/IE4/IE5 formulation.
/// Original point status is retained (candidate geometry is kSingular plus
/// its typed failure). Optional failed_point is reset to-1, then receives the
/// failing0-based IP when known. All ng candidate F buffers are recomputed on
/// success; any failure invalidates all of them and all requested outputs.
/// IE5-R correction uses candidate geometry. The optional selector uses only
/// committed Cauchy stress in geometric-stiffness coefficients; default false
/// uses candidate stress. It never changes f, candidate stress, or exact K.
StatusCode EvaluateIncrementalSolidPoints(const Matrix* const* gradients,
    const Vector* weights, const Vector* committed_u, const Vector* candidate_u,
    const Matrix* committed_stress, const ElasticMaterial* material,
    IncrementalElasticWorkspace* workspace, ElementEvaluationRequest request,
    Vector* force, Matrix* exact_tangent, Matrix* symmetric_tangent,
    Matrix* candidate_stress, TrialEvaluationFailure* failure, Diagnostic* d,
    bool use_committed_geometric_stress = false, int* failed_point = nullptr);
}  // namespace fesa
#endif  // FESA_ELEMENTS_INCREMENTAL_SOLID_OPS_H_
