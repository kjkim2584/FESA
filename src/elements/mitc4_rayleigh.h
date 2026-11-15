/// @file
/// @brief Current-plane MITC4 stiffness-proportional damping, RD1--RD9.
#ifndef FESA_ELEMENTS_MITC4_RAYLEIGH_H_
#define FESA_ELEMENTS_MITC4_RAYLEIGH_H_

#include "fesa/model/element_dynamic.h"
#include "fesa/model/node.h"
#include "fesa/model/shell_section.h"

namespace fesa {
/// Computes a x b in the inputs' common Cartesian basis, without allocation.
/// Borrows nonnull length-3 inputs/output and exclusive 3x3 scratch for this
/// call; never retains pointers. Caller owns scratch and must finish using its
/// previous value before borrowing it. Output must differ from both inputs.
/// Overwrites all nine scratch coefficients before the external matrix-vector
/// multiply, so a previous nonfinite scratch value is allowed. Inputs stay
/// unchanged. Null/shape/alias/nonfinite-input preflight preserves scratch and
/// output; calculation failure invalidates scratch and output. No state is
/// retained between calls or shared across workers.
/// @see docs/HHT/mitc4_rayleigh.md, approved df/dp lifetime boundaries.
StatusCode CrossMitc4RayleighVectors(Matrix* scratch, const Vector* a,
                                    const Vector* b, Vector* out,
                                    Diagnostic* d);

/// Borrows four validated nodes/section and finite E0 motion through this call.
/// Overwrites global node-major force24 [F,FL] and, for the linearization
/// request, GD/CD24x24 in spatial pose/velocity coordinates. Unrequested
/// matrices must be null. This is only the beta_R part; excludes transverse
/// shear, drilling, mass drag, HHT weighting and official stress/results.
/// Caller first applies ValidateElementDynamicEvaluation and mass capability.
/// Required failure is cleared; invalid initial geometry is fatal, current
/// nonpositive/collapsed volume is typed geometry with the Gauss location.
/// Preflight preserves outputs; calculation failure invalidates all requested
/// outputs. Inputs remain unchanged; scratch/factors are local, never retained.
/// @see docs/HHT/mitc4_rayleigh.md, RD1--RD9.
StatusCode EvaluateMitc4RayleighBeta(const Node* const* nodes,
                                     const ShellSection* section,
                                     const ElementMotionInput* input,
                                     DynamicEvaluationRequest request,
                                     Vector* force, Matrix* gd, Matrix* cd,
                                     TrialEvaluationFailure* failure,
                                     Diagnostic* d);
}  // namespace fesa
#endif  // FESA_ELEMENTS_MITC4_RAYLEIGH_H_
