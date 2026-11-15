/// @file
/// @brief Immutable initial mapping and adopted convergence length contracts.
#ifndef FESA_ELEMENTS_INCREMENTAL_SOLID_GEOMETRY_H_
#define FESA_ELEMENTS_INCREMENTAL_SOLID_GEOMETRY_H_
#include "fesa/elements/incremental_solid_element.h"
namespace fesa {
/// Borrows element/reference and ng owned output buffers for the call. Builds
/// H_g (3xn,1/L), q_g (L^3) at stiffness points after initial geometry guards.
/// Preflight preserves outputs; computation failure invalidates all outputs.
/// No candidate/history mutation; only straight initial Tetra10 is supported.
StatusCode BuildIncrementalInitialGeometry(const Element* element,
    const SolidReferenceData* stiffness, Matrix* const* gradients,
    Vector* weights, Diagnostic* d);
/// Initial ordered-node empirical Lc in L, per geometry_contract.md. Failure
/// preserves out; no geometry/history/result cache is read or changed.
StatusCode ComputeIncrementalCharacteristicLength(
    const IncrementalSolidElement* element, double* out, Diagnostic* d);
}  // namespace fesa
#endif  // FESA_ELEMENTS_INCREMENTAL_SOLID_GEOMETRY_H_
