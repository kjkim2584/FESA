/// @file
/// @brief Internal scalar coefficient placement for Tetra4/10 interpolation.
#ifndef FESA_ELEMENTS_TETRA_INTERPOLATION_H_
#define FESA_ELEMENTS_TETRA_INTERPOLATION_H_
#include "fesa/math/matrix.h"
namespace fesa {
/// Fills natural shape row 1 by nodes and gradients 3 by nodes at finite
/// xi/eta/zeta, using Abaqus ordered corners/edges. nodes must be 4 or 10.
/// Inputs/outputs are borrowed for the call. Invalid arguments preserve
/// outputs; later backend failure invalidates them. No geometry/material/result
/// mutation.
StatusCode FillTetraInterpolation(int nodes, double xi, double eta, double zeta,
                                  Matrix* shape, Matrix* gradient,
                                  Diagnostic* d);
}  // namespace fesa
#endif  // FESA_ELEMENTS_TETRA_INTERPOLATION_H_
