/// @file
/// @brief Natural Wedge6 scalar interpolation from formulation F4.
#ifndef FESA_ELEMENTS_WEDGE_INTERPOLATION_H_
#define FESA_ELEMENTS_WEDGE_INTERPOLATION_H_
#include "fesa/math/matrix.h"
namespace fesa {
/// Fills borrowed shape1x6 and gradient3x6 in xi,eta>=0,xi+eta<=1,|zeta|<=1.
/// Node order is the lower triangle1..3 followed by upper4..6. Invalid inputs
/// preserve outputs; a later backend failure invalidates them. No geometry,
/// material, or stored-result mutation.
StatusCode FillWedgeInterpolation(double xi, double eta, double zeta,
                                  Matrix* shape, Matrix* gradient,
                                  Diagnostic* d);
}  // namespace fesa
#endif  // FESA_ELEMENTS_WEDGE_INTERPOLATION_H_
