/// @file
/// @brief Rational Pyramid5 natural interpolation from formulation F5.
#ifndef FESA_ELEMENTS_PYRAMID_INTERPOLATION_H_
#define FESA_ELEMENTS_PYRAMID_INTERPOLATION_H_
#include "fesa/math/matrix.h"
namespace fesa {
/// Fills borrowed shape1x5 and optional gradient3x5 in the natural pyramid.
/// Base order is (--),(+-),(++),(-+), followed by the apex. At (0,0,1)
/// only shape with null gradient is defined: [0,0,0,0,1]. Requesting its
/// path-dependent gradient is an error, with no epsilon substitution.
/// Validation preserves outputs; backend failure invalidates them. Does not
/// evaluate geometry, material, or modify official results.
StatusCode FillPyramidInterpolation(double xi, double eta, double zeta,
                                    Matrix* shape, Matrix* gradient,
                                    Diagnostic* d);
}  // namespace fesa
#endif  // FESA_ELEMENTS_PYRAMID_INTERPOLATION_H_
