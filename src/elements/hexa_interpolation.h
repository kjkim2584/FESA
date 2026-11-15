/// @file
/// @brief Natural Hexa8 scalar interpolation coefficients from formulation F3.
#ifndef FESA_ELEMENTS_HEXA_INTERPOLATION_H_
#define FESA_ELEMENTS_HEXA_INTERPOLATION_H_
#include "fesa/math/matrix.h"
namespace fesa {
/// Fills shape1x8 and natural gradient3x8 in the natural cube [-1,1]^3.
/// Borrows all objects for this call, changes only outputs. Invalid arguments
/// preserve outputs; later backend failure invalidates them. No geometry or
/// material evaluation. Node order is the eight F3/Abaqus ordered signs.
StatusCode FillHexaInterpolation(double xi, double eta, double zeta,
                                 Matrix* shape, Matrix* gradient,
                                 Diagnostic* d);
}  // namespace fesa
#endif  // FESA_ELEMENTS_HEXA_INTERPOLATION_H_
