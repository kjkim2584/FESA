/// @file
/// @brief Initial flat MITC4 geometry validation, without persistent data.
#ifndef FESA_ELEMENTS_MITC4_GEOMETRY_H_
#define FESA_ELEMENTS_MITC4_GEOMETRY_H_
#include "fesa/math/matrix.h"
#include "fesa/model/node.h"
namespace fesa {
/// Borrows four distinct six-component nodes. Writes initial 4x3 local
/// coordinates relative to X1 in model length units and 3x3 global-to-local
/// row basis [a,b,n]: a=(X2-X1)/norm, n=normalized edge12 cross edge14,
/// b=n cross a. diameter is maximum pair distance, length mean edge.
/// All outputs required; matrices must be distinct and correctly sized.
/// diameter and length must not alias. d is optional and borrowed.
/// Validation/invalid geometry preserves outputs; backend failure invalidates
/// matrices but preserves scalar outputs. No input mutation or retained data.
/// See docs/MITC4/geometry_and_convergence_contract.md for scaled guards.
StatusCode ValidateMitc4Geometry(const Node* const* nodes, int count,
                                 Matrix* local_coordinates, Matrix* frame,
                                 double* diameter, double* length,
                                 Diagnostic* d);
}  // namespace fesa
#endif  // FESA_ELEMENTS_MITC4_GEOMETRY_H_
