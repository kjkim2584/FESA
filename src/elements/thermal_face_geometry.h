/// @file
/// @brief Shared planar thermal quad mapping and nodal face-area weights.
#ifndef FESA_SRC_ELEMENTS_THERMAL_FACE_GEOMETRY_H_
#define FESA_SRC_ELEMENTS_THERMAL_FACE_GEOMETRY_H_

#include "fesa/math/matrix.h"
#include "fesa/model/node.h"

namespace fesa {
class ThermalElement;
namespace thermal_face_geometry {
/// Borrows four valid nodes. Caller owns returned normalized4x2 coordinates;
/// *scale is the physical diameter [L]. Same DS4/formulation3.3 normalized
/// height64epsilon guard. Failure nulls *out; no input mutation or retention.
StatusCode CreateCoordinates(const Node* const* nodes, Matrix** out,
                             double* scale, Diagnostic* d);
/// Maps one 2x2 Gauss station of normalized coordinates. Optional gradient
/// receives dN/dxhat, optional shape receives four N. Positive detJ is
/// dimensionless. Failure may invalidate outputs; inputs remain unchanged.
StatusCode MapPoint(const Matrix* xy, int point, Matrix* gradient,
                    double* shape, double* determinant, Diagnostic* d);
/// Overwrites full element DOF area weights[L²], zero off this face. Borrows
/// valid local node indices(count3/4), using this section_point at each node.
/// Triangles use A/3; planar quads use2x2 Gauss row sums(formulation8.4).
/// No retained pointers/model mutation; failure may invalidate caller's out.
StatusCode AssignWeights(const ThermalElement* element, const int* local_nodes,
                         int count, int section_point, Vector* out,
                         Diagnostic* d);
}  // namespace thermal_face_geometry
}  // namespace fesa
#endif  // FESA_SRC_ELEMENTS_THERMAL_FACE_GEOMETRY_H_
