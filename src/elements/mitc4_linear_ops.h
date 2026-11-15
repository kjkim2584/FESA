/// @file
/// @brief Call-local linear MITC4 operators shared by force and recovery.
#ifndef FESA_ELEMENTS_MITC4_LINEAR_OPS_H_
#define FESA_ELEMENTS_MITC4_LINEAR_OPS_H_
#include "fesa/math/matrix.h"
#include "fesa/model/elastic_material.h"
namespace fesa {
/// Initial local coordinates4x3 -> four natural tying rows4x24.
/// Order r(0,-1),r(0,1),s(-1,0),s(1,0). Borrowed inputs unchanged;
/// outputs fully overwritten. Preflight preserves out, backend invalidates it.
StatusCode BuildMitc4LinearTying(const Matrix* coordinates, Matrix* out,
                                 Diagnostic* d);
/// Plane stress D5x5, engineering order xx,yy,xy,yz,xz, kappa5/6 once.
/// Same ownership/failure contract. d optional for all functions here.
StatusCode BuildMitc4LinearConstitutive(const ElasticMaterial* material,
                                        Matrix* out, Diagnostic* d);
/// At finite natural r,s in[-1,1], computes grad2x4 (x,y rows), transverse
/// engineering shears2x24 (xz,yz rows), positive detJ in length squared.
/// coordinates4x3 and ties4x24 are initial validated local inputs. All matrices
/// distinct. Full overwrite, no retained factors; scalar only on success.
StatusCode ComputeMitc4LinearPointOperators(const Matrix* coordinates,
                                            const Matrix* ties, double r,
                                            double s, Matrix* gradients,
                                            Matrix* shears, double* determinant,
                                            Diagnostic* d);
/// Builds full5x24 B from grad2x4/shears2x24 and finite physical z. Engineering
/// row order xx,yy,xy,yz,xz; global-independent local node-major6 components.
/// Separate input/output matrices, preflight preserves out, backend
/// invalidates.
StatusCode FillMitc4LinearB(const Matrix* gradients, const Matrix* shears,
                            double z, Matrix* out, Diagnostic* d);
}  // namespace fesa
#endif  // FESA_ELEMENTS_MITC4_LINEAR_OPS_H_
