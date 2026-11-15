/// @file
/// @brief Initial linear MITC4 operators, independent of finite trial state.
#ifndef FESA_ELEMENTS_MITC4_INITIAL_H_
#define FESA_ELEMENTS_MITC4_INITIAL_H_
#include "fesa/math/matrix.h"
#include "fesa/model/node.h"
#include "fesa/model/shell_section.h"
namespace fesa {
/// All global operators use node-major [U1,U2,U3,UR1,UR2,UR3] columns
/// and [F1,F2,F3,M1,M2,M3] rows. Matrix outputs are fully overwritten.
/// Formulas, units and independent cases:
/// docs/MITC4/initial_kernel_contract.md. Borrows four distinct six-DOF nodes
/// and section. Writes full global24x24 K and kd: min positive local
/// tangential-rotation diagonal times1e-3. K includes the four local drilling
/// diagonal springs. Initial 2x2x2 rule; no state/mass/results evaluation. No
/// cross-call data are retained. Output pointers required; preflight/geometry
/// errors preserve outputs; backend failure invalidates matrix but preserves
/// kd. d is optional.
StatusCode ComputeMitc4InitialStiffness(const Node* const* nodes, int count,
                                        const ShellSection* section,
                                        Matrix* stiffness, double* drilling,
                                        Diagnostic* d);
/// RD12 initial elastic membrane/bending operator, global24x24 [F,FL]/[L,rad].
/// Same integration/geometry as initial K; excludes transverse shear and all
/// drilling springs. Borrows nodes/section for this call, retains no data.
/// Preflight/geometry errors preserve out; later failure invalidates it.
StatusCode ComputeMitc4InitialElasticStiffness(const Node* const* nodes,
                                               int count,
                                               const ShellSection* section,
                                               Matrix* out, Diagnostic* d);
/// Writes initial global24x24 consistent mass by2x2x2 integration. Density is
/// required only here; missing density is kInputError. No constitutive matrix,
/// stiffness, drilling inertia or official results. Same validity contract.
StatusCode ComputeMitc4InitialMass(const Node* const* nodes, int count,
                                   const ShellSection* section, Matrix* mass,
                                   Diagnostic* d);
/// Converts borrowed global24x24 consistent mass to local same-component
/// node row sums, then transforms the lumped matrix back to global axes.
/// Reconstructs initial axes from borrowed nodes, without new integration,
/// density or K evaluation. Both matrices required, matching and distinct.
/// Preflight/geometry failure preserves inputs/outputs; backend failure
/// invalidates lumped output only. No input mutation or retained data.
StatusCode LumpMitc4InitialMass(const Node* const* nodes, int count,
                                const Matrix* consistent, Matrix* lumped,
                                Diagnostic* d);
/// Integrates linear B^T D B*u directly, without constructing K or mass.
/// u is borrowed global node-major [U1,U2,U3,UR1,UR2,UR3], size24;
/// force is separate size24, rows [F1,F2,F3,M1,M2,M3]. drilling is the fixed
/// initial positive spring coefficient supplied by the element factory.
/// Preflight/geometry failure preserves force; backend failure invalidates it.
/// No official result mutation or retained data; d is optional.
StatusCode ComputeMitc4LinearInternalForce(const Node* const* nodes, int count,
                                           const ShellSection* section,
                                           const Vector* u, double drilling,
                                           Vector* force, Diagnostic* d);
/// Computes immutable output-model initial axes3x3 and initial kd while
/// deriving and discarding temporary K. Shares the same geometry/integration
/// as ComputeMitc4InitialStiffness without a second frame reconstruction.
/// Both outputs required, kd only on full success; backend invalidates axes.
StatusCode ComputeMitc4InitialProperties(const Node* const* nodes, int count,
                                        const ShellSection* section,
                                        Matrix* axes, double* drilling,
                                        Diagnostic* d);
}  // namespace fesa
#endif  // FESA_ELEMENTS_MITC4_INITIAL_H_
