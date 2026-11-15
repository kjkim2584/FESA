/// @file
/// @brief P17 node-local finite-rotation physical initial compatibility.
#ifndef FESA_ANALYSIS_FINITE_SHELL_INITIALIZATION_H_
#define FESA_ANALYSIS_FINITE_SHELL_INITIALIZATION_H_
#include "fesa/math/vector.h"
namespace fesa {
/// Evaluates N6/N7 normalized errors for one allowed current mass-null node.
/// Caller supplies validated unit initial normal/current director, positive
/// kappa[FL/rad], body_spin[rad] in its existing material components, spatial
/// omega[rad/T], physical applied M[FL] and right Mdot[FL/T]; all vectors3.
/// Inputs borrowed/read-only and may alias each other. Distinct scalar outputs
/// required, both preserved on every failure. Zero denominator/numerator
/// gives0. Exponent-scaled terms avoid overflow of kappa*spin and M*omega. No
/// projection, model/state/result changes, decision or retained scratch. Caller
/// owns geometry/rank/mapping validation. @see
/// docs/HHT/initialization.md.
StatusCode EvaluateFiniteShellNodeCompatibility(
    const Vector* initial_normal, const Vector* director, double kappa,
    const Vector* body_spin, const Vector* omega, const Vector* moment,
    const Vector* moment_rate, double* configuration_error,
    double* velocity_error, Diagnostic* d);
/// Applies P17 e_q,e_v<=1e-10 to the same borrowed data, with original node ID.
/// Incompatibility returns kEquilibriumFailure and node/phase diagnostic;
/// invalid inputs and numerical/backend failures retain their actual status.
/// No inputs/state mutated or retained. This is physical initialization only,
/// never an increment convergence/reference decision or motion correction.
StatusCode CheckFiniteShellNodeCompatibility(
    int node_id, const Vector* initial_normal, const Vector* director,
    double kappa, const Vector* body_spin, const Vector* omega,
    const Vector* moment, const Vector* moment_rate, Diagnostic* d);
/// Computes N8 gamma=(bg-A.dot(eta_range))/kappa [rad/T^2], where
/// A=kappa*d-d cross M and bg=(omega cross(omega cross d)).M
/// +2*(omega cross d).Mdot+d.Mddot. Caller applies P17 first and supplies
/// validated unit current director, right physical loads/rates, and the
/// angular range solution of current mass. Vector inputs3, read-only and
/// borrowed; scalar gamma required and preserved on failure. No mass solve,
/// projection, motion mutation, retained scratch or artificial inertia.
/// Factored exponents avoid overflowing intermediate physical products; an
/// unrepresentable nonzero result is a fatal numerical failure.
StatusCode ComputeFiniteShellNodeNullAcceleration(
    const Vector* director, double kappa, const Vector* omega,
    const Vector* moment, const Vector* moment_rate,
    const Vector* moment_second_rate, const Vector* range_acceleration,
    double* gamma, Diagnostic* d);
}  // namespace fesa
#endif  // FESA_ANALYSIS_FINITE_SHELL_INITIALIZATION_H_
