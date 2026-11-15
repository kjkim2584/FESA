/// @file
/// @brief Borrowed physical dynamic input/results for the common Element API.
#ifndef FESA_MODEL_ELEMENT_DYNAMIC_H_
#define FESA_MODEL_ELEMENT_DYNAMIC_H_

#include "fesa/model/element_trial.h"

namespace fesa {
/// Linear Cartesian DOFs or finite spatial rotations; not an analysis policy.
enum class DynamicKinematics { kLinear, kFiniteRotation };
/// Both physical forces are required; unrequested derivative buffers are null.
enum class DynamicEvaluationRequest { kForces, kForcesAndLinearization };
/// Caller owns all buffers and keeps them alive through evaluation. No object
/// or numerical buffer is retained by an element. Inputs may alias each other.
/// d=node_count*dofs_per_node; vectors use connectivity then component order.
struct ElementMotionInput {
  DynamicKinematics kinematics = DynamicKinematics::kLinear;  ///< Pose kind.
  /// Borrowed d-vector [u] or [u,small rotation], L/rad; linear mode only.
  const Vector* linear_displacement = nullptr;
  /// Borrowed trial pose, finite mode only. translations 3n in global L,
  /// orientations n-by-4 unit wxyz initial-to-current. Optional trial history
  /// keeps its existing meaning and lifetime; it is not a new commit source.
  const ElementTrialState* finite_state = nullptr;
  /// Borrowed d-vector [v] or [v,omega], global L/T and spatial rad/T.
  const Vector* velocity = nullptr;
  /// Borrowed d-vector [a] or [a,eta], global L/T^2 and spatial rad/T^2.
  const Vector* acceleration = nullptr;
  /// Optional Beam internal eta velocity (2), connectivity a/b, m/s.
  /// Caller-owned call-lifetime borrow, never mutated or retained. Beam finite
  /// role only; other families leave null. Read-only inputs may alias.
  const Vector* internal_velocity = nullptr;
  /// Optional Beam internal eta acceleration (2), connectivity a/b, m/s^2.
  /// Same ownership, lifetime, alias and finite-role contract as above.
  const Vector* internal_acceleration = nullptr;
};
/// Caller-owned local d-vectors and d-by-d matrices, borrowed for this call.
/// d=physical+internal: generic node_count*dofs_per_node, finite Beam 12+2=14.
/// Rows are global physical force/moment (F/FL); Beam appends internal force,
/// giving [F,FL,F].
/// Columns are spatial q corrections (L/rad/m), V (L/T,rad/T,m/s), or
/// A (L/T^2,rad/T^2,m/s^2); internal derivatives use eta in m. Row-major.
/// Every output is disjoint from inputs and other outputs of the same type.
/// No HHT weighting, external/structural force or official result storage.
/// A supported evaluation overwrites all requested blocks, including zeros.
/// Preflight errors preserve values; failures after calculation invalidate all
/// requested outputs. No partial success or cached derivative blocks.
struct ElementDynamicResult {
  Vector* inertia_force = nullptr;  ///< fI(z,V,A).
  Vector* damping_force = nullptr;  ///< fD(z,V).
  Matrix* inertia_displacement = nullptr;  ///< GI: dfI/dz, V/A fixed.
  Matrix* inertia_velocity = nullptr;      ///< CI: dfI/dV, z/A fixed.
  Matrix* inertia_acceleration = nullptr;  ///< MI: dfI/dA, z/V fixed.
  Matrix* damping_displacement = nullptr;  ///< GD: dfD/dz, V fixed.
  Matrix* damping_velocity = nullptr;      ///< CD: dfD/dV, z fixed.
};
/// Checks positive node count, 3/6 DOFs, supported pose representation,
/// finite input coefficients, buffer shapes, requests and output aliasing.
/// Input-only aliasing is allowed. Optional finite trial history is checked
/// when supplied; quaternion normalization/geometry are the element's job.
/// Finite rotation requires 6 DOFs; finite Solid dynamics is unsupported.
/// Generic role rejects all nonnull Beam internal pose/motion inputs.
/// No allocations, numerical algebra, writes to inputs/outputs, capability or
/// mass-mode decisions. Diagnostic is optional; success preserves it.
/// The actual Element API also validates MassMode and const workspace; see
/// docs/HHT/implementation_contract.md for final signature and ownership.
StatusCode ValidateElementDynamicEvaluation(
    int node_count, int dofs_per_node, const ElementMotionInput* input,
    DynamicEvaluationRequest request, const ElementDynamicResult* out,
    Diagnostic* d);
/// Read-only, allocation/backend-free finite Beam preflight: physical 2x6
/// motion plus eta/velocity/acceleration (2 each), local outputs 14/14x14.
/// Requires finite rotation; optional history retains the generic contract.
/// Input-only aliasing is allowed, outputs must be disjoint. No geometry,
/// mass or capability decision. Failure preserves buffers and reports the
/// original validation operation; success preserves optional Diagnostic.
StatusCode ValidateBeamElementDynamicEvaluation(
    const ElementMotionInput* input, DynamicEvaluationRequest request,
    const ElementDynamicResult* out, Diagnostic* d);
}  // namespace fesa
#endif  // FESA_MODEL_ELEMENT_DYNAMIC_H_
