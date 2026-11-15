/// @file
/// @brief Private shared selected-D geometry graph; no numerical state
/// retained.
#ifndef FESA_ELEMENTS_BEAM3D_KINEMATIC_GRAPH_H_
#define FESA_ELEMENTS_BEAM3D_KINEMATIC_GRAPH_H_

#include <cmath>
#include <limits>

#include "elements/beam3d_kinematics.h"

namespace fesa {
namespace beam3d_graph {

/// Roundoff allowance for unit triads and unresolved rotation branches.
constexpr double kRoundoff = 64 * std::numeric_limits<double>::epsilon();
/// Principal twist branch limit.
constexpr double kPi = 3.141592653589793238462643383279502884;
/// Number of independent spatial coordinates, including both eta variables.
constexpr int kVariables = 14;

/// Reports invalid borrowed input without repairing its geometry.
inline StatusCode Invalid(Diagnostic* d) {
  return ReportFailure(StatusCode::kInvalidArgument,
                       DiagnosticStage::kValidation,
                       "BeamKinematics input shape/finite/unit triad", d);
}
/// Sets the required typed failure and preserves the original operation label.
inline StatusCode Geometry(TrialEvaluationFailure kind, const char* operation,
                           TrialEvaluationFailure* failure, Diagnostic* d) {
  *failure = kind;
  return ReportFailure(StatusCode::kEquilibriumFailure,
                       DiagnosticStage::kComputation, operation, d);
}

/// Tests borrowed coefficients for finiteness without modifying them.
inline bool Finite(const double* data, int count) {
  for (int i = 0; i < count; ++i)
    if (!std::isfinite(data[i])) return false;
  return true;
}
/// Requires finite endpoint rows with shape2x3.
inline bool Shape(const Matrix* matrix) {
  return matrix && matrix->Rows() == 2 && matrix->Columns() == 3 &&
         Finite(matrix->Data(), 6);
}
/// K1 argument validation shared with the time coefficient evaluator.
inline bool Valid(const BeamKinematicInput* input, double xi) {
  return input && std::isfinite(xi) && xi >= 0 && xi <= 1 &&
         std::isfinite(input->length) && input->length > 0 &&
         Shape(input->positions) && Shape(input->tangents) &&
         (!input->relative_positions || Shape(input->relative_positions)) &&
         Shape(input->normals) && input->internal_position &&
         input->internal_position->Size() == 2 &&
         Finite(input->internal_position->Data(), 2) && input->direction &&
         input->direction->Size() == kVariables &&
         Finite(input->direction->Data(), kVariables);
}

/// Minimal rotation transport from the current reference-a. Both a and normal
/// remain differential quantities, including in the S derivative below.
/// Inputs/ops are borrowed during this call. base receives an arena pointer;
/// base_prime is required and overwritten only when target_prime is nonnull.
/// Required failure receives kRotationBranch for unresolved antiparallel
/// transport; optional d preserves original backend diagnostics. Discard output
/// on any failure. No input quantity or model state is changed.
template <typename Ops, typename Quantity>
StatusCode Transport(Ops* ops, const Quantity* a, const Quantity* normal,
                     const Quantity* target, const Quantity* target_prime,
                     Quantity** base, Quantity** base_prime,
                     TrialEvaluationFailure* failure, Diagnostic* d) {
  const Quantity* denominator = ops->Offset(ops->Dot(a, target), 1);
  if (ops->status() != StatusCode::kOk) return ops->status();
  if (ops->Value(denominator) <= kRoundoff)
    return Geometry(TrialEvaluationFailure::kRotationBranch,
                    "BeamKinematics antiparallel transport", failure, d);
  const Quantity* inverse = ops->Reciprocal(denominator);
  const Quantity* factor = ops->Product(ops->Dot(normal, target), inverse);
  const Quantity* together = ops->Add(a, target);
  *base = ops->Add(normal, ops->Product(together, factor), -1);
  if (target_prime) {
    const Quantity* derivative = ops->Product(
        ops->Add(ops->Dot(normal, target_prime),
                 ops->Product(ops->Dot(a, target_prime), factor), -1),
        inverse);
    *base_prime = ops->Scale(ops->Add(ops->Product(together, derivative),
                                      ops->Product(target_prime, factor)),
                             -1);
  }
  return ops->status();
}

/// Evaluates the original K1 Hermite/transport/normal/S-derivative/strain
/// graph. All quantity pointers are call-lifetime borrows of the caller's
/// arena. Ops supplies exact algebra and sticky common-math failure
/// propagation. fields receives six arena pointers on success; failure
/// invalidates them. length and xi use initial S, metres and dimensionless
/// reference position.
/// Inputs/fields/failure are required; d is optional. Unit/orthogonality errors
/// return kInvalidArgument, zero stretch/rotation branches return typed
/// kEquilibriumFailure, and backend errors preserve their original diagnostics.
/// @see docs/Linear static/b33_kinematic_kernel_contract.md.
template <typename Ops, typename Quantity>
StatusCode Evaluate(Ops* ops, const Quantity* xa, const Quantity* xb,
                    const Quantity* relative_xa, const Quantity* relative_xb,
                    const Quantity* ta, const Quantity* tb, const Quantity* na,
                    const Quantity* nb, const Quantity* la, const Quantity* lb,
                    double length, double xi, Quantity** fields,
                    TrialEvaluationFailure* failure, Diagnostic* d) {
  const Quantity* reference_squared = nullptr;
  for (int i = 0; i < 2; ++i) {
    const Quantity* t = i == 0 ? ta : tb;
    const Quantity* n = i == 0 ? na : nb;
    const Quantity* tt = ops->Dot(t, t);
    if (i == 0) reference_squared = tt;
    const Quantity* nn = ops->Dot(n, n);
    const Quantity* tn = ops->Dot(t, n);
    if (ops->status() != StatusCode::kOk) return ops->status();
    if (std::abs(ops->Value(tt) - 1) > kRoundoff ||
        std::abs(ops->Value(nn) - 1) > kRoundoff ||
        std::abs(ops->Value(tn)) > kRoundoff)
      return Invalid(d);
  }
  const Quantity* pa = ops->Product(ta, la);
  const Quantity* pb = ops->Product(tb, lb);
  const double x2 = xi * xi, x3 = x2 * xi;
  // Hermite derivatives are with respect to initial S, hence 1/L and 1/L^2.
  Quantity* position =
      ops->Combination(xa, pa, xb, pb, 1 - 3 * x2 + 2 * x3, xi - 2 * x2 + x3,
                       3 * x2 - 2 * x3, -x2 + x3);
  const Quantity* first =
      ops->Combination(relative_xa, pa, relative_xb, pb,
                       (-6 * xi + 6 * x2) / length,
                       (1 - 4 * xi + 3 * x2) / length,
                       (6 * xi - 6 * x2) / length, (-2 * xi + 3 * x2) / length);
  const Quantity* second = ops->Combination(
      relative_xa, pa, relative_xb, pb, (-6 + 12 * xi) / (length * length),
      (-4 + 6 * xi) / (length * length), (6 - 12 * xi) / (length * length),
      (-2 + 6 * xi) / (length * length));
  const Quantity* squared = ops->Dot(first, first);
  if (ops->status() != StatusCode::kOk) return ops->status();
  if (!(ops->Value(squared) > 0))
    return Geometry(TrialEvaluationFailure::kGeometry,
                    "BeamKinematics zero centerline derivative", failure, d);
  Quantity* stretch = ops->SquareRoot(squared);
  const Quantity* inverse = ops->Reciprocal(stretch);
  Quantity* tangent = ops->Product(first, inverse);
  const Quantity* tangent_prime = ops->Product(
      ops->Add(second, ops->Product(tangent, ops->Dot(tangent, second)), -1),
      inverse);
  Quantity *base = nullptr, *base_prime = nullptr, *endpoint = nullptr;
  StatusCode status = Transport<Ops, Quantity>(
      ops, ta, na, tangent, tangent_prime, &base, &base_prime, failure, d);
  if (status != StatusCode::kOk) return status;
  status = Transport<Ops, Quantity>(ops, ta, na, tb, nullptr, &endpoint,
                                    nullptr, failure, d);
  if (status != StatusCode::kOk) return status;
  const Quantity* angle = ops->Angle(ops->Dot(tb, ops->Cross(endpoint, nb)),
                                     ops->Dot(endpoint, nb));
  if (ops->status() != StatusCode::kOk) return ops->status();
  if (kPi - std::abs(ops->Value(angle)) <= kRoundoff)
    return Geometry(TrialEvaluationFailure::kRotationBranch,
                    "BeamKinematics ambiguous pi endpoint twist", failure, d);
  const Quantity* phase = ops->Scale(angle, xi);
  const Quantity* cosine = ops->Cos(phase);
  const Quantity* sine = ops->Sin(phase);
  const Quantity* twist_rate = ops->Scale(angle, 1 / length);
  const Quantity* across = ops->Cross(tangent, base);
  Quantity* normal1 =
      ops->Add(ops->Product(base, cosine), ops->Product(across, sine));
  // Exact S derivative of cos(xi*phi)*base + sin(xi*phi)*(t x base).
  const Quantity* normal1_prime =
      ops->Combination(ops->Product(base, ops->Product(sine, twist_rate)),
                       ops->Product(base_prime, cosine),
                       ops->Product(across, ops->Product(cosine, twist_rate)),
                       ops->Product(ops->Add(ops->Cross(tangent_prime, base),
                                             ops->Cross(tangent, base_prime)),
                                    sine),
                       -1, 1, 1, 1);
  Quantity* normal2 = ops->Cross(tangent, normal1);
  // Exact Green factorization, including the moving reference tangent's
  // variations and its roundoff-unit defect. No strain clipping/repair.
  const Quantity* difference = ops->Add(first, ta, -1);
  const Quantity* together = ops->Add(first, ta);
  const Quantity* green = ops->Add(ops->Dot(difference, together),
                                  ops->Offset(reference_squared, -1));
  ops->Pack(ops->Scale(green, .5), 0);
  ops->Pack(ops->Scale(ops->Dot(normal2, second), -1), 1);
  ops->Pack(ops->Dot(normal1, second), 2);
  ops->Pack(ops->Dot(normal2, normal1_prime), 3);
  if (ops->status() != StatusCode::kOk) return ops->status();
  fields[0] = position;
  fields[1] = tangent;
  fields[2] = normal1;
  fields[3] = normal2;
  fields[4] = ops->Strain();
  fields[5] = stretch;
  return ops->status();
}

}  // namespace beam3d_graph
}  // namespace fesa

#endif  // FESA_ELEMENTS_BEAM3D_KINEMATIC_GRAPH_H_
