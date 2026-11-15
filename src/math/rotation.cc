#include "fesa/math/rotation.h"

#include <cmath>
#include <limits>
#include <memory>

namespace fesa {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kUnitTolerance =
    4096.0 * std::numeric_limits<double>::epsilon();
constexpr double kBranchMinimum =
    256.0 * std::numeric_limits<double>::epsilon();

StatusCode Invalid(const char* operation, Diagnostic* d) {
  return ReportFailure(StatusCode::kInvalidArgument,
                       DiagnosticStage::kValidation, operation, d);
}

bool FiniteVector(const Vector* value, int size) {
  if (!value || value->Size() != size) return false;
  for (int i = 0; i < size; ++i)
    if (!std::isfinite(value->Data()[i])) return false;
  return true;
}

bool MatrixShape(const Matrix* value) {
  return value && value->Rows() == 3 && value->Columns() == 3;
}

StatusCode UnitNorm(const Vector* q, double* norm, Diagnostic* d) {
  if (!FiniteVector(q, 4)) return Invalid("Rotation finite quaternion", d);
  const auto status = q->Norm2(norm, d);
  if (status != StatusCode::kOk) return status;
  if (std::abs(*norm - 1.0) > kUnitTolerance)
    return Invalid("Rotation quaternion unit tolerance", d);
  return StatusCode::kOk;
}

StatusCode CreateVector(int size, std::unique_ptr<Vector>* owner,
                        Diagnostic* d) {
  Vector* raw = nullptr;
  const auto status = Vector::Create(size, &raw, d);
  owner->reset(raw);
  return status;
}

StatusCode CreateMatrix(int size, std::unique_ptr<Matrix>* owner,
                        Diagnostic* d) {
  Matrix* raw = nullptr;
  const auto status = Matrix::Create(size, size, &raw, d);
  owner->reset(raw);
  return status;
}

StatusCode HamiltonMatrix(const Vector* left, Matrix* out, Diagnostic* d) {
  const double* q = left->Data();
  // Hamilton left-product operator in wxyz order. These are scalar
  // coefficient definitions; the quaternion product is a library matvec.
  const double coefficients[16] = {q[0],  -q[1], -q[2], -q[3], q[1], q[0],
                                   -q[3], q[2],  q[2],  q[3],  q[0], -q[1],
                                   q[3],  -q[2], q[1],  q[0]};
  for (int i = 0; i < 16; ++i) {
    const auto status = out->Set(i / 4, i % 4, coefficients[i], d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}

StatusCode QuaternionLog(const Vector* q, bool principal, Vector* out,
                         Diagnostic* d) {
  std::unique_ptr<Vector> imaginary;
  auto status = CreateVector(3, &imaginary, d);
  if (status != StatusCode::kOk) return status;
  const int indices[3] = {1, 2, 3};
  status = imaginary->GatherFrom(q, indices, 3, 0, d);
  if (status != StatusCode::kOk) return status;
  double sine = 0.0;
  status = imaginary->Norm2(&sine, d);
  if (status != StatusCode::kOk) return status;
  const double w = q->Data()[0];
  const double angle = 2.0 * std::atan2(sine, principal ? std::abs(w) : w);
  double sign = principal && w < 0.0 ? -1.0 : 1.0;
  if (principal && angle == kPi) {
    // The first largest absolute component resolves exact-pi ties X,Y,Z.
    // Canonicalization belongs only to official principal output.
    double maximum = 0.0;
    int axis = 0;
    status = imaginary->MaximumAbs(&maximum, &axis, d);
    if (status != StatusCode::kOk) return status;
    sign = imaginary->Data()[axis] < 0.0 ? -1.0 : 1.0;
  }
  status = out->CopyFrom(imaginary.get(), d);
  if (status != StatusCode::kOk) return status;
  return out->Scale(sine == 0.0 ? 2.0 : sign * angle / sine, d);
}

StatusCode RelativeIncrement(const Vector* reference, const Vector* current,
                             Vector* out, bool* branch_exceeded, bool body,
                             Diagnostic* d) {
  if (!branch_exceeded) return Invalid("Rotation increment branch flag", d);
  *branch_exceeded = false;
  if (!out || out->Size() != 3 || out == reference || out == current)
    return Invalid("Rotation increment output/alias", d);
  double norm = 0.0;
  auto status = UnitNorm(reference, &norm, d);
  if (status != StatusCode::kOk) return status;
  status = UnitNorm(current, &norm, d);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<Matrix> conjugation;
  std::unique_ptr<Vector> conjugate;
  std::unique_ptr<Vector> relative;
  status = CreateMatrix(4, &conjugation, d);
  if (status != StatusCode::kOk) return status;
  status = CreateVector(4, &conjugate, d);
  if (status != StatusCode::kOk) return status;
  status = CreateVector(4, &relative, d);
  if (status != StatusCode::kOk) return status;
  for (int i = 0; i < 4; ++i) {
    status = conjugation->Set(i, i, i == 0 ? 1.0 : -1.0, d);
    if (status != StatusCode::kOk) return status;
  }
  status = Matrix::MultiplyVector(conjugation.get(), false, reference,
                                  conjugate.get(), d);
  if (status != StatusCode::kOk) return status;
  status = body ? ComposeRotations(conjugate.get(), current, relative.get(), d)
                : ComposeRotations(current, conjugate.get(), relative.get(), d);
  if (status != StatusCode::kOk) return status;
  // state_and_rotation_contract.md: never flip the cover to rescue an
  // increment outside the continuous log chart of the committed reference.
  if (!(relative->Data()[0] > kBranchMinimum)) {
    *branch_exceeded = true;
    return ReportFailure(StatusCode::kEquilibriumFailure,
                         DiagnosticStage::kComputation,
                         "Rotation increment pi branch exceeded", d);
  }
  return QuaternionLog(relative.get(), false, out, d);
}

StatusCode SkewMatrix(const Vector* b, Matrix* out, Diagnostic* d) {
  const double* v = b->Data();
  const double coefficients[9] = {0.0,   -v[2], v[1], v[2], 0.0,
                                  -v[0], -v[1], v[0], 0.0};
  for (int i = 0; i < 9; ++i) {
    const auto status = out->Set(i / 3, i % 3, coefficients[i], d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}

}  // namespace

StatusCode NormalizeQuaternion(const Vector* q, Vector* out, Diagnostic* d) {
  if (!out || out->Size() != 4 || q == out)
    return Invalid("NormalizeQuaternion output/alias", d);
  double norm = 0.0;
  auto status = UnitNorm(q, &norm, d);
  if (status != StatusCode::kOk) return status;
  status = out->CopyFrom(q, d);
  if (status != StatusCode::kOk) return status;
  return out->Scale(1.0 / norm, d);
}

StatusCode RotationExp(const Vector* v, Vector* out, Diagnostic* d) {
  if (!FiniteVector(v, 3) || !out || out->Size() != 4 || out == v)
    return Invalid("RotationExp shape/alias/nonfinite", d);
  double angle = 0.0;
  auto status = v->Norm2(&angle, d);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<Vector> quaternion;
  status = CreateVector(4, &quaternion, d);
  if (status != StatusCode::kOk) return status;
  const int indices[3] = {0, 1, 2};
  status = quaternion->GatherFrom(v, indices, 3, 1, d);
  if (status != StatusCode::kOk) return status;
  // sin(theta/2)/theta is regular at zero. No principal-sign rule is
  // applied: Exp beyond pi must preserve the continuous quaternion cover.
  const double coefficient = angle == 0.0 ? 0.5 : std::sin(angle / 2.0) / angle;
  status = quaternion->Scale(coefficient, d);
  if (status != StatusCode::kOk) return status;
  status = quaternion->Set(0, std::cos(angle / 2.0), d);
  if (status != StatusCode::kOk) return status;
  return NormalizeQuaternion(quaternion.get(), out, d);
}

StatusCode ComposeRotations(const Vector* left, const Vector* right,
                            Vector* out, Diagnostic* d) {
  if (!out || out->Size() != 4 || out == left || out == right)
    return Invalid("ComposeRotations output/alias", d);
  double left_norm = 0.0;
  double right_norm = 0.0;
  auto status = UnitNorm(left, &left_norm, d);
  if (status != StatusCode::kOk) return status;
  status = UnitNorm(right, &right_norm, d);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<Matrix> hamilton;
  std::unique_ptr<Vector> product;
  status = CreateMatrix(4, &hamilton, d);
  if (status != StatusCode::kOk) return status;
  status = CreateVector(4, &product, d);
  if (status != StatusCode::kOk) return status;
  status = HamiltonMatrix(left, hamilton.get(), d);
  if (status != StatusCode::kOk) return status;
  status =
      Matrix::MultiplyVector(hamilton.get(), false, right, product.get(), d);
  if (status != StatusCode::kOk) return status;
  // Each accepted input may carry its own rounding-scale norm error. Remove
  // both positive scales before checking and normalizing the computed product.
  status = product->Scale(1.0 / left_norm, d);
  if (status != StatusCode::kOk) return status;
  status = product->Scale(1.0 / right_norm, d);
  if (status != StatusCode::kOk) return status;
  return NormalizeQuaternion(product.get(), out, d);
}

StatusCode RotationMatrix(const Vector* q, Matrix* out, Diagnostic* d) {
  if (!MatrixShape(out)) return Invalid("RotationMatrix output shape", d);
  double norm = 0.0;
  auto status = UnitNorm(q, &norm, d);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<Vector> unit;
  status = CreateVector(4, &unit, d);
  if (status != StatusCode::kOk) return status;
  status = unit->CopyFrom(q, d);
  if (status != StatusCode::kOk) return status;
  status = unit->Scale(1.0 / norm, d);
  if (status != StatusCode::kOk) return status;
  const double w = unit->Data()[0];
  const double x = unit->Data()[1];
  const double y = unit->Data()[2];
  const double z = unit->Data()[3];
  // Standard unit-quaternion rotation coefficients, initial to current axes.
  const double coefficients[9] = {
      1.0 - 2.0 * (y * y + z * z), 2.0 * (x * y - w * z),
      2.0 * (x * z + w * y),       2.0 * (x * y + w * z),
      1.0 - 2.0 * (x * x + z * z), 2.0 * (y * z - w * x),
      2.0 * (x * z - w * y),       2.0 * (y * z + w * x),
      1.0 - 2.0 * (x * x + y * y)};
  for (int i = 0; i < 9; ++i) {
    status = out->Set(i / 3, i % 3, coefficients[i], d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}

StatusCode PrincipalRotationVector(const Vector* q, Vector* out,
                                   Diagnostic* d) {
  if (!out || out->Size() != 3 || out == q)
    return Invalid("PrincipalRotationVector output/alias", d);
  double norm = 0.0;
  const auto status = UnitNorm(q, &norm, d);
  if (status != StatusCode::kOk) return status;
  return QuaternionLog(q, true, out, d);
}

StatusCode BodyRotationIncrement(const Vector* reference, const Vector* current,
                                 Vector* out, bool* branch_exceeded,
                                 Diagnostic* d) {
  return RelativeIncrement(reference, current, out, branch_exceeded, true, d);
}

StatusCode SpatialRotationIncrement(const Vector* reference,
                                    const Vector* current, Vector* out,
                                    bool* branch_exceeded, Diagnostic* d) {
  return RelativeIncrement(reference, current, out, branch_exceeded, false, d);
}

StatusCode InverseLeftRotationJacobian(const Vector* b, Matrix* out,
                                       Diagnostic* d) {
  if (!FiniteVector(b, 3) || !MatrixShape(out))
    return Invalid("InverseLeftRotationJacobian shape/nonfinite", d);
  double angle = 0.0;
  auto status = b->Norm2(&angle, d);
  if (status != StatusCode::kOk) return status;
  if (!(angle < kPi)) return Invalid("InverseLeftRotationJacobian pi chart", d);
  std::unique_ptr<Matrix> skew;
  std::unique_ptr<Matrix> square;
  std::unique_ptr<Matrix> identity;
  status = CreateMatrix(3, &skew, d);
  if (status != StatusCode::kOk) return status;
  status = CreateMatrix(3, &square, d);
  if (status != StatusCode::kOk) return status;
  status = CreateMatrix(3, &identity, d);
  if (status != StatusCode::kOk) return status;
  status = SkewMatrix(b, skew.get(), d);
  if (status != StatusCode::kOk) return status;
  status =
      Matrix::Multiply(skew.get(), false, skew.get(), false, square.get(), d);
  if (status != StatusCode::kOk) return status;
  for (int i = 0; i < 3; ++i) {
    status = identity->Set(i, i, 1.0, d);
    if (status != StatusCode::kOk) return status;
  }
  // Solà B145-146 and the R6 contract. The small-angle branch avoids the
  // cancellation in 1/theta^2-1/(2 theta tan(theta/2)).
  const double squared = angle * angle;
  const double coefficient =
      angle <= 0.01
          ? 1.0 / 12.0 + squared / 720.0 + squared * squared / 30240.0
          : 1.0 / squared - 1.0 / (2.0 * angle * std::tan(angle / 2.0));
  status = out->CopyFrom(identity.get(), d);
  if (status != StatusCode::kOk) return status;
  status = out->AddScaled(-0.5, skew.get(), d);
  if (status != StatusCode::kOk) return status;
  return out->AddScaled(coefficient, square.get(), d);
}

}  // namespace fesa
