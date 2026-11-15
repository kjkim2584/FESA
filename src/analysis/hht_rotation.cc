#include "analysis/hht_rotation.h"

#include <cmath>
#include <limits>
#include <memory>

#include "fesa/math/rotation.h"
#include "model/model_storage.h"

namespace fesa {
namespace {
StatusCode Invalid(Diagnostic* d) {
  return model_storage::Invalid("HHT rotational kinematics inputs", d);
}
StatusCode NewVector(int count, std::unique_ptr<Vector>* owner, Diagnostic* d) {
  Vector* raw = nullptr;
  const auto status = Vector::Create(count, &raw, d);
  owner->reset(raw);
  return status;
}
StatusCode NewMatrix(int rows, int columns, std::unique_ptr<Matrix>* owner,
                     Diagnostic* d) {
  Matrix* raw = nullptr;
  const auto status = Matrix::Create(rows, columns, &raw, d);
  owner->reset(raw);
  return status;
}
bool Shape(const Matrix* matrix, int rows, int columns) {
  return matrix && matrix->Rows() == rows && matrix->Columns() == columns;
}
StatusCode Validate(const RotationalStateView* state,
                    const RotationalCandidate* candidate,
                    HhtRotationOperation operation,
                    const RotationalLinearization* out, Diagnostic* d) {
  if (!state || !candidate || !state->orientations ||
      state->orientations->Rows() <= 0 ||
      state->orientations->Rows() > std::numeric_limits<int>::max() / 3)
    return Invalid(d);
  const int count = state->orientations->Rows();
  if (!Shape(state->orientations, count, 4) ||
      !Shape(candidate->orientations, count, 4) ||
      state->orientations == candidate->orientations)
    return Invalid(d);
  const Vector* values[] = {state->velocity, state->acceleration,
                            candidate->velocity, candidate->acceleration};
  for (int i = 0; i < 4; ++i) {
    if (!values[i] || values[i]->Size() != 3 * count) return Invalid(d);
    for (int j = 0; j < i; ++j)
      if (values[i] == values[j]) return Invalid(d);
    if (i < 2)
      for (int entry = 0; entry < 3 * count; ++entry)
        if (!std::isfinite(values[i]->Data()[entry])) return Invalid(d);
  }
  if (operation == HhtRotationOperation::kLinearize) {
    if (!out || !Shape(out->velocity_per_rotation, count, 9) ||
        !Shape(out->acceleration_per_rotation, count, 9) ||
        out->velocity_per_rotation == out->acceleration_per_rotation)
      return Invalid(d);
  } else if (out || (operation != HhtRotationOperation::kPredict &&
                     operation != HhtRotationOperation::kUpdate)) {
    return Invalid(d);
  }
  return StatusCode::kOk;
}
StatusCode Row(const Matrix* matrix, int row, Vector* out, Diagnostic* d) {
  int indices[4] = {4 * row, 4 * row + 1, 4 * row + 2, 4 * row + 3};
  return matrix->GatherTo(indices, 4, out, 0, d);
}
// Unit/finite quaternion errors are argument failures: validate every node
// before any candidate or derivative block can be written.
StatusCode ValidateOrientations(const RotationalStateView* state,
                                const RotationalCandidate* candidate,
                                HhtRotationOperation operation, Diagnostic* d) {
  std::unique_ptr<Vector> q, normalized;
  auto status = NewVector(4, &q, d);
  if (status != StatusCode::kOk) return status;
  status = NewVector(4, &normalized, d);
  if (status != StatusCode::kOk) return status;
  const Matrix* sources[] = {
      state->orientations, operation == HhtRotationOperation::kPredict
                               ? nullptr
                               : candidate->orientations};
  for (const Matrix* source : sources) {
    if (!source) continue;
    for (int node = 0; node < source->Rows(); ++node) {
      status = Row(source, node, q.get(), d);
      if (status != StatusCode::kOk) return status;
      status = NormalizeQuaternion(q.get(), normalized.get(), d);
      if (status != StatusCode::kOk) return status;
    }
  }
  return StatusCode::kOk;
}
StatusCode Skew(const Vector* v, Matrix* out, Diagnostic* d) {
  const double* x = v->Data();
  const double data[] = {0, -x[2], x[1], x[2], 0, -x[0], -x[1], x[0], 0};
  for (int i = 0; i < 9; ++i) {
    const auto status = out->Set(i / 3, i % 3, data[i], d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}
StatusCode CheckPredictor(const Vector* theta, bool* branch, Diagnostic* d) {
  double angle = 0;
  const auto status = theta->Norm2(&angle, d);
  if (status != StatusCode::kOk) return status;
  constexpr double kPi = 3.141592653589793238462643383279502884;
  if (angle >= kPi ||
      std::cos(angle / 2) <= 256 * std::numeric_limits<double>::epsilon()) {
    *branch = true;
    return ReportFailure(StatusCode::kEquilibriumFailure,
                         DiagnosticStage::kComputation,
                         "HHT predictor rotation crosses pi", d);
  }
  return StatusCode::kOk;
}
}  // namespace
StatusCode EvaluateHhtRotations(const HhtParameters* parameters,
                                const RotationalStateView* state,
                                const RotationalCandidate* candidate, double dt,
                                HhtRotationOperation operation,
                                RotationalLinearization* out, bool* branch,
                                Diagnostic* d) {
  if (branch) *branch = false;
  if (!branch || !parameters || !std::isfinite(parameters->alpha) ||
      parameters->alpha < -1. / 3 || parameters->alpha > 0 ||
      !std::isfinite(dt) || dt <= 0)
    return Invalid(d);
  auto status = Validate(state, candidate, operation, out, d);
  if (status != StatusCode::kOk) return status;
  status = ValidateOrientations(state, candidate, operation, d);
  if (status != StatusCode::kOk) return status;
  const double beta = (1 - parameters->alpha) * (1 - parameters->alpha) / 4;
  const double gamma = .5 - parameters->alpha;
  const double c0 = (1 / beta / dt) / dt, c1 = gamma / beta / dt;
  const double h2 = dt * dt;
  if (!std::isfinite(c0) || !(c0 > 0) || !std::isfinite(c1) || !(c1 > 0) ||
      !std::isfinite(h2) || !(h2 > 0))
    return ReportFailure(StatusCode::kEquilibriumFailure,
                         DiagnosticStage::kComputation,
                         "HHT rotational coefficients not representable", d);
  std::unique_ptr<Vector> old_q, target_q, increment_q;
  for (auto* owner : {&old_q, &target_q, &increment_q}) {
    status = NewVector(4, owner, d);
    if (status != StatusCode::kOk) return status;
  }
  std::unique_ptr<Vector> omega, eta, b, c, qb, qc, theta, trial_v, trial_a;
  for (auto* owner :
       {&omega, &eta, &b, &c, &qb, &qc, &theta, &trial_v, &trial_a}) {
    status = NewVector(3, owner, d);
    if (status != StatusCode::kOk) return status;
  }
  std::unique_ptr<Matrix> rotation, inverse, skew, derivative;
  status = NewMatrix(3, 3, &rotation, d);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<Vector> block;
  if (operation == HhtRotationOperation::kLinearize) {
    for (auto* owner : {&inverse, &skew, &derivative}) {
      status = NewMatrix(3, 3, owner, d);
      if (status != StatusCode::kOk) return status;
    }
    status = NewVector(9, &block, d);
    if (status != StatusCode::kOk) return status;
  }
  const int identity[] = {0, 1, 2}, all[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  for (int node = 0; node < state->orientations->Rows(); ++node) {
    status = Row(state->orientations, node, old_q.get(), d);
    if (status != StatusCode::kOk) return status;
    const int indices[] = {3 * node, 3 * node + 1, 3 * node + 2};
    status = omega->GatherFrom(state->velocity, indices, 3, 0, d);
    if (status != StatusCode::kOk) return status;
    status = eta->GatherFrom(state->acceleration, indices, 3, 0, d);
    if (status != StatusCode::kOk) return status;
    if (operation == HhtRotationOperation::kPredict) {
      status = theta->CopyFrom(omega.get(), d);
      if (status != StatusCode::kOk) return status;
      status = theta->Scale(dt, d);
      if (status != StatusCode::kOk) return status;
      status = theta->AddScaled(h2 * (.5 - beta), eta.get(), d);
      if (status != StatusCode::kOk) return status;
      status = CheckPredictor(theta.get(), branch, d);
      if (status != StatusCode::kOk) return status;
    } else {
      status = Row(candidate->orientations, node, target_q.get(), d);
      if (status != StatusCode::kOk) return status;
      status = SpatialRotationIncrement(old_q.get(), target_q.get(),
                                        theta.get(), branch, d);
      if (status != StatusCode::kOk) return status;
    }
    status = RotationExp(theta.get(), increment_q.get(), d);
    if (status != StatusCode::kOk) return status;
    if (operation == HhtRotationOperation::kPredict) {
      status =
          ComposeRotations(increment_q.get(), old_q.get(), target_q.get(), d);
      if (status != StatusCode::kOk) return status;
    }
    status = RotationMatrix(increment_q.get(), rotation.get(), d);
    if (status != StatusCode::kOk) return status;
    status = b->CopyFrom(omega.get(), d);
    if (status != StatusCode::kOk) return status;
    status = b->Scale(1 - gamma / beta, d);
    if (status != StatusCode::kOk) return status;
    status = b->AddScaled(dt * (1 - gamma / (2 * beta)), eta.get(), d);
    if (status != StatusCode::kOk) return status;
    status = c->CopyFrom(omega.get(), d);
    if (status != StatusCode::kOk) return status;
    status = c->Scale(1 / beta / dt, d);
    if (status != StatusCode::kOk) return status;
    status = c->AddScaled(1 / (2 * beta) - 1, eta.get(), d);
    if (status != StatusCode::kOk) return status;
    status =
        Matrix::MultiplyVector(rotation.get(), false, b.get(), qb.get(), d);
    if (status != StatusCode::kOk) return status;
    status =
        Matrix::MultiplyVector(rotation.get(), false, c.get(), qc.get(), d);
    if (status != StatusCode::kOk) return status;
    if (operation == HhtRotationOperation::kLinearize) {
      status = InverseLeftRotationJacobian(theta.get(), inverse.get(), d);
      if (status != StatusCode::kOk) return status;
      const Vector* rotated[] = {qb.get(), qc.get()};
      Matrix* outputs[] = {out->velocity_per_rotation,
                           out->acceleration_per_rotation};
      const double factors[] = {c1, c0}, signs[] = {-1, 1};
      for (int quantity = 0; quantity < 2; ++quantity) {
        status = Skew(rotated[quantity], skew.get(), d);
        if (status != StatusCode::kOk) return status;
        status = derivative->CopyFrom(inverse.get(), d);
        if (status != StatusCode::kOk) return status;
        status = derivative->Scale(factors[quantity], d);
        if (status != StatusCode::kOk) return status;
        status = derivative->AddScaled(signs[quantity], skew.get(), d);
        if (status != StatusCode::kOk) return status;
        status = derivative->GatherTo(all, 9, block.get(), 0, d);
        if (status != StatusCode::kOk) return status;
        status = outputs[quantity]->SetRow(node, block.get(), d);
        if (status != StatusCode::kOk) return status;
      }
    } else {
      status = trial_v->CopyFrom(theta.get(), d);
      if (status != StatusCode::kOk) return status;
      status = trial_v->Scale(c1, d);
      if (status != StatusCode::kOk) return status;
      status = trial_v->AddScaled(1, qb.get(), d);
      if (status != StatusCode::kOk) return status;
      status = trial_a->CopyFrom(theta.get(), d);
      if (status != StatusCode::kOk) return status;
      status = trial_a->Scale(c0, d);
      if (status != StatusCode::kOk) return status;
      status = trial_a->AddScaled(-1, qc.get(), d);
      if (status != StatusCode::kOk) return status;
      if (operation == HhtRotationOperation::kPredict) {
        status = candidate->orientations->SetRow(node, target_q.get(), d);
        if (status != StatusCode::kOk) return status;
      }
      status = candidate->velocity->GatherFrom(trial_v.get(), identity, 3,
                                               3 * node, d);
      if (status != StatusCode::kOk) return status;
      status = candidate->acceleration->GatherFrom(trial_a.get(), identity, 3,
                                                   3 * node, d);
      if (status != StatusCode::kOk) return status;
    }
  }
  return StatusCode::kOk;
}
}  // namespace fesa
