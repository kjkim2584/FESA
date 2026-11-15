#include "analysis/finite_shell_initialization.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>

#include "fesa/math/matrix.h"
#include "model/model_storage.h"

namespace fesa {
namespace {
// Coefficient magnitudes, not physical vector storage. Products stay factored
// until their common binary exponent is removed from a dimensionless ratio.
struct Magnitude {
  double mantissa = 0;
  int exponent = 0;
};
void MultiplyMagnitude(const Magnitude* a, const Magnitude* b, Magnitude* out) {
  out->mantissa = a->mantissa * b->mantissa;
  out->exponent = a->exponent + b->exponent;
}
StatusCode NewVector(std::unique_ptr<Vector>* owner, Diagnostic* d) {
  Vector* raw = nullptr;
  const auto status = Vector::Create(3, &raw, d);
  owner->reset(raw);
  return status;
}
StatusCode Normalize(const Vector* source, Vector* normalized, Magnitude* scale,
                     double* norm, Diagnostic* d) {
  auto status = normalized->CopyFrom(source, d);
  if (status != StatusCode::kOk) return status;
  double maximum = 0;
  int index = 0;
  status = source->MaximumAbs(&maximum, &index, d);
  if (status != StatusCode::kOk) return status;
  *norm = 0;
  if (maximum == 0) return StatusCode::kOk;
  scale->mantissa = std::frexp(maximum, &scale->exponent);
  int remaining = -scale->exponent;
  while (remaining != 0) {
    const int shift = std::max(-512, std::min(512, remaining));
    status = normalized->Scale(std::scalbn(1., shift), d);
    if (status != StatusCode::kOk) return status;
    remaining -= shift;
  }
  status = normalized->Scale(1 / scale->mantissa, d);
  return status == StatusCode::kOk ? normalized->Norm2(norm, d) : status;
}
StatusCode Weights(const Magnitude* scales, int count, double* out,
                   Diagnostic* d) {
  int common = std::numeric_limits<int>::min();
  for (int i = 0; i < count; ++i)
    if (scales[i].mantissa != 0) common = std::max(common, scales[i].exponent);
  for (int i = 0; i < count; ++i) {
    out[i] = 0;
    if (scales[i].mantissa == 0) continue;
    out[i] = std::scalbn(scales[i].mantissa, scales[i].exponent - common);
    if (out[i] == 0)
      return ReportFailure(StatusCode::kLibraryFailure,
                           DiagnosticStage::kComputation,
                           "Finite shell compatibility scale range", d);
  }
  return StatusCode::kOk;
}
StatusCode Cross(const Vector* director, const Vector* omega, Vector* out,
                 Diagnostic* d) {
  Matrix* raw = nullptr;
  auto status = Matrix::Create(3, 3, &raw, d);
  std::unique_ptr<Matrix> skew(raw);
  if (status != StatusCode::kOk) return status;
  const double* x = director->Data();
  const double data[] = {0, -x[2], x[1], x[2], 0, -x[0], -x[1], x[0], 0};
  for (int i = 0; i < 9; ++i) {
    status = raw->Set(i / 3, i % 3, data[i], d);
    if (status != StatusCode::kOk) return status;
  }
  return Matrix::MultiplyVector(raw, false, omega, out, d);
}
}  // namespace
StatusCode EvaluateFiniteShellNodeCompatibility(
    const Vector* initial_normal, const Vector* director, double kappa,
    const Vector* body_spin, const Vector* omega, const Vector* moment,
    const Vector* moment_rate, double* configuration_error,
    double* velocity_error, Diagnostic* d) {
  if (!configuration_error || !velocity_error ||
      configuration_error == velocity_error || !std::isfinite(kappa) ||
      kappa <= 0)
    return model_storage::Invalid("Finite shell compatibility scalar arguments",
                                  d);
  for (const Vector* value :
       {initial_normal, director, body_spin, omega, moment, moment_rate}) {
    if (!value || value->Size() != 3)
      return model_storage::Invalid("Finite shell compatibility vector shape",
                                    d);
    for (int i = 0; i < 3; ++i)
      if (!std::isfinite(value->Data()[i]))
        return model_storage::Invalid("Finite shell compatibility finite data",
                                      d);
  }
  std::unique_ptr<Vector> normalized[4];
  const Vector* fields[] = {body_spin, omega, moment, moment_rate};
  Magnitude scales[4];
  double norms[4] = {};
  for (int i = 0; i < 4; ++i) {
    auto status = NewVector(&normalized[i], d);
    if (status != StatusCode::kOk) return status;
    status =
        Normalize(fields[i], normalized[i].get(), &scales[i], &norms[i], d);
    if (status != StatusCode::kOk) return status;
  }
  Magnitude spring;
  spring.mantissa = std::frexp(kappa, &spring.exponent);
  Magnitude configuration[2], velocity[3];
  MultiplyMagnitude(&spring, &scales[0], &configuration[0]);
  configuration[1].mantissa = scales[2].mantissa;
  configuration[1].exponent = scales[2].exponent;
  MultiplyMagnitude(&spring, &scales[1], &velocity[0]);
  MultiplyMagnitude(&scales[2], &scales[1], &velocity[1]);
  velocity[2].mantissa = scales[3].mantissa;
  velocity[2].exponent = scales[3].exponent;
  double q_weights[2] = {}, v_weights[3] = {};
  auto status = Weights(configuration, 2, q_weights, d);
  if (status == StatusCode::kOk) status = Weights(velocity, 3, v_weights, d);
  if (status != StatusCode::kOk) return status;
  double spin_projection = 0, moment_projection = 0, omega_projection = 0,
         rate_projection = 0, geometric = 0;
  status = normalized[0]->Dot(initial_normal, &spin_projection, d);
  if (status == StatusCode::kOk)
    status = normalized[2]->Dot(director, &moment_projection, d);
  if (status == StatusCode::kOk)
    status = normalized[1]->Dot(director, &omega_projection, d);
  if (status == StatusCode::kOk)
    status = normalized[3]->Dot(director, &rate_projection, d);
  std::unique_ptr<Vector> cross;
  if (status == StatusCode::kOk) status = NewVector(&cross, d);
  if (status == StatusCode::kOk)
    status = Cross(director, normalized[1].get(), cross.get(), d);
  if (status == StatusCode::kOk)
    status = normalized[2]->Dot(cross.get(), &geometric, d);
  if (status != StatusCode::kOk) return status;
  const double q_den = q_weights[0] * norms[0] + q_weights[1] * norms[2];
  const double v_den = v_weights[0] * norms[1] +
                       v_weights[1] * norms[2] * norms[1] +
                       v_weights[2] * norms[3];
  const double q = q_den == 0 ? 0
                              : std::abs(q_weights[0] * spin_projection -
                                         q_weights[1] * moment_projection) /
                                    q_den;
  const double v = v_den == 0 ? 0
                              : std::abs(v_weights[0] * omega_projection +
                                         v_weights[1] * geometric -
                                         v_weights[2] * rate_projection) /
                                    v_den;
  if (!std::isfinite(q) || !std::isfinite(v))
    return ReportFailure(StatusCode::kLibraryFailure,
                         DiagnosticStage::kComputation,
                         "Finite shell compatibility ratio range", d);
  *configuration_error = q;
  *velocity_error = v;
  return StatusCode::kOk;
}
StatusCode CheckFiniteShellNodeCompatibility(
    int node_id, const Vector* initial_normal, const Vector* director,
    double kappa, const Vector* body_spin, const Vector* omega,
    const Vector* moment, const Vector* moment_rate, Diagnostic* d) {
  if (node_id <= 0)
    return model_storage::Invalid("Finite shell compatibility node", d);
  double q = 0, v = 0;
  auto status = EvaluateFiniteShellNodeCompatibility(
      initial_normal, director, kappa, body_spin, omega, moment, moment_rate,
      &q, &v, d);
  if (status != StatusCode::kOk) return status;
  if (q <= 1e-10 && v <= 1e-10) return StatusCode::kOk;
  status = ReportFailure(
      StatusCode::kEquilibriumFailure, DiagnosticStage::kComputation,
      q > 1e-10 ? "Finite shell initial configuration compatibility"
                : "Finite shell initial velocity compatibility",
      d);
  if (d) d->node_id = node_id;
  return status;
}

StatusCode ComputeFiniteShellNodeNullAcceleration(
    const Vector* director, double kappa, const Vector* omega,
    const Vector* moment, const Vector* moment_rate,
    const Vector* moment_second_rate, const Vector* range_acceleration,
    double* gamma, Diagnostic* d) {
  if (!gamma || !std::isfinite(kappa) || kappa <= 0)
    return model_storage::Invalid("Finite shell N8 scalar arguments", d);
  for (const Vector* value : {director, omega, moment, moment_rate,
                              moment_second_rate, range_acceleration}) {
    if (!value || value->Size() != 3)
      return model_storage::Invalid("Finite shell N8 vector shape", d);
    for (int i = 0; i < 3; ++i)
      if (!std::isfinite(value->Data()[i]))
        return model_storage::Invalid("Finite shell N8 finite data", d);
  }
  const Vector* fields[] = {omega, moment, moment_rate, moment_second_rate,
                            range_acceleration};
  std::unique_ptr<Vector> normalized[5];
  Magnitude scales[5];
  for (int i = 0; i < 5; ++i) {
    auto status = NewVector(&normalized[i], d);
    if (status != StatusCode::kOk) return status;
    double norm = 0;
    status = Normalize(fields[i], normalized[i].get(), &scales[i], &norm, d);
    if (status != StatusCode::kOk) return status;
  }
  Magnitude spring, terms[5], squared;
  spring.mantissa = std::frexp(kappa, &spring.exponent);
  MultiplyMagnitude(&scales[0], &scales[0], &squared);
  MultiplyMagnitude(&squared, &scales[1], &terms[0]);
  MultiplyMagnitude(&scales[0], &scales[2], &terms[1]);
  terms[1].mantissa *= 2;
  terms[2].mantissa = scales[3].mantissa;
  terms[2].exponent = scales[3].exponent;
  MultiplyMagnitude(&spring, &scales[4], &terms[3]);
  MultiplyMagnitude(&scales[1], &scales[4], &terms[4]);
  double weights[5] = {}, values[5] = {};
  auto status = Weights(terms, 5, weights, d);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<Vector> rate, second, moment_cross;
  for (auto* owner : {&rate, &second, &moment_cross}) {
    status = NewVector(owner, d);
    if (status != StatusCode::kOk) return status;
  }
  status = Cross(normalized[0].get(), director, rate.get(), d);
  if (status != StatusCode::kOk) return status;
  status = Cross(normalized[0].get(), rate.get(), second.get(), d);
  if (status != StatusCode::kOk) return status;
  status = normalized[1]->Dot(second.get(), &values[0], d);
  if (status == StatusCode::kOk)
    status = normalized[2]->Dot(rate.get(), &values[1], d);
  if (status == StatusCode::kOk)
    status = normalized[3]->Dot(director, &values[2], d);
  if (status == StatusCode::kOk)
    status = normalized[4]->Dot(director, &values[3], d);
  if (status == StatusCode::kOk)
    status = Cross(director, normalized[1].get(), moment_cross.get(), d);
  if (status == StatusCode::kOk)
    status = normalized[4]->Dot(moment_cross.get(), &values[4], d);
  if (status != StatusCode::kOk) return status;
  const double numerator = weights[0] * values[0] + weights[1] * values[1] +
                           weights[2] * values[2] - weights[3] * values[3] +
                           weights[4] * values[4];
  double value = 0;
  if (numerator != 0) {
    int common = std::numeric_limits<int>::min();
    for (int i = 0; i < 5; ++i)
      if (terms[i].mantissa != 0) common = std::max(common, terms[i].exponent);
    value = std::scalbn(numerator / spring.mantissa, common - spring.exponent);
    if (!std::isfinite(value) || value == 0)
      return ReportFailure(StatusCode::kLibraryFailure,
                           DiagnosticStage::kComputation,
                           "Finite shell N8 acceleration range", d);
  }
  *gamma = value;
  return StatusCode::kOk;
}
}  // namespace fesa
