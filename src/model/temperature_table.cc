#include "fesa/model/temperature_table.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <new>

#include "model/model_storage.h"

namespace fesa {
TemperatureTable::TemperatureTable(Vector* values, Vector* temperatures)
    : values_(values), temperatures_(temperatures) {}
TemperatureTable::~TemperatureTable() {
  delete values_;
  delete temperatures_;
}
StatusCode TemperatureTable::Create(const Vector* values,
                                    const Vector* temperatures,
                                    TemperatureTable** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("TemperatureTable output", d);
  *out = nullptr;
  if (!values || !temperatures || values->Size() != temperatures->Size())
    return model_storage::Invalid("TemperatureTable dimensions", d);
  for (int i = 0; i < values->Size(); ++i) {
    if (!std::isfinite(values->Data()[i]) ||
        !std::isfinite(temperatures->Data()[i]) ||
        (i > 0 && temperatures->Data()[i] <= temperatures->Data()[i - 1]))
      return model_storage::Invalid("TemperatureTable finite increasing knots",
                                    d);
  }
  std::unique_ptr<TemperatureTable> result(
      new (std::nothrow) TemperatureTable(nullptr, nullptr));
  if (!result) return model_storage::AllocationFailure("TemperatureTable", d);
  auto status = Vector::Create(values->Size(), &result->values_, d);
  if (status == StatusCode::kOk)
    status = Vector::Create(values->Size(), &result->temperatures_, d);
  if (status == StatusCode::kOk) status = result->values_->CopyFrom(values, d);
  if (status == StatusCode::kOk)
    status = result->temperatures_->CopyFrom(temperatures, d);
  if (status != StatusCode::kOk) return status;
  *out = result.release();
  return StatusCode::kOk;
}
StatusCode TemperatureTable::Clone(TemperatureTable** out,
                                   Diagnostic* d) const {
  return Create(values_, temperatures_, out, d);
}
StatusCode TemperatureTable::Evaluate(double theta, double* value,
                                      double* slope, Diagnostic* d) const {
  if (!value || !slope || value == slope || !std::isfinite(theta))
    return model_storage::Invalid("TemperatureTable evaluation", d);
  const int count = values_->Size();
  const double* x = temperatures_->Data();
  const double* y = values_->Data();
  double evaluated = y[0], derivative = 0;
  if (theta >= x[count - 1]) {
    evaluated = y[count - 1];
  } else if (count > 1 && theta >= x[0]) {
    const int right =
        static_cast<int>(std::upper_bound(x, x + count, theta) - x);
    const int left = right - 1;
    const double span = x[right] - x[left];
    derivative = (y[right] - y[left]) / span;
    evaluated = y[left] + (y[right] - y[left]) * ((theta - x[left]) / span);
    if (!std::isfinite(span) || !std::isfinite(derivative) ||
        !std::isfinite(evaluated))
      return ReportFailure(StatusCode::kLibraryFailure,
                           DiagnosticStage::kComputation,
                           "TemperatureTable interpolation range", d);
  }
  *value = evaluated;
  *slope = derivative;
  return StatusCode::kOk;
}
StatusCode TemperatureTable::Integrate(double from, double to, double* out,
                                       Diagnostic* d) const {
  if (!out || !std::isfinite(from) || !std::isfinite(to))
    return model_storage::Invalid("TemperatureTable integral bounds", d);
  if (from == to) {
    *out = 0;
    return StatusCode::kOk;
  }
  const double lo = std::min(from, to), hi = std::max(from, to);
  const double sign = to > from ? 1 : -1;
  const int count = values_->Size();
  const double* x = temperatures_->Data();
  Vector* raw = nullptr;
  auto status = Vector::Create(count, &raw, d);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<Vector> weights(raw);
  for (int i = 0; i < count; ++i) {
    // Scalar integrals of the piecewise-linear interpolation basis. Each
    // coefficient is defined once; the weighted sum is delegated to BLAS.
    double weight = 0;
    if (i == 0 && lo < x[0]) weight += std::min(hi, x[0]) - lo;
    if (i == count - 1 && hi > x[i]) weight += hi - std::max(lo, x[i]);
    if (i > 0) {
      const double left = std::max(lo, x[i - 1]);
      const double right = std::min(hi, x[i]);
      if (right > left)
        weight += .5 * (right - left) *
                  ((left - x[i - 1]) / (x[i] - x[i - 1]) +
                   (right - x[i - 1]) / (x[i] - x[i - 1]));
    }
    if (i + 1 < count) {
      const double left = std::max(lo, x[i]);
      const double right = std::min(hi, x[i + 1]);
      if (right > left)
        weight += .5 * (right - left) *
                  ((x[i + 1] - left) / (x[i + 1] - x[i]) +
                   (x[i + 1] - right) / (x[i + 1] - x[i]));
    }
    if (!std::isfinite(weight))
      return ReportFailure(StatusCode::kLibraryFailure,
                           DiagnosticStage::kComputation,
                           "TemperatureTable integral coefficient range", d);
    status = weights->Set(i, sign * weight, d);
    if (status != StatusCode::kOk) return status;
  }
  return weights->Dot(values_, out, d);
}
bool TemperatureTable::IsConstant() const { return values_->Size() == 1; }
bool TemperatureTable::HasPositiveValues() const {
  for (int i = 0; i < values_->Size(); ++i)
    if (values_->Data()[i] <= 0) return false;
  return true;
}
bool TemperatureTable::HasNonnegativeValues() const {
  for (int i = 0; i < values_->Size(); ++i)
    if (values_->Data()[i] < 0) return false;
  return true;
}
bool TemperatureTable::IsZero() const {
  for (int i = 0; i < values_->Size(); ++i)
    if (values_->Data()[i] != 0) return false;
  return true;
}
}  // namespace fesa
