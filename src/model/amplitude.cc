#include "fesa/model/amplitude.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <new>

#include "model/model_storage.h"

namespace fesa {
namespace {
StatusCode CheckAmplitudeContext(const AmplitudeEvaluationContext* context,
                                 double* value, Diagnostic* d) {
  if (!context || !value || value == &context->definition_period ||
      value == &context->step_time || value == &context->total_time ||
      context->definition_step < 1 ||
      context->current_step < context->definition_step ||
      !std::isfinite(context->definition_period) ||
      context->definition_period <= 0 || !std::isfinite(context->step_time) ||
      context->step_time < 0 || !std::isfinite(context->total_time) ||
      context->total_time < context->step_time ||
      (context->current_step == context->definition_step &&
       context->step_time > context->definition_period))
    return model_storage::Invalid("Amplitude step context", d);
  return StatusCode::kOk;
}
int PairExponent(double a, double b) {
  int ea = 0, eb = 0;
  std::frexp(a, &ea);
  std::frexp(b, &eb);
  return a == 0 ? eb : b == 0 ? ea : std::max(ea, eb);
}
}  // namespace
Amplitude::~Amplitude() {
  delete[] name_;
  delete times_;
  delete values_;
}
StatusCode Amplitude::Create(int id, const char* name, AmplitudeTimeBasis basis,
                             const Vector* times, const Vector* values,
                             Amplitude** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("Amplitude output", d);
  *out = nullptr;
  if (!model_storage::ValidId(id) || !name || !*name || !times || !values ||
      times->Size() != values->Size() ||
      (basis != AmplitudeTimeBasis::kStepTime &&
       basis != AmplitudeTimeBasis::kTotalTime))
    return model_storage::Invalid("Amplitude definition", d);
  for (int i = 0; i < times->Size(); ++i) {
    if (!std::isfinite(times->Data()[i]) || !std::isfinite(values->Data()[i]) ||
        (i > 0 && times->Data()[i] <= times->Data()[i - 1]))
      return model_storage::Invalid("Amplitude finite increasing table", d);
  }
  std::unique_ptr<Amplitude> result(new (std::nothrow) Amplitude);
  if (!result) return model_storage::AllocationFailure("Amplitude", d);
  result->id_ = id;
  result->basis_ = basis;
  result->name_ = model_storage::CopyText(name);
  if (!result->name_)
    return model_storage::AllocationFailure("Amplitude name", d);
  auto status = Vector::Create(times->Size(), &result->times_, d);
  if (status == StatusCode::kOk)
    status = Vector::Create(values->Size(), &result->values_, d);
  if (status == StatusCode::kOk) status = result->times_->CopyFrom(times, d);
  if (status == StatusCode::kOk) status = result->values_->CopyFrom(values, d);
  if (status != StatusCode::kOk) return status;
  *out = result.release();
  return StatusCode::kOk;
}
int Amplitude::Id() const { return id_; }
const char* Amplitude::Name() const { return name_; }
AmplitudeTimeBasis Amplitude::TimeBasis() const { return basis_; }
int Amplitude::PointCount() const { return times_->Size(); }
StatusCode Amplitude::EvaluateInStep(const AmplitudeEvaluationContext* context,
                                     double* value, Diagnostic* d) const {
  const auto status = CheckAmplitudeContext(context, value, d);
  if (status != StatusCode::kOk) return status;
  const double time = basis_ == AmplitudeTimeBasis::kTotalTime
                          ? context->total_time
                          : (context->current_step == context->definition_step
                                 ? context->step_time
                                 : context->definition_period);
  return Evaluate(time, value, d);
}
StatusCode Amplitude::EvaluateRightDerivative(double time, double* rate,
                                              Diagnostic* d) const {
  if (!rate || !std::isfinite(time))
    return model_storage::Invalid("Amplitude right derivative", d);
  const double* times = times_->Data();
  const double* values = values_->Data();
  const int count = times_->Size();
  if (time < times[0] || time >= times[count - 1]) {
    *rate = 0;
    return StatusCode::kOk;
  }
  const int right =
      static_cast<int>(std::upper_bound(times, times + count, time) - times);
  const int left = right - 1;
  const int ey = PairExponent(values[left], values[right]);
  const int et = PairExponent(times[left], times[right]);
  // A power-of-two normalization keeps both endpoint differences finite,
  // including subnormal endpoints. Apply the exponent ratio only once.
  const double difference =
      std::ldexp(values[right], -ey) - std::ldexp(values[left], -ey);
  const double span =
      std::ldexp(times[right], -et) - std::ldexp(times[left], -et);
  const double slope = std::ldexp(difference / span, ey - et);
  if (!std::isfinite(slope) || (difference != 0 && slope == 0))
    return ReportFailure(StatusCode::kLibraryFailure,
                         DiagnosticStage::kComputation,
                         "Amplitude derivative representable range", d);
  *rate = slope;
  return StatusCode::kOk;
}
StatusCode Amplitude::EvaluateRightDerivativeInStep(
    const AmplitudeEvaluationContext* context, double* rate,
    Diagnostic* d) const {
  const auto status = CheckAmplitudeContext(context, rate, d);
  if (status != StatusCode::kOk) return status;
  if (basis_ == AmplitudeTimeBasis::kStepTime &&
      (context->current_step > context->definition_step ||
       context->step_time == context->definition_period)) {
    *rate = 0;
    return StatusCode::kOk;
  }
  return EvaluateRightDerivative(basis_ == AmplitudeTimeBasis::kTotalTime
                                     ? context->total_time
                                     : context->step_time,
                                 rate, d);
}
StatusCode Amplitude::EvaluateRange(double start, double end, double* minimum,
                                    double* maximum, Diagnostic* d) const {
  if (!minimum || !maximum || minimum == maximum || !std::isfinite(start) ||
      !std::isfinite(end) || start > end)
    return model_storage::Invalid("Amplitude evaluation range", d);
  double first = 0, last = 0;
  auto status = Evaluate(start, &first, d);
  if (status == StatusCode::kOk) status = Evaluate(end, &last, d);
  if (status != StatusCode::kOk) return status;
  double low = std::min(first, last), high = std::max(first, last);
  for (int i = 0; i < times_->Size(); ++i) {
    if (times_->Data()[i] <= start || times_->Data()[i] >= end) continue;
    low = std::min(low, values_->Data()[i]);
    high = std::max(high, values_->Data()[i]);
  }
  *minimum = low;
  *maximum = high;
  return StatusCode::kOk;
}
StatusCode Amplitude::Evaluate(double time, double* value,
                               Diagnostic* d) const {
  if (!value || !std::isfinite(time))
    return model_storage::Invalid("Amplitude evaluation", d);
  const double* times = times_->Data();
  const double* values = values_->Data();
  const int last = times_->Size() - 1;
  if (time <= times[0]) {
    *value = values[0];
    return StatusCode::kOk;
  }
  if (time >= times[last]) {
    *value = values[last];
    return StatusCode::kOk;
  }
  const int right =
      static_cast<int>(std::lower_bound(times, times + last + 1, time) - times);
  if (time == times[right]) {
    *value = values[right];
    return StatusCode::kOk;
  }
  const int left = right - 1;
  const double span = times[right] - times[left];
  // Compute each endpoint distance directly: 1-weight can lose a small
  // distance when the other weight rounds to one near a knot.
  const double left_weight = std::isfinite(span)
      ? (times[right] - time) / span
      : (times[right] * .5 - time * .5) /
            (times[right] * .5 - times[left] * .5);
  const double right_weight = std::isfinite(span)
      ? (time - times[left]) / span
      : (time * .5 - times[left] * .5) /
            (times[right] * .5 - times[left] * .5);
  const double a = values[left], b = values[right];
  if (std::signbit(a) != std::signbit(b)) {
    // Opposite signs need a convex combination to avoid difference overflow.
    *value = left_weight * a + right_weight * b;
  } else if (right_weight <= .5) {
    *value = a + right_weight * (b - a);
  } else {
    *value = b + left_weight * (a - b);
  }
  return StatusCode::kOk;
}
}  // namespace fesa
