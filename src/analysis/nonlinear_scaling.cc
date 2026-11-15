#include "analysis/nonlinear_scaling.h"

#include <cmath>

namespace fesa {
namespace {

bool Nonnegative(double value) {
  return std::isfinite(value) && value >= 0;
}

bool ValidField(const FieldResidual* field) {
  return field && std::isfinite(field->raw_mean) && field->raw_mean > 0 &&
         Nonnegative(field->spatial_mean);
}

// Implements nonlinear_scaling_contract.md without materializing an irrelevant
// candidate. The finite positive inputs have bounded binary exponents; products
// and ratios of their mantissas stay safely within the normal double range.
bool ComputeMean(const FieldResidual* target, double other_raw,
                 double coefficient, double length, bool divide_length,
                 double* used) {
  *used = target->raw_mean;
  if (coefficient == 0 || target->spatial_mean > target->raw_mean) return true;

  int coefficient_exponent = 0;
  int raw_exponent = 0;
  int length_exponent = 0;
  const double coefficient_mantissa =
      std::frexp(coefficient, &coefficient_exponent);
  const double raw_mantissa = std::frexp(other_raw, &raw_exponent);
  const double length_mantissa = std::frexp(length, &length_exponent);
  const double product = coefficient_mantissa * raw_mantissa;
  int normalization_exponent = 0;
  const double candidate_mantissa = std::frexp(
      divide_length ? product / length_mantissa : product * length_mantissa,
      &normalization_exponent);
  const int candidate_exponent =
      coefficient_exponent + raw_exponent + normalization_exponent +
      (divide_length ? -length_exponent : length_exponent);

  int target_exponent = 0;
  const double target_mantissa =
      std::frexp(target->raw_mean, &target_exponent);
  if (candidate_exponent < target_exponent ||
      (candidate_exponent == target_exponent &&
       candidate_mantissa <= target_mantissa))
    return true;

  const double candidate = std::scalbn(candidate_mantissa, candidate_exponent);
  if (!std::isfinite(candidate) || !(candidate > 0)) return false;
  // A larger normalized value can round back to raw in the subnormal range.
  if (candidate > target->raw_mean) *used = candidate;
  return true;
}

}  // namespace

StatusCode ConvertNonlinearFieldMeans(const ControlSettings* controls,
                                       double characteristic_length,
                                       FieldResidual* force,
                                       FieldResidual* moment, Diagnostic* d) {
  if (!controls || !controls->displacement || !controls->rotation ||
      controls->displacement == controls->rotation || !ValidField(force) ||
      !ValidField(moment) || force == moment ||
      !Nonnegative(controls->displacement->cf) ||
      !Nonnegative(controls->rotation->cf) ||
      !std::isfinite(characteristic_length) || !(characteristic_length > 0))
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation,
                         "Nonlinear field mean conversion inputs", d);

  double force_used = 0;
  double moment_used = 0;
  const double force_cf =
      controls->displacement->qu_explicit ? 0 : controls->displacement->cf;
  const double moment_cf =
      controls->rotation->qu_explicit ? 0 : controls->rotation->cf;
  if (!ComputeMean(force, moment->raw_mean, force_cf,
                   characteristic_length, true, &force_used))
    return ReportFailure(StatusCode::kLibraryFailure,
                         DiagnosticStage::kComputation,
                         "Nonlinear force mean conversion", d);
  if (!ComputeMean(moment, force->raw_mean, moment_cf,
                   characteristic_length, false, &moment_used))
    return ReportFailure(StatusCode::kLibraryFailure,
                         DiagnosticStage::kComputation,
                         "Nonlinear moment mean conversion", d);

  // Publish only after both simultaneous RAW-based candidates are valid.
  force->used_mean = force_used;
  moment->used_mean = moment_used;
  force->converted = force_used > force->raw_mean;
  moment->converted = moment_used > moment->raw_mean;
  return StatusCode::kOk;
}

}  // namespace fesa
