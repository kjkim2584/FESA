#include "fesa/analysis/dynamic_step_settings.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <new>

#include "analysis/direct_increment.h"

namespace fesa {
DynamicStepSettings::~DynamicStepSettings() {
  delete parameters_;
  delete hht_;
  delete source_file_;
}
StatusCode DynamicStepSettings::CreateDirect(
    const DynamicStepParameters* parameters, const char* source_file,
    int source_line, int step, DynamicStepSettings** out, Diagnostic* d) {
  return Create(parameters, source_file, source_line, step, false, out, d);
}
StatusCode DynamicStepSettings::CreateAutomatic(
    const DynamicStepParameters* parameters, const char* source_file,
    int source_line, int step, DynamicStepSettings** out, Diagnostic* d) {
  return Create(parameters, source_file, source_line, step, true, out, d);
}
StatusCode DynamicStepSettings::Create(
    const DynamicStepParameters* parameters, const char* source_file,
    int source_line, int step, bool automatic, DynamicStepSettings** out,
    Diagnostic* d) {
  if (out) *out = nullptr;
  if (!out || !parameters || !source_file || source_line <= 0 || step <= 0 ||
      parameters->increment_limit <= 0 || !std::isfinite(parameters->alpha) ||
      parameters->alpha < -1.0 / 3.0 || parameters->alpha > 0)
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation, "Dynamic step settings",
                         d);
  double nominal = 0;
  const auto status = ResolveDirectIncrement(
      parameters->initial_increment, parameters->period, parameters->minimum,
      parameters->maximum, &nominal, d);
  if (status != StatusCode::kOk) return status;
  double minimum = parameters->minimum, maximum = parameters->maximum;
  if (automatic) {
    if (minimum == 0)
      minimum = std::min(parameters->initial_increment,
                         parameters->period * 1e-5);
    if (maximum == 0)
      maximum = std::max(parameters->initial_increment,
                         parameters->period / 100);
    if (minimum <= 0 || minimum > maximum)
      return ReportFailure(StatusCode::kInputError,
                           DiagnosticStage::kValidation,
                           "Automatic unrepresentable/inconsistent limits", d);
    // Preserve the logical request when T<minimum. Only its endpoint duration
    // is clipped by the controller, not the minimum-valid nominal metadata.
    nominal = std::min(parameters->initial_increment, maximum);
  }
  try {
    std::unique_ptr<DynamicStepSettings> result(new DynamicStepSettings);
    result->parameters_ = new DynamicStepParameters;
    auto* copy = result->parameters_;
    copy->initial_increment = parameters->initial_increment;
    copy->period = parameters->period;
    copy->minimum = parameters->minimum;
    copy->maximum = parameters->maximum;
    copy->alpha = parameters->alpha;
    copy->increment_limit = parameters->increment_limit;
    copy->minimum_specified = parameters->minimum_specified;
    copy->maximum_specified = parameters->maximum_specified;
    copy->alpha_specified = parameters->alpha_specified;
    copy->increment_limit_specified = parameters->increment_limit_specified;
    result->hht_ = new HhtParameters;
    result->hht_->alpha = parameters->alpha;
    result->source_file_ = new std::string(source_file);
    result->source_line_ = source_line;
    result->step_ = step;
    result->nominal_ = nominal;
    result->minimum_ = minimum;
    result->maximum_ = maximum;
    result->automatic_ = automatic;
    *out = result.release();
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return ReportFailure(StatusCode::kAllocationFailure,
                         DiagnosticStage::kAllocation, "Dynamic step settings",
                         d);
  }
}
const DynamicStepParameters* DynamicStepSettings::Parameters() const {
  return parameters_;
}
const HhtParameters* DynamicStepSettings::IntegratorParameters() const {
  return hht_;
}
double DynamicStepSettings::NominalIncrement() const { return nominal_; }
bool DynamicStepSettings::IsAutomatic() const { return automatic_; }
double DynamicStepSettings::MinimumIncrement() const { return minimum_; }
double DynamicStepSettings::MaximumIncrement() const { return maximum_; }
const char* DynamicStepSettings::SourceFile() const {
  return source_file_->c_str();
}
int DynamicStepSettings::SourceLine() const { return source_line_; }
int DynamicStepSettings::Step() const { return step_; }
}  // namespace fesa
