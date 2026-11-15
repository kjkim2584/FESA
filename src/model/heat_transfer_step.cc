#include "fesa/model/heat_transfer_step.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <new>

#include "model/model_storage.h"

namespace fesa {
double HeatTransferStepParameters::EffectiveMinimum() const {
  return std::min(.8 * initial_increment,
                  has_minimum ? minimum : 1e-5 * period);
}
bool HeatTransferStepParameters::HasMaximum() const { return has_maximum; }
ThermalStepAmplitude HeatTransferStepParameters::DefaultAmplitude() const {
  return has_default_amplitude ? default_amplitude
                               : (steady ? ThermalStepAmplitude::kRamp
                                         : ThermalStepAmplitude::kStep);
}
StatusCode ValidateHeatTransferStepParameters(
    const HeatTransferStepParameters* p, Diagnostic* d) {
  if (!p) return model_storage::Invalid("Heat transfer parameters", d);
  if (!std::isfinite(p->period) || p->period <= 0 ||
      !std::isfinite(p->initial_increment) || p->initial_increment <= 0 ||
      (p->has_minimum && (!std::isfinite(p->minimum) || p->minimum <= 0)) ||
      (p->has_maximum && (!std::isfinite(p->maximum) || p->maximum <= 0)) ||
      (p->has_deltmx &&
       (p->steady || !std::isfinite(p->deltmx) || p->deltmx <= 0)) ||
      p->max_increments <= 0 || !(p->EffectiveMinimum() > 0) ||
      (p->has_maximum && p->maximum < p->EffectiveMinimum()) ||
      (p->has_default_amplitude &&
       p->default_amplitude != ThermalStepAmplitude::kStep &&
       p->default_amplitude != ThermalStepAmplitude::kRamp))
    return ReportFailure(StatusCode::kInputError, DiagnosticStage::kValidation,
                         "Invalid heat transfer time/limit settings", d);
  return StatusCode::kOk;
}
StatusCode CloneHeatTransferStepParameters(const HeatTransferStepParameters* p,
                                           HeatTransferStepParameters** out,
                                           Diagnostic* d) {
  if (!out) return model_storage::Invalid("Heat transfer parameter output", d);
  *out = nullptr;
  const auto status = ValidateHeatTransferStepParameters(p, d);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<HeatTransferStepParameters> copy(
      new (std::nothrow) HeatTransferStepParameters);
  if (!copy)
    return model_storage::AllocationFailure("Heat transfer settings", d);
  copy->steady = p->steady;
  copy->initial_increment = p->initial_increment;
  copy->period = p->period;
  copy->has_minimum = p->has_minimum;
  copy->minimum = p->minimum;
  copy->has_maximum = p->has_maximum;
  copy->maximum = p->maximum;
  copy->has_deltmx = p->has_deltmx;
  copy->deltmx = p->deltmx;
  copy->has_default_amplitude = p->has_default_amplitude;
  copy->default_amplitude = p->default_amplitude;
  copy->max_increments = p->max_increments;
  *out = copy.release();
  return StatusCode::kOk;
}
}  // namespace fesa
