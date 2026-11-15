#include "fesa/model/thermal_physical_constants.h"

#include <cmath>
#include <memory>
#include <new>

#include "model/model_storage.h"

namespace fesa {
StatusCode ValidateThermalPhysicalConstants(const ThermalPhysicalConstants* p,
                                            Diagnostic* d) {
  if (!p || p->source_line <= 0 ||
      (p->has_absolute_zero && !std::isfinite(p->absolute_zero)) ||
      (p->has_stefan_boltzmann &&
       (!std::isfinite(p->stefan_boltzmann) || p->stefan_boltzmann <= 0)))
    return model_storage::Invalid("Thermal physical constants", d);
  return StatusCode::kOk;
}
StatusCode CloneThermalPhysicalConstants(const ThermalPhysicalConstants* p,
                                         ThermalPhysicalConstants** out,
                                         Diagnostic* d) {
  if (!out) return model_storage::Invalid("Thermal constants output", d);
  *out = nullptr;
  const auto status = ValidateThermalPhysicalConstants(p, d);
  if (status != StatusCode::kOk) return status;
  auto* copy = new (std::nothrow) ThermalPhysicalConstants;
  if (!copy) return model_storage::AllocationFailure("Thermal constants", d);
  copy->has_absolute_zero = p->has_absolute_zero;
  copy->absolute_zero = p->absolute_zero;
  copy->has_stefan_boltzmann = p->has_stefan_boltzmann;
  copy->stefan_boltzmann = p->stefan_boltzmann;
  copy->source_line = p->source_line;
  *out = copy;
  return StatusCode::kOk;
}
}  // namespace fesa
