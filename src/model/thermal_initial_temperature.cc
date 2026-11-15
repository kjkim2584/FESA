#include "fesa/model/thermal_initial_temperature.h"

#include <cmath>
#include <memory>
#include <new>

#include "model/model_storage.h"

namespace fesa {
ThermalInitialTemperature::ThermalInitialTemperature(int node_id,
                                                     Vector* values,
                                                     char* source_file,
                                                     int source_line)
    : node_id_(node_id),
      values_(values),
      source_file_(source_file),
      source_line_(source_line) {}
ThermalInitialTemperature::~ThermalInitialTemperature() {
  delete values_;
  delete[] source_file_;
}
StatusCode ThermalInitialTemperature::Create(int node_id, const Vector* values,
                                             const char* source_file,
                                             int source_line,
                                             ThermalInitialTemperature** out,
                                             Diagnostic* d) {
  if (!out) return model_storage::Invalid("Thermal IC output", d);
  *out = nullptr;
  if (!model_storage::ValidId(node_id) || !values || values->Size() > 20 ||
      !source_file || !*source_file || source_line <= 0)
    return model_storage::Invalid("Thermal IC definition", d);
  for (int i = 0; i < values->Size(); ++i)
    if (!std::isfinite(values->Data()[i]))
      return model_storage::Invalid("Thermal IC nonfinite component", d);
  Vector* raw = nullptr;
  auto status = Vector::Create(values->Size(), &raw, d);
  std::unique_ptr<Vector> copy(raw);
  if (status != StatusCode::kOk) return status;
  status = copy->CopyFrom(values, d);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<char[]> file(model_storage::CopyText(source_file));
  if (!file) return model_storage::AllocationFailure("Thermal IC source", d);
  auto* result = new (std::nothrow)
      ThermalInitialTemperature(node_id, copy.get(), file.get(), source_line);
  if (!result) return model_storage::AllocationFailure("Thermal IC", d);
  copy.release();
  file.release();
  *out = result;
  return StatusCode::kOk;
}
int ThermalInitialTemperature::Id() const { return node_id_; }
const Vector* ThermalInitialTemperature::values() const { return values_; }
const char* ThermalInitialTemperature::source_file() const {
  return source_file_;
}
int ThermalInitialTemperature::source_line() const { return source_line_; }
}  // namespace fesa
