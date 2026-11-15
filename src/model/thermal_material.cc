#include "fesa/model/thermal_material.h"

#include <cmath>
#include <memory>
#include <new>

#include "model/model_storage.h"

namespace fesa {
ThermalMaterial::ThermalMaterial(int id, char* name,
                                 TemperatureTable* conductivity,
                                 bool has_density, double density,
                                 TemperatureTable* heat)
    : Material(id, name),
      conductivity_(conductivity),
      specific_heat_(heat),
      has_density_(has_density),
      density_(density) {}
ThermalMaterial::~ThermalMaterial() {
  delete conductivity_;
  delete specific_heat_;
}
StatusCode ThermalMaterial::Create(int id, const char* name,
                                   const TemperatureTable* conductivity,
                                   bool has_density, double density,
                                   const TemperatureTable* specific_heat,
                                   ThermalMaterial** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("ThermalMaterial output", d);
  *out = nullptr;
  if (!model_storage::ValidId(id) || !name || !*name || !conductivity ||
      !conductivity->HasPositiveValues() ||
      (specific_heat && !specific_heat->HasPositiveValues()) ||
      (has_density && (!std::isfinite(density) || density <= 0)))
    return model_storage::Invalid("ThermalMaterial properties", d);
  std::unique_ptr<char[]> owned_name(model_storage::CopyText(name));
  if (!owned_name)
    return model_storage::AllocationFailure("ThermalMaterial name", d);
  std::unique_ptr<ThermalMaterial> result(new (std::nothrow) ThermalMaterial(
      id, owned_name.get(), nullptr, has_density, has_density ? density : 0,
      nullptr));
  if (!result) return model_storage::AllocationFailure("ThermalMaterial", d);
  owned_name.release();
  auto status = conductivity->Clone(&result->conductivity_, d);
  if (status == StatusCode::kOk && specific_heat)
    status = specific_heat->Clone(&result->specific_heat_, d);
  if (status != StatusCode::kOk) return status;
  *out = result.release();
  return StatusCode::kOk;
}
const ThermalMaterial* ThermalMaterial::AsThermal() const { return this; }
StatusCode ThermalMaterial::Conductivity(double theta, double* out,
                                         Diagnostic* d) const {
  double slope = 0;
  return conductivity_->Evaluate(theta, out, &slope, d);
}
StatusCode ThermalMaterial::Capacity(double theta, double* out,
                                     Diagnostic* d) const {
  if (!out) return model_storage::Invalid("ThermalMaterial capacity output", d);
  if (!HasCapacity())
    return ReportFailure(StatusCode::kInputError, DiagnosticStage::kValidation,
                         "ThermalMaterial density or specific heat missing", d);
  double heat = 0, slope = 0;
  const auto status = specific_heat_->Evaluate(theta, &heat, &slope, d);
  if (status != StatusCode::kOk) return status;
  const double result = density_ * heat;
  if (!std::isfinite(result) || result <= 0)
    return ReportFailure(StatusCode::kLibraryFailure,
                         DiagnosticStage::kComputation,
                         "ThermalMaterial capacity range", d);
  *out = result;
  return StatusCode::kOk;
}
StatusCode ThermalMaterial::EnergyDifference(double from, double to,
                                             double* out, Diagnostic* d) const {
  if (!out) return model_storage::Invalid("ThermalMaterial energy output", d);
  if (!HasCapacity())
    return ReportFailure(StatusCode::kInputError, DiagnosticStage::kValidation,
                         "ThermalMaterial density or specific heat missing", d);
  double integral = 0;
  const auto status = specific_heat_->Integrate(from, to, &integral, d);
  if (status != StatusCode::kOk) return status;
  const double result = density_ * integral;
  if (!std::isfinite(result) || (integral != 0 && result == 0))
    return ReportFailure(StatusCode::kLibraryFailure,
                         DiagnosticStage::kComputation,
                         "ThermalMaterial energy range", d);
  *out = result;
  return StatusCode::kOk;
}
bool ThermalMaterial::HasCapacity() const {
  return has_density_ && specific_heat_;
}
bool ThermalMaterial::IsTemperatureDependent() const {
  return !conductivity_->IsConstant() ||
         (specific_heat_ && !specific_heat_->IsConstant());
}
bool ThermalMaterial::IsConductivityTemperatureDependent() const {
  return !conductivity_->IsConstant();
}
}  // namespace fesa
