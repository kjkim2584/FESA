#include "fesa/model/nodal_temperature_result.h"

#include <cmath>
#include <memory>
#include <new>

#include "model/model_storage.h"

namespace fesa {
NodalTemperatureResult::~NodalTemperatureResult() {
  delete temperatures_;
  delete reaction_fluxes_;
}
StatusCode NodalTemperatureResult::Create(int count,
                                          NodalTemperatureResult** out,
                                          Diagnostic* d) {
  if (!out) return model_storage::Invalid("NodalTemperatureResult output", d);
  *out = nullptr;
  if (count < 1 || count > 20)
    return ReportFailure(StatusCode::kInputError, DiagnosticStage::kValidation,
                         "Temperature component count must be 1..20", d);
  std::unique_ptr<NodalTemperatureResult> result(new (std::nothrow)
                                                     NodalTemperatureResult);
  if (!result)
    return model_storage::AllocationFailure("NodalTemperatureResult", d);
  auto status = Vector::Create(count, &result->temperatures_, d);
  if (status == StatusCode::kOk)
    status = Vector::Create(count, &result->reaction_fluxes_, d);
  if (status != StatusCode::kOk) return status;
  *out = result.release();
  return StatusCode::kOk;
}
int NodalTemperatureResult::DofCount() const { return temperatures_->Size(); }
const Vector* NodalTemperatureResult::Temperatures() const {
  return temperatures_;
}
const Vector* NodalTemperatureResult::ReactionFluxes() const {
  return reaction_fluxes_;
}
StatusCode NodalTemperatureResult::Update(const Vector* nt, const Vector* rfl,
                                          Diagnostic* d) {
  if (!nt || !rfl || nt->Size() != DofCount() || rfl->Size() != DofCount())
    return model_storage::Invalid("NodalTemperatureResult input sizes", d);
  if (nt == reaction_fluxes_ || rfl == temperatures_)
    return model_storage::Invalid("NodalTemperatureResult cross-alias input",
                                  d);
  for (int i = 0; i < DofCount(); ++i)
    if (!std::isfinite(nt->Data()[i]) || !std::isfinite(rfl->Data()[i]))
      return model_storage::Invalid("NodalTemperatureResult finite inputs", d);
  auto status = temperatures_->CopyFrom(nt, d);
  if (status == StatusCode::kOk) status = reaction_fluxes_->CopyFrom(rfl, d);
  return status;
}
}  // namespace fesa
