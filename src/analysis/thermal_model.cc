#include "fesa/analysis/thermal_model.h"

#include <memory>
#include <new>

namespace fesa {
ThermalModel::~ThermalModel() {
  delete map_;
  delete elements_;
}
StatusCode ThermalModel::Create(ThermalModel** out, Diagnostic* d) {
  if (!out)
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation,
                         "Null thermal view output", d);
  *out = nullptr;
  try {
    std::unique_ptr<ThermalModel> model(new ThermalModel);
    model->elements_ = new std::vector<const ThermalElement*>(0);
    *out = model.release();
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return ReportFailure(StatusCode::kAllocationFailure,
                         DiagnosticStage::kValidation,
                         "Thermal view allocation", d);
  }
}
StatusCode ThermalModel::BindStep(const Domain* domain,
                                  const ThermalConditions* conditions,
                                  Diagnostic* d) {
  ThermalDofMap* raw = nullptr;
  auto status = ThermalDofMap::Create(domain, conditions, &raw, d);
  std::unique_ptr<ThermalDofMap> candidate_map(raw);
  if (status != StatusCode::kOk) return status;
  try {
    std::vector<const ThermalElement*> candidate(0);
    candidate.reserve(domain->ThermalElements()->size());
    for (const auto entry : *domain->ThermalElements())
      candidate.push_back(entry.second);
    // Nothing below can allocate: publish list, map and borrowed inputs
    // together.
    elements_->swap(candidate);
    delete map_;
    map_ = candidate_map.release();
    domain_ = domain;
    conditions_ = conditions;
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return ReportFailure(StatusCode::kAllocationFailure,
                         DiagnosticStage::kValidation,
                         "Thermal view binding allocation", d);
  }
}
void ThermalModel::Clear() noexcept {
  delete map_;
  map_ = nullptr;
  domain_ = nullptr;
  conditions_ = nullptr;
  elements_->clear();
}
const std::vector<const ThermalElement*>* ThermalModel::Elements() const {
  return elements_;
}
const ThermalConditions* ThermalModel::Conditions() const {
  return conditions_;
}
const ThermalDofMap* ThermalModel::DofMap() const { return map_; }
const Domain* ThermalModel::SourceDomain() const { return domain_; }
}  // namespace fesa
