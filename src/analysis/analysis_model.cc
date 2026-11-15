#include "fesa/analysis/analysis_model.h"

#include <new>

#include "model/model_storage.h"

namespace fesa {

AnalysisModel::~AnalysisModel() { Clear(); }
StatusCode AnalysisModel::Create(AnalysisModel** out, Diagnostic* d) {
  if (out == nullptr) return model_storage::Invalid("AnalysisModel::Create", d);
  *out = nullptr;
  try {
    *out = new (std::nothrow) AnalysisModel;
    if (*out == nullptr) {
      return model_storage::AllocationFailure("AnalysisModel::Create", d);
    }
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return model_storage::AllocationFailure("AnalysisModel::Create", d);
  }
}
StatusCode AnalysisModel::Bind(const Domain* domain, Diagnostic* d) {
  if (domain == nullptr) {
    return model_storage::Invalid("AnalysisModel::Bind", d);
  }
  try {
    std::map<int, Node*> nodes;
    std::map<int, Element*> elements;
    for (auto entry : *domain->Nodes()) {
      nodes.emplace(entry.first, entry.second);
    }
    for (auto entry : *domain->Elements()) {
      elements.emplace(entry.first, entry.second);
    }
    nodes_.swap(nodes);
    elements_.swap(elements);
    domain_ = domain;
    conditions_ = nullptr;
  } catch (const std::bad_alloc&) {
    return model_storage::AllocationFailure("AnalysisModel::Bind", d);
  }
  return StatusCode::kOk;
}
StatusCode AnalysisModel::BindStep(const Domain* domain,
                                   const EffectiveConditions* conditions,
                                   Diagnostic* d) {
  Clear();
  if (!domain || !conditions || conditions->SourceDomain() != domain)
    return model_storage::Invalid("AnalysisModel::BindStep", d);
  StatusCode status = Bind(domain, d);
  if (status != StatusCode::kOk) return status;
  conditions_ = conditions;
  return StatusCode::kOk;
}
const EffectiveConditions* AnalysisModel::Conditions() const {
  return conditions_;
}
void AnalysisModel::Clear() {
  conditions_ = nullptr;
  nodes_.clear();
  elements_.clear();
  domain_ = nullptr;
}
const std::map<int, Node*>* AnalysisModel::Nodes() const { return &nodes_; }
const std::map<int, Element*>* AnalysisModel::Elements() const {
  return &elements_;
}

}  // namespace fesa
