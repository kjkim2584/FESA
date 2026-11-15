#include "fesa/analysis/incremental_solid_history.h"

#include <cassert>
#include <cmath>
#include <memory>
#include <new>
#include <utility>

#include "fesa/elements/incremental_solid_element.h"
#include "model/model_storage.h"
namespace fesa {
struct IncrementalSolidHistory::Slot {
  Slot() = default;
  Slot(const Slot&) = delete;
  Slot& operator=(const Slot&) = delete;
  Slot(Slot&&) = delete;
  Slot& operator=(Slot&&) = delete;
  ~Slot() {
    delete committed;
    delete candidate;
  }
  Matrix* committed = nullptr;
  Matrix* candidate = nullptr;
  bool staged = false;
  bool validated = false;
};
StatusCode IncrementalSolidHistory::Create(const Domain* domain,
    IncrementalSolidHistory** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("Incremental history output", d);
  *out = nullptr;
  if (!domain || domain->Elements()->empty())
    return model_storage::Invalid("Incremental history Domain", d);
  try {
    std::unique_ptr<IncrementalSolidHistory> history(
        new IncrementalSolidHistory);
    history->slots_ = new std::map<int, Slot*>;
    for (const auto& entry : *domain->Elements()) {
      const auto* element =
          dynamic_cast<const IncrementalSolidElement*>(entry.second);
      if (!element) {
        if (d) d->element_id = entry.first;
        return ReportFailure(StatusCode::kUnsupported,
            DiagnosticStage::kValidation,
            "Incremental history requires exclusive incremental solids", d);
      }
      std::unique_ptr<Slot> slot(new Slot);
      const int points = element->integration_point_count();
      auto status = Matrix::Create(points, 6, &slot->committed, d);
      if (status == StatusCode::kOk)
        status = Matrix::Create(points, 6, &slot->candidate, d);
      if (status != StatusCode::kOk) return status;
      const auto inserted = history->slots_->emplace(entry.first, slot.get());
      if (!inserted.second)
        return model_storage::Invalid("Incremental history duplicate ID", d);
      slot.release();
    }
    *out = history.release();
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return model_storage::AllocationFailure("Incremental history storage", d);
  }
}
IncrementalSolidHistory::~IncrementalSolidHistory() {
  if (slots_) {
    for (const auto& entry : *slots_) delete entry.second;
    delete slots_;
  }
}
const Matrix* IncrementalSolidHistory::CommittedStress(int element_id) const {
  const auto found = slots_->find(element_id);
  return found == slots_->end() ? nullptr : found->second->committed;
}
const Matrix* IncrementalSolidHistory::CandidateStress(int element_id) const {
  const auto found = slots_->find(element_id);
  return !prepared_ || found == slots_->end() || !found->second->staged
             ? nullptr : found->second->candidate;
}
StatusCode IncrementalSolidHistory::PrepareCandidate(Diagnostic* d) {
  if (!slots_) return model_storage::Invalid("Incremental history layout", d);
  DiscardCandidate();
  prepared_ = true;
  return StatusCode::kOk;
}
StatusCode IncrementalSolidHistory::StageStress(int element_id,
                                               const Matrix* stress,
                                               Diagnostic* d) {
  const auto found = slots_->find(element_id);
  if (!prepared_ || found == slots_->end())
    return model_storage::Invalid("Incremental stress stage identity", d);
  auto* slot = found->second;
  slot->staged = false;
  slot->validated = false;
  if (!stress || stress->Rows() != slot->candidate->Rows() ||
      stress->Columns() != 6)
    return model_storage::Invalid("Incremental stress stage dimensions", d);
  for (int i = 0; i < stress->Rows() * 6; ++i)
    if (!std::isfinite(stress->Data()[i]))
      return model_storage::Invalid("Incremental stress stage finite", d);
  const auto status = slot->candidate->CopyFrom(stress, d);
  slot->staged = status == StatusCode::kOk;
  return status;
}
StatusCode IncrementalSolidHistory::ValidateCandidate(Diagnostic* d) {
  if (!prepared_)
    return model_storage::Invalid("Incremental candidate not prepared", d);
  for (const auto& entry : *slots_)
    if (!entry.second->staged) {
      if (d) d->element_id = entry.first;
      return model_storage::Invalid("Incremental candidate missing IP", d);
    }
  for (const auto& entry : *slots_) entry.second->validated = true;
  return StatusCode::kOk;
}
bool IncrementalSolidHistory::CanCommit() const noexcept {
  if (!prepared_) return false;
  for (const auto& entry : *slots_)
    if (!entry.second->staged || !entry.second->validated) return false;
  return true;
}
void IncrementalSolidHistory::DiscardCandidate() noexcept {
  prepared_ = false;
  for (const auto& entry : *slots_) {
    entry.second->staged = false;
    entry.second->validated = false;
  }
}
void IncrementalSolidHistory::CommitCandidate() noexcept {
  assert(CanCommit());
  for (const auto& entry : *slots_)
    std::swap(entry.second->committed, entry.second->candidate);
  DiscardCandidate();
}
}  // namespace fesa
