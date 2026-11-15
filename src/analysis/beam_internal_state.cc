#include "fesa/analysis/beam_internal_state.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <memory>
#include <new>
#include <utility>

#include "fesa/elements/beam3d.h"
#include "fesa/model/domain.h"
#include "model/model_storage.h"

namespace fesa {
StatusCode BeamInternalState::Create(const Domain* domain,
                                      BeamInternalState** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("Beam internal output", d);
  *out = nullptr;
  if (!domain || domain->Elements()->empty() ||
      domain->Elements()->size() > std::numeric_limits<int>::max() / 2)
    return model_storage::Invalid("Beam internal layout", d);
  for (const auto& entry : *domain->Elements()) {
    if (!dynamic_cast<const Beam3D*>(entry.second)) {
      if (d) d->element_id = entry.first;
      return ReportFailure(StatusCode::kUnsupported,
                           DiagnosticStage::kValidation,
                           "Beam internal state requires exclusive Beam3D", d);
    }
  }
  std::unique_ptr<BeamInternalState> state(
      new (std::nothrow) BeamInternalState);
  if (!state)
    return model_storage::AllocationFailure("Beam internal state", d);
  state->element_count_ = static_cast<int>(domain->Elements()->size());
  state->element_ids_ = new (std::nothrow) int[state->element_count_];
  if (!state->element_ids_)
    return model_storage::AllocationFailure("Beam internal IDs", d);
  int row = 0;
  for (const auto& entry : *domain->Elements())
    state->element_ids_[row++] = entry.first;
  auto status = Vector::Create(2 * state->element_count_,
                                &state->committed_position_, d);
  if (status == StatusCode::kOk)
    status = Vector::Create(2 * state->element_count_,
                            &state->candidate_position_, d);
  if (status != StatusCode::kOk) return status;
  *out = state.release();
  return StatusCode::kOk;
}
BeamInternalState::~BeamInternalState() {
  delete candidate_acceleration_;
  delete candidate_velocity_;
  delete committed_acceleration_;
  delete committed_velocity_;
  delete candidate_position_;
  delete committed_position_;
  delete[] element_ids_;
}
int BeamInternalState::ElementCount() const { return element_count_; }
int BeamInternalState::ElementId(int row) const {
  return row >= 0 && row < element_count_ ? element_ids_[row] : -1;
}
int BeamInternalState::Offset(int element_id) const {
  const int* begin = element_ids_;
  const int* end = begin + element_count_;
  const int* found = std::lower_bound(begin, end, element_id);
  return found != end && *found == element_id
             ? 2 * static_cast<int>(found - element_ids_)
             : -1;
}
const Vector* BeamInternalState::CommittedPosition() const {
  return committed_position_;
}
const Vector* BeamInternalState::CandidatePosition() const {
  return prepared_ ? candidate_position_ : nullptr;
}
const Vector* BeamInternalState::CommittedVelocity() const {
  return committed_motion_valid_ ? committed_velocity_ : nullptr;
}
const Vector* BeamInternalState::CommittedAcceleration() const {
  return committed_motion_valid_ ? committed_acceleration_ : nullptr;
}
const Vector* BeamInternalState::CandidateVelocity() const {
  return candidate_motion_valid_ ? candidate_velocity_ : nullptr;
}
const Vector* BeamInternalState::CandidateAcceleration() const {
  return candidate_motion_valid_ ? candidate_acceleration_ : nullptr;
}
StatusCode BeamInternalState::EnableDynamics(Diagnostic* d) {
  if (prepared_)
    return model_storage::Invalid("Beam motion requires inactive candidate", d);
  for (auto** store : {&committed_velocity_, &committed_acceleration_,
                       &candidate_velocity_, &candidate_acceleration_}) {
    if (*store) continue;
    const auto status = Vector::Create(2 * element_count_, store, d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}
StatusCode BeamInternalState::StageMotion(const Vector* velocity,
                                         const Vector* acceleration,
                                         Diagnostic* d) {
  DiscardMotion();
  if (!prepared_ || !candidate_velocity_ || !candidate_acceleration_ ||
      velocity == acceleration)
    return model_storage::Invalid("Beam internal motion state/inputs", d);
  for (const auto* input : {velocity, acceleration}) {
    if (!input || input->Size() != 2 * element_count_ ||
        input == committed_position_ || input == candidate_position_ ||
        input == committed_velocity_ || input == committed_acceleration_ ||
        input == candidate_velocity_ || input == candidate_acceleration_)
      return model_storage::Invalid("Beam internal motion input", d);
    for (int i = 0; i < input->Size(); ++i)
      if (!std::isfinite(input->Data()[i]))
        return model_storage::Invalid("Beam internal nonfinite motion", d);
  }
  auto status = candidate_velocity_->CopyFrom(velocity, d);
  if (status == StatusCode::kOk)
    status = candidate_acceleration_->CopyFrom(acceleration, d);
  candidate_motion_valid_ = status == StatusCode::kOk;
  return status;
}
void BeamInternalState::DiscardMotion() noexcept {
  candidate_motion_valid_ = false;
}
bool BeamInternalState::CanCommit(bool require_motion) const {
  return prepared_ && (!require_motion || candidate_motion_valid_);
}
StatusCode BeamInternalState::PrepareCandidate(Diagnostic* d) {
  DiscardMotion();
  prepared_ = false;
  const auto status = candidate_position_->CopyFrom(committed_position_, d);
  prepared_ = status == StatusCode::kOk;
  return status;
}
StatusCode BeamInternalState::ApplyCorrection(const Vector* correction,
                                              Diagnostic* d) {
  DiscardMotion();
  if (!prepared_ || !correction || correction == committed_position_ ||
      correction == candidate_position_ || correction == committed_velocity_ ||
      correction == committed_acceleration_ ||
      correction == candidate_velocity_ ||
      correction == candidate_acceleration_) {
    DiscardCandidate();
    return model_storage::Invalid("Beam internal correction", d);
  }
  const auto status = candidate_position_->AddScaled(1, correction, d);
  if (status != StatusCode::kOk) DiscardCandidate();
  return status;
}
void BeamInternalState::DiscardCandidate() noexcept {
  prepared_ = false;
  DiscardMotion();
}
void BeamInternalState::CommitCandidate(bool require_motion) noexcept {
  assert(CanCommit(require_motion));
  std::swap(candidate_position_, committed_position_);
  if (require_motion) {
    std::swap(candidate_velocity_, committed_velocity_);
    std::swap(candidate_acceleration_, committed_acceleration_);
    committed_motion_valid_ = true;
  }
  DiscardCandidate();
}
}  // namespace fesa
