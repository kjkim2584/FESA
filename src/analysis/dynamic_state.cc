#include "fesa/analysis/dynamic_state.h"

#include <cassert>
#include <cmath>
#include <memory>
#include <new>
#include <utility>

#include "adapters/mkl/dense_ops.h"
#include "fesa/analysis/nonlinear_field_history.h"
#include "model/model_storage.h"

namespace fesa {
DynamicState::~DynamicState() {
  if (owns_history_) delete convergence_history_;
  delete committed_velocity_;
  delete committed_acceleration_;
  delete candidate_velocity_;
  delete candidate_acceleration_;
}
StatusCode DynamicState::Create(int size, DynamicState** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("DynamicState output", d);
  *out = nullptr;
  std::unique_ptr<DynamicState> state(new (std::nothrow) DynamicState);
  if (!state) return model_storage::AllocationFailure("DynamicState", d);
  for (auto** target :
       {&state->committed_velocity_, &state->committed_acceleration_,
        &state->candidate_velocity_, &state->candidate_acceleration_}) {
    const auto status = Vector::Create(size, target, d);
    if (status != StatusCode::kOk) return status;
  }
  *out = state.release();
  return StatusCode::kOk;
}
StatusCode DynamicState::Stage(const Vector* velocity,
                               const Vector* acceleration, double step_time,
                               double total_time, Diagnostic* d) {
  Discard();
  if (!std::isfinite(step_time) || !std::isfinite(total_time) ||
      step_time < 0 || total_time < step_time ||
      (committed_valid_ && total_time < total_time_))
    return model_storage::Invalid("DynamicState time", d);
  for (const auto* input : {velocity, acceleration}) {
    if (!input || input->Size() != candidate_velocity_->Size() ||
        input == candidate_velocity_ || input == candidate_acceleration_ ||
        !dense_ops::AllFinite(input->Size(), input->Data()))
      return model_storage::Invalid("DynamicState kinematic input", d);
  }
  auto status = candidate_velocity_->CopyFrom(velocity, d);
  if (status == StatusCode::kOk)
    status = candidate_acceleration_->CopyFrom(acceleration, d);
  if (status != StatusCode::kOk) return status;
  candidate_step_time_ = step_time;
  candidate_total_time_ = total_time;
  candidate_valid_ = true;
  return StatusCode::kOk;
}
void DynamicState::Discard() noexcept {
  candidate_valid_ = false;
  history_step_ = 0;
  history_initial_ = false;
  if (convergence_history_) convergence_history_->DiscardCandidate();
}
const NonlinearFieldHistory* DynamicState::ConvergenceHistory() const {
  return convergence_history_;
}
bool DynamicState::HistoryMatches(int step, int increment) const {
  return !convergence_history_ ||
         (convergence_history_->CanCommit() && history_step_ == step &&
          (history_initial_ ? increment == 0 : increment > 0));
}
void DynamicState::Commit() noexcept {
  assert(candidate_valid_);
  if (convergence_history_) convergence_history_->Commit();
  std::swap(committed_velocity_, candidate_velocity_);
  std::swap(committed_acceleration_, candidate_acceleration_);
  step_time_ = candidate_step_time_;
  total_time_ = candidate_total_time_;
  committed_valid_ = true;
  Discard();
}
const Vector* DynamicState::CommittedVelocity() const {
  return committed_valid_ ? committed_velocity_ : nullptr;
}
const Vector* DynamicState::CommittedAcceleration() const {
  return committed_valid_ ? committed_acceleration_ : nullptr;
}
const Vector* DynamicState::CandidateVelocity() const {
  return candidate_valid_ ? candidate_velocity_ : nullptr;
}
const Vector* DynamicState::CandidateAcceleration() const {
  return candidate_valid_ ? candidate_acceleration_ : nullptr;
}
double DynamicState::StepTime() const { return step_time_; }
double DynamicState::TotalTime() const { return total_time_; }
double DynamicState::CandidateStepTime() const { return candidate_step_time_; }
double DynamicState::CandidateTotalTime() const {
  return candidate_total_time_;
}
}  // namespace fesa
