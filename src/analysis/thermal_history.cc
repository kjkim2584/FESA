#include "fesa/analysis/thermal_history.h"

#include <cassert>
#include <cmath>
#include <limits>
#include <new>

#include "fesa/io/control_settings.h"

namespace fesa {
namespace {
StatusCode Error(const char* message, Diagnostic* d,
                 StatusCode code = StatusCode::kInvalidArgument) {
  return ReportFailure(code, DiagnosticStage::kValidation, message, d);
}
}  // namespace
StatusCode HeatFluxHistory::Create(HeatFluxHistory** out, Diagnostic* d) {
  if (!out) return Error("Null thermal history output", d);
  *out = new (std::nothrow) HeatFluxHistory;
  return *out ? StatusCode::kOk
              : Error("Thermal history allocation", d,
                      StatusCode::kAllocationFailure);
}
StatusCode HeatFluxHistory::BeginStep(int step, const FieldControls* c,
                                      Diagnostic* d) {
  if (!c || staged_ || step <= 0 || step_ == std::numeric_limits<int>::max() ||
      step != step_ + 1)
    return Error("Invalid thermal history step", d);
  if (!std::isfinite(c->epsilon) || c->epsilon < 0 ||
      (c->qu_explicit && (!std::isfinite(c->qu) || c->qu <= 0)) ||
      (!c->qu_explicit && c->q0_explicit &&
       (!std::isfinite(c->q0) || c->q0 <= 0)))
    return Error("Invalid thermal history controls", d,
                 StatusCode::kInputError);
  // Capture the ended step using its own mode/epsilon/h0, never new controls.
  double carry = carry_, seed = maximum_seed_;
  if (step_ > 0) {
    const bool extra = user_ ? active_ && last_ > 0 : last_ > epsilon_ * h0_;
    if (extra && included_ == std::numeric_limits<int>::max())
      return Error("Thermal history count overflow", d,
                   StatusCode::kLibraryFailure);
    const int count = included_ + (extra ? 1 : 0);
    const double sum = sum_ + (extra ? last_ : 0);
    carry = count > 0 ? sum / count : h0_;
    seed = maximum_count_ > 0 ? maximum_sum_ / maximum_count_ : 0;
    if (!std::isfinite(sum) || !std::isfinite(carry) || carry <= 0 ||
        !std::isfinite(seed))
      return Error("Thermal history carry overflow", d,
                   StatusCode::kLibraryFailure);
  }
  const double initial = c->qu_explicit ? .01 : c->q0_explicit ? c->q0 : carry;
  if (!std::isfinite(c->epsilon * initial))
    return Error("Thermal history gate overflow", d,
                 StatusCode::kLibraryFailure);
  carry_ = carry;
  maximum_seed_ = seed;
  step_ = step;
  h0_ = initial;
  epsilon_ = c->epsilon;
  user_ = c->qu_explicit;
  qu_ = c->qu;
  included_ = maximum_count_ = 0;
  sum_ = maximum_sum_ = last_ = 0;
  active_ = staged_ = false;
  return StatusCode::kOk;
}
StatusCode HeatFluxHistory::Preview(double b, double maximum, double* used,
                                    Diagnostic* d) {
  DiscardCandidate();
  if (step_ == 0 || !used || !std::isfinite(b) || b < 0 ||
      !std::isfinite(maximum) || maximum < 0)
    return Error("Invalid thermal history preview", d);
  const double prior = included_ > 0 ? sum_ / included_ : h0_;
  const double gate = epsilon_ * (user_ ? h0_ : prior);
  const bool active = active_ || (user_ && b > gate);
  const bool include = user_ ? active && b > 0 : b > gate;
  if (!std::isfinite(gate) ||
      maximum_count_ == std::numeric_limits<int>::max() ||
      (include && included_ == std::numeric_limits<int>::max()))
    return Error("Thermal history count/gate overflow", d,
                 StatusCode::kLibraryFailure);
  const double sum = sum_ + (include ? b : 0);
  const double maxima = maximum_sum_ + maximum;
  const int count = included_ + (include ? 1 : 0);
  const double value = user_ ? qu_ : count > 0 ? sum / count : h0_;
  if (!std::isfinite(sum) || !std::isfinite(maxima) || !std::isfinite(value) ||
      value <= 0)
    return Error("Thermal history sum overflow", d,
                 StatusCode::kLibraryFailure);
  candidate_sum_ = sum;
  candidate_included_ = count;
  candidate_maximum_sum_ = maxima;
  candidate_maximum_count_ = maximum_count_ + 1;
  candidate_last_ = b;
  candidate_used_ = value;
  candidate_active_ = active;
  staged_ = true;
  *used = value;
  return StatusCode::kOk;
}
void HeatFluxHistory::Commit() noexcept {
  assert(staged_);
  sum_ = candidate_sum_;
  included_ = candidate_included_;
  maximum_sum_ = candidate_maximum_sum_;
  maximum_count_ = candidate_maximum_count_;
  last_ = candidate_last_;
  committed_used_ = candidate_used_;
  active_ = candidate_active_;
  staged_ = false;
}
void HeatFluxHistory::DiscardCandidate() noexcept { staged_ = false; }
bool HeatFluxHistory::CanCommit() const { return staged_; }
int HeatFluxHistory::StepNumber() const { return step_; }
int HeatFluxHistory::AcceptedCount() const { return maximum_count_; }
int HeatFluxHistory::IncludedCount() const { return included_; }
double HeatFluxHistory::PreviousMaximumMean() const {
  return maximum_count_ > 0 ? maximum_sum_ / maximum_count_ : maximum_seed_;
}
double HeatFluxHistory::DefaultCarry() const { return carry_; }
double HeatFluxHistory::StepInitialQ0() const { return h0_; }
double HeatFluxHistory::CommittedTimeMean() const { return committed_used_; }
StatusCode RateHistory::Create(RateHistory** out, Diagnostic* d) {
  if (!out) return Error("Null thermal rate history output", d);
  *out = new (std::nothrow) RateHistory;
  return *out ? StatusCode::kOk
              : Error("Thermal rate history allocation", d,
                      StatusCode::kAllocationFailure);
}
void RateHistory::BeginStep() noexcept {
  count_ = candidate_count_ = 0;
  for (int i = 0; i < 3; ++i) rates_[i] = candidate_rates_[i] = 0;
  staged_ = false;
}
StatusCode RateHistory::Preview(double maximum, double dt, Diagnostic* d) {
  DiscardCandidate();
  if (!std::isfinite(maximum) || maximum < 0 || !std::isfinite(dt) || dt <= 0)
    return Error("Invalid thermal rate history input", d);
  const double rate = maximum / dt;
  if (!std::isfinite(rate) || count_ == std::numeric_limits<int>::max())
    return Error("Thermal rate history overflow", d,
                 StatusCode::kLibraryFailure);
  for (int i = 0; i < 3; ++i) candidate_rates_[i] = rates_[i];
  if (count_ < 3)
    candidate_rates_[count_] = rate;
  else {
    candidate_rates_[0] = rates_[1];
    candidate_rates_[1] = rates_[2];
    candidate_rates_[2] = rate;
  }
  candidate_count_ = count_ + 1;
  staged_ = true;
  return StatusCode::kOk;
}
void RateHistory::Commit() noexcept {
  assert(staged_);
  for (int i = 0; i < 3; ++i) rates_[i] = candidate_rates_[i];
  count_ = candidate_count_;
  staged_ = false;
}
void RateHistory::DiscardCandidate() noexcept { staged_ = false; }
bool RateHistory::CanCommit() const { return staged_; }
int RateHistory::AcceptedCount() const { return count_; }
StatusCode RateHistory::AllRecentBelow(double dt, double threshold, bool* out,
                                       Diagnostic* d) const {
  if (!out || !std::isfinite(dt) || dt <= 0 || !std::isfinite(threshold) ||
      threshold < 0)
    return Error("Invalid thermal rate comparison", d);
  bool below = count_ >= 3;
  if (below)
    for (int i = 0; i < 3; ++i) {
      const double value = dt * rates_[i];
      if (!std::isfinite(value))
        return Error("Thermal rate comparison overflow", d,
                     StatusCode::kLibraryFailure);
      below = below && value < threshold;
    }
  *out = below;
  return StatusCode::kOk;
}
}  // namespace fesa
