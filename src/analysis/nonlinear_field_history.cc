#include "fesa/analysis/nonlinear_field_history.h"

#include <cassert>
#include <cmath>
#include <limits>
#include <new>

#include "model/model_storage.h"
namespace fesa {
namespace {
bool Nonnegative(double value) { return std::isfinite(value) && value >= 0; }
}  // namespace
StatusCode NonlinearFieldHistory::Create(NonlinearFieldHistory** out,
                                         Diagnostic* d) {
  return Create(true, out, d);
}
StatusCode NonlinearFieldHistory::Create(bool moment_active,
                                         NonlinearFieldHistory** out,
                                         Diagnostic* d) {
  if (!out) return model_storage::Invalid("Nonlinear field history output", d);
  *out = new (std::nothrow) NonlinearFieldHistory;
  if (*out) {
    (*out)->moment_active_ = moment_active;
    if (!moment_active) (*out)->committed_mean_[1] = 0;
  }
  return *out ? StatusCode::kOk
              : model_storage::AllocationFailure("Nonlinear field history", d);
}
StatusCode NonlinearFieldHistory::PrepareCandidate(
    int step, const ControlSettings* controls, Diagnostic* d) {
  DiscardCandidate();
  const bool same = step == committed_step_;
  const bool next = committed_step_ < std::numeric_limits<int>::max() &&
                    step == committed_step_ + 1;
  if (step <= 0 || (!same && !next) || !controls || !controls->displacement ||
      (moment_active_ && !controls->rotation))
    return model_storage::Invalid("Field history step/controls", d);
  const FieldControls* fields[2] = {controls->displacement, controls->rotation};
  for (int field = 0; field < 2; ++field) {
    if (field == 1 && !moment_active_) continue;
    const auto* c = fields[field];
    const double seed = c->q0_explicit ? c->q0 : committed_mean_[field];
    if (c->qu_explicit) {
      if (!std::isfinite(c->qu) || !(c->qu > 0))
        return model_storage::Invalid("Field history constant raw scale", d);
    } else if (!Nonnegative(c->epsilon) || c->epsilon >= 1 ||
               !std::isfinite(seed) || !(seed > 0))
      return model_storage::Invalid("Field history automatic seed/epsilon", d);
    constant_[field] = c->qu_explicit;
    constant_value_[field] = c->qu;
    epsilon_[field] = c->epsilon;
    seed_[field] = seed;
    candidate_initial_maximum_[field] =
        same ? committed_initial_maximum_[field] : 0;
  }
  candidate_step_ = step;
  prepared_ = true;
  return StatusCode::kOk;
}
StatusCode NonlinearFieldHistory::PrepareDynamicInitial(
    int step, const ControlSettings* controls, double initial_maximum,
    Diagnostic* d) {
  DiscardCandidate();
  if (moment_active_)
    return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                         "Dynamic initial history requires force-only field",
                         d);
  return PrepareDynamicInitial(step, controls, initial_maximum, 0, d);
}
StatusCode NonlinearFieldHistory::PrepareDynamicInitial(
    int step, const ControlSettings* controls, double force_maximum,
    double moment_maximum, Diagnostic* d) {
  DiscardCandidate();
  if (!Nonnegative(force_maximum) ||
      (moment_active_ && !Nonnegative(moment_maximum)) ||
      committed_step_ == std::numeric_limits<int>::max() ||
      step != committed_step_ + 1)
    return model_storage::Invalid("Dynamic initial history step/maximum", d);
  const auto status = PrepareCandidate(step, controls, d);
  if (status != StatusCode::kOk) return status;
  const double maxima[2] = {force_maximum, moment_maximum};
  candidate_increments_ = 0;
  for (int field = 0; field < 2; ++field) {
    if (field == 1 && !moment_active_) continue;
    // (M*H+q)/(M+1), without overflowing the weighted numerator.
    candidate_initial_maximum_[field] =
        committed_increments_ == 0
            ? maxima[field]
            : committed_qmax_[field] +
                  (maxima[field] - committed_qmax_[field]) /
                      (static_cast<double>(committed_increments_) + 1);
    candidate_mean_[field] =
        constant_[field] ? constant_value_[field] : seed_[field];
    candidate_n_[field] = 0;
    candidate_qmax_[field] = 0;
  }
  initial_candidate_ = staged_ = true;
  return StatusCode::kOk;
}
StatusCode NonlinearFieldHistory::PrepareDynamicInitialAfterLinearStatic(
    int previous_step, const ControlSettings* controls, double previous_raw,
    double previous_maximum, double initial_maximum, Diagnostic* d) {
  if (moment_active_) {
    DiscardCandidate();
    return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                         "Dynamic static handoff requires force-only field", d);
  }
  const double raw[2] = {previous_raw, 0};
  const double maximum[2] = {previous_maximum, 0};
  const double initial[2] = {initial_maximum, 0};
  return PrepareDynamicInitialAfterLinearStatic(previous_step, controls, raw,
                                                maximum, initial, d);
}
StatusCode NonlinearFieldHistory::PrepareDynamicInitialAfterLinearStatic(
    int previous_step, const ControlSettings* controls,
    const double* previous_raw, const double* previous_maximum,
    const double* initial_maximum, Diagnostic* d) {
  DiscardCandidate();
  if (!previous_raw || !previous_maximum || !initial_maximum ||
      committed_step_ != 0 || previous_step <= 0 ||
      previous_step == std::numeric_limits<int>::max())
    return model_storage::Invalid("Linear static dynamic history handoff", d);
  for (int field = 0; field < (moment_active_ ? 2 : 1); ++field)
    if (!std::isfinite(previous_raw[field]) || previous_raw[field] <= 0 ||
        !Nonnegative(previous_maximum[field]))
      return model_storage::Invalid("Linear static dynamic history field", d);
  // Reuse initial validation/controls/count staging on this fresh history.
  // Only the prepared transaction is retagged; committed history stays empty.
  const auto status = PrepareDynamicInitial(1, controls, initial_maximum[0],
                                            initial_maximum[1], d);
  if (status != StatusCode::kOk) return status;
  candidate_step_ = previous_step + 1;
  const FieldControls* fields[2] = {controls->displacement, controls->rotation};
  for (int field = 0; field < (moment_active_ ? 2 : 1); ++field) {
    candidate_initial_maximum_[field] =
        previous_maximum[field] +
        (initial_maximum[field] - previous_maximum[field]) / 2;
    if (!fields[field]->q0_explicit && !fields[field]->qu_explicit)
      candidate_mean_[field] = previous_raw[field];
  }
  return StatusCode::kOk;
}
StatusCode NonlinearFieldHistory::Preview(double force_mean, double moment_mean,
                                          double force_max, double moment_max,
                                          double* force_raw, double* moment_raw,
                                          Diagnostic* d) {
  const bool ready = prepared_;
  staged_ = false;
  if (!ready || initial_candidate_ || !force_raw || !moment_raw ||
      force_raw == moment_raw || !Nonnegative(force_mean) ||
      !Nonnegative(force_max) ||
      (moment_active_ &&
       (!Nonnegative(moment_mean) || !Nonnegative(moment_max)))) {
    DiscardCandidate();
    return model_storage::Invalid("Field history preview inputs", d);
  }
  const bool same = candidate_step_ == committed_step_;
  const int increments = same ? committed_increments_ : 0;
  if (increments == std::numeric_limits<int>::max()) {
    DiscardCandidate();
    return model_storage::Invalid("Field history increment count exhausted", d);
  }
  const double means[2] = {force_mean, moment_mean};
  const double maxima[2] = {force_max, moment_max};
  double raw[2] = {}, max_mean[2] = {};
  int included[2] = {};
  for (int field = 0; field < 2; ++field) {
    if (field == 1 && !moment_active_) continue;
    const double previous_max = same ? committed_qmax_[field] : 0;
    max_mean[field] =
        increments == 0
            ? maxima[field]
            : previous_max + (maxima[field] - previous_max) /
                                 (static_cast<double>(increments) + 1);
    const int n = same ? committed_n_[field] : 0;
    const double mu = same ? committed_mean_[field] : seed_[field];
    if (constant_[field]) {
      raw[field] = constant_value_[field];
      included[field] = 0;
    } else {
      const double threshold =
          n == 0 ? epsilon_[field] * seed_[field]
                 : epsilon_[field] *
                       (static_cast<double>(n) /
                        (static_cast<double>(n) + 1 - epsilon_[field])) *
                       mu;
      const bool include = means[field] > threshold;
      included[field] = n + (include ? 1 : 0);
      raw[field] =
          !include ? mu
          : n == 0 ? means[field]
                   : mu + (means[field] - mu) / (static_cast<double>(n) + 1);
    }
    if (!std::isfinite(raw[field]) || !(raw[field] > 0) ||
        !Nonnegative(max_mean[field])) {
      DiscardCandidate();
      return ReportFailure(StatusCode::kLibraryFailure,
                           DiagnosticStage::kComputation,
                           "Field history scalar mean range", d);
    }
  }
  for (int field = 0; field < 2; ++field) {
    candidate_mean_[field] = raw[field];
    candidate_qmax_[field] = max_mean[field];
    candidate_n_[field] = included[field];
  }
  candidate_increments_ = increments + 1;
  *force_raw = raw[0];
  *moment_raw = raw[1];
  staged_ = true;
  return StatusCode::kOk;
}
void NonlinearFieldHistory::Commit() noexcept {
  assert(CanCommit());
  committed_step_ = candidate_step_;
  committed_increments_ = candidate_increments_;
  for (int field = 0; field < 2; ++field) {
    committed_initial_maximum_[field] = candidate_initial_maximum_[field];
    committed_mean_[field] = candidate_mean_[field];
    committed_qmax_[field] = candidate_qmax_[field];
    committed_n_[field] = candidate_n_[field];
  }
  DiscardCandidate();
}
void NonlinearFieldHistory::DiscardCandidate() noexcept {
  prepared_ = staged_ = false;
  initial_candidate_ = false;
}
bool NonlinearFieldHistory::CanCommit() const { return prepared_ && staged_; }
double NonlinearFieldHistory::CommittedRawMean(int field) const {
  assert(field >= 0 && field < 2);
  return committed_mean_[field];
}
double NonlinearFieldHistory::CommittedMaximumMean(int field) const {
  assert(field >= 0 && field < 2);
  return committed_qmax_[field];
}
int NonlinearFieldHistory::CommittedIncludedCount(int field) const {
  assert(field >= 0 && field < 2);
  return committed_n_[field];
}
int NonlinearFieldHistory::CommittedIncrementCount() const {
  return committed_increments_;
}
int NonlinearFieldHistory::CommittedStepNumber() const {
  return committed_step_;
}
double NonlinearFieldHistory::PreviousMaximumMean(int field) const {
  assert(prepared_ && field >= 0 && field < 2);
  if (initial_candidate_) return candidate_initial_maximum_[field];
  if (candidate_step_ != committed_step_) return 0;
  if (committed_increments_ == 0) return committed_initial_maximum_[field];
  return committed_qmax_[field];
}
int NonlinearFieldHistory::CandidateIncludedCount(int field) const {
  assert(CanCommit() && field >= 0 && field < 2);
  return candidate_n_[field];
}
}  // namespace fesa
