/// @file
/// @brief Transactional thermal norm statistics and three accepted rate
/// samples.
#ifndef FESA_ANALYSIS_THERMAL_HISTORY_H_
#define FESA_ANALYSIS_THERMAL_HISTORY_H_
#include "fesa/base/diagnostic.h"
namespace fesa {
struct FieldControls;
/// Owns scalar SUM/COUNT sufficient statistics, never model objects.
/// Powers use E/T; qbar/qmax follow formulation7. Preview never changes the
/// committed baseline. Only a final accepted preview may be committed.
class HeatFluxHistory final {
 public:
  static StatusCode Create(HeatFluxHistory** out, Diagnostic* d);
  ~HeatFluxHistory() = default;
  HeatFluxHistory(const HeatFluxHistory&) = delete;
  HeatFluxHistory& operator=(const HeatFluxHistory&) = delete;
  HeatFluxHistory(HeatFluxHistory&&) = delete;
  HeatFluxHistory& operator=(HeatFluxHistory&&) = delete;
  /// Borrows controls for this call, retains only relevant scalar values.
  /// Requires next positive step and no staged preview. Captures ended-step
  /// carry/maximum seed before clearing step statistics. Preserves the last
  /// committed output norm. Invalid controls/sum overflow preserve all state.
  StatusCode BeginStep(int step, const FieldControls* controls, Diagnostic* d);
  /// Stages one final-iteration sample from the committed baseline. Nonnegative
  /// finite qbar/qmax required. Success writes the norm used for convergence
  /// (qu in user mode); error preserves *used and all committed data, but
  /// invalidates the staged preview. No allocation or retained argument.
  StatusCode Preview(double qbar, double qmax, double* used, Diagnostic* d);
  /// Requires CanCommit(); scalar assignments only, no allocation/backend.
  void Commit() noexcept;
  void DiscardCandidate() noexcept;
  bool CanCommit() const;
  int StepNumber() const;
  int AcceptedCount() const;
  int IncludedCount() const;
  double PreviousMaximumMean() const;
  double DefaultCarry() const;
  double StepInitialQ0() const;
  double CommittedTimeMean()
      const;  ///< Last committed output norm, including qu.
 private:
  HeatFluxHistory() = default;
  int step_ = 0;
  int included_ = 0, maximum_count_ = 0;
  double sum_ = 0, maximum_sum_ = 0, last_ = 0;
  double carry_ = .01, maximum_seed_ = 0, h0_ = .01, epsilon_ = 1e-5;
  double qu_ = 0, committed_used_ = .01;
  bool user_ = false, active_ = false, staged_ = false;
  int candidate_included_ = 0, candidate_maximum_count_ = 0;
  double candidate_sum_ = 0, candidate_maximum_sum_ = 0;
  double candidate_last_ = 0, candidate_used_ = .01;
  bool candidate_active_ = false;
};
/// Fixed IT=3 ring of accepted S/dt [theta/T], with a separate staged sample.
/// No configurable history size; TIME INCREMENTATION data are unsupported.
class RateHistory final {
 public:
  static StatusCode Create(RateHistory** out, Diagnostic* d);
  ~RateHistory() = default;
  RateHistory(const RateHistory&) = delete;
  RateHistory& operator=(const RateHistory&) = delete;
  RateHistory(RateHistory&&) = delete;
  RateHistory& operator=(RateHistory&&) = delete;
  void BeginStep() noexcept;  ///< Clears committed/candidate step rates only.
  /// Stages S/dt from accepted-step baseline. Invalid/nonfinite input preserves
  /// committed ring and invalidates candidate; no allocation or retention.
  StatusCode Preview(double increment_maximum, double dt, Diagnostic* d);
  void Commit() noexcept;  ///< Requires CanCommit(); scalar-only.
  void DiscardCandidate() noexcept;
  bool CanCommit() const;
  int AcceptedCount() const;
  /// True only with3 accepted rates and current_dt*r_i < threshold for all3.
  /// Both arguments finite, dt>0, threshold>=0. Error preserves *out.
  StatusCode AllRecentBelow(double current_dt, double threshold, bool* out,
                            Diagnostic* d) const;

 private:
  RateHistory() = default;
  int count_ = 0, candidate_count_ = 0;
  double rates_[3] = {}, candidate_rates_[3] = {};
  bool staged_ = false;
};
}  // namespace fesa
#endif  // FESA_ANALYSIS_THERMAL_HISTORY_H_
