/// @file
/// @brief Transactional per-step time-averaged force/moment field statistics.
#ifndef FESA_ANALYSIS_NONLINEAR_FIELD_HISTORY_H_
#define FESA_ANALYSIS_NONLINEAR_FIELD_HISTORY_H_
#include "fesa/io/control_settings.h"
namespace fesa {
/// AnalysisState owns this optional nonlinear history. A linear dynamic owner
/// may also use a separate force-only instance without changing state kind.
/// Statistics are physical
/// scalars, not element scratch or a matrix cache. Single control-thread
/// writer.
/// @see docs/MITC4/field_history_contract.md and abaqus_increment_control.md.
class NonlinearFieldHistory final {
 public:
  /// Initial raw seeds .01, step/counts0. Caller owns success; failure
  /// nulls*out.
  static StatusCode Create(NonlinearFieldHistory** out, Diagnostic* d);
  /// Creates optional moment statistics. Inactive moment values/counts stay
  /// zero and its controls/inputs are not consumed. Caller owns success;
  /// failure nulls *out. No borrowed pointers are retained.
  static StatusCode Create(bool moment_active, NonlinearFieldHistory** out,
                           Diagnostic* d);
  ~NonlinearFieldHistory() = default;
  NonlinearFieldHistory(const NonlinearFieldHistory&) = delete;
  NonlinearFieldHistory& operator=(const NonlinearFieldHistory&) = delete;
  NonlinearFieldHistory(NonlinearFieldHistory&&) = delete;
  NonlinearFieldHistory& operator=(NonlinearFieldHistory&&) = delete;
  /// Same committed step or immediately next positive step, controls borrowed
  /// only for this call. Copies scalar seed/epsilon/qu settings. Any failure
  /// invalidates candidate and preserves committed state. No backend needed.
  StatusCode PrepareCandidate(int step, const ControlSettings* controls,
                              Diagnostic* d);
  /// P12 force-only dynamic initial frame, immediately next step. Caller
  /// establishes direct/continuous start; static-preload transfer is not
  /// supported by this entry point. Borrows controls only during this call.
  /// Stages raw seed/qu and active baseline (M*H+initial_maximum)/(M+1),
  /// without counting an increment or including the initial mean. Moment
  /// histories are unsupported. Failure discards candidate, preserves commit.
  /// Caller commits only after initial recovery/write/flush; Preview is not
  /// allowed on this initial candidate. No allocation/backend work.
  StatusCode PrepareDynamicInitial(int step, const ControlSettings* controls,
                                   double initial_maximum, Diagnostic* d);
  /// Field-wise P12 initial staging, force [F] and moment [FL]. Moment input
  /// is ignored when inactive. Same next-step/controls/failure and no-count
  /// contract as the force-only overload; no cross-field sums or conversion.
  StatusCode PrepareDynamicInitial(int step, const ControlSettings* controls,
                                   double force_maximum, double moment_maximum,
                                   Diagnostic* d);
  /// Fresh force-only history after one successful linear static increment in
  /// previous_step. Stages next-step initial seed from previous_raw [F] and
  /// baseline (previous_maximum+initial_maximum)/2 [F]; q0/qu override seed.
  /// Does not fabricate committed static statistics or increment counts.
  /// Requires no committed history. Controls borrowed for call; failure
  /// discards candidate and preserves committed data. Joint output commit only.
  StatusCode PrepareDynamicInitialAfterLinearStatic(
      int previous_step, const ControlSettings* controls, double previous_raw,
      double previous_maximum, double initial_maximum, Diagnostic* d);
  /// Two-field variant; required arrays contain two entries, force then moment
  /// [F,FL], and are borrowed only through this call. Inactive moment ignored.
  /// Same single-increment preload, ownership and failure contract as above.
  StatusCode PrepareDynamicInitialAfterLinearStatic(
      int previous_step, const ControlSettings* controls,
      const double* previous_raw, const double* previous_maximum,
      const double* initial_maximum, Diagnostic* d);
  /// Recomputes current candidate from the SAME committed baseline each call.
  /// qbar/qmax finite nonnegative, force F, moment F*L in consistent model
  /// units. Two required distinct scalar outputs receive raw decision means,
  /// before field conversion. Any failure preserves outputs/committed,
  /// invalidates candidate. Caller uses only the final converged preview for
  /// commit.
  StatusCode Preview(double force_mean, double moment_mean, double force_max,
                     double moment_max, double* force_raw, double* moment_raw,
                     Diagnostic* d);
  /// Requires CanCommit(); scalar-only/noexcept, no
  /// allocation/callback/backend.
  void Commit() noexcept;
  void DiscardCandidate() noexcept;  ///< Candidate access invalidated only.
  bool CanCommit() const;  ///< A successful preview exists for this candidate.
  /// field index0 force,1 moment required; valid for object lifetime.
  double CommittedRawMean(int field) const;
  double CommittedMaximumMean(int field) const;
  int CommittedIncludedCount(int field) const;
  int CommittedIncrementCount() const;
  int CommittedStepNumber() const;
  /// Requires prepared candidate and field0/1. Previous qmax mean of THIS step;
  /// zero at a new static step, excludes the current preview. A
  /// dynamic initial frame stages a P12 baseline used until first success.
  double PreviousMaximumMean(int field) const;
  /// Requires CanCommit() and field0/1; preview count, not iteration count.
  int CandidateIncludedCount(int field) const;

 private:
  NonlinearFieldHistory() = default;
  bool moment_active_ = true;  ///< Immutable physical field presence.
  int committed_step_ = 0, candidate_step_ = 0;
  int committed_increments_ = 0, candidate_increments_ = 0;
  int committed_n_[2] = {}, candidate_n_[2] = {};
  double committed_mean_[2] = {.01, .01};  ///< Raw positive mean/seed.
  double committed_qmax_[2] = {}, candidate_qmax_[2] = {};
  double candidate_mean_[2] = {};
  double epsilon_[2] = {}, seed_[2] = {}, constant_value_[2] = {};
  bool constant_[2] = {};
  bool prepared_ = false, staged_ = false;
  bool initial_candidate_ = false;  ///< Initial frame, never an increment.
  double committed_initial_maximum_[2] = {};  ///< Per-field P12 baseline.
  double candidate_initial_maximum_[2] = {};  ///< Staged until joint commit.
};
}  // namespace fesa
#endif  // FESA_ANALYSIS_NONLINEAR_FIELD_HISTORY_H_
