/// @file
/// @brief Read-only integration motion owned and committed by AnalysisState.
#ifndef FESA_ANALYSIS_DYNAMIC_STATE_H_
#define FESA_ANALYSIS_DYNAMIC_STATE_H_

#include "fesa/math/vector.h"

namespace fesa {
class NonlinearFieldHistory;
/// Global Cartesian v/a and times [T]. Three-DOF rows store translation
/// [L/T,L/T^2]; six-DOF rows append angular motion [rad/T,rad/T^2], in the
/// owning AnalysisState's physical node/component order. Angular values are
/// small-rotation derivatives or finite spatial motion, selected at activation.
/// HHT acceleration entries are algorithmic integration variables, not a
/// separately reconstructed physical end-equilibrium acceleration. Initial
/// physical acceleration seeds them; a continuous step inherits their values.
/// No finite rotation configuration, equation, force, or integration scratch.
/// Callers cannot
/// construct, stage, discard, or independently commit this child object.
class DynamicState final {
 public:
  ~DynamicState();
  DynamicState(const DynamicState&) = delete;
  DynamicState& operator=(const DynamicState&) = delete;
  DynamicState(DynamicState&&) = delete;
  DynamicState& operator=(DynamicState&&) = delete;
  /// Null until the first joint physical commit; borrowed until owner deletion.
  /// Values may change at the next joint commit. Initial zero storage is
  /// invalid.
  const Vector* CommittedVelocity() const;
  const Vector* CommittedAcceleration() const;
  /// Null unless full kinematics staging succeeded. Borrow expires on restage,
  /// discard, joint commit, or owner deletion; caller never retains it.
  const Vector* CandidateVelocity() const;
  const Vector* CandidateAcceleration() const;
  /// Optional force/moment R11 history, borrowed until owner deletion. Null
  /// before explicit activation. Candidate access follows joint staging.
  const NonlinearFieldHistory* ConvergenceHistory() const;
  /// Committed times; meaningful only when committed getters are nonnull.
  double StepTime() const;
  double TotalTime() const;
  /// Staged times; meaningful only when candidate getters are nonnull.
  double CandidateStepTime() const;
  double CandidateTotalTime() const;

 private:
  friend class AnalysisState;
  DynamicState() = default;
  /// Allocates four fixed stores; failure nulls output and frees partial
  /// stores.
  static StatusCode Create(int size, DynamicState** out, Diagnostic* d);
  /// Copies finite same-size input vectors, distinct from internal candidates.
  /// Failure invalidates candidate only; committed values are never modified.
  StatusCode Stage(const Vector* velocity, const Vector* acceleration,
                   double step_time, double total_time, Diagnostic* d);
  /// Allocation/backend-free validity change and joint promotion respectively.
  void Discard() noexcept;
  void Commit() noexcept;
  /// True when optional staged history matches the output transaction.
  bool HistoryMatches(int step, int increment) const;
  /// Owned for linear dynamics; borrowed from parent for finite dynamics.
  /// Parent destroys this child first, so the borrowed history remains alive.
  NonlinearFieldHistory* convergence_history_ = nullptr;
  bool owns_history_ = true;  ///< Delete history only in the owning role.
  int history_step_ = 0;          ///< Zero when history candidate is invalid.
  bool history_initial_ = false;  ///< Initial frame, not an increment.
  Vector* committed_velocity_ = nullptr;      ///< Owned.
  Vector* committed_acceleration_ = nullptr;  ///< Owned.
  Vector* candidate_velocity_ = nullptr;      ///< Owned.
  Vector* candidate_acceleration_ = nullptr;  ///< Owned.
  double step_time_ = 0, total_time_ = 0;
  double candidate_step_time_ = 0, candidate_total_time_ = 0;
  bool committed_valid_ = false, candidate_valid_ = false;
};
}  // namespace fesa
#endif  // FESA_ANALYSIS_DYNAMIC_STATE_H_
