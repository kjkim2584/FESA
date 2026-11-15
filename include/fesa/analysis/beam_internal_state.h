/// @file
/// @brief AnalysisState-owned transactional B33 internal length coordinates.
#ifndef FESA_ANALYSIS_BEAM_INTERNAL_STATE_H_
#define FESA_ANALYSIS_BEAM_INTERNAL_STATE_H_

#include "fesa/math/vector.h"

namespace fesa {
class AnalysisState;
class Domain;

/// Owns eta=L*(stretch-1) [L], two per element in increasing original ID,
/// connectivity-end a then b. Initial committed values are zero. Physical
/// node layout/output is separate. The parent alone controls transactions;
/// no element/worker may mutate or independently commit this store.
/// Optional internal v/a are allocated on finite dynamic activation; no Domain,
/// solver or clock is retained and no initial motion is inferred.
/// @see docs/Linear static/b33_nonlinear_contracts.md.
class BeamInternalState final {
 public:
  ~BeamInternalState();
  BeamInternalState(const BeamInternalState&) = delete;
  BeamInternalState& operator=(const BeamInternalState&) = delete;
  BeamInternalState(BeamInternalState&&) = delete;
  BeamInternalState& operator=(BeamInternalState&&) = delete;
  /// Immutable number of elements, each with two independent coordinates.
  int ElementCount() const;
  /// Original element ID for zero-based row; -1 for invalid row.
  int ElementId(int row) const;
  /// First coordinate for an original element ID; -1 if absent.
  int Offset(int element_id) const;
  /// Read-only owned eta, valid through the next parent commit/destruction.
  const Vector* CommittedPosition() const;
  /// Read-only prepared eta; null after discard/commit or failed preparation.
  /// Borrow expires at the next parent mutation; caller never retains it.
  const Vector* CandidatePosition() const;
  /// Internal eta rate [m/s], same element-ID/end order as position. Null
  /// until first joint dynamic commit. Borrow expires on commit/destruction.
  const Vector* CommittedVelocity() const;
  /// Internal eta acceleration [m/s^2], same validity/ownership as velocity.
  const Vector* CommittedAcceleration() const;
  /// Read-only staged motion; null until complete parent-coordinated staging.
  /// Borrow expires on parent mutation, commit/discard, or destruction.
  const Vector* CandidateVelocity() const;
  /// Internal acceleration candidate [m/s^2], same validity as velocity.
  const Vector* CandidateAcceleration() const;

 private:
  friend class AnalysisState;
  BeamInternalState() = default;
  /// Copies IDs from a nonempty exclusive Beam3D Domain; retains no pointers.
  /// Parent owns success; failure nulls out and cleans partial storage.
  static StatusCode Create(const Domain* domain, BeamInternalState** out,
                           Diagnostic* d);
  /// Lazily allocates four motion stores with no active candidate. Failure
  /// preserves positions/committed motion and retains private stores for retry.
  StatusCode EnableDynamics(Diagnostic* d);
  /// Copies finite distinct nonalias 2E inputs into motion candidates. Parent
  /// requires staged physical motion. Failure invalidates internal motion only.
  StatusCode StageMotion(const Vector* velocity, const Vector* acceleration,
                         Diagnostic* d);
  /// Clears candidate motion only, without touching eta or committed stores.
  void DiscardMotion() noexcept;
  /// Checks prepared eta and, for dynamic parents, complete staged motion.
  bool CanCommit(bool require_motion) const;
  /// Copies committed eta to candidate and invalidates motion candidates.
  /// Failure invalidates candidate only; committed values are preserved.
  StatusCode PrepareCandidate(Diagnostic* d);
  /// Adds borrowed finite same-size nonalias correction [L]. Failure
  /// invalidates candidate only; parent then discards its whole transaction.
  StatusCode ApplyCorrection(const Vector* correction, Diagnostic* d);
  /// No allocation/backend calls; preserves all committed coordinates.
  void DiscardCandidate() noexcept;
  /// Requires CanCommit(require_motion). Parent calls after output/flush;
  /// swaps eta and required motion buffers only, then invalidates candidates.
  void CommitCandidate(bool require_motion) noexcept;
  int element_count_ = 0;
  int* element_ids_ = nullptr;             ///< Owned increasing original IDs.
  Vector* committed_position_ = nullptr;  ///< Owned eta [L], size 2*count.
  Vector* candidate_position_ = nullptr;  ///< Owned eta [L], size 2*count.
  Vector* committed_velocity_ = nullptr;      ///< Owned eta rate [m/s], 2E.
  Vector* committed_acceleration_ = nullptr;  ///< Owned eta rate [m/s^2], 2E.
  Vector* candidate_velocity_ = nullptr;      ///< Owned eta rate [m/s], 2E.
  Vector* candidate_acceleration_ = nullptr;  ///< Owned eta rate [m/s^2], 2E.
  bool prepared_ = false;                 ///< Parent coordinator writes only.
  bool committed_motion_valid_ = false;  ///< First joint dynamic commit only.
  bool candidate_motion_valid_ = false;  ///< Both motion copies succeeded.
};
}  // namespace fesa
#endif  // FESA_ANALYSIS_BEAM_INTERNAL_STATE_H_
