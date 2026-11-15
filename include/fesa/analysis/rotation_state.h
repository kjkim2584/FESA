/// @file
/// @brief Candidate/committed finite orientations and accepted body-spin
/// history.
#ifndef FESA_ANALYSIS_ROTATION_STATE_H_
#define FESA_ANALYSIS_ROTATION_STATE_H_

#include "fesa/math/matrix.h"

namespace fesa {
/// AnalysisState owns this optional nonlinear store; linear states omit it.
/// Rows follow immutable physical node order, all nodes sixDOF in this scope.
/// q is unit wxyz initial-global to current-global; a is accumulated body spin
/// in radians. No directors, element-specific axes or numerical caches stored.
class RotationState final {
 public:
  /// Positive count. Initializes committed identity q and zero a; no candidate.
  /// Caller owns success; failure nulls *out and releases partial storage.
  static StatusCode Create(int node_count, RotationState** out, Diagnostic* d);
  ~RotationState();
  RotationState(const RotationState&) = delete;
  RotationState& operator=(const RotationState&) = delete;
  RotationState(RotationState&&) = delete;
  RotationState& operator=(RotationState&&) = delete;
  /// Copies committed q/a to candidate. Any failure invalidates candidate;
  /// committed buffers are unchanged. No input retained, d optional.
  StatusCode PrepareCandidate(Diagnostic* d);
  /// Left-multiplies one candidate q by Exp(spatial_delta[3]), then recomputes
  /// a_trial=a_committed+Log(q_committed^-1*q_trial), never sums Newton steps.
  /// Guards applied path and committed-relative chart, both angle<pi with
  /// 256epsilon relative-scalar guard. Branch sets required flag and returns
  /// kEquilibriumFailure. Any failure invalidates candidate, preserves
  /// committed. Inputs borrowed/distinct from owned stores. Node index
  /// zero-based.
  StatusCode ApplyCorrection(int node, const Vector* spatial_delta,
                             bool* branch_exceeded, Diagnostic* d);
  /// Sets one candidate unit q, preserving cover sign, updates a from
  /// committed. Same failure/branch contract. Caller MUST separately guard
  /// known unwrapped prescribed increment path: endpoint q alone cannot detect
  /// a4pi command.
  StatusCode SetCandidateOrientation(int node, const Vector* q,
                                     bool* branch_exceeded, Diagnostic* d);
  /// Writes principal candidate UR to caller-owned3N vector in row node order.
  /// Does not mutate state; not-ready/shape/alias error preserves out. Backend
  /// failure may partially change out; discard output, committed stays intact.
  StatusCode PrincipalOutput(Vector* out, Diagnostic* d) const;
  /// Requires prepared candidate; exchanges q/a pointers only, no allocation,
  /// library operation or callback. Caller commits together with other state.
  void Commit() noexcept;
  /// Invalidates trial access; no numerical operation, committed unaffected.
  void DiscardCandidate() noexcept;
  int NodeCount() const;
  /// Borrowed until next commit/destruction. Never normalize/canonicalize
  /// caller.
  const Matrix* CommittedOrientations() const;
  const Vector* CommittedBodySpin() const;
  /// Null unless prepared; lifetime ends on failed
  /// update/prepare/discard/commit.
  const Matrix* CandidateOrientations() const;
  const Vector* CandidateBodySpin() const;

 private:
  RotationState() = default;
  int node_count_ = 0;
  bool candidate_ready_ = false;
  Matrix* committed_q_ = nullptr;  ///< OwnedN×4.
  Matrix* candidate_q_ = nullptr;  ///< OwnedN×4.
  Vector* committed_a_ = nullptr;  ///< Owned3N, node-major body components.
  Vector* candidate_a_ = nullptr;  ///< Owned3N.
};
}  // namespace fesa
#endif  // FESA_ANALYSIS_ROTATION_STATE_H_
