/// @file
/// @brief Transactional global Cauchy stress history for isolated Solid IPs.
#ifndef FESA_ANALYSIS_INCREMENTAL_SOLID_HISTORY_H_
#define FESA_ANALYSIS_INCREMENTAL_SOLID_HISTORY_H_
#include <map>
#include "fesa/math/matrix.h"
#include "fesa/model/domain.h"
namespace fesa {
/// AnalysisState-owned ng×6 committed/candidate tensors, global tensor order
/// [11,22,33,12,13,23], F/L². Immutable element/IP mapping, initially zero.
/// Each worker stages a distinct element slot; coordinator prepares/validates/
/// discards/commits only after all workers join. No shared worker flag writes.
class IncrementalSolidHistory {
 public:
  /// Borrows Domain during creation. Requires exclusively incremental solids;
  /// copies layout. Caller owns success. Failure cleans and nulls out.
  static StatusCode Create(const Domain* domain, IncrementalSolidHistory** out,
                           Diagnostic* d);
  ~IncrementalSolidHistory();
  IncrementalSolidHistory(const IncrementalSolidHistory&) = delete;
  IncrementalSolidHistory& operator=(const IncrementalSolidHistory&) = delete;
  IncrementalSolidHistory(IncrementalSolidHistory&&) = delete;
  IncrementalSolidHistory& operator=(IncrementalSolidHistory&&) = delete;
  /// Invalidates all candidate slots, preserving committed tensors.
  StatusCode PrepareCandidate(Diagnostic* d);
  /// Invalidates candidates, preserving committed values; after worker join.
  void DiscardCandidate() noexcept;
  /// Copies finite matching ng×6 into its distinct prepared slot. Failure
  /// invalidates that slot; coordinator discards all on any worker failure.
  /// Restaging invalidates this slot's validation, never alters committed data.
  StatusCode StageStress(int element_id, const Matrix* stress, Diagnostic* d);
  /// Requires all slots staged/finite and marks each validated. No allocation
  /// or numeric copy. Failure forbids commit and preserves committed data.
  StatusCode ValidateCandidate(Diagnostic* d);
  /// Checks validation flags without allocation, modification or diagnostics.
  bool CanCommit() const noexcept;
  /// Requires CanCommit; pointer swap only, after successful official flushes.
  /// No allocation/fallible operations. Invalidates candidate views.
  void CommitCandidate() noexcept;
  /// Borrowed history-owned tensor until next commit; unknown ID returns null.
  const Matrix* CommittedStress(int element_id) const;
  /// Borrowed staged candidate until prepare/discard/commit; absent returns null.
  const Matrix* CandidateStress(int element_id) const;
 private:
  struct Slot;
  IncrementalSolidHistory() = default;
  std::map<int, Slot*>* slots_ = nullptr;  ///< Owned mapping and slot resources.
  bool prepared_ = false;  ///< Coordinator writes only outside dispatch.
};
}  // namespace fesa
#endif  // FESA_ANALYSIS_INCREMENTAL_SOLID_HISTORY_H_
