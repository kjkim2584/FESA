/// @file
/// @brief Owned nominal thermal layout, IC buffers and transactional histories.
#ifndef FESA_ANALYSIS_THERMAL_STATE_H_
#define FESA_ANALYSIS_THERMAL_STATE_H_
#include "fesa/analysis/result_state.h"
#include "fesa/analysis/thermal_history.h"
#include "fesa/math/vector.h"
namespace fesa {
class Domain;
class ThermalDofMap;
/// Owns5 fixed-size Vectors and an optional first-accepted film-temperature
/// Vector, copied node layout and both histories. No Domain
/// or mapping pointers are retained. Global ordering is sorted original
/// ID/label.
class ThermalState final {
 public:
  /// Copies immutable IC to committed T, zero elsewhere including unused nodes.
  /// Explicit IC on a nonparticipant is an input error. Borrows Domain for this
  /// call only; caller owns success. Failure nulls *out and releases
  /// candidates.
  static StatusCode CreateForDomain(const Domain* domain, ThermalState** out,
                                    Diagnostic* d);
  ~ThermalState();
  ThermalState(const ThermalState&) = delete;
  ThermalState& operator=(const ThermalState&) = delete;
  ThermalState(ThermalState&&) = delete;
  ThermalState& operator=(ThermalState&&) = delete;
  /// Compares nominal IDs/components only; changed constraints are allowed.
  /// Borrows mapping, retains none; failure leaves every state buffer intact.
  StatusCode ValidateMapping(const ThermalDofMap* mapping, Diagnostic* d) const;
  /// Nominal layout identity, includes non-participating nodes; constraints
  /// irrelevant. Borrows Domain only for this call, retaining none. Compares
  /// sorted IDs and per-node counts; mismatch kInvalidArgument, no allocation,
  /// numerical library work or state mutation.
  StatusCode ValidateDomainLayout(const Domain* domain, Diagnostic* d) const;
  int GlobalCount() const;
  /// No allocation; copies committed T and immutable zero RFL using Vector API.
  /// Error invalidates candidates/histories, preserves committed/results.
  StatusCode PrepareCandidate(Diagnostic* d);
  bool HasValidCandidate() const;
  /// Candidate/committed T/RFL pointers below expire at next PrepareCandidate,
  /// DiscardCandidate or Commit. Do not cache them across those calls. Mutable
  /// candidate getters return null unless prepared; committed buffers always
  /// remain readable.
  Vector* CandidateTemperature();
  const Vector* CandidateTemperature() const;
  Vector* CandidateReactionFlux();
  const Vector* CandidateReactionFlux() const;
  const Vector* CommittedTemperature() const;
  const Vector* CommittedReactionFlux() const;
  /// Stages candidate T for the first accepted film-table increment. Requires
  /// a prepared candidate. No-op after freezing or without a film property.
  /// Copies before output; failure/unconfirmed output never publishes it.
  /// Borrows d; retains no new inputs. Existing committed freeze is preserved.
  StatusCode StageFirstAcceptedTemperature(Diagnostic* d);
  /// Borrowed first accepted T until next successful BeginStep/destruction.
  /// Null before its successful Commit or without film tables. Used only to
  /// evaluate immutable row h(T); no Domain/condition pointer is retained.
  const Vector* FirstAcceptedTemperature() const;
  HeatFluxHistory* FluxHistory();  ///< Borrowed until state destruction.
  const HeatFluxHistory* FluxHistory() const;
  RateHistory* IncrementRates();  ///< Same lifetime; rates optional in steady.
  const RateHistory* IncrementRates() const;
  /// Requires no prepared attempt. Strong failure guarantee across histories;
  /// success captures carry/seed then resets rate ring and current increment.
  StatusCode BeginStep(int step, const FieldControls* controls, Diagnostic* d);
  void DiscardCandidate() noexcept;  ///< Never restores Node result objects.
  /// Requires prepared candidate, current step and next positive increment.
  /// If first film-table T is unset, StageFirstAcceptedTemperature must
  /// succeed. Call immediately before Node/result writes; no restoration on
  /// flush error.
  void MarkResultsUnconfirmed(int step, int increment);
  /// Requires matching unconfirmed tags, prepared T/RFL, staged heat preview
  /// and successful output flush. Pointer swaps/scalar history commit only.
  /// Commits rate if staged; steady callers need no rate. No
  /// allocation/backend.
  void Commit(int step, int increment) noexcept;
  ResultState Results() const;
  int ResultStepNumber() const;
  int ResultIncrementNumber() const;

 private:
  ThermalState() = default;
  int node_count_ = 0, global_count_ = 0;
  int* nodes_ = nullptr;         ///< Owned sorted node IDs.
  int* offsets_ = nullptr;       ///< Owned node_count+1 component offsets.
  Vector* committed_ = nullptr;  ///< Owned temperature [theta].
  Vector* candidate_ = nullptr;  ///< Owned temperature [theta].
  Vector* committed_reaction_ = nullptr;  ///< Owned RFL [E/T].
  Vector* candidate_reaction_ = nullptr;  ///< Owned RFL [E/T].
  Vector* first_accepted_ = nullptr;  ///< Owned optional film-only T buffer.
  bool first_accepted_set_ = false;
  bool first_accepted_staged_ = false;
  Vector* zero_ = nullptr;           ///< Owned immutable global zero vector.
  HeatFluxHistory* flux_ = nullptr;  ///< Owned.
  RateHistory* rates_ = nullptr;     ///< Owned.
  bool candidate_ready_ = false;
  int step_ = 0, accepted_increment_ = 0;
  ResultState result_state_ = ResultState::kNoResult;
  int result_step_ = 0, result_increment_ = 0;
};
}  // namespace fesa
#endif  // FESA_ANALYSIS_THERMAL_STATE_H_
