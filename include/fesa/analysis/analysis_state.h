/// @file
/// @brief Committed/candidate static state with optional finite rotations.
#ifndef FESA_ANALYSIS_ANALYSIS_STATE_H_
#define FESA_ANALYSIS_ANALYSIS_STATE_H_

#include "fesa/analysis/dynamic_state.h"
#include "fesa/analysis/result_state.h"
#include "fesa/math/vector.h"
namespace fesa {
class DofManager;
class RotationState;
class NonlinearFieldHistory;
class IncrementalSolidHistory;
class BeamInternalState;
struct ControlSettings;
class Domain;
/// Owns preallocated candidate/committed total U/UR and reaction stores (E4).
/// Each vector has GlobalCount() physical components, original node-ID order;
/// global axes,
/// translations length, rotations radians, forces force and moments
/// force*length. Object results remain a single set: failure never promises
/// their restoration.
class AnalysisState {
 public:
  /// Positive node_count required. Initial committed totals/reactions are zero
  /// storage, step/inc=0 and no_result. Caller owns *out; failure nulls it.
  static StatusCode Create(int node_count, AnalysisState** out, Diagnostic* d);
  /// Copies physical layout from the borrowed step mapping; does not retain it.
  /// Caller owns success; failure nulls output and cleans partial storage.
  static StatusCode Create(const DofManager* dofs, AnalysisState** out,
                           Diagnostic* d);
  /// Copies immutable physical node layout before step execution begins.
  /// Borrows Domain only for this call, retains no mapping/definition pointer.
  /// Caller owns success; failure nulls output and cleans partial storage.
  static StatusCode CreateForDomain(const Domain* domain, AnalysisState** out,
                                    Diagnostic* d);
  /// Nonlinear sixDOF layout in immutable node order. Also owns finite q/a;
  /// linear factories above allocate no rotation store. Same creation
  /// ownership.
  static StatusCode CreateNonlinear(int node_count, AnalysisState** out,
                                    Diagnostic* d);
  /// Copies immutable Domain IDs/offsets, borrows Domain only for this call.
  /// Owns history and optional rotations; caller deletes success. Failure
  /// nulls *out and releases all partial storage.
  static StatusCode CreateNonlinearForDomain(const Domain* domain,
                                            AnalysisState** out, Diagnostic* d);
  /// Candidate validity, independent of rotations and output readiness.
  bool HasValidCandidate() const;
  /// Checks a borrowed step mapping against owned Domain IDs/physical offsets.
  /// Ignores equation/constraint changes. Requires a Domain-created state;
  /// mismatch/null returns invalid argument without changing any state.
  StatusCode ValidateMapping(const DofManager* dofs, Diagnostic* d) const;
  ~AnalysisState();
  AnalysisState(const AnalysisState&) = delete;
  AnalysisState& operator=(const AnalysisState&) = delete;
  AnalysisState(AnalysisState&&) = delete;
  AnalysisState& operator=(AnalysisState&&) = delete;
  /// Number of physical displacement/reaction components, immutable.
  int GlobalCount() const;
  /// Copies committed totals to candidate and zeros candidate reactions via
  /// common math API. Failure invalidates candidate, preserves committed
  /// buffers, step/inc and object-result validity. Optional Beam eta is copied
  /// while both physical/internal motion candidates are invalidated.
  /// No Node/Element result is restored.
  StatusCode PrepareCandidate(Diagnostic* d);
  /// Actual nonlinear attempts also prepare transactional per-step field
  /// time history. Static role only: rejects enabled dynamic storage; dynamic
  /// callers use PrepareCandidate/PrepareDynamicConvergenceHistory. Requires
  /// nonlinear field history and matching committed step tags. Same/next
  /// positive step allowed; controls borrowed for this call.
  /// Any failure discards all candidates, preserves every committed store.
  StatusCode PrepareNonlinearCandidate(int step,
                                       const ControlSettings* controls,
                                       Diagnostic* d);
  /// Read-only optional field time history; null for linear factories.
  /// PreviousMaximumMean is a per-step average; CommittedMaximumFlux below
  /// retains its separate legacy meaning of last successful current maximum.
  const NonlinearFieldHistory* FieldHistory() const;
  /// Borrowed owned IP history for exclusively incremental Solid models;
  /// null otherwise. Survives steps, owned until state destruction. Workers
  /// stage distinct candidate slots only; coordinator controls transactions.
  IncrementalSolidHistory* SolidHistory() const;
  /// Optional read-only owned B33 internal eta and motion; null without Beam3D.
  /// Borrowed until owner destruction. Candidate views follow parent validity.
  const BeamInternalState* BeamInternals() const;
  /// Adds B33 internal eta corrections [L] in element-ID/end order. Requires
  /// prepared nonlinear Beam state and a distinct finite matching vector.
  /// Invalidates dynamic staging/output projection. Any failure discards all
  /// candidates and preserves committed state; input is borrowed for the call.
  StatusCode ApplyBeamInternalCorrection(const Vector* correction,
                                         Diagnostic* d);
  /// Stages both raw time means and current maxima from the same committed
  /// baseline; preview iterations/retries never increment committed counts.
  /// Requires PrepareNonlinearCandidate. Scalar output pointers required and
  /// distinct; errors preserve outputs/committed and discard all candidates.
  StatusCode PreviewNonlinearFluxHistory(double force_mean, double moment_mean,
                                         double force_max, double moment_max,
                                         double* force_raw, double* moment_raw,
                                         Diagnostic* d);
  /// Enables optional linear 3N motion storage before preparing a candidate.
  /// Preserves static preload/results on success or failure. Existing storage
  /// is retained. Rotation/nonlinear layouts are explicitly unsupported here.
  /// Allocation failure leaves Dynamics null; no physical initialization
  /// occurs.
  StatusCode EnableTranslationalDynamics(Diagnostic* d);
  /// Enables linear motion storage in this state's physical 3/6DOF mapping.
  /// Six-DOF rows are [v1,v2,v3,omega1,omega2,omega3] and analogous a/eta,
  /// global axes, units L/T,rad/T and L/T^2,rad/T^2. Small rotations only:
  /// omega=UR_dot; no finite-rotation interpretation. Rejects nonlinear state
  /// or prepared candidate. Same allocation/preload preservation as above.
  StatusCode EnableLinearDynamics(Diagnostic* d);
  /// Enables spatial v/omega and a/eta in a nonlinear all-sixDOF state.
  /// Units/layout follow EnableLinearDynamics, but angular values are spatial
  /// finite-rotation motion, never outputUR derivatives. Preserves preload,
  /// rotations/results on activation failure; active candidate is rejected.
  /// Allocation failure leaves Dynamics null. No pose/physical initialization.
  /// Beam state also allocates internal motion stores; activation alone does
  /// not enable Beam solver execution or supply initial internal motion.
  StatusCode EnableFiniteDynamics(Diagnostic* d);
  /// Activates force/moment R11 history for the physical layout, before
  /// any candidate or dynamic commit. A preceding linear static step must
  /// have exactly one successful increment. Finite dynamics instead borrows
  /// this state's existing nonlinear history, retaining every preload
  /// increment; it is neither copied nor independently owned by DynamicState.
  /// Failure preserves all stores.
  StatusCode EnableDynamicConvergenceHistory(Diagnostic* d);
  /// Prepares R11 statistics AFTER kinematics staging. Initial requires next
  /// step, zero step time and unchanged total time (first total time zero).
  /// After linear static preload the caller supplies its cumulative total
  /// time; the single previous increment seeds the new history. Nonlinear
  /// preload uses its existing history and all successful increments instead.
  /// End requires current step and strictly increasing step/total time.
  /// Initial maximum is physical unweighted force [F]; ignored for end.
  /// Controls are borrowed only for this call. Failure discards all candidates
  /// and preserves committed state. Restaging kinematics invalidates history.
  /// This force-only overload rejects rotational layouts.
  StatusCode PrepareDynamicConvergenceHistory(int step,
                                              const ControlSettings* controls,
                                              bool initial,
                                              double initial_maximum,
                                              Diagnostic* d);
  /// Two-field variant, initial maxima force [F] and moment [FL]. Separate
  /// P12 baselines/seeds, never mixed. Beam state requires staged internal
  /// motion too. Same staging/time/failure contract.
  StatusCode PrepareDynamicConvergenceHistory(int step,
                                              const ControlSettings* controls,
                                              bool initial,
                                              double initial_force_maximum,
                                              double initial_moment_maximum,
                                              Diagnostic* d);
  /// Previews end-state force history from the same committed baseline on
  /// every call; no iteration counting. Mean/max are physical force [F].
  /// Required raw output is written only on success. Failure discards all
  /// candidates, preserves output/commit. Initial candidates cannot preview.
  /// This force-only overload rejects rotational layouts.
  StatusCode PreviewDynamicFluxHistory(double force_mean, double force_max,
                                       double* force_raw, Diagnostic* d);
  /// Two-field preview [F,FL]. Required distinct raw outputs, unchanged on
  /// failure; same joint candidate discard/committed preservation as above.
  StatusCode PreviewDynamicFluxHistory(double force_mean, double moment_mean,
                                       double force_max, double moment_max,
                                       double* force_raw, double* moment_raw,
                                       Diagnostic* d);
  /// Borrowed read-only child, or null before successful activation.
  const DynamicState* Dynamics() const;
  /// Stages computed global v/a and times [T] with prepared u/RF. Units/layout
  /// follow the chosen EnableLinearDynamics/EnableFiniteDynamics role,
  /// including angular components on 6DOF rows.
  /// Inputs must be finite, same physical size, and not alias internal dynamic
  /// candidates. Times are finite, 0<=step_time<=total_time, total time cannot
  /// decrease from committed motion. Failure discards ALL candidates and
  /// preserves committed state. Does not solve initial equilibrium or commit.
  /// Restaging also invalidates staged Beam internal motion and flux history.
  StatusCode StageDynamicKinematics(const Vector* velocity,
                                    const Vector* acceleration,
                                    double step_time, double total_time,
                                    Diagnostic* d);
  /// Copies internal eta velocity [m/s] and acceleration [m/s^2], size 2E in
  /// element-ID order, connectivity end a then b. Requires prepared Beam state
  /// and staged physical dynamic motion. Inputs are borrowed for this call,
  /// mutually distinct and nonalias with every child-owned eta/motion store.
  /// Restaging invalidates dynamic convergence-history readiness. Success
  /// writes candidate motion only; failure discards all parent/child candidates
  /// and preserves committed buffers, times and result tags. No time/solver
  /// pointer is retained and no initial velocity/acceleration is synthesized.
  StatusCode StageBeamInternalKinematics(const Vector* velocity,
                                         const Vector* acceleration,
                                         Diagnostic* d);
  /// Invalidates candidate access only; preserves committed state and current
  /// result validity/tags. No allocation/copy. Next Prepare restores committed.
  void DiscardCandidate() noexcept;
  /// Sequential single-increment convenience: marks last-success-step+1, inc1.
  /// Requires prepared candidate and last-success-step < INT_MAX. S8 uses the
  /// explicit overload below to identify the current output attempt.
  void MarkResultsUnconfirmed();
  /// Requires prepared candidate and positive step/increment identifiers.
  /// Dynamic initial frame may use increment=0 after kinematics staging.
  /// Dynamic Beam state requires complete physical and internal motion.
  /// Dedicated nonlinear attempts also require the prepared history step.
  /// Call immediately before result mutation. Tags the common object result set
  /// unconfirmed without changing committed buffers or last-success metadata.
  void MarkResultsUnconfirmed(int step_number, int increment_number);
  /// Requires prepared candidate, unconfirmed results, positive step/inc
  /// matching result tags and the prepared history transaction step (if used),
  /// including validated complete Solid IP history when present,
  /// and complete Beam eta/physical/internal motion when dynamics is active,
  /// with increment=0 permitted only for staged dynamic initial frames,
  /// and successful output/flush. Exchanges
  /// preallocated stores/metadata only; no allocation/copy/backend call.
  /// Invalid call is a programmer error.
  void Commit(int step_number, int increment_number) noexcept;
  /// Borrowed mutable linear trial buffer only while candidate is prepared.
  /// Nonlinear states always return null; use translation/rotation methods.
  /// Null after creation, failed preparation, discard or commit; do not retain
  /// across PrepareCandidate/DiscardCandidate/Commit. Callers must discard
  /// candidate after any failure.
  Vector* CandidateDisplacement();
  /// Optional finite state, borrowed read-only until destruction. Nonlinear
  /// trial access obeys its readiness; null for all linear factories.
  const RotationState* Rotations() const;
  /// Nonlinear correction in physical node-major order, length GlobalCount().
  /// Three-DOF rows contain du; six-DOF rows contain [du,dw].
  /// Physical translations add; rotations compose, never add to outputUR.
  /// Pose updates invalidate physical/internal motion, history and Solid IPs.
  /// Inputs are borrowed/distinct. Any failure discards candidate
  /// and preserves committed state; required branch flag is true only for
  /// recoverable chart/path failure. Success invalidates projected outputUR.
  StatusCode ApplyNonlinearCorrection(const Vector* correction,
                                      bool* branch_exceeded, Diagnostic* d);
  /// Nonlinear translation3 in global length units, zero-based node index.
  /// Invalidates physical/internal motion, history and Solid stress. Failure
  /// discards candidate; gather failure preserves state and follows
  /// common math output failure contract. No argument is retained.
  StatusCode SetCandidateTranslation(int node, const Vector* translation,
                                     Diagnostic* d);
  StatusCode GatherCandidateTranslation(int node, Vector* out,
                                        Diagnostic* d) const;
  /// Nonlinear unit wxyz target, recomputes body history from committed state.
  /// Requires a six-DOF row. Invalidates physical/internal motion and history.
  /// Same failure/branch contract; caller guards
  /// unwrapped command path.
  StatusCode SetCandidateOrientation(int node, const Vector* q,
                                     bool* branch_exceeded, Diagnostic* d);
  /// Prepares the physical GlobalCount() output vector. Six-DOF rows project
  /// candidate q to principalUR; three-DOF rows retain translations only.
  /// Call after numerical convergence and before writer. Failure discards
  /// trial; success makes const CandidateDisplacement available, never mutable
  /// access.
  StatusCode PrepareNonlinearOutput(Diagnostic* d);
  /// Nonlinear mutable access above is always null. This read-only output is
  /// also null until PrepareNonlinearOutput succeeds; any state update makes
  /// the projection unavailable again. Linear access retains its old contract.
  /// Read-only trial totals for the output writer; same validity/lifetime.
  const Vector* CandidateDisplacement() const;
  /// Current trial RF/RM; zeroed on preparation, never accumulated from
  /// history.
  Vector* CandidateReaction();
  /// Read-only trial reactions for the output writer; same validity/lifetime.
  const Vector* CandidateReaction() const;
  /// Last committed totals, borrowed until next commit/destruction.
  const Vector* CommittedDisplacement() const;
  /// Last successful reaction, same lifetime as committed displacement.
  const Vector* CommittedReaction() const;
  /// Stages paired force/moment raw scales and current maxima for the
  /// prepared candidate. Three-component layouts ignore moment arguments and
  /// retain zero inactive moment history, without averaging fake contributions.
  /// Failure preserves trial and committed history.
  StatusCode StageFluxHistory(double force_raw, double moment_raw,
                              double force_max, double moment_max,
                              Diagnostic* d);
  /// Last committed raw scale; field 0 is force, 1 is moment.
  double CommittedRawMean(int field) const;
  /// Last committed current maximum; field 0 is force, 1 is moment.
  double CommittedMaximumFlux(int field) const;
  /// Common Node/Element result validity.
  ResultState Results() const;
  /// Current common result-set step; zero with no_result. Identifies the
  /// attempted result even after failure, independently of last success.
  int ResultStepNumber() const;
  /// Current common result-set increment; same lifetime as its step tag.
  int ResultIncrementNumber() const;
  /// Last successful output step and increment; initially zero.
  int StepNumber() const;
  /// Last successful increment number.
  int IncrementNumber() const;

 private:
  AnalysisState() = default;
  /// Offset for a validated node row; legacy factories have six components.
  int PhysicalOffset(int node) const;
  int node_count_ = 0;  ///< Immutable node count.
  int* node_ids_ = nullptr;  ///< Owned sorted Domain IDs; legacy has none.
  int* node_offsets_ = nullptr;  ///< Owned N+1 physical offsets.
  /// Allocates fixed global stores; failure nulls output and cleans resources.
  static StatusCode CreateLayout(int global_count, bool has_rotation,
                                 AnalysisState** out, Diagnostic* d);
  int nonlinear_history_step_ =
      0;  ///< Prepared history/output transaction step.
  NonlinearFieldHistory* field_history_ =
      nullptr;  ///< Owned optional statistics.
  IncrementalSolidHistory* solid_history_ = nullptr;  ///< Owned optional IPs.
  BeamInternalState* beam_internals_ = nullptr;  ///< Owned B33 eta/motion.
  bool nonlinear_history_required_ =
      false;                            ///< Real attempt needs staged stats.
  RotationState* rotations_ = nullptr;  ///< Owned optional nonlinear q/a state.
  bool nonlinear_output_ready_ = false;  ///< OutputUR matches candidate q.
  bool has_rotation_ = true;             ///< Immutable physical field presence.
  DynamicState* dynamics_ = nullptr;     ///< Owned optional physical motion.
  Vector* committed_u_ = nullptr;        ///< Owned.
  Vector* candidate_u_ = nullptr;        ///< Owned.
  Vector* committed_reaction_ = nullptr;  ///< Owned.
  Vector* candidate_reaction_ = nullptr;  ///< Owned.
  Vector* zero_ = nullptr;  ///< Owned immutable zero source for retry safety.
  double committed_raw_[2] = {.01, .01};  ///< Force/moment unconverted scales.
  double committed_max_[2] = {0, 0};      ///< Most recent force/moment maxima.
  double candidate_raw_[2] = {.01, .01};  ///< Staged pair, same units.
  double candidate_max_[2] = {0, 0};      ///< Staged current maxima.
  bool candidate_ready_ = false;
  ResultState result_state_ = ResultState::kNoResult;
  int result_step_number_ = 0;       ///< Current common result-set identity.
  int result_increment_number_ = 0;  ///< Current common result-set identity.
  int step_number_ = 0;
  int increment_number_ = 0;
};
}  // namespace fesa

#endif  // FESA_ANALYSIS_ANALYSIS_STATE_H_
