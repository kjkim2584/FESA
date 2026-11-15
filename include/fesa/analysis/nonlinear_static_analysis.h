/// @file
/// @brief Load-increment nonlinear static procedure with output/commit control.
#ifndef FESA_ANALYSIS_NONLINEAR_STATIC_ANALYSIS_H_
#define FESA_ANALYSIS_NONLINEAR_STATIC_ANALYSIS_H_

#include "fesa/analysis/analysis.h"
#include "fesa/analysis/newton_raphson.h"

namespace fesa {
class ElementWorkspace;
class SolidReferenceData;
struct IncrementClock;
struct IncrementTarget;
struct NonlinearStepConditionInput;
/// Owns one step's increment controller and numerical scratch. Common Newton
/// owns iterations; this procedure owns load/prescription evaluation, rollback,
/// official recovery and successful writer-flush/physical commit ordering.
class NonlinearStaticAnalysis final : public Analysis,
                                      private NonlinearProblem,
                                      private NewtonProgress {
 public:
  /// Caller deletes success; failure nulls output and cleans partial storage.
  /// Optional pure B33 static policy uses B^T C B for the first correction
  /// (iteration==0, attempted_iteration==1), then the consistent tangent.
  /// Selection is recomputed for each assembly; no tangent/factorization or
  /// selector history is cached between Newton corrections or attempts.
  static StatusCode Create(NonlinearStaticAnalysis** out, Diagnostic* d,
                           bool b33_material_first = false);
  ~NonlinearStaticAnalysis() override;
  NonlinearStaticAnalysis(const NonlinearStaticAnalysis&) = delete;
  NonlinearStaticAnalysis& operator=(const NonlinearStaticAnalysis&) = delete;
  NonlinearStaticAnalysis(NonlinearStaticAnalysis&&) = delete;
  NonlinearStaticAnalysis& operator=(NonlinearStaticAnalysis&&) = delete;
  /// Executes the next nonlinear step in physical three/six-DOF layout.
  /// Nonlinear Solid requires the incremental elastic family. Deprecated StVK
  /// elements are rejected before iteration/output; committed state is kept.
  /// Material-first rejects non-B33 models with kUnsupported and empty models
  /// with kInvalidArgument before iteration/recovery/commit.
  /// Requires preflighted same-Domain conditions, nonlinear state and
  /// Newton/increment controls. Inputs borrowed only for this call. report
  /// fields borrow owned Newton report until next SolveStep/destruction. Every
  /// failed attempt discards candidate, preserves committed state; successful
  /// partial increments stay committed on later failure. No resume of a
  /// partially completed step. Detailed progress is mandatory. Numerical
  /// failure outranks logging error. Output flush immediately precedes no-fail
  /// memory commit; writer is never owned/finalized here. See
  /// nonlinear_analysis_execution_contract.md.
  StatusCode SolveStep(const AnalysisModel* model, const DofManager* dofs,
                       const EffectiveConditions* conditions,
                       const ControlSettings* controls, AnalysisState* state,
                       ResultWriter* writer, std::int64_t first_frame_number,
                       double total_time_before_step, AnalysisReport* report,
                       AnalysisProgress* progress, Diagnostic* d) override;

 private:
  NonlinearStaticAnalysis() = default;
  bool b33_material_first_ = false;  ///< Immutable policy, not attempt state.
  /// Callback contracts follow NonlinearProblem, valid only during SolveStep.
  StatusCode EvaluateResidual(int completed_corrections, FieldResidual* force,
                              FieldResidual* moment, double* length,
                              NonlinearFailure* failure,
                              Diagnostic* d) override;
  /// Owns the published EquationSystem (Kff and RHS) after this callback until
  /// ReleaseEquationSystem(), next assembly, or ClearExecution(). Newton
  /// borrows it through its solve. Selection uses the current Newton report
  /// because the NonlinearProblem callback has no iteration argument.
  StatusCode AssembleTangent(const EquationSystem** system, Vector** correction,
                             NonlinearFailure* failure, Diagnostic* d) override;
  StatusCode ApplyCorrection(const Vector* correction,
                             NonlinearFailure* failure, Diagnostic* d) override;
  /// Idempotent numerical cleanup; physical state and reports remain valid.
  void ReleaseEquationSystem() noexcept override;
  void ReleaseFreeCorrection() noexcept override;
  void ReleaseCorrectionAssessment() noexcept override;
  void ReleaseAttemptResources() noexcept override;
  /// Synchronous Newton observation; primary diagnostic stays separate from d.
  StatusCode OnIteration(const NewtonReport* report, Diagnostic* d) override;
  /// Releases per-call owned scratch and clears borrowers; Newton report lives.
  void ClearExecution() noexcept;
  StatusCode PrepareStep(Diagnostic* d);
  /// Owns immutable Solid tables and Solid/B33 worker contexts through
  /// recovery. First validates material-first policy: empty model is
  /// kInvalidArgument, any non-B33 element is kUnsupported; no numerical/model
  /// state mutation. Failure leaves partial ownership for ClearExecution; never
  /// mutates model.
  StatusCode PrepareSolidWorkspaces(Diagnostic* d);
  StatusCode PrepareAttempt(NonlinearFailure* failure, Diagnostic* d);
  /// Inspects a valid failed candidate before discard; no physical mutation.
  /// Updates failed classification and synchronous attempt metadata; d is the
  /// separate inspection diagnostic, borrowed only until this attempt ends.
  /// Returns fatal inspection status or kOk (including retryable geometry).
  StatusCode InspectFailedCandidate(NewtonIterationInfo* failed, Diagnostic* d);
  StatusCode FinishAttempt(ResultWriter* writer, std::int64_t frame,
                           double total_time_before_step, bool had_cutback,
                           NonlinearFailure* failure, Diagnostic* d);
  StatusCode RunIncrements(ResultWriter* writer, std::int64_t first_frame,
                           double total_time_before_step, Diagnostic* d);
  StatusCode NotifyIncrement(IncrementEvent event,
                             const NewtonIterationInfo* iteration,
                             const Diagnostic* primary, Diagnostic* d);
  /// Reports outside the write/commit gate and retains the first log failure.
  StatusCode NotifyProgress(AnalysisEvent event, Diagnostic* d);
  void UpdateReport() noexcept;

  NewtonRaphson* newton_ = nullptr;        ///< Owned, stateless iteration loop.
  NewtonReport* newton_report_ = nullptr;  ///< Owned, retained report fields.
  LinearSolver* solver_ = nullptr;         ///< Owned general linear solver.
  EquationSystem* equations_ = nullptr;    ///< Owned current Kff/RHS.
  Vector* free_correction_ = nullptr;  ///< Owned current reduced correction.
  Vector* correction_ = nullptr;       ///< Owned actual global correction.
  Vector* residual_ = nullptr;         ///< Owned current global p-fint.
  Vector* external_ = nullptr;         ///< Owned current attempt dead load.
  Vector* force_origins_ = nullptr;    ///< Owned nonzero per-origin flux.
  Vector* moment_origins_ = nullptr;   ///< Owned nonzero per-origin flux.
  Vector* start_u_ = nullptr;  ///< Owned fixed actual step-start u/output UR.
  Vector* start_reaction_ = nullptr;  ///< Owned fixed step-start RF/RM.
  Matrix* start_q_ = nullptr;         ///< Owned fixed actual step-start unit q.
  Vector* predictor_ = nullptr;  ///< Owned last successful physical increment.
  Vector* next_predictor_ = nullptr;  ///< Owned precomputed commit candidate.
  Vector* internal_residual_ = nullptr;  ///< Owned eta residual [F], 2E.
  Vector* internal_correction_ = nullptr;  ///< Owned actual eta correction [L].
  Vector* internal_predictor_ = nullptr;  ///< Owned last committed eta increment.
  Vector* next_internal_predictor_ = nullptr;  ///< Owned proposed eta increment.
  IncrementClock* clock_ = nullptr;   ///< Owned accepted controller metadata.
  IncrementClock* next_clock_ = nullptr;  ///< Owned precomputed next metadata.
  IncrementTarget* target_ = nullptr;     ///< Owned current proposal.
  IncrementAttemptInfo* attempt_ = nullptr;  ///< Owned synchronous log context.
  NonlinearStepConditionInput* path_ =
      nullptr;                            ///< Owned borrowed snapshot view.
  ElementWorkspace** workers_ =
      nullptr;  ///< Owned entries/array, release after join.
  int worker_count_ = 0;
  SolidReferenceData** solid_references_ =
      nullptr;  ///< Owned; outlive workers.
  int solid_reference_count_ = 0;
  const AnalysisModel* model_ = nullptr;  ///< Borrowed during SolveStep only.
  const DofManager* dofs_ = nullptr;      ///< Same call lifetime.
  const EffectiveConditions* conditions_ = nullptr;  ///< Same call lifetime.
  const ControlSettings* controls_ = nullptr;        ///< Same call lifetime.
  AnalysisState* state_ = nullptr;    ///< Borrowed, committed only after flush.
  AnalysisReport* report_ = nullptr;  ///< Borrowed during SolveStep only.
  AnalysisProgress* progress_ = nullptr;  ///< Borrowed synchronous receiver.
  Diagnostic* attempt_diagnostic_ = nullptr;  ///< Borrowed local attempt scope.
  double length_ = 0;           ///< Last residual's characteristic length.
  double element_started_ = 0;  ///< Dispatch counter baseline for this step.
};
}  // namespace fesa
#endif  // FESA_ANALYSIS_NONLINEAR_STATIC_ANALYSIS_H_
