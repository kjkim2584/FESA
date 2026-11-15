/// @file
/// @brief Single-temperature-field Newton roles and synchronous attempt report.
#ifndef FESA_ANALYSIS_THERMAL_NEWTON_H_
#define FESA_ANALYSIS_THERMAL_NEWTON_H_
#include "fesa/analysis/heat_flux_residual.h"
#include "fesa/analysis/linear_solver.h"
#include "fesa/analysis/newton_controls.h"
namespace fesa {
struct FieldControls;
/// Analysis owns physical candidate and all equation/correction storage.
/// Newton borrows this role for Solve only; no output or commit is hidden here.
class ThermalNonlinearProblem {
 public:
  virtual ~ThermalNonlinearProblem() = default;
  /// Recomputes at current candidate; completed_corrections==0 has no actual
  /// correction. Stages the norm from committed history, without committing.
  /// Borrows out; failure invalidates measurements, never committed state.
  virtual StatusCode EvaluateResidual(int completed_corrections,
                                      HeatFluxResidual* out, Diagnostic* d) = 0;
  /// Builds current J_ff*correction=r_f. Outputs borrow owned distinct system
  /// and writable correction until their release callbacks. Failure nulls both
  /// outputs and leaves partial storage available to attempt cleanup.
  virtual StatusCode AssembleTangent(const EquationSystem** system,
                                     Vector** correction, Diagnostic* d) = 0;
  /// Applies free correction to candidate only; failure invalidates candidate.
  /// Input is borrowed through this call; retains correction assessment until
  /// the following residual measurement consumes it.
  virtual StatusCode ApplyCorrection(const Vector* correction,
                                     Diagnostic* d) = 0;
  /// Idempotent no-allocation/backend cleanup, preserves physical state.
  virtual void ReleaseEquationSystem() noexcept = 0;
  /// Releases the free correction after ApplyCorrection consumes it.
  virtual void ReleaseFreeCorrection() noexcept = 0;
  /// Releases assessment storage after EvaluateResidual consumes it.
  virtual void ReleaseCorrectionAssessment() noexcept = 0;
  /// Releases all partial attempt scratch on every Solve exit, preserving the
  /// candidate and committed state. No allocation/backend work is allowed.
  virtual void ReleaseAttemptResources() noexcept = 0;
};
/// Only the analysis decides retry size and whether an attempt is committed.
enum class ThermalNewtonOutcome {
  kNotStarted,
  kInProgress,
  kConverged,
  kDivergence,
  kSlowConvergence,
  kIterationLimit,
  kFatal
};
/// Scalar metadata at the last evaluated or failed correction, model units.
struct ThermalIterationInfo {
  int iteration = 0;
  int linear_solver_calls = 0;
  int residual_evaluations = 0;
  double phi = 0;
  double correction_max = 0;
  double increment_max = 0;
  ThermalResidualBranch branch = ThermalResidualBranch::kNoDofs;
  bool residual_valid = false;
  bool estimate_available = false;
  double estimate_max = 0;
  bool prediction_available = false;
  double predicted_total = 0;  ///< May be infinity; no value if unavailable.
  StatusCode status = StatusCode::kOk;
  StatusCode numerical_status =
      StatusCode::kOk;  ///< Independent numerical result.
  StatusCode logging_status = StatusCode::kOk;  ///< First observer failure.
  bool logging_failed = false;
  ThermalNewtonOutcome outcome = ThermalNewtonOutcome::kNotStarted;
};
/// Caller owns report and its scalar child records, no model/state buffers.
class ThermalIterationReport final {
 public:
  /// Creation-only output; failure nulls it and frees partial children.
  static StatusCode Create(ThermalIterationReport** out, Diagnostic* d);
  ~ThermalIterationReport();
  ThermalIterationReport(const ThermalIterationReport&) = delete;
  ThermalIterationReport& operator=(const ThermalIterationReport&) = delete;
  ThermalIterationReport(ThermalIterationReport&&) = delete;
  ThermalIterationReport& operator=(ThermalIterationReport&&) = delete;
  void Reset() noexcept;  ///< Scalar-only, invalidates borrowed report data.
  const ThermalIterationInfo* Info()
      const;  ///< Borrowed until Reset/destruction.
  const HeatFluxResidual* Residual()
      const;  ///< Same lifetime, check residual_valid.
  /// First logging error, borrowed until Reset/destruction; null if absent.
  const Diagnostic* LoggingDiagnostic() const;
  /// Copies the first failure without allocation. Stable diagnostic string
  /// lifetimes apply. Numerical status and physical state remain unchanged.
  void RecordLoggingFailure(StatusCode status, const Diagnostic* d) noexcept;

 private:
  friend class ThermalNewton;
  ThermalIterationReport() = default;
  ThermalIterationInfo* info_ = nullptr;      ///< Owned.
  HeatFluxResidual* residual_ = nullptr;      ///< Owned.
  Diagnostic* logging_diagnostic_ = nullptr;  ///< Owned first error copy.
};
/// Optional synchronous observer, borrowing reports for the call only.
class ThermalNewtonProgress {
 public:
  virtual ~ThermalNewtonProgress() = default;
  /// Called after an actual decision, or once on a failed numerical stage.
  /// Never mutate candidate/state/history, Preview/Commit or retain pointers.
  /// Failure aborts; a simultaneous numerical primary stays in report/d.
  virtual StatusCode OnIteration(const ThermalIterationReport* report,
                                 Diagnostic* d) = 0;
};
/// Stateless loop, with formulation7's single-field criteria and stop order.
class ThermalNewton final {
 public:
  /// Caller owns success. Failure nulls *out; allocates no numerical buffers.
  static StatusCode Create(ThermalNewton** out, Diagnostic* d);
  ~ThermalNewton() = default;
  ThermalNewton(const ThermalNewton&) = delete;
  ThermalNewton& operator=(const ThermalNewton&) = delete;
  ThermalNewton(ThermalNewton&&) = delete;
  ThermalNewton& operator=(ThermalNewton&&) = delete;
  /// Borrows all inputs for this call. Resets report/history per attempt; free
  /// DOFs require at least one actual correction, no-DOF trials may accept i0.
  /// Releases attempt scratch on every exit, retains converged candidate and
  /// scalar reports. kEquilibriumFailure is retryable; other failures are
  /// fatal. IC is prediction threshold, ISc protection cap. No
  /// output/commit/discard, increment selection, or mechanical small-increment
  /// exception occurs.
  StatusCode Solve(ThermalNonlinearProblem* problem, LinearSolver* solver,
                   const FieldControls* controls, const NewtonControls* counts,
                   ThermalIterationReport* report, Diagnostic* d,
                   ThermalNewtonProgress* progress = nullptr);

 private:
  ThermalNewton() = default;
};
}  // namespace fesa
#endif  // FESA_ANALYSIS_THERMAL_NEWTON_H_
