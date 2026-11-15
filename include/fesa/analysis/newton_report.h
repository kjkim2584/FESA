/// @file
/// @brief Owned Newton field diagnostics and synchronous attempt metadata.
#ifndef FESA_ANALYSIS_NEWTON_REPORT_H_
#define FESA_ANALYSIS_NEWTON_REPORT_H_
#include <cstdint>

#include "fesa/analysis/residual_report.h"
namespace fesa {
/// Final attempt category; only Analysis decides cutback and increment size.
enum class NewtonOutcome {
  kNotStarted,
  kInProgress,
  kConverged,
  kDivergence,
  kSlowConvergence,
  kIterationLimit,
  kGeometry,
  kRotationBranch,
  kFatal
};
/// Last completed/failed stage. Timings below are inclusive attempt totals.
enum class NewtonStage {
  kIdle,
  kResidual,
  kCriteria,
  kTangent,
  kLinearSolve,
  kUpdate,
  kFinished
};
/// Read-only borrowed through NewtonReport, never retained by progress sinks.
struct NewtonIterationInfo {
  StatusCode status =
      StatusCode::kOk;  ///< Return status including logging error.
  StatusCode numerical_status = StatusCode::kOk;  ///< Primary numeric status.
  NewtonOutcome outcome = NewtonOutcome::kNotStarted;
  NewtonStage stage = NewtonStage::kIdle;
  int iteration = 0;  ///< Completed corrections at last residual evaluation.
  int attempted_iteration =
      0;  ///< Target correction during tangent/solve/update.
  int linear_solver_calls = 0;
  std::int64_t residual_evaluations = 0;
  bool residual_valid = false;
  bool correction_valid = false;
  bool logging_failed = false;  ///< Abort even if numerical primary retained.
  bool estimate_available[2] = {};
  double estimate_max[2] = {};  ///< Length/radian units; unavailable means NA.
  bool prediction_available[2] = {};
  double predicted_total[2] = {};  ///< May be infinity, unavailable means NA.
  double residual_seconds = 0, tangent_seconds = 0, solve_seconds = 0;
  double update_seconds = 0, criteria_seconds = 0;
};
/// Caller owns report. Children are owned; no physical state/reaction storage.
/// Borrowed getters expire on reset/next Solve or destruction.
class NewtonReport final {
 public:
  /// Creation-only out; failure nulls and releases partial owned metadata.
  static StatusCode Create(NewtonReport** out, Diagnostic* d);
  ~NewtonReport();
  NewtonReport(const NewtonReport&) = delete;
  NewtonReport& operator=(const NewtonReport&) = delete;
  NewtonReport(NewtonReport&&) = delete;
  NewtonReport& operator=(NewtonReport&&) = delete;
  /// Resets scalar metadata only, no allocation/backend/physical mutation.
  void Reset() noexcept;
  /// Borrowed attempt metadata, current until reset/next Solve/destruction.
  const NewtonIterationInfo* Info() const;
  /// Borrowed force-field summary; physical data valid only if residual_valid.
  const FieldResidual* Displacement() const;
  /// Borrowed moment-field summary, same validity and lifetime as force.
  const FieldResidual* Rotation() const;
  /// Borrows the first progress failure independently of numerical status.
  /// Null until a logging failure; expires on Reset/next Solve/destruction.
  /// String fields retain Diagnostic's static/stable-owner lifetime contract.
  const Diagnostic* LoggingDiagnostic() const;
  /// Retains the first failed progress status/metadata without allocation.
  /// Success is a no-op; never changes numerical status or physical state.
  /// Diagnostic strings must outlive reporting through this report.
  void RecordLoggingFailure(StatusCode status,
                            const Diagnostic* diagnostic) noexcept;

 private:
  friend class NewtonRaphson;
  NewtonReport() = default;
  NewtonIterationInfo* info_ = nullptr;  ///< Owned metadata.
  FieldResidual* force_ = nullptr;   ///< Owned physical force field summary.
  FieldResidual* moment_ = nullptr;  ///< Owned moment field summary.
  Diagnostic* logging_diagnostic_ = nullptr;  ///< Owned first progress error.
};
/// Optional standalone Newton observer; actual Analysis requires its detailed
/// log path. Calls synchronous; may not retain report or mutate physical state.
class NewtonProgress {
 public:
  virtual ~NewtonProgress() = default;
  /// Called after each residual decision, or once on a failed numerical stage.
  /// Failure aborts attempt; existing numerical primary remains in report/d.
  virtual StatusCode OnIteration(const NewtonReport* report, Diagnostic* d) = 0;
};
}  // namespace fesa
#endif  // FESA_ANALYSIS_NEWTON_REPORT_H_
