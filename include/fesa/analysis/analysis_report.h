/// @file
/// @brief Caller-owned execution metadata with explicitly borrowed fields.
#ifndef FESA_ANALYSIS_ANALYSIS_REPORT_H_
#define FESA_ANALYSIS_ANALYSIS_REPORT_H_

#include <cstdint>

#include "fesa/analysis/residual_report.h"
#include "fesa/base/status_code.h"

namespace fesa {
struct HeatFluxResidual;

/// Scalar report, never an owner of physical state or field result storage.
/// The caller owns this object. A valid residual pointer borrows the current
/// Analysis and expires before its next calculation or destruction. Progress
/// receivers may inspect this report only during their synchronous callback.
struct AnalysisReport {
  /// Invalidates borrowed fields and resets metadata without allocation.
  /// Does not delete residual data or modify any analysis/model state.
  void Reset() noexcept;

  StatusCode status = StatusCode::kOk;  ///< Final step execution status.
  StatusCode calculation_status = StatusCode::kOk;  ///< Numerical status.
  bool calculation_finished = false;  ///< Numerical work has returned.
  bool residual_valid = false;  ///< Permits access to common field pointers.
  bool step_completed = false;  ///< Final increment output and commit done.
  const ResidualReport* residual = nullptr;  ///< Borrowed, never deleted here.
  /// Borrowed force/moment summaries, valid only when residual_valid is true;
  /// expire before the next calculation/destruction. May alias NewtonReport
  /// without a ResidualReport or another reaction vector allocation.
  const FieldResidual* displacement_field = nullptr;
  const FieldResidual* rotation_field = nullptr;
  const HeatFluxResidual* heat_flux_field = nullptr;  ///< Borrowed thermal field.
  /// Optional first progress error, independent of residual_valid and the
  /// primary numerical status. Borrowed until next SolveStep/destruction;
  /// Reset clears this reference only. Diagnostic string lifetimes apply.
  const Diagnostic* logging_diagnostic = nullptr;
  double characteristic_length = 0;  ///< Initial model length, L units.
  int equation_count = 0;  ///< Actual free equations, independent of solves.
  int step_number = 0;  ///< Current one-based step, zero before execution.
  int increment_number = 0;  ///< Current one-based increment.
  int attempt_count = 0;  ///< Actual increment attempts begun.
  int iteration_count = 0;  ///< Newton iterations; zero for direct linear work.
  int linear_solver_calls = 0;  ///< Actual reduced solve calls.
  int last_committed_step = 0;  ///< Step of the latest committed increment.
  int last_committed_increment = 0;  ///< Latest committed increment ID.
  std::int64_t last_committed_frame = 0;  ///< Run-wide committed frame ID.
  double last_committed_total_time = 0;  ///< Committed physical model time.
  double assembly_seconds = 0;  ///< Inclusive assembly wall seconds.
  double solve_seconds = 0;  ///< Reduced solve wall seconds.
  double equilibrium_recovery_seconds = 0;  ///< Inclusive check/recovery time.
  double element_dispatch_seconds = 0;  ///< Nested element-dispatch wall time.
  double calculation_seconds = 0;  ///< Numerical calculation wall time.
  double output_seconds = 0;  ///< Frame write/flush/commit wall time.
  int observed_element_team_max = 0;  ///< Observed element team maximum.
};

}  // namespace fesa

#endif  // FESA_ANALYSIS_ANALYSIS_REPORT_H_
