/// @file
/// @brief Synchronous increment log metadata; never owns analysis state.
#ifndef FESA_ANALYSIS_INCREMENT_REPORT_H_
#define FESA_ANALYSIS_INCREMENT_REPORT_H_

#include "fesa/analysis/newton_report.h"

namespace fesa {
/// All events occur outside the successful writer/physical commit gate.
enum class IncrementEvent {
  kAttemptStarted,
  kIteration,
  kAttemptFailed,
  kCutback,
  kCommitted
};

/// Extra failed-trial sampling result; skipped/unavailable never means passed.
enum class GeometryInspection {
  kNotRequested,
  kUnavailableCandidate,
  kUnsupported,
  kPassed,
  kGeometry,
  kFailed
};

/// Analysis supplies scalar context for one synchronous callback. No ownership
/// or retained numerical buffers. Times/logical sizes use model time units;
/// characteristic_length uses L. Proposed size is valid only with its flag.
struct IncrementAttemptInfo {
  int step_number = 0, increment_number = 0, attempt_number = 0;
  double step_time = 0, total_time = 0;
  double increment_size = 0, nominal_increment = 0;
  double characteristic_length = 0;
  bool completes_step = false;
  bool predictor_used = false;
  double predictor_ratio = 0;
  /// Thermal-only S=max_free|theta-theta_n| [theta], formulation6.2, available
  /// only after successful free-change measurement. Includes0 for no free DOFs;
  /// never substitutes the prescribed-including convergence increment_max.
  bool thermal_change_available = false;
  double thermal_free_change = 0;
  bool proposed_available = false;
  double proposed_increment = 0;
  GeometryInspection geometry_inspection = GeometryInspection::kNotRequested;
  NewtonOutcome original_outcome = NewtonOutcome::kNotStarted;
  /// Borrowed extra diagnostic through the synchronous attempt/cutback
  /// callback. Separate from the original primary diagnostic; never retained by
  /// sinks.
  const Diagnostic* geometry_diagnostic = nullptr;
};
}  // namespace fesa
#endif  // FESA_ANALYSIS_INCREMENT_REPORT_H_
