/// @file
/// @brief Vendor-independent synchronous step progress reporting boundary.
#ifndef FESA_ANALYSIS_ANALYSIS_PROGRESS_H_
#define FESA_ANALYSIS_ANALYSIS_PROGRESS_H_

#include "fesa/analysis/analysis_report.h"
#include "fesa/analysis/increment_report.h"
#include "fesa/base/diagnostic.h"

namespace fesa {
struct ThermalIterationInfo;

/// Events occur outside the successful WriteFrame/Commit memory gate.
enum class AnalysisEvent { kCalculationFinished, kIncrementCommitted };

/// Receives borrowed reports synchronously on the analysis control thread.
/// The receiver owns neither the report nor its field pointers and must not
/// retain them. It may report logging failure but never commits physical state.
class AnalysisProgress {
 public:
  virtual ~AnalysisProgress() = default;
  /// Thermal-only detailed receiver. Borrows metadata only for this call;
  /// never retains pointers or changes physical state. A numerical primary
  /// remains independent of logging output d. Default is explicitly unsupported.
  virtual StatusCode OnThermalIncrementProgress(
      IncrementEvent, const IncrementAttemptInfo*, const ThermalIterationInfo*,
      const HeatFluxResidual*, const Diagnostic* primary, Diagnostic* d) {
    if (primary && primary == d) return StatusCode::kInvalidArgument;
    return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                         "Detailed thermal progress receiver unavailable", d);
  }
  /// Reads report only during this call. A failed receiver may leave partial
  /// log bytes but must not modify the report or numerical state. Failure
  /// writes d when provided; success preserves d. The caller preserves an
  /// already established numerical failure over a subsequent receiver error.
  virtual StatusCode OnProgress(AnalysisEvent event,
                                const AnalysisReport* report,
                                Diagnostic* d) = 0;
  /// Optional to linear receivers; required for actual nonlinear execution.
  /// Reads borrowed increment/iteration/field metadata only during this call.
  /// Never called inside the successful write/commit gate. primary is read-only
  /// numerical error context, separate from logging-error output d. No physical
  /// mutation or retained pointers. The default fails explicitly so an absent
  /// detailed log receiver cannot silently drop required nonlinear reporting.
  virtual StatusCode OnIncrementProgress(
      IncrementEvent, const IncrementAttemptInfo*, const NewtonIterationInfo*,
      const FieldResidual*, const FieldResidual*, const Diagnostic* primary,
      Diagnostic* d) {
    if (primary && primary == d) return StatusCode::kInvalidArgument;
    return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                         "Detailed nonlinear progress receiver unavailable", d);
  }
};

}  // namespace fesa

#endif  // FESA_ANALYSIS_ANALYSIS_PROGRESS_H_
