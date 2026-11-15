/// @file
/// @brief Owned UTF-8 analysis log with checked writes/flush/close.
#ifndef FESA_IO_ANALYSIS_LOG_H_
#define FESA_IO_ANALYSIS_LOG_H_

#include <string>

#include "fesa/analysis/increment_report.h"
#include "fesa/base/diagnostic.h"
#include "fesa/io/control_settings.h"

namespace fesa {
/// Owns its native stream; destruction closes it. Single-threaded controller.
class AnalysisLog {
 public:
  ~AnalysisLog();
  AnalysisLog(const AnalysisLog&) = delete;
  AnalysisLog& operator=(const AnalysisLog&) = delete;
  AnalysisLog(AnalysisLog&&) = delete;
  AnalysisLog& operator=(AnalysisLog&&) = delete;
  /// Appends borrowed UTF-8 text and flushes. Failure may leave partial bytes;
  /// no solver state is committed. Success preserves d. Write/close diagnostics
  /// borrow this object's owned input path until its destruction.
  StatusCode Write(const char* text, Diagnostic* d);
  /// Records source, all scalar values, automatic flags and unused Rp.
  /// Borrows settings; step 0 means initial values. Same failure contract.
  StatusCode WriteControls(const ControlSettings* settings, int step,
                           Diagnostic* d);
  /// Logs borrowed increment/iteration metadata synchronously, no numerical
  /// evaluation or physical mutation. iteration may be null outside iteration
  /// events; fields are required only for a valid iteration residual. Invalid
  /// fields/estimates print NA. primary is optional read-only numerical failure
  /// context; it must differ from d and outlive this call. Write/flush failure
  /// changes only d and may leave partial log bytes, preserving primary.
  StatusCode WriteIncrementProgress(
      IncrementEvent event, const IncrementAttemptInfo* attempt,
      const NewtonIterationInfo* iteration, const FieldResidual* displacement,
      const FieldResidual* rotation, const Diagnostic* primary, Diagnostic* d);
  /// Checked final close; stream becomes closed even on failure. Idempotent.
  StatusCode Finish(Diagnostic* d);

 private:
  friend class LogStdioAccess;  ///< Internal native-stream verification seam.
  friend StatusCode OpenAnalysisLog(const char*, AnalysisLog**, Diagnostic*);
  explicit AnalysisLog(std::string* input_path);
  void* stream_;             ///< Owned FILE hidden from the public interface.
  std::string* input_path_;  ///< Owned source for write/close diagnostics.
};
/// Opens <input-directory>/<input-base>_analysis.log, replacing that log only.
/// Called after CLI/input path validation. Caller deletes *out; failure nulls
/// it and prevents the run. On failure d.source_file borrows caller input_path
/// with kAnalysisInputFile kind; FormatDiagnostic composes the log filename.
StatusCode OpenAnalysisLog(const char* input_path, AnalysisLog** out,
                           Diagnostic* d);
/// Formats borrowed diagnostic source/context into caller-owned text. Resolves
/// derived cfg/log paths using path_kind without changing d. Allocation errors
/// return kAllocationFailure; output text is then unspecified.
StatusCode FormatDiagnostic(const Diagnostic* d, std::string* text);
}  // namespace fesa
#endif  // FESA_IO_ANALYSIS_LOG_H_
