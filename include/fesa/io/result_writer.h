/// @file
/// @brief Vendor-independent completed-frame serialization boundary.
#ifndef FESA_IO_RESULT_WRITER_H_
#define FESA_IO_RESULT_WRITER_H_
#include "fesa/analysis/analysis_state.h"
#include "fesa/base/diagnostic.h"
#include "fesa/io/frame_info.h"
#include "fesa/model/domain.h"
namespace fesa {
class ThermalState;
/// Deterministic schema preflight, to be called BEFORE OpenResultWriter by the
/// controller so invalid model data cannot truncate an existing result. Checks
/// counts/byte bounds, IDs/connectivity, UTF-8 identity and supported stored
/// element shapes, without numerical allocation, result recovery or file
/// access. Borrows domain/d for this call; Domain remains unchanged and d
/// records failure only.
/// WriteModel repeats exactly these checks. Programmatic producers supply valid
/// UTF-8 names.
StatusCode ValidateResultModel(const Domain* domain, Diagnostic* d);
/// Caller-owned sequential writer; delete through this virtual interface.
/// Never computes results or changes Domain/AnalysisState. Destruction only
/// cleans up; Finish must succeed explicitly for a successful execution.
class ResultWriter {
 public:
  virtual ~ResultWriter() = default;
  ResultWriter(const ResultWriter&) = delete;
  ResultWriter& operator=(const ResultWriter&) = delete;
  ResultWriter(ResultWriter&&) = delete;
  ResultWriter& operator=(ResultWriter&&) = delete;
  /// Writes model once, in ascending original ID order. Borrows domain until
  /// Finish/destruction; definitions must remain unchanged and alive.
  /// Failure latches this writer failed; later frames are forbidden.
  virtual StatusCode WriteModel(const Domain* domain, Diagnostic* d) = 0;
  /// Writes prepared candidate global U/UR/RF/RM and already recovered element
  /// results under unconfirmed current-result tags matching frame.
  /// Common FESA_RESULT follows docs/HHT/result_h5_structure.md.
  /// Linear Solid dynamic frames also require staged V/A and exactly matching
  /// times. An initial frame has increment/time zero and cannot complete a
  /// step. No physical initialization or reaction calculation occurs here.
  /// Finite Shell stores spatial angular motion with nonlinear geometry;
  /// dynamic Solid remains linear. No new schema or numerical conversion.
  /// Solid writes U/RF; Shell
  /// also writes UR/RM, current axes and thickness-resolved tensor payload.
  /// Shell engineering shear is converted on output; producer data is
  /// unchanged. Frame completion and whole-step completion are distinct. Both
  /// flushes must succeed. Caller alone commits memory afterwards; failure
  /// preserves all memory and latches failure. File failure markers are best
  /// effort only.
  virtual StatusCode WriteFrame(const FrameInfo* frame, const Domain* domain,
                                const AnalysisState* state, Diagnostic* d) = 0;
  /// Writes prepared candidate NT/RFL in the approved thermal schema. Borrows
  /// inputs for the call; Domain matches WriteModel and nominal state layout,
  /// with matching unconfirmed tags. No state mutation/commit. Actual adapter
  /// failures latch the writer; default unsupported preserves all inputs so
  /// existing mechanical-only implementations need no new override.
  virtual StatusCode WriteThermalFrame(const FrameInfo*, const Domain*,
                                       const ThermalState*, Diagnostic* d) {
    return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                         "Thermal frames unsupported by writer", d);
  }
  /// Checks last_completed_step against successful writes, records run state,
  /// flushes and closes explicitly. Failure (including close) is observable;
  /// completed memory/frame history is not retroactively undone. Finish(false)
  /// closes a failed run without manufacturing overall success. After a
  /// successful Finish(true), exactly one Finish(false, same last step) is
  /// permitted for a later controller/log failure: best-effort reopen, verify
  /// same file identity/schema/input/history, mark failed, flush and close.
  /// Only run_state changes; completed frames/history are preserved. Failed or
  /// unknown Finish/close does not authorize reopening. Other closed calls
  /// fail. Borrowed paths/Domain must remain alive through this optional call
  /// and its diagnostics. Reopen/verification/mark/close failure is observable;
  /// it does not guarantee recovery or atomicity against concurrent in-place
  /// edits.
  virtual StatusCode Finish(bool success, int last_completed_step,
                            Diagnostic* d) = 0;

 protected:
  ResultWriter() = default;
};
/// Creates/truncates a file after UTF-8/path/environment preflight. input_path
/// must be absolute. output naming/location is the caller's responsibility.
/// Caller owns *out; failure nulls it and releases partial resources. Both
/// paths are borrowed unchanged for writer lifetime AND subsequent diagnostic
/// reporting, including factory failure; caller must keep their storage alive
/// accordingly. No HDF5 API/type is exposed here. d is optional and borrowed
/// for each call.
StatusCode OpenResultWriter(const char* output_path, const char* input_path,
                            ResultWriter** out, Diagnostic* d);
}  // namespace fesa
#endif  // FESA_IO_RESULT_WRITER_H_
