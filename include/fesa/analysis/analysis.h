/// @file
/// @brief Vendor-independent step calculation interface.
#ifndef FESA_ANALYSIS_ANALYSIS_H_
#define FESA_ANALYSIS_ANALYSIS_H_
#include "fesa/analysis/analysis_model.h"
#include "fesa/analysis/analysis_progress.h"
#include "fesa/analysis/analysis_state.h"
#include "fesa/analysis/dof_manager.h"
#include "fesa/io/control_settings.h"
#include "fesa/io/result_writer.h"
namespace fesa {
/// Borrows solver-owned definitions/mapping/state for each call. Derived
/// procedures own scratch and coordinate output/flush before committing state.
class Analysis {
 public:
  virtual ~Analysis() = default;
  /// Executes a complete step, including recovered results, frame output and
  /// successful flush followed immediately by memory commit. Required pointers
  /// are borrowed for this call; d is optional. first_frame_number is positive
  /// and matches the writer's next frame; total_time_before_step is finite and
  /// nonnegative. The caller owns/reset-lifetimes of report; progress may only
  /// inspect reports synchronously. Scalar report data survives return, while
  /// valid field pointers expire before the next calculation or destruction.
  /// Failure discards candidate access and preserves the latest committed
  /// state, including any commit preceding a late progress failure. Partial
  /// object results are never exposed as confirmed. Numerical failure outranks
  /// secondary progress errors. This method never owns or finalizes writer.
  virtual StatusCode SolveStep(
      const AnalysisModel* model, const DofManager* dofs,
      const EffectiveConditions* conditions, const ControlSettings* controls,
      AnalysisState* state, ResultWriter* writer,
      std::int64_t first_frame_number, double total_time_before_step,
      AnalysisReport* report, AnalysisProgress* progress, Diagnostic* d) = 0;
};
}  // namespace fesa
#endif  // FESA_ANALYSIS_ANALYSIS_H_
