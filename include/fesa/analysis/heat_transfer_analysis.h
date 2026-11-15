/// @file
/// @brief Typed uncoupled thermal step orchestration and output/commit gate.
#ifndef FESA_ANALYSIS_HEAT_TRANSFER_ANALYSIS_H_
#define FESA_ANALYSIS_HEAT_TRANSFER_ANALYSIS_H_
#include <cstdint>

#include "fesa/analysis/analysis_progress.h"
#include "fesa/analysis/thermal_model.h"
#include "fesa/analysis/thermal_newton.h"
#include "fesa/analysis/thermal_state.h"
#include "fesa/io/control_settings.h"
#include "fesa/io/result_writer.h"
namespace fesa {
/// Owns immutable step settings and scalar iteration reports. No model/map or
/// condition pointer survives SolveStep. Only successful output is committed.
class HeatTransferAnalysis final {
 public:
  /// Clones step; previous is optional borrowed validation context, not
  /// retained. Caller owns success. Failure nulls *out and releases partial
  /// resources.
  static StatusCode Create(const HeatTransferStepParameters* step,
                           const ThermalConditions* previous,
                           HeatTransferAnalysis** out, Diagnostic* d);
  ~HeatTransferAnalysis();
  HeatTransferAnalysis(const HeatTransferAnalysis&) = delete;
  HeatTransferAnalysis& operator=(const HeatTransferAnalysis&) = delete;
  HeatTransferAnalysis(HeatTransferAnalysis&&) = delete;
  HeatTransferAnalysis& operator=(HeatTransferAnalysis&&) = delete;
  /// Geometry/capacity/connectivity preflight before file truncation. Borrows
  /// all inputs, changes no state; allocates indices but no equation storage.
  static StatusCode ValidateStep(const ThermalModel* model,
                                 const ThermalDofMap* dofs,
                                 const ThermalConditions* conditions,
                                 const HeatTransferStepParameters* step,
                                 Diagnostic* d);
  /// Borrows all objects for this call only. Field pointers in report borrow
  /// this analysis until its next SolveStep/destruction. Each accepted
  /// increment writes unconfirmed results, flushes, then commits without
  /// allocation or callbacks between write and commit. Failure discards
  /// candidate access and keeps the last commit, including commits preceding
  /// late log failure. Does not finalize/own writer. Numerical errors outrank
  /// secondary log errors; any log failure is fatal. progress may be null for
  /// standalone tests.
  /// Same-step LINEAR extrapolation uses two call-local temperature-increment
  /// buffers. The previous increment/dt are promoted only after output/commit;
  /// rejected trials leave them intact and each SolveStep starts with carry.
  StatusCode SolveStep(const ThermalModel* model, const ThermalDofMap* dofs,
                       const ThermalConditions* conditions,
                       const ControlSettings* controls, ThermalState* state,
                       ResultWriter* writer, std::int64_t first_frame_number,
                       double total_time_before_step, AnalysisReport* report,
                       AnalysisProgress* progress, Diagnostic* d);

 private:
  HeatTransferAnalysis() = default;
  HeatTransferStepParameters* settings_ = nullptr;  ///< Owned immutable clone.
  ThermalIterationReport* iterations_ = nullptr;    ///< Owned report storage.
  ThermalNewton* newton_ = nullptr;                 ///< Owned stateless loop.
  Diagnostic* logging_error_ = nullptr;             ///< Owned first log error.
};
}  // namespace fesa
#endif  // FESA_ANALYSIS_HEAT_TRANSFER_ANALYSIS_H_
