/// @file
/// @brief Sequential static-step execution owner.
#ifndef FESA_SOLVER_FESA_SOLVER_H_
#define FESA_SOLVER_FESA_SOLVER_H_
#include <string>

#include "fesa/analysis/analysis.h"
#include "fesa/analysis/analysis_state.h"
#include "fesa/analysis/dynamic_step_settings.h"
#include "fesa/io/analysis_log.h"
#include "fesa/io/input_deck.h"
#include "fesa/io/result_writer.h"
#include "fesa/io/run_options.h"
#include "fesa/model/domain.h"
namespace fesa {
class ThermalConditions;
class ThermalModel;
class ThermalState;
class HeatTransferAnalysis;
/// Owns execution objects and diagnostic path storage until next Run/deletion.
class FESASolver : public AnalysisProgress {
 public:
  /// Caller deletes successful *out; failure nulls output.
  static StatusCode Create(FESASolver** out, Diagnostic* d);
  ~FESASolver() override;
  FESASolver(const FESASolver&) = delete;
  FESASolver& operator=(const FESASolver&) = delete;
  FESASolver(FESASolver&&) = delete;
  FESASolver& operator=(FESASolver&&) = delete;
  /// Borrows options for call. Returned diagnostic survives until next Run or
  /// destruction. Direct callers supply executable_directory only if neither
  /// INP controls nor explicit config_path select the source. Single-threaded;
  /// calls may change process-local numerical thread settings. New Run releases
  /// all prior execution objects. Input/active-cfg collisions with destinations
  /// fail before replacement (absolute lexical/existing native identity checks,
  /// not a concurrent-filesystem guarantee). Deck is parsed once before log
  /// activation; valid INP controls bypass all cfg path inspection. Malformed
  /// input gets an error log only if conservative source guards establish
  /// safety. Errors print to stderr and are logged once safe log activation
  /// succeeds. Failure never commits an incomplete frame; postcommit I/O
  /// failure retains completed state. d is optional; no caller-owned source
  /// backing is retained.
  StatusCode Run(const RunOptions* options, Diagnostic* d);
  /// Synchronous receiver used only by the active Analysis during Run. Borrows
  /// its owned report for this call, writes checked logs and promotes completed
  /// step conditions before postcommit logging. Never commits numerical state.
  /// Failure may leave partial log bytes; already committed state/conditions
  /// remain promoted. No report/field pointer is retained from the callback.
  StatusCode OnProgress(AnalysisEvent event, const AnalysisReport* report,
                        Diagnostic* d) override;
  /// Forwards synchronous borrowed nonlinear metadata to the owned checked
  /// log. Numerical primary stays separate from logging output d. Never
  /// retains fields or changes physical state/conditions. Run-only callback.
  StatusCode OnIncrementProgress(IncrementEvent event,
                                 const IncrementAttemptInfo* attempt,
                                 const NewtonIterationInfo* iteration,
                                 const FieldResidual* force,
                                 const FieldResidual* moment,
                                 const Diagnostic* primary,
                                 Diagnostic* d) override;

  /// Synchronous thermal-only detailed log; borrows metadata for the call,
  /// retains none and never commits/changes physical state.
  StatusCode OnThermalIncrementProgress(IncrementEvent event,
                                        const IncrementAttemptInfo* attempt,
                                        const ThermalIterationInfo* iteration,
                                        const HeatFluxResidual* heat,
                                        const Diagnostic* primary,
                                        Diagnostic* d) override;

 private:
  friend class SolverAccess;
  FESASolver() = default;
  /// Deletes borrowers before model/path owners; ends diagnostic lifetime.
  void Reset();
  /// Copies request and validates definitions before HDF5 truncate.
  StatusCode Prepare(const RunOptions* options, Diagnostic* d);
  /// Numerical/output gate; only successful WriteFrame permits Commit.
  StatusCode Execute(Diagnostic* d);
  /// Logs completed/partial calculation data before output may begin.
  StatusCode LogCalculation(const AnalysisReport* report, Diagnostic* d);
  /// Discards trial access, closes output and reports without replacing d.
  StatusCode Fail(StatusCode status, Diagnostic* d);
  /// Thermal preparation/execute helpers preserve the common output lifetime.
  StatusCode PrepareThermal(Diagnostic* d);
  StatusCode ExecuteThermal(Diagnostic* d);
  StatusCode LogThermalCalculation(const AnalysisReport* report, Diagnostic* d);
  /// Finalizes output after a fully committed step sequence.
  StatusCode FinishOutput(Diagnostic* d);
  bool heat_transfer_ =
      false;  ///< Thermal run dispatch, independent of geometry.
  ThermalConditions* thermal_conditions_ =
      nullptr;  ///< Owned current snapshot.
  ThermalConditions* trial_thermal_conditions_ = nullptr;  ///< Owned trial.
  ThermalModel* thermal_model_ = nullptr;  ///< Owns the sole thermal DOF map.
  ThermalState* thermal_state_ =
      nullptr;  ///< Owned transactional temperatures.
  HeatTransferAnalysis* thermal_analysis_ = nullptr;  ///< Owned step procedure.
  bool writer_finalized_ = false;  ///< Explicit Finish has been attempted.
  bool writer_succeeded_ =
      false;                    ///< Successful close allows late failed marker.
  bool log_finalized_ = false;  ///< Explicit log close has been attempted.
  RunOptions* request_ = nullptr;        ///< Owned copied request/context.
  std::string* input_ = nullptr;         ///< Owned stable absolute path.
  std::string* output_ = nullptr;        ///< Owned stable result path.
  InputDeck* deck_ = nullptr;            ///< Owned diagnostic source.
  ControlSettings* controls_ = nullptr;  ///< Owned.
  EffectiveConditions* conditions_ =
      nullptr;  ///< Owned last success; borrows Domain.
  EffectiveConditions* trial_conditions_ =
      nullptr;                      ///< Owned trial; borrows Domain.
  DofManager* dofs_ = nullptr;      ///< Owned.
  AnalysisModel* model_ = nullptr;  ///< Owned borrowed view.
  Analysis* analysis_ = nullptr;    ///< Owned step procedure.
  DynamicStepSettings* dynamic_settings_ =
      nullptr;                          ///< Owned; outlives Analysis.
  Vector* initial_velocity_ = nullptr;  ///< Owned; outlives direct Analysis.
  AnalysisReport* analysis_report_ =
      nullptr;                      ///< Owned; fields borrow analysis.
  ResultWriter* writer_ = nullptr;  ///< Owned; borrows paths/Domain.
  AnalysisLog* log_ = nullptr;      ///< Owned, retained for diagnostic backing.
  AnalysisState* state_ = nullptr;  ///< Owned.
  Domain* domain_ = nullptr;        ///< Owned.
  int last_completed_step_ = 0;     ///< Whole-step completion, not partial inc.
  std::int64_t next_frame_number_ = 1;    ///< Next run-wide frame identity.
  double completed_step_total_time_ = 0;  ///< Physical time at last whole step.
  double calculation_started_ = 0;  ///< Monotonic wall seconds, includes logs.
  double output_started_ = 0;  ///< Monotonic wall seconds at the output gate.
  double output_seconds_ = 0;  ///< Accumulated frame/commit/finalization time.
  /// Last logged cumulative step scopes; log lines expose interval differences
  /// so existing measurement aggregation counts each completed scope once.
  double logged_assembly_ = 0, logged_solve_ = 0;
  double logged_recovery_ = 0, logged_dispatch_ = 0;
};
}  // namespace fesa
#endif  // FESA_SOLVER_FESA_SOLVER_H_
