/// @file
/// @brief Common HHT initialization, retries, output and joint commit.
#ifndef FESA_ANALYSIS_DYNAMIC_ANALYSIS_H_
#define FESA_ANALYSIS_DYNAMIC_ANALYSIS_H_

#include "fesa/analysis/analysis.h"
#include "fesa/analysis/dynamic_step_settings.h"
#include "fesa/analysis/increment_report.h"

namespace fesa {
/// Call-scoped borrowed inputs; no member is owned. The analysis and callbacks
/// must not retain this local context past SolveStep or its pointed-to owners.
struct DynamicStepContext {
  const AnalysisModel* model = nullptr;  ///< Borrowed active model.
  const DofManager* dofs = nullptr;      ///< Borrowed current mapping.
  const EffectiveConditions* conditions =
      nullptr;                                ///< Borrowed effective loads.
  const ControlSettings* controls = nullptr;  ///< Borrowed field controls.
  AnalysisState* state = nullptr;             ///< Borrowed mutable candidate.
  ResultWriter* writer = nullptr;             ///< Borrowed output boundary.
  AnalysisReport* report = nullptr;           ///< Borrowed output metadata.
  AnalysisProgress* progress = nullptr;       ///< Borrowed synchronous log.
  IncrementAttemptInfo* attempt = nullptr;  ///< Borrowed runner metadata.
  NewtonIterationInfo* newton = nullptr;  ///< Borrowed per-attempt output.
  double total_time_before_step = 0;          ///< Physical model time [T].
};
/// Owns the integrator; borrows immutable solver-owned settings for its
/// lifetime. Derived classes own numerical resources and implement candidate
/// operations. The common runner alone publishes/commits. ResultWriter::Finish
/// stays with the caller. DIRECT stops on failure; AGGRESSIVE retries only
/// typed preconvergence numerical failures from the same committed state.
class DynamicAnalysis : public Analysis {
 public:
  ~DynamicAnalysis() override;
  DynamicAnalysis(const DynamicAnalysis&) = delete;
  DynamicAnalysis& operator=(const DynamicAnalysis&) = delete;
  DynamicAnalysis(DynamicAnalysis&&) = delete;
  DynamicAnalysis& operator=(DynamicAnalysis&&) = delete;
  /// Executes the next sequential dynamic step, with initial increment0 then
  /// advancing frames. Numerical/recovery/write failures discard only the
  /// current candidate; previous commits survive. Logs occur outside commit.
  /// A terminal post-commit log failure preserves step_completed=true.
  /// All arguments borrowed for this call. Report fields borrow this analysis
  /// until its next SolveStep/destruction. No output Finish. Automatic retries
  /// recompute kinematics at new h; recovery/output/log failures terminate.
  StatusCode SolveStep(const AnalysisModel* model, const DofManager* dofs,
                       const EffectiveConditions* conditions,
                       const ControlSettings* controls, AnalysisState* state,
                       ResultWriter* writer, std::int64_t first_frame_number,
                       double total_time_before_step, AnalysisReport* report,
                       AnalysisProgress* progress, Diagnostic* d) final;

 protected:
  DynamicAnalysis() = default;
  /// Factory-only preparation; borrows settings until destruction. Failure
  /// releases partial allocations and leaves this uninitialized.
  StatusCode InitializeIntegration(const DynamicStepSettings* settings,
                                   Diagnostic* d);
  /// Borrowed immutable integration/control settings.
  const DynamicStepSettings* Settings() const;
  /// Borrowed owned integrator, valid until destruction.
  const TimeIntegrator* Integrator() const;
  /// Prepares resources and validates supported physics without committing.
  virtual StatusCode BeginStep(const DynamicStepContext* context,
                               Diagnostic* d) = 0;
  /// Stages initial physical state/history at time0, no output/result recovery.
  virtual StatusCode InitializeCandidate(const DynamicStepContext* context,
                                         Diagnostic* d) = 0;
  /// Stages one accepted numerical candidate/history, preserving old commit.
  /// Times are physical T units. No output, result recovery, commit or retry.
  /// A Newton implementation fills context.newton with actual scalar outcome;
  /// a direct linear implementation leaves it kNotStarted (iterations NA).
  virtual StatusCode AdvanceCandidate(const DynamicStepContext* context,
                                      double increment_size, double step_time,
                                      double total_time, Diagnostic* d) = 0;
  /// Recovers official results and marks matching step/inc unconfirmed.
  /// No write/commit; failure can leave partial unconfirmed object results.
  virtual StatusCode RecoverCandidate(const DynamicStepContext* context,
                                      Diagnostic* d) = 0;
  /// Idempotent cleanup on every exit, including failed validation. Must retain
  /// valid report fields until next SolveStep/destruction; never mutate state.
  virtual void ReleaseStepResources() noexcept = 0;

 private:
  /// Runs numerical/recovery callbacks, records time and notifies calculation.
  /// Optional numerical_failure marks Advance failure only, never recovery.
  StatusCode CalculateFrame(const DynamicStepContext* context, bool initial,
                            double increment_size, double step_time,
                            double total_time, bool* numerical_failure,
                            Diagnostic* d);
  /// Writes/flushes and commits without an intervening fallible callback.
  StatusCode PublishFrame(const DynamicStepContext* context,
                          const FrameInfo* frame, Diagnostic* d);
  /// Preserves an earlier numerical error over logging failure; owns log error.
  StatusCode Notify(const DynamicStepContext* context, AnalysisEvent event,
                    StatusCode primary, Diagnostic* d);
  /// Nonlinear detailed log only, outside commit; retains primary failure.
  StatusCode NotifyIncrement(const DynamicStepContext* context,
                             IncrementEvent event, StatusCode primary,
                             Diagnostic* d);
  const DynamicStepSettings* settings_ = nullptr;  ///< Borrowed from solver.
  TimeIntegrator* integrator_ = nullptr;           ///< Owned.
  Diagnostic* logging_error_ = nullptr;            ///< Owned report backing.
};
}  // namespace fesa
#endif  // FESA_ANALYSIS_DYNAMIC_ANALYSIS_H_
