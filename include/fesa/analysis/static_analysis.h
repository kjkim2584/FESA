/// @file
/// @brief Complete direct linear-static step with output/commit coordination.
#ifndef FESA_ANALYSIS_STATIC_ANALYSIS_H_
#define FESA_ANALYSIS_STATIC_ANALYSIS_H_
#include "fesa/analysis/analysis.h"
#include "fesa/analysis/equation_system.h"
#include "fesa/analysis/linear_algorithm.h"
#include "fesa/analysis/linear_solver.h"
#include "fesa/analysis/residual_report.h"
#include "fesa/analysis/static_workspace.h"
namespace fesa {
class SolidReferenceData;
/// Owns workspace, reduced equations, linear solver and residual report.
/// Each SolveStep rebuilds scratch; Domain is never owned or cached here.
class StaticAnalysis final : public Analysis, private LinearProblem {
 public:
  /// Caller deletes success; failure nulls output and releases partial state.
  static StatusCode Create(StaticAnalysis** out, Diagnostic* d);
  ~StaticAnalysis() override;
  StaticAnalysis(const StaticAnalysis&) = delete;
  StaticAnalysis& operator=(const StaticAnalysis&) = delete;
  StaticAnalysis(StaticAnalysis&&) = delete;
  StaticAnalysis& operator=(StaticAnalysis&&) = delete;
  /// Requires the step immediately after the last committed step. Reanalysis
  /// of the same Domain uses a fresh AnalysisState. Never reads stored results
  /// as numerical input. All-constrained/free-node cases skip reduced solve.
  /// Success means increment 1 was written, flushed and committed. All required
  /// pointers are borrowed. Caller owns the report; valid residual fields
  /// borrow this Analysis until its next SolveStep/destruction. Failure clears
  /// candidate access and scratch, preserving completed commits. Only complete
  /// equilibrium-failure fields survive numerical failure cleanup. A late
  /// progress error can return failure with report->step_completed true.
  StatusCode SolveStep(const AnalysisModel* model, const DofManager* dofs,
                       const EffectiveConditions* conditions,
                       const ControlSettings* controls, AnalysisState* state,
                       ResultWriter* writer, std::int64_t first_frame_number,
                       double total_time_before_step, AnalysisReport* report,
                       AnalysisProgress* progress, Diagnostic* d) override;

 private:
  StaticAnalysis() = default;
  /// LinearProblem callbacks borrow the current Calculate inputs only.
  StatusCode AssembleSystem(const EquationSystem** system, Vector** solution,
                            bool* solve_required, Diagnostic* d) override;
  StatusCode ApplySolution(const Vector* solution, Diagnostic* d) override;
  StatusCode EvaluateResidual(Diagnostic* d) override;
  void ReleaseEquationSystem() noexcept override;
  void ReleaseSolution() noexcept override;
  void ReleaseAttemptResources() noexcept override;
  /// Releases joined worker contexts before their borrowed immutable catalog.
  void ReleaseWorkers() noexcept;
  double assembly_seconds_ = 0;  ///< AssembleReducedSystem inclusive wall time.
  double solve_seconds_ = 0;  ///< LinearSolver::Solve inclusive wall time.
  double recovery_seconds_ = 0;  ///< Equilibrium through result recovery.
  /// Releases dependent trial workspace/report before another calculation.
  void ClearScratch();
  /// Records actual completed/partial calculation scopes, including allocation
  /// exceptions, without publishing or committing results. Updates only trial
  /// state/scratch and caller report; preserves committed state on failure.
  StatusCode CalculateWithReport(const AnalysisModel* model,
                                 const DofManager* dofs,
                                 const EffectiveConditions* conditions,
                                 const ControlSettings* controls,
                                 AnalysisState* state, AnalysisReport* report,
                                 Diagnostic* d);
  /// Computes the candidate; caller discards it on every failure.
  StatusCode Calculate(const AnalysisModel* model, const DofManager* dofs,
                       const EffectiveConditions* conditions,
                       const ControlSettings* controls, AnalysisState* state,
                       Diagnostic* d);
  ElementWorkspace** workers_ = nullptr;  ///< Owns entries and array.
  SolidReferenceData** solid_references_ =
      nullptr;                     ///< Owned step-only tables/array.
  int solid_reference_count_ = 0;  ///< All entries released after workers.
  int worker_count_ = 0;  ///< Prepared workers, released by ClearScratch.
  StaticWorkspace* workspace_ = nullptr;  ///< Owned.
  EquationSystem* equations_ = nullptr;   ///< Owned.
  LinearSolver* solver_ = nullptr;        ///< Owned.
  LinearAlgorithm* algorithm_ = nullptr;  ///< Owned; no numerical cache.
  const AnalysisModel* model_ = nullptr;  ///< Borrowed during Calculate only.
  const DofManager* dofs_ = nullptr;       ///< Same call lifetime.
  const EffectiveConditions* conditions_ = nullptr;  ///< Same call lifetime.
  const ControlSettings* controls_ = nullptr;        ///< Same call lifetime.
  AnalysisState* state_ = nullptr;  ///< Borrowed candidate, never committed here.
  ResidualReport* report_ =
      nullptr;  ///< Owned; complete fields borrowed by AnalysisReport.
};
}  // namespace fesa
#endif  // FESA_ANALYSIS_STATIC_ANALYSIS_H_
