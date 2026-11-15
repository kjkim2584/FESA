/// @file
/// @brief Finite Shell/retained-internal Beam HHT through common Newton.
#ifndef FESA_ANALYSIS_NONLINEAR_DYNAMIC_ANALYSIS_H_
#define FESA_ANALYSIS_NONLINEAR_DYNAMIC_ANALYSIS_H_
#include "fesa/analysis/dynamic_analysis.h"
#include "fesa/analysis/newton_raphson.h"
namespace fesa {
class ElementWorkspace;
/// Owns per-step Newton/solver/report and B33 worker numerical scratch. Common
/// DynamicAnalysis alone owns the frame/recovery/write/joint-commit sequence.
/// @see docs/HHT/implementation_contract.md for units/start/failure contracts.
class NonlinearDynamicAnalysis final : public DynamicAnalysis,
                                       private NonlinearProblem,
                                       private NewtonProgress {
 public:
  /// Borrows immutable settings/previous conditions/optional direct-start
  /// global6N initial velocity until destruction. Domain/settings outlive this
  /// owner. Caller deletes success; all failures null output/clean resources.
  static StatusCode Create(const DynamicStepSettings* settings, MassMode mode,
                           const EffectiveConditions* previous,
                           const Vector* initial_velocity,
                           NonlinearDynamicAnalysis** out, Diagnostic* d);
  /// Metadata/capability/density/unchanged-command preflight only. Borrows
  /// matching bound model/mapping and previous conditions; no state/evaluation.
  static StatusCode ValidateStep(const AnalysisModel* model,
                                 const DofManager* dofs, MassMode mode,
                                 const EffectiveConditions* previous,
                                 Diagnostic* d);
  ~NonlinearDynamicAnalysis() override;
  NonlinearDynamicAnalysis(const NonlinearDynamicAnalysis&) = delete;
  NonlinearDynamicAnalysis& operator=(const NonlinearDynamicAnalysis&) = delete;
  NonlinearDynamicAnalysis(NonlinearDynamicAnalysis&&) = delete;
  NonlinearDynamicAnalysis& operator=(NonlinearDynamicAnalysis&&) = delete;

 private:
  NonlinearDynamicAnalysis() = default;
  /// Dynamic hooks change candidate only; base runner handles output/commit.
  StatusCode BeginStep(const DynamicStepContext* context,
                       Diagnostic* d) override;
  StatusCode InitializeCandidate(const DynamicStepContext* context,
                                 Diagnostic* d) override;
  StatusCode AdvanceCandidate(const DynamicStepContext* context, double size,
                              double time, double total,
                              Diagnostic* d) override;
  StatusCode RecoverCandidate(const DynamicStepContext* context,
                              Diagnostic* d) override;
  /// Idempotent cleanup; preserves state and retained report backing.
  void ReleaseStepResources() noexcept override;
  /// R1 predictor is not an equilibrium correction; no-DOF is handled by
  /// Newton.
  bool RequiresInitialCorrection() const noexcept override { return true; }
  /// Physical current R11 residual/history and H6 callbacks, borrowed only
  /// within Advance. Same transactional failure contract as NonlinearProblem.
  StatusCode EvaluateResidual(int completed, FieldResidual* force,
                              FieldResidual* moment, double* length,
                              NonlinearFailure* failure,
                              Diagnostic* d) override;
  StatusCode AssembleTangent(const EquationSystem** system, Vector** correction,
                             NonlinearFailure* failure, Diagnostic* d) override;
  StatusCode ApplyCorrection(const Vector* correction,
                             NonlinearFailure* failure, Diagnostic* d) override;
  void ReleaseEquationSystem() noexcept override;
  void ReleaseFreeCorrection() noexcept override;
  void ReleaseCorrectionAssessment() noexcept override;
  void ReleaseAttemptResources() noexcept override;
  /// Synchronous iteration metadata; logging failure never hides a prior
  /// numerical error (Newton owns that precedence).
  StatusCode OnIteration(const NewtonReport* report, Diagnostic* d) override;
  /// Copies report scalars/borrows retained Newton fields; no physical writes.
  void UpdateReport() noexcept;
  const EffectiveConditions* previous_ = nullptr;  ///< Borrowed for lifetime.
  const Vector* initial_velocity_ = nullptr;       ///< Borrowed direct input.
  const DynamicStepContext* context_ = nullptr;    ///< Borrowed inside Advance.
  const Diagnostic* numerical_diagnostic_ = nullptr;  ///< Advance call only.
  NewtonRaphson* newton_ = nullptr;                   ///< Owned iteration role.
  NewtonReport* newton_report_ = nullptr;             ///< Owned report backing.
  LinearSolver* solver_ = nullptr;           ///< Owned PARDISO adapter.
  EquationSystem* equations_ = nullptr;      ///< Owned current reduced system.
  Vector* free_correction_ = nullptr;        ///< Owned reduced update.
  Vector* correction_ = nullptr;  ///< Owned actual global update assessment.
  Vector* residual_ = nullptr;  ///< Owned current physical/HHT global residual.
  Vector* internal_residual_ = nullptr;  ///< Owned retained eta HHT residual2E.
  Vector* internal_correction_ = nullptr;  ///< Owned actual eta correction2E.
  Vector* external_ = nullptr;  ///< Owned current physical external.
  Vector* old_external_ = nullptr;    ///< Owned previous physical external.
  Vector* force_origins_ = nullptr;   ///< Owned individual current F origins.
  Vector* moment_origins_ = nullptr;  ///< Owned individual current FL origins.
  /// Owned B33 worker-exclusive scratch, released with this step. No state cache.
  ElementWorkspace** workers_ = nullptr;
  int worker_count_ = 0;  ///< Owned context array length.
  MassMode mode_ = MassMode::kConsistent;
  double size_ = 0, time_ = 0, total_ = 0, length_ = 0, element_started_ = 0;
};
}  // namespace fesa
#endif  // FESA_ANALYSIS_NONLINEAR_DYNAMIC_ANALYSIS_H_
