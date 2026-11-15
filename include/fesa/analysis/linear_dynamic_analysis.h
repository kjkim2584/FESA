/// @file
/// @brief Linear Solid/MITC4 HHT with Rayleigh and the common linear solve.
#ifndef FESA_ANALYSIS_LINEAR_DYNAMIC_ANALYSIS_H_
#define FESA_ANALYSIS_LINEAR_DYNAMIC_ANALYSIS_H_

#include "fesa/analysis/dynamic_analysis.h"
#include "fesa/analysis/linear_algorithm.h"
#include "fesa/analysis/linear_solver.h"
#include "fesa/analysis/static_workspace.h"
#include "fesa/elements/solid_reference_data.h"

namespace fesa {
/// Owns numerical work and workers; DynamicAnalysis owns output/commit flow.
/// Supports direct, static preload, dynamic continuation and load jumps.
/// All data use consistent units and global3N or6N small-rotation coordinates.
class LinearDynamicAnalysis final : public DynamicAnalysis,
                                    private LinearProblem {
 public:
  /// Borrows immutable settings, previous conditions (null for first step) and
  /// optional global initial velocity through destruction. Caller owns success;
  /// failure nulls *out and cleans partial resources. Initial velocity is only
  /// valid for direct start. Same Domain must own current/previous conditions.
  /// kDefault resolves each Solid's existing mass policy independently; a
  /// homogeneous MITC4 model uses its existing row-sum default.
  static StatusCode Create(const DynamicStepSettings* settings, MassMode mode,
                           const EffectiveConditions* previous,
                           const Vector* initial_velocity,
                           LinearDynamicAnalysis** out, Diagnostic* d);
  /// Read-only preflight of the bound model, mass capability/density and fixed
  /// condition transition. Borrows inputs for the call; changes only d.
  static StatusCode ValidateStep(const AnalysisModel* model, MassMode mode,
                                 const EffectiveConditions* previous,
                                 Diagnostic* d);
  ~LinearDynamicAnalysis() override;
  LinearDynamicAnalysis(const LinearDynamicAnalysis&) = delete;
  LinearDynamicAnalysis& operator=(const LinearDynamicAnalysis&) = delete;
  LinearDynamicAnalysis(LinearDynamicAnalysis&&) = delete;
  LinearDynamicAnalysis& operator=(LinearDynamicAnalysis&&) = delete;

 private:
  LinearDynamicAnalysis() = default;
  /// Numerical hooks preserve committed state; common runner discards failures.
  StatusCode BeginStep(const DynamicStepContext* context,
                       Diagnostic* d) override;
  StatusCode InitializeCandidate(const DynamicStepContext* context,
                                 Diagnostic* d) override;
  StatusCode AdvanceCandidate(const DynamicStepContext* context, double size,
                              double time, double total,
                              Diagnostic* d) override;
  StatusCode RecoverCandidate(const DynamicStepContext* context,
                              Diagnostic* d) override;
  void ReleaseStepResources() noexcept override;
  /// LinearProblem owns the H5 system, total solution and real H1 assessment.
  StatusCode AssembleSystem(const EquationSystem** system, Vector** solution,
                            bool* solve_required, Diagnostic* d) override;
  StatusCode ApplySolution(const Vector* solution, Diagnostic* d) override;
  StatusCode EvaluateResidual(Diagnostic* d) override;
  void ReleaseEquationSystem() noexcept override;
  void ReleaseSolution() noexcept override;
  void ReleaseAttemptResources() noexcept override;
  const EffectiveConditions* previous_ =
      nullptr;                                ///< Borrowed until destruction.
  const Vector* initial_velocity_ = nullptr;  ///< Borrowed until destruction.
  const DynamicStepContext* context_ = nullptr;  ///< Borrowed within Advance.
  MassMode mode_ = MassMode::kConsistent;
  LinearAlgorithm* algorithm_ = nullptr;  ///< Owned.
  LinearSolver* solver_ = nullptr;        ///< Owned PARDISO adapter.
  ElementWorkspace** workers_ = nullptr;  ///< Owned array/entries.
  int worker_count_ = 0;
  bool shell_ = false;  ///< Bound homogeneous MITC4 model, selected per step.
  SolidReferenceData** references_ = nullptr;  ///< Owned array/tables.
  int reference_count_ = 0;  ///< Tables outlive borrowing worker catalogs.
  StaticWorkspace* workspace_ = nullptr;     ///< Owned attempt scratch.
  EquationSystem* equations_ = nullptr;      ///< Owned until solve returns.
  Vector* predictor_ = nullptr;              ///< Owned attempt kinematics.
  Vector* velocity_ = nullptr;               ///< Owned attempt kinematics.
  Vector* acceleration_ = nullptr;           ///< Owned attempt kinematics.
  Vector* external_ = nullptr;               ///< Owned current physical load.
  Vector* external_origins_ = nullptr;       ///< Owned individual CLOADs.
  Vector* weighted_ = nullptr;               ///< Owned HHT weighted load.
  FieldResidual* force_ = nullptr;           ///< Owned report backing.
  FieldResidual* moment_ = nullptr;          ///< Owned moment field (absent3N).
  double size_ = 0, time_ = 0, total_ = 0;
};
}  // namespace fesa
#endif  // FESA_ANALYSIS_LINEAR_DYNAMIC_ANALYSIS_H_
