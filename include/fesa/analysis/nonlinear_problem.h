/// @file
/// @brief Analysis-provided residual/tangent/update roles for one Newton trial.
#ifndef FESA_ANALYSIS_NONLINEAR_PROBLEM_H_
#define FESA_ANALYSIS_NONLINEAR_PROBLEM_H_
#include "fesa/analysis/equation_system.h"
#include "fesa/analysis/residual_report.h"
namespace fesa {
/// Only valid-reference trial geometry/chart failures permit element cutback.
/// No element family, static/dynamic formulation or vendor type is exposed.
enum class NonlinearFailure { kNone, kGeometry, kRotationBranch };
/// Analysis implements this directly, owns physical state and equation buffers.
/// Newton borrows only for the attempt; no solve/commit/output hidden here.
class NonlinearProblem {
 public:
  virtual ~NonlinearProblem() = default;
  /// Execution metadata only: active trial equations require one actual
  /// correction before accepting convergence. Default preserves static i0
  /// exceptions; no-DOF trials remain correction-free. No state mutation.
  virtual bool RequiresInitialCorrection() const noexcept { return false; }
  /// Recomputes residual and field measurements at current trial. Fills all
  /// field scalars (no-DOF fields excluded), positive used means for active
  /// fields, and finite nonnegative characteristic_length. Raw/used means are
  /// supplied from Analysis-owned increment history. correction unavailable
  /// when completed_corrections==0; caller must not claim an actual zero.
  /// Inputs/outputs borrowed; failure invalidates measured fields. Clears
  /// required failure flag unless returning a typed recoverable failure.
  virtual StatusCode EvaluateResidual(int completed_corrections,
                                      FieldResidual* force,
                                      FieldResidual* moment,
                                      double* characteristic_length,
                                      NonlinearFailure* failure,
                                      Diagnostic* d) = 0;
  /// Builds CURRENT correction system K*delta=residual, never total U/UR.
  /// Outputs borrow Analysis-owned immutable system and separate writable
  /// free-correction vector, matching positive size. Borrow lasts through the
  /// immediately following solve/update; next assembly may replace old stores.
  /// Failure nulls both outputs and owns partial cleanup; no physical commit.
  virtual StatusCode AssembleTangent(const EquationSystem** system,
                                     Vector** correction,
                                     NonlinearFailure* failure,
                                     Diagnostic* d) = 0;
  /// Applies solved free correction using Analysis's mapping/kinematics.
  /// Input borrowed/read-only. Failure invalidates candidate, never committed
  /// state; typed geometry/chart errors distinguished from resource/library.
  virtual StatusCode ApplyCorrection(const Vector* correction,
                                     NonlinearFailure* failure,
                                     Diagnostic* d) = 0;
  /// Releases equation storage after the linear solve, including failure.
  /// Cleanup callbacks are idempotent and perform no allocation/backend work.
  virtual void ReleaseEquationSystem() noexcept = 0;
  /// Releases the free correction after ApplyCorrection consumes it.
  virtual void ReleaseFreeCorrection() noexcept = 0;
  /// Releases actual correction storage after residual measurements consume it.
  virtual void ReleaseCorrectionAssessment() noexcept = 0;
  /// Releases all attempt numerical scratch on every Solve exit. Preserves
  /// physical candidate/committed state, copied reports and recovery contexts.
  virtual void ReleaseAttemptResources() noexcept = 0;
};
}  // namespace fesa
#endif  // FESA_ANALYSIS_NONLINEAR_PROBLEM_H_
