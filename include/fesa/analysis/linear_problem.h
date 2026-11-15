/// @file
/// @brief Analysis-owned linear candidate assembly, update and assessment.
#ifndef FESA_ANALYSIS_LINEAR_PROBLEM_H_
#define FESA_ANALYSIS_LINEAR_PROBLEM_H_
#include "fesa/analysis/equation_system.h"
namespace fesa {
/// Borrowed for one solve. Owns its numerical resources and physical candidate;
/// the algorithm never owns them, publishes output, or commits physical state.
class LinearProblem {
 public:
  virtual ~LinearProblem() = default;
  /// Returns borrowed immutable equations and distinct writable solution.
  /// Positive matching sizes are required when solve_required is true.
  /// Fully constrained problems return false and two null pointers, with the
  /// prescribed candidate already prepared. Failure nulls both pointers;
  /// the owner releases partial resources. No arguments are retained.
  virtual StatusCode AssembleSystem(const EquationSystem** system,
                                    Vector** solution, bool* solve_required,
                                    Diagnostic* d) = 0;
  /// Applies the total solved unknown using the problem's physical mapping.
  /// Borrowed solution remains valid through this call only. Failure leaves
  /// candidate invalid for commit; committed state is never modified.
  virtual StatusCode ApplySolution(const Vector* solution, Diagnostic* d) = 0;
  /// Recomputes candidate residual/assessment and retains the problem's report.
  /// Does not solve, publish or commit. Distinguishes numerical failure from
  /// equilibrium failure; no residual or correction is replaced by zero.
  virtual StatusCode EvaluateResidual(Diagnostic* d) = 0;
  /// Called after solver returns; releases only owned equation resources.
  /// Idempotent, no allocation/backend call/exception/diagnostic mutation.
  virtual void ReleaseEquationSystem() noexcept = 0;
  /// Called after solution consumption; same cleanup guarantees as above.
  virtual void ReleaseSolution() noexcept = 0;
  /// Called on every exit for remaining algorithm calculation resources.
  /// Preserves physical state, reports and worker/catalog contexts still
  /// borrowed by later official recovery; their owner handles last use.
  /// Idempotent; no allocation/backend call/exception/diagnostic mutation.
  virtual void ReleaseAttemptResources() noexcept = 0;
};
}  // namespace fesa
#endif  // FESA_ANALYSIS_LINEAR_PROBLEM_H_
