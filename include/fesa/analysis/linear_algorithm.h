/// @file
/// @brief Common direct linear solve ordering and last-consumer cleanup.
#ifndef FESA_ANALYSIS_LINEAR_ALGORITHM_H_
#define FESA_ANALYSIS_LINEAR_ALGORITHM_H_
#include "fesa/analysis/linear_problem.h"
#include "fesa/analysis/linear_solver.h"
namespace fesa {
/// Owned by an Analysis. Holds timing only, no numerical buffers or borrowers.
class LinearAlgorithm final {
 public:
  /// Caller deletes success; failure nulls output. Output must be empty.
  static StatusCode Create(LinearAlgorithm** out, Diagnostic* d);
  ~LinearAlgorithm() = default;
  LinearAlgorithm(const LinearAlgorithm&) = delete;
  LinearAlgorithm& operator=(const LinearAlgorithm&) = delete;
  LinearAlgorithm(LinearAlgorithm&&) = delete;
  LinearAlgorithm& operator=(LinearAlgorithm&&) = delete;
  /// Borrows both roles for this call only. Assemble -> solve -> apply ->
  /// residual; constrained problems skip solve/apply. Releases equations after
  /// solve, solution after apply, remaining calculation resources on all exits.
  /// Converts bad_alloc to allocation failure; other exceptions unwind cleanup.
  /// Preserves failure diagnostics; no output/commit or retry. Null problem
  /// is rejected, a nonnull problem is cleaned even when solver is null.
  StatusCode Solve(LinearProblem* problem, LinearSolver* solver, Diagnostic* d);
  /// Most recent Solve's inclusive solver wall time in seconds; zero if skipped.
  double SolveSeconds() const;
 private:
  LinearAlgorithm() = default;
  double solve_seconds_ = 0;
};
}  // namespace fesa
#endif  // FESA_ANALYSIS_LINEAR_ALGORITHM_H_
