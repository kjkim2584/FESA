/// @file
/// @brief Replaceable reduced linear solve contract.
#ifndef FESA_ANALYSIS_LINEAR_SOLVER_H_
#define FESA_ANALYSIS_LINEAR_SOLVER_H_
#include "fesa/analysis/equation_system.h"
namespace fesa {
/// Borrowed inputs; implementations own their factorization resources.
class LinearSolver {
 public:
  virtual ~LinearSolver() = default;
  /// Solves the caller-defined reduced system (total U/UR, correction or
  /// physical initial acceleration with a mass operator).
  /// Output is distinct from system RHS and has
  /// matching size. Validation preserves output; later failure invalidates it.
  /// No automatic matrix-type fallback, regularization or retry is permitted.
  virtual StatusCode Solve(const EquationSystem* system, Vector* free_u,
                           Diagnostic* d) = 0;
};
}  // namespace fesa
#endif  // FESA_ANALYSIS_LINEAR_SOLVER_H_
