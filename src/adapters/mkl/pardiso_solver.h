/// @file
/// @brief oneMKL typed PARDISO adapter; vendor types stay private.
#ifndef FESA_ADAPTERS_MKL_PARDISO_SOLVER_H_
#define FESA_ADAPTERS_MKL_PARDISO_SOLVER_H_
#include "fesa/analysis/linear_solver.h"
namespace fesa {
/// Internal per-thread Solve entry count, including failed solves.
unsigned long long PardisoSolveCallCount();
/// Owns opaque factorization storage, released after each solve and on
/// deletion.
class PardisoSolver final : public LinearSolver {
 public:
  PardisoSolver() = default;
  ~PardisoSolver() override;
  PardisoSolver(const PardisoSolver&) = delete;
  PardisoSolver& operator=(const PardisoSolver&) = delete;
  PardisoSolver(PardisoSolver&&) = delete;
  PardisoSolver& operator=(PardisoSolver&&) = delete;
  /// Upper SPD mtype=2, upper symmetric indefinite -2 or full general11.
  /// Phases11/22/33/-1. Symmetric indefinite enables scaling and weighted
  /// matching before phase11 (iparm[10]=iparm[12]=1), per user P15.
  /// No fallback. Indefinite/general with iparm[13]>0 (perturbed pivots)
  /// fails kSingular; its solution is never accepted. See
  /// docs/Integration/hht_initial_mass_pardiso_handoff.md and
  /// docs/MITC4/sparse_tangent_contract.md for source and numerical scope.
  /// Coordinator only: rejects active outer regions and serialized FESA
  /// dispatch. Applies the configured MKL limit for solve/release and restores
  /// the prior local setting on all returns. Failure invalidates free_u.
  /// Resolves the coordinator affinity ceiling before backend allocation;
  /// iparm[33]=ceiling, [59]=0 enable in-core CNR before phase 11. Requested N
  /// remains independent. See docs/Linear static/parallel_reproducibility.md
  /// for stable-affinity, sequential-analysis and RNG interference conditions.
  StatusCode Solve(const EquationSystem* system, Vector* free_u,
                   Diagnostic* d) override;

 private:
  /// Releases the current backend context, including after phase failure.
  StatusCode Release(Diagnostic* d);
  void* context_ = nullptr;  ///< Owned internal state; no vendor public types.
};
}  // namespace fesa
#endif  // FESA_ADAPTERS_MKL_PARDISO_SOLVER_H_
