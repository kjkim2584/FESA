#include "fesa/analysis/linear_algorithm.h"

#include <chrono>
#include <new>

#include "model/model_storage.h"

namespace fesa {
namespace {
class Cleanup final {
 public:
  explicit Cleanup(LinearProblem* problem) : problem_(problem) {}
  ~Cleanup() {
    if (problem_) problem_->ReleaseAttemptResources();
  }
  Cleanup(const Cleanup&) = delete;
  Cleanup& operator=(const Cleanup&) = delete;
  Cleanup(Cleanup&&) = delete;
  Cleanup& operator=(Cleanup&&) = delete;

 private:
  LinearProblem* problem_;  // Borrowed for this scope only.
};
}  // namespace
StatusCode LinearAlgorithm::Create(LinearAlgorithm** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("LinearAlgorithm output", d);
  *out = new (std::nothrow) LinearAlgorithm;
  return *out ? StatusCode::kOk
              : model_storage::AllocationFailure("LinearAlgorithm", d);
}
double LinearAlgorithm::SolveSeconds() const { return solve_seconds_; }
StatusCode LinearAlgorithm::Solve(LinearProblem* problem, LinearSolver* solver,
                                   Diagnostic* d) {
  solve_seconds_ = 0;
  Cleanup cleanup(problem);
  if (!problem || !solver)
    return model_storage::Invalid("LinearAlgorithm input", d);
  try {
    const EquationSystem* system = nullptr;
    Vector* solution = nullptr;
    bool required = false;
    auto status = problem->AssembleSystem(&system, &solution, &required, d);
    if (status != StatusCode::kOk) return status;
    if ((!required && (system || solution)) ||
        (required && (!system || !solution || solution == system->Rhs() ||
                      solution->Size() != system->Rhs()->Size())))
      return model_storage::Invalid("LinearAlgorithm equation contract", d);
    if (required) {
      const auto started = std::chrono::steady_clock::now();
      try {
        status = solver->Solve(system, solution, d);
      } catch (const std::bad_alloc&) {
        solve_seconds_ = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - started).count();
        throw;
      }
      solve_seconds_ = std::chrono::duration<double>(
          std::chrono::steady_clock::now() - started).count();
      problem->ReleaseEquationSystem();
      if (status != StatusCode::kOk) return status;
      status = problem->ApplySolution(solution, d);
      problem->ReleaseSolution();
      if (status != StatusCode::kOk) return status;
    }
    return problem->EvaluateResidual(d);
  } catch (const std::bad_alloc&) {
    return model_storage::AllocationFailure("LinearAlgorithm allocation", d);
  }
}
}  // namespace fesa
