#include "adapters/mkl/pardiso_solver.h"

#include <new>

#include "adapters/mkl/dense_ops.h"
#include "adapters/openmp/thread_runtime.h"
#include "mkl.h"
#include "model/model_storage.h"
namespace fesa {
namespace {
thread_local unsigned long long solve_calls = 0;
struct LocalThreads {
  explicit LocalThreads(int count)
      : previous(mkl_set_num_threads_local(count)) {}
  LocalThreads(const LocalThreads&) = delete;
  LocalThreads& operator=(const LocalThreads&) = delete;
  LocalThreads(LocalThreads&&) = delete;
  LocalThreads& operator=(LocalThreads&&) = delete;
  ~LocalThreads() { mkl_set_num_threads_local(previous); }
  int previous;
};
struct Context {
  void* handles[64] = {};
  int parameters[64] = {};
  int size = 0;
  int type = 2;
};
StatusCode NativeFailure(int error, Diagnostic* d) {
  const StatusCode code = error == -4 || error == -7 ? StatusCode::kSingular
                          : error == -2 || error == -9
                              ? StatusCode::kAllocationFailure
                              : StatusCode::kLibraryFailure;
  if (d) d->Record(code, DiagnosticStage::kComputation, "PARDISO", error);
  return code;
}
}  // namespace
unsigned long long PardisoSolveCallCount() { return solve_calls; }
PardisoSolver::~PardisoSolver() { Release(nullptr); }
StatusCode PardisoSolver::Release(Diagnostic* d) {
  if (!context_) return StatusCode::kOk;
  auto* state = static_cast<Context*>(context_);
  const int one = 1, phase = -1, message = 0;
  int error = 0, dummy = 0;
  double scalar = 0;
  pardiso(state->handles, &one, &one, &state->type, &phase, &state->size,
          &scalar, &dummy, &dummy, &dummy, &one, state->parameters, &message,
          &scalar, &scalar, &error);
  delete state;
  context_ = nullptr;
  return error == 0 ? StatusCode::kOk : NativeFailure(error, d);
}
StatusCode PardisoSolver::Solve(const EquationSystem* system, Vector* free_u,
                                Diagnostic* d) {
  ++solve_calls;
  if (InElementExecution() || !system || !free_u || free_u == system->Rhs() ||
      free_u->Size() != system->Stiffness()->Size())
    return model_storage::Invalid("PardisoSolver::Solve", d);
  int cnr_ceiling = 0;
  StatusCode status = PardisoCnrThreadCeiling(&cnr_ceiling, d);
  if (status != StatusCode::kOk) return status;
  LocalThreads threads(ConfiguredThreadCount());
  const SparseMatrix* matrix = system->Stiffness();
  if (system->Type() == MatrixType::kSymmetricPositiveDefinite) {
    // Positive energy requires every diagonal to be strictly positive.
    for (int row = 0; row < matrix->Size(); ++row) {
      const int start = matrix->RowOffsets()[row];
      if (start == matrix->RowOffsets()[row + 1] ||
          matrix->Columns()[start] != row || matrix->Values()[start] <= 0)
        return NativeFailure(-7, d);
    }
  }
  status = Release(d);
  if (status != StatusCode::kOk) return status;
  auto* state = new (std::nothrow) Context;
  if (!state) return model_storage::AllocationFailure("PARDISO", d);
  context_ = state;
  state->size = matrix->Size();
  switch (system->Type()) {
    case MatrixType::kSymmetricPositiveDefinite: state->type = 2; break;
    case MatrixType::kSymmetricIndefinite: state->type = -2; break;
    case MatrixType::kGeneral: state->type = 11; break;
  }
  const int one = 1, message = 0;
  pardisoinit(state->handles, &state->type, state->parameters);
  state->parameters[0] = 1;
  state->parameters[1] = 2;   // METIS nested dissection.
  state->parameters[5] = 0;   // Preserve RHS; solution goes to x.
  state->parameters[7] = 2;   // Direct-solver iterative refinement only.
  if (state->type == -2) {
    // Use weighted matching/scaling for symmetric saddle-point systems.
    // Intel oneMKL iparm[10]/[12]; docs/Integration/
    // hht_initial_mass_pardiso_handoff.md. Pivot substitution stays forbidden.
    state->parameters[10] = 1;
    state->parameters[12] = 1;
  }
  state->parameters[26] = 1;  // CSR index checker.
  state->parameters[27] = 0;  // Double precision.
  state->parameters[33] = cnr_ceiling;  // Common CNR scale, independent of N.
  state->parameters[34] = 1;            // Zero-based indices.
  state->parameters[59] = 0;            // In-core CNR; no out-of-core fallback.
  int dummy = 0, error = 0;
  const int phases[] = {11, 22, 33};
  for (const int phase : phases) {
    if (!dense_ops::BeginOperation()) {
      error = -5;
      break;
    }
    free_u->finite_verified_ = false;
    pardiso(state->handles, &one, &one, &state->type, &phase, &state->size,
            matrix->Values(), matrix->RowOffsets(), matrix->Columns(), &dummy,
            &one, state->parameters, &message,
            const_cast<double*>(system->Rhs()->Data()), free_u->data_, &error);
    if (error != 0) break;
    // Indefinite/general factorizations may replace small pivots.
    // Do not accept a different operator, even for a consistent singular RHS.
    if (phase == 22 && state->type != 2 && state->parameters[13] > 0) {
      const int perturbed = state->parameters[13];
      Release(nullptr);
      if (d)
        d->Record(StatusCode::kSingular, DiagnosticStage::kComputation,
                  "PARDISO perturbed pivots", perturbed);
      return StatusCode::kSingular;
    }
  }
  if (error != 0) {
    Release(nullptr);
    return NativeFailure(error, d);
  }
  status = Release(d);
  if (status != StatusCode::kOk) return status;
  return dense_ops::AllFinite(free_u->Size(), free_u->Data())
             ? StatusCode::kOk
             : NativeFailure(-5, d);
}
}  // namespace fesa
