#include "fesa/analysis/static_workspace.h"

#include <limits>
#include <memory>
#include <new>

#include "fesa/analysis/dof_manager.h"
#include "model/model_storage.h"

namespace fesa {
StatusCode PrepareWorkerWorkspaces(int count, int element_dofs,
                                   ElementWorkspace*** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("Worker workspace output", d);
  *out = nullptr;
  if (count <= 0 || element_dofs != 12)
    return model_storage::Invalid("Worker workspace shape", d);
  auto** workers = new (std::nothrow) ElementWorkspace*[count]();
  if (!workers) return model_storage::AllocationFailure("Worker array", d);
  for (int i = 0; i < count; ++i) {
    const auto status = ElementWorkspace::Create(element_dofs, &workers[i], d);
    if (status != StatusCode::kOk) {
      for (int j = 0; j < i; ++j) delete workers[j];
      delete[] workers;
      return status;
    }
  }
  *out = workers;
  return StatusCode::kOk;
}
StaticWorkspace::~StaticWorkspace() {
  delete external_force_;
  delete internal_force_;
  delete increment_;
  delete reduced_rhs_;
  delete reduced_solution_;
}
StatusCode StaticWorkspace::Create(const DofManager* dofs,
                                   StaticWorkspace** out, Diagnostic* d) {
  if (!dofs) {
    if (out) *out = nullptr;
    return model_storage::Invalid("StaticWorkspace mapping", d);
  }
  return CreateLayout(dofs->GlobalCount(), dofs->FreeCount(), out, d);
}
StatusCode StaticWorkspace::Create(int node_count, int free_count,
                                   StaticWorkspace** out, Diagnostic* d) {
  if (node_count <= 0 || node_count > std::numeric_limits<int>::max() / 6) {
    if (out) *out = nullptr;
    return model_storage::Invalid("StaticWorkspace::Create", d);
  }
  return CreateLayout(6 * node_count, free_count, out, d);
}
StatusCode StaticWorkspace::CreateLayout(int global_count, int free_count,
                                         StaticWorkspace** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("StaticWorkspace::Create", d);
  *out = nullptr;
  if (global_count <= 0 || free_count < 0 || free_count > global_count)
    return model_storage::Invalid("StaticWorkspace::Create", d);
  std::unique_ptr<StaticWorkspace> result(new (std::nothrow) StaticWorkspace);
  if (!result)
    return model_storage::AllocationFailure("StaticWorkspace::Create", d);
  Vector** stores[] = {&result->external_force_, &result->internal_force_,
                       &result->increment_};
  for (auto* store : stores) {
    const StatusCode status = Vector::Create(global_count, store, d);
    if (status != StatusCode::kOk) return status;
  }
  if (free_count > 0) {
    StatusCode status = Vector::Create(free_count, &result->reduced_rhs_, d);
    if (status != StatusCode::kOk) return status;
    status = Vector::Create(free_count, &result->reduced_solution_, d);
    if (status != StatusCode::kOk) return status;
  }
  *out = result.release();
  return StatusCode::kOk;
}
Vector* StaticWorkspace::ExternalForce() { return external_force_; }
Vector* StaticWorkspace::InternalForce() { return internal_force_; }
Vector* StaticWorkspace::Increment() { return increment_; }
Vector* StaticWorkspace::ReducedRhs() { return reduced_rhs_; }
Vector* StaticWorkspace::ReducedSolution() { return reduced_solution_; }
void StaticWorkspace::ReleaseReducedBuffers() noexcept {
  delete reduced_rhs_;
  reduced_rhs_ = nullptr;
  delete reduced_solution_;
  reduced_solution_ = nullptr;
}
}  // namespace fesa
