/// @file
/// @brief Linear-static static workspace.
#ifndef FESA_ANALYSIS_STATIC_WORKSPACE_H_
#define FESA_ANALYSIS_STATIC_WORKSPACE_H_

#include "fesa/elements/element_workspace.h"
#include "fesa/math/vector.h"
namespace fesa {
class DofManager;
/// Creates count>0 independent 12-DOF element workspaces before dispatch.
/// On success caller owns the array AND every entry: delete each, then
/// delete[] the array. Failure cleans all partial storage and nulls *out.
StatusCode PrepareWorkerWorkspaces(int count, int element_dofs,
                                   ElementWorkspace*** out, Diagnostic* d);
/// Owns iteration scratch only; no committed state or DOF mappings.
/// Global vectors follow the immutable physical mapping in original-ID order,
/// [U] for three components and [U,UR] for six components, with
/// external/internal forces [F,M] and displacement increments [length,radian].
/// S6 owns sparse stiffness assembly; no dense global stiffness is allocated.
class StaticWorkspace {
 public:
  /// node_count>0, 0<=free_count<=6*N. Reduced vectors are null for zero free
  /// DOFs: no zero-dimensional math factory call. Caller owns *out on success;
  /// failure nulls it and releases partial allocations. d is optional.
  static StatusCode Create(int node_count, int free_count,
                           StaticWorkspace** out, Diagnostic* d);
  /// Copies physical layout from the borrowed step mapping; does not retain it.
  /// Caller owns success; failure nulls output and cleans partial storage.
  static StatusCode Create(const DofManager* dofs, StaticWorkspace** out,
                           Diagnostic* d);
  ~StaticWorkspace();
  StaticWorkspace(const StaticWorkspace&) = delete;
  StaticWorkspace& operator=(const StaticWorkspace&) = delete;
  StaticWorkspace(StaticWorkspace&&) = delete;
  StaticWorkspace& operator=(StaticWorkspace&&) = delete;
  /// Borrowed writable global external force scratch, force/moment units.
  Vector* ExternalForce();
  /// Borrowed writable global internal force scratch, force/moment units.
  Vector* InternalForce();
  /// Borrowed writable global displacement correction/increment scratch.
  Vector* Increment();
  /// Borrowed reduced force RHS, null if free_count==0.
  Vector* ReducedRhs();
  /// Borrowed reduced total displacement, null if free_count==0.
  Vector* ReducedSolution();
  /// Releases reduced RHS/solution after their last consumers. Idempotent and
  /// allocation/backend-free. Global residual buffers stay valid; another
  /// assembly requires a new workspace. Accessors then return null.
  void ReleaseReducedBuffers() noexcept;

 private:
  StaticWorkspace() = default;
  /// Allocates global/reduced stores; zero free count skips reduced arrays.
  static StatusCode CreateLayout(int global_count, int free_count,
                                 StaticWorkspace** out, Diagnostic* d);
  Vector* external_force_ = nullptr;    ///< Owned.
  Vector* internal_force_ = nullptr;    ///< Owned.
  Vector* increment_ = nullptr;         ///< Owned.
  Vector* reduced_rhs_ = nullptr;       ///< Owned or null.
  Vector* reduced_solution_ = nullptr;  ///< Owned or null.
};
}  // namespace fesa

#endif  // FESA_ANALYSIS_STATIC_WORKSPACE_H_
