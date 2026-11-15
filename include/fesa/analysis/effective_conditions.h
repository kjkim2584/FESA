/// @file
/// @brief Effective step conditions and optional nonlinear physical transitions.
#ifndef FESA_ANALYSIS_EFFECTIVE_CONDITIONS_H_
#define FESA_ANALYSIS_EFFECTIVE_CONDITIONS_H_

#include <vector>

#include "fesa/loads/beam_line_load.h"
#include "fesa/model/domain.h"
namespace fesa {
/// A load origin (node_id=0 before expansion) or one expanded target.
/// All definition pointers borrow immutable Domain storage, never a snapshot.
struct EffectiveLoad {
  /// Current external-flux inclusion; false for current zeros/removals.
  bool contributes_flux = false;
  bool active = false;  ///< Effective history key; false for removal/new zero.
  int node_id = 0;      ///< Original target node, never an equation index.
  double value = 0;     ///< Current endpoint force/moment, zero on removal.
  double change = 0;    ///< Signed per-origin current minus previous value.
  /// Current input definition; on removal, the removed definition (old value).
  /// Consumers must use value/flags, not this source's Value(), for removal.
  const LoadDefinition* definition = nullptr;
  const LoadDefinition* origin_definition = nullptr;  ///< First effective row.
  /// Immediately preceding active source; null for a newly introduced origin.
  const LoadDefinition* previous_definition = nullptr;
};
/// An effective prescribed component and its original definition.
struct EffectiveBoundary {
  int node_id = 0;                                 ///< Original target node.
  const BoundaryDefinition* definition = nullptr;  ///< Borrowed from Domain.
  const BoundaryDefinition* origin_definition = nullptr;    ///< First row.
  const BoundaryDefinition* previous_definition = nullptr;  ///< Prior source.
  /// Nonlinear-only prior active definition for this physical node/DOF,
  /// independent of original input identity or MOD/NEW. Null for first
  /// addition and all linear snapshots. Borrows Domain, not previous snapshot.
  const BoundaryDefinition* physical_previous_definition = nullptr;
};
/// Summed endpoint load for a node/component, separate from its origins.
struct EffectiveLoadTotal {
  int node_id = 0;   ///< Original node ID, not global/equation index.
  int dof = 0;       ///< Global component 1..6.
  double value = 0;  ///< Total force (1..3) or moment (4..6).
};
/// Immutable effective step conditions; Domain must outlive this object.
class EffectiveConditions {
 public:
  ~EffectiveConditions();
  EffectiveConditions(const EffectiveConditions&) = delete;
  EffectiveConditions& operator=(const EffectiveConditions&) = delete;
  EffectiveConditions(EffectiveConditions&&) = delete;
  EffectiveConditions& operator=(EffectiveConditions&&) = delete;
  /// Borrowed expanded entries; owns entries, never their definitions.
  const std::vector<const EffectiveLoad*>* Loads() const;
  /// Borrowed expanded constraints, unique per node/DOF.
  const std::vector<const EffectiveBoundary*>* Boundaries() const;
  /// Nonlinear-only removed physical constraints, with old definition and
  /// original node ID. Not active constraints. Null for linear snapshots;
  /// nonnull (possibly empty) for nonlinear. Owns entries/list; definitions
  /// borrow Domain. Remains valid even after previous snapshot is destroyed.
  const std::vector<const EffectiveBoundary*>* ReleasedBoundaries() const;
  /// Borrowed summed values; cancellations retain their individual origins.
  const std::vector<const EffectiveLoadTotal*>* LoadTotals() const;
  /// Borrowed source model for binding/mapping consistency checks.
  const Domain* SourceDomain() const;
  /// One-based step number.
  int StepNumber() const;
  /// Borrowed immutable effective line loads for this step, or null when the
  /// Domain has no program. Same Domain/step as this snapshot; expires here.
  const BeamLineLoadState* LineLoadState() const;

 private:
  friend StatusCode BuildEffectiveConditions(const Domain*, int,
                                             const EffectiveConditions*,
                                             EffectiveConditions**,
                                             Diagnostic*);
  EffectiveConditions() = default;
  const Domain* domain_ = nullptr;  ///< Borrowed.
  BeamLineLoadState* line_load_state_ = nullptr;  ///< Owned; borrows Domain.
  int step_number_ = 0;
  std::vector<const EffectiveLoadTotal*>* totals_ = nullptr;  ///< Owned.
  std::vector<const EffectiveLoad*>* loads_ = nullptr;  ///< Owned entries/list.
  std::vector<const EffectiveBoundary*>* boundaries_ = nullptr;  ///< Owned.
  std::vector<const EffectiveBoundary*>* released_boundaries_ =
      nullptr;  ///< Owned optional nonlinear transition entries/list.
  /// Owned unexpanded original-row records; definitions borrow Domain only.
  std::vector<const EffectiveLoad*>* load_rows_ = nullptr;
  std::vector<const EffectiveBoundary*>* boundary_rows_ = nullptr;  ///< Owned.
};
/// Resolves MOD/NEW at original scope/target/DOF before expanding targets.
/// Mechanical-only: heat-transfer steps return kUnsupported before allocation,
/// with *out null and Domain/previous unchanged. Thermal conditions resolve
/// separately. Invalid history retains its existing validation precedence.
/// Step1 requires null previous; later steps require the same Domain's
/// preceding step. Runtime previous is the last successful snapshot; upfront
/// validation may instead chain immediately preceding resolved snapshots. This
/// function does not certify calculation success. Inputs are borrowed for this
/// call. Caller owns success and may then destroy previous; Domain must outlive
/// *out. Failure nulls *out, preserves Domain/previous, and diagnostic borrows
/// Domain.
StatusCode BuildEffectiveConditions(const Domain* domain, int step_number,
                                    const EffectiveConditions* previous,
                                    EffectiveConditions** out, Diagnostic* d);
}  // namespace fesa

#endif  // FESA_ANALYSIS_EFFECTIVE_CONDITIONS_H_
