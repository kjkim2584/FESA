/// @file
/// @brief Physical DOF mapping with an explicit nonlinear Beam internal role.
#ifndef FESA_ANALYSIS_DOF_MANAGER_H_
#define FESA_ANALYSIS_DOF_MANAGER_H_

#include "fesa/analysis/effective_conditions.h"
namespace fesa {
/// Sole owner of original-ID/global/reduced mapping, rebuilt transactionally.
/// Global layout includes every node in increasing original ID with its
/// immutable three or six physical components. Reduced equations exclude
/// constraints/free nodes; their physical global entries remain.
class DofManager {
 public:
  /// Caller owns *out only on successful Build; failure nulls it.
  static StatusCode Create(const Domain* domain,
                           const EffectiveConditions* conditions,
                           DofManager** out, Diagnostic* d);
  /// Creates the explicit nonlinear Beam role: physical mappings stay intact
  /// and two unconstrained eta equations per original Beam3D follow them.
  /// Requires a nonempty exclusive Beam3D Domain. Inputs borrowed only during
  /// creation/Build; caller owns success. Failure cleans and nulls *out.
  /// @see docs/Linear static/b33_algebraic_contract.md.
  static StatusCode CreateForNonlinearBeam(
      const Domain* domain, const EffectiveConditions* conditions,
      DofManager** out, Diagnostic* d);
  ~DofManager();
  DofManager(const DofManager&) = delete;
  DofManager& operator=(const DofManager&) = delete;
  DofManager(DofManager&&) = delete;
  DofManager& operator=(DofManager&&) = delete;
  /// Borrows inputs for call; failure preserves all previous mappings.
  /// Includes connectivity, current BCs, individual nonzero load flux or delta;
  /// typed solid face interpolation must be conforming
  /// (connectivity_contract.md). current zero totals/cancellation never erase
  /// per-origin participation.
  /// A CreateForNonlinearBeam object retains its Beam role across rebuilds;
  /// unsupported elements or any failure preserve the entire previous layout.
  StatusCode Build(const Domain* domain, const EffectiveConditions* conditions,
                   Diagnostic* d);
  /// Total nodes, including excluded free nodes.
  int NodeCount() const;
  /// Total physical components, including excluded free nodes.
  int GlobalCount() const;
  /// Physical components for an original node; zero when the ID is missing.
  int NodeDofCount(int node_id) const;
  /// Number of reduced PHYSICAL unknowns; zero is valid. Excludes eta.
  int FreeCount() const;
  /// Physical free count plus internal eta count, used for algebraic systems.
  int EquationCount() const;
  /// Two eta coordinates per Beam in the explicit role, otherwise zero.
  int InternalDofCount() const;
  /// Original internal element ID in sorted order; -1 for invalid row.
  int InternalElementId(int row) const;
  /// Internal offset for original element/end=0,1; -1 for invalid identity.
  int InternalDof(int element_id, int end) const;
  /// Algebraic equation for eta; -1 for invalid identity. No physical global
  /// index/node is invented for these coordinates.
  int InternalEquation(int element_id, int end) const;
  /// Original node ID for zero-based row; -1 for invalid row.
  int NodeId(int row) const;
  /// Global component index for original node and DOF 1..6; -1 if invalid.
  int GlobalDof(int node_id, int dof) const;
  /// Reduced equation for global index; -1 constrained, -2 excluded/invalid.
  int Equation(int global_dof) const;
  /// Global index for a physical reduced equation; -1 for internal/invalid.
  int GlobalDofForEquation(int equation) const;
  /// True for included nodes, even if all their components are constrained.
  bool Participates(int node_id) const;

 private:
  DofManager() = default;
  /// Binary-searches immutable node metadata; missing ID returns -1.
  int NodeRow(int node_id) const;
  int node_count_ = 0;
  int free_count_ = 0;
  int global_count_ = 0;
  bool includes_beam_internal_ = false;  ///< Immutable factory-selected role.
  int internal_element_count_ = 0;       ///< Each has two unconstrained eta.
  int* internal_element_ids_ = nullptr;  ///< Owned sorted original IDs.
  int* node_offsets_ = nullptr;   ///< Owned prefix offsets, node_count+1.
  int* node_ids_ = nullptr;       ///< Owned sorted original IDs.
  int* global_dofs_ = nullptr;    ///< Owned reduced-to-global mapping, or null.
  int* equations_ = nullptr;      ///< Owned global-to-reduced mapping.
  bool* participates_ = nullptr;  ///< Owned node participation flags.
};
}  // namespace fesa

#endif  // FESA_ANALYSIS_DOF_MANAGER_H_
