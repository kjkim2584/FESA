/// @file
/// @brief Original-node thermal component and reduced-equation mapping.
#ifndef FESA_ANALYSIS_THERMAL_DOF_MAP_H_
#define FESA_ANALYSIS_THERMAL_DOF_MAP_H_
#include "fesa/analysis/thermal_conditions.h"
namespace fesa {
/// Immutable owned index arrays; no Matrix/Vector or reduced numerical arrays.
/// Domain is borrowed and must outlive this map. All nodes have nominal global
/// slots; only element-connected nodes participate. Conditions add constraints.
class ThermalDofMap final {
 public:
  /// Caller owns success; failure nulls *out, cleans candidates and preserves
  /// Domain/conditions. Requires a snapshot of this Domain. Rejects mixed
  /// DC/DS, shared component-count mismatch and invalid labels before
  /// truncation.
  static StatusCode Create(const Domain* domain,
                           const ThermalConditions* conditions,
                           ThermalDofMap** out, Diagnostic* d);
  ~ThermalDofMap();
  ThermalDofMap(const ThermalDofMap&) = delete;
  ThermalDofMap& operator=(const ThermalDofMap&) = delete;
  ThermalDofMap(ThermalDofMap&&) = delete;
  ThermalDofMap& operator=(ThermalDofMap&&) = delete;
  /// Sorted original-node rows, including nonparticipants.
  int NodeCount() const;
  int NodeId(int row) const;         ///< -1 for an invalid row.
  int NodeDofCount(int node) const;  ///< Nominal count; zero for missing ID.
  int GlobalCount() const;
  int FreeCount() const;  ///< Zero is valid; no zero-sized Matrix is created.
  /// Label11..10+ns -> physical global slot; -1 for missing node/label.
  int GlobalDof(int node, int label) const;
  /// Reduced equation; -1 constrained, -2 nonparticipating or invalid index.
  int Equation(int global) const;
  int GlobalDofForEquation(int equation) const;  ///< -1 for invalid equation.
  bool Participates(int node) const;
  const Domain* SourceDomain() const;  ///< Borrowed until destruction.
  int StepNumber() const;

 private:
  ThermalDofMap() = default;
  int NodeRow(int node) const;
  const Domain* domain_ = nullptr;  ///< Borrowed.
  int step_ = 0;
  int node_count_ = 0;
  int global_count_ = 0;
  int free_count_ = 0;
  int* nodes_ = nullptr;      ///< Owned sorted IDs.
  int* offsets_ = nullptr;    ///< Owned node_count+1 prefixes.
  int* equations_ = nullptr;  ///< Owned global-to-reduced slots.
  int* globals_ = nullptr;  ///< Owned reduced-to-global, null if no free DOFs.
  bool* participates_ = nullptr;  ///< Owned node participation flags.
};
}  // namespace fesa
#endif  // FESA_ANALYSIS_THERMAL_DOF_MAP_H_
