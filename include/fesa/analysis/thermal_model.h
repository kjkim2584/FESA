/// @file
/// @brief Transactionally bound thermal element/condition view and owned map.
#ifndef FESA_ANALYSIS_THERMAL_MODEL_H_
#define FESA_ANALYSIS_THERMAL_MODEL_H_
#include "fesa/analysis/thermal_dof_map.h"
namespace fesa {
/// Owns list storage/map, borrows Domain elements and a condition snapshot.
/// Domain and bound conditions must outlive the view or be unbound via Clear.
class ThermalModel final {
 public:
  /// Creates an empty reusable view; caller owns success, failure nulls *out.
  static StatusCode Create(ThermalModel** out, Diagnostic* d);
  ~ThermalModel();
  ThermalModel(const ThermalModel&) = delete;
  ThermalModel& operator=(const ThermalModel&) = delete;
  ThermalModel(ThermalModel&&) = delete;
  ThermalModel& operator=(ThermalModel&&) = delete;
  /// Replaces borrowed views and owned map only on success. Failure preserves
  /// the previous binding and all inputs; caller owns supplied conditions.
  StatusCode BindStep(const Domain* domain, const ThermalConditions* conditions,
                      Diagnostic* d);
  /// Releases the map and borrowed views without changing Domain/conditions;
  /// no allocation or failure. Empty element list remains valid for reuse.
  void Clear() noexcept;
  const std::vector<const ThermalElement*>* Elements() const;
  const ThermalConditions* Conditions() const;
  const ThermalDofMap* DofMap() const;
  const Domain* SourceDomain() const;

 private:
  ThermalModel() = default;
  const Domain* domain_ = nullptr;                 ///< Borrowed while bound.
  const ThermalConditions* conditions_ = nullptr;  ///< Borrowed while bound.
  ThermalDofMap* map_ = nullptr;                   ///< Owned.
  std::vector<const ThermalElement*>* elements_ = nullptr;  ///< Owned list.
};
}  // namespace fesa
#endif  // FESA_ANALYSIS_THERMAL_MODEL_H_
