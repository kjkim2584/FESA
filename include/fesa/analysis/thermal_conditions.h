/// @file
/// @brief Immutable thermal step conditions with original-row load history.
#ifndef FESA_ANALYSIS_THERMAL_CONDITIONS_H_
#define FESA_ANALYSIS_THERMAL_CONDITIONS_H_

#include <vector>

#include "fesa/model/domain.h"

namespace fesa {
/// Transient effective power evaluation, distinct from endpoint CFL output.
enum class ThermalFluxEvaluation { kEndpoint, kEndpointAverage };
/// Original row's role in this step; removed rows are not carried onward.
enum class ThermalFluxHistory { kNew, kRedefined, kCarried, kRemoved };
/// Expanded temperature component. All pointers borrow immutable Domain data.
struct EffectiveTemperature {
  int node_id = 0;
  int label = 11;
  double end_value = 0;  ///< Effective prescribed temperature at step end.
  const Amplitude* amplitude = nullptr;
  bool ramp = false;
  bool start_from_committed_temperature = false;
  int definition_step = 0;
  double definition_period = 0;
  const ThermalConditionDefinition* origin = nullptr;
};
/// One original power row, optionally expanded to a node (zero before
/// expansion). Power units E/T; origin/amplitude borrow Domain, never another
/// snapshot.
struct EffectiveConcentratedFlux {
  int node_id = 0;
  int label = 11;
  double end_value = 0;  ///< Effective row power at this step's end.
  double previous_end_value = 0;
  const Amplitude* amplitude = nullptr;
  bool ramp = false;
  int definition_step = 0;
  double definition_period = 0;
  ThermalFluxEvaluation evaluation = ThermalFluxEvaluation::kEndpointAverage;
  ThermalFluxHistory history = ThermalFluxHistory::kNew;
  const ThermalConditionDefinition* origin = nullptr;
};
/// One ended film/radiation definition. Domain origin/table are borrowed;
/// values are frozen at the previous step end. No prior-snapshot pointer is
/// retained.
struct PreviousThermalExchange {
  const ThermalConditionDefinition* origin = nullptr;
  double coefficient = 0;  ///< Effective h for film; epsilon for radiation.
  double reference_temperature = 0;  ///< Effective sink/ambient [theta].
};
/// One distributed original row expanded to one element. Multiple faces of
/// this same row/element are summed before norm masking; face expansion is
/// borrowed from immutable Domain origin. Units E/T/L^3 or E/T/L^2.
struct EffectiveElementLoad {
  int element_id = 0;
  int previous_exchange_begin = 0;
  int previous_exchange_count = 0;
  double end_value = 0;
  double previous_end_value = 0;
  const Amplitude* amplitude = nullptr;
  bool ramp = false;
  int definition_step = 0;
  double definition_period = 0;
  ThermalFluxHistory history = ThermalFluxHistory::kNew;
  const ThermalConditionDefinition* origin = nullptr;
};
/// Owns resolved entries/rows. Domain must outlive this and all borrowed views.
/// Pure time evaluation never changes Domain, previous snapshots or history.
class ThermalConditions final {
 public:
  ~ThermalConditions();
  ThermalConditions(const ThermalConditions&) = delete;
  ThermalConditions& operator=(const ThermalConditions&) = delete;
  ThermalConditions(ThermalConditions&&) = delete;
  ThermalConditions& operator=(ThermalConditions&&) = delete;
  /// Borrowed expanded entries, unique temperature node/label components.
  const std::vector<const EffectiveTemperature*>* Temperatures() const;
  /// Borrowed expanded original power rows, including zero/removal histories.
  /// Assembly sums rows; norm counts each nonzero effective row separately.
  const std::vector<const EffectiveConcentratedFlux*>* Fluxes() const;
  /// Borrowed original-row/element contributions, including ramp removals.
  const std::vector<const EffectiveElementLoad*>* ElementLoads() const;
  /// Borrowed source model and one-based step identity.
  const Domain* SourceDomain() const;
  int StepNumber() const;
  /// Endpoint temperature at step_time [T]. committed_start is this DOF's
  /// actual step-start committed temperature, including IC on the first step.
  /// Borrows entry from this snapshot. Finite valid inputs required; failure
  /// preserves *value. No allocations, model mutation or retained arguments.
  StatusCode EvaluateTemperature(const EffectiveTemperature* entry,
                                 double step_time, double committed_start,
                                 double* value, Diagnostic* d) const;
  /// Effective and endpoint power over [start,end] in current step time [T].
  /// Requires 0<=start<end<=period and distinct caller-owned output scalars.
  /// Steady uses endpoint; transient follows formulation5.1 row semantics.
  /// Borrows entry from this snapshot; failure preserves both outputs.
  /// Pure evaluation: rejected attempts never become a new load baseline.
  StatusCode EvaluateFlux(const EffectiveConcentratedFlux* entry, double start,
                          double end, double* effective, double* endpoint,
                          Diagnostic* d) const;
  /// Pure distributed magnitude evaluation over this attempt. Transient
  /// always averages endpoints, including MOD/NEW redefinition; steady uses
  /// endpoint. Same output preservation/lifetime contract as EvaluateFlux.
  StatusCode EvaluateElementFlux(const EffectiveElementLoad* entry,
                                 double start, double end, double* effective,
                                 double* endpoint, Diagnostic* d) const;

  /// Populates film/radiation scalar/table/time fields, retaining no out
  /// pointer. Borrows entry from this snapshot; finite time in[0,period].
  /// Failure preserves out. Does not change face/local theta_s/p pointers;
  /// caller supplies those.
  StatusCode PopulateExchangeContext(const EffectiveElementLoad* entry,
                                     double time, ThermalLoadContext* out,
                                     Diagnostic* d) const;
  /// Borrowed ended definition, null for an invalid index; Domain outlives it.
  const PreviousThermalExchange* PreviousExchangeAt(int index) const;

 private:
  friend StatusCode BuildThermalConditions(const Domain*, int,
                                           const ThermalConditions*,
                                           ThermalConditions**, Diagnostic*);
  ThermalConditions() = default;
  const Domain* domain_ = nullptr;  ///< Borrowed, immutable while bound.
  int step_number_ = 0;
  double period_ = 0;
  double total_start_ = 0;
  bool steady_ = false;
  /// Owned lists and entries; origin pointers only borrow Domain.
  std::vector<const EffectiveTemperature*>* temperatures_ = nullptr;
  std::vector<const EffectiveConcentratedFlux*>* fluxes_ = nullptr;
  /// Owned unexpanded original rows, independent of any previous snapshot.
  std::vector<const EffectiveConcentratedFlux*>* rows_ = nullptr;
  std::vector<const EffectiveElementLoad*>* element_rows_ =
      nullptr;  ///< Owned.
  std::vector<const PreviousThermalExchange*>* previous_exchanges_ =
      nullptr;  ///< Owned,lazy.
  std::vector<const EffectiveElementLoad*>* element_loads_ =
      nullptr;  ///< Owned.
};
/// Resolves one heat-transfer step (model temperature definitions at step0).
/// First step requires null previous; later steps require this Domain's
/// immediately preceding snapshot. Owns rows/expansions, borrows Domain only;
/// previous may be destroyed after success. Caller owns *out on success.
/// Failure nulls *out, releases candidates and preserves Domain/previous.
/// BC/CFLUX do not create participation; invalid direct orphan/empty filtered
/// set/label fails before expansion can drop an explicitly requested component.
StatusCode BuildThermalConditions(const Domain* domain, int step,
                                  const ThermalConditions* previous,
                                  ThermalConditions** out, Diagnostic* d);
}  // namespace fesa
#endif  // FESA_ANALYSIS_THERMAL_CONDITIONS_H_
