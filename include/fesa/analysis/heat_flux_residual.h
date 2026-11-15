/// @file
/// @brief Single temperature-field residual in original thermal node/label
/// units.
#ifndef FESA_ANALYSIS_HEAT_FLUX_RESIDUAL_H_
#define FESA_ANALYSIS_HEAT_FLUX_RESIDUAL_H_
#include "fesa/base/diagnostic.h"
namespace fesa {
class ThermalAssembly;
class HeatFluxHistory;
struct FieldControls;
/// Thermal convergence branch; no mechanical force/moment field identity.
enum class ThermalResidualBranch { kNoDofs, kZeroFlux, kLinear, kGeneral };
/// Caller-owned scalar report; fluxes E/T, temperature increments theta.
struct HeatFluxResidual {
  int dof_count = 0;
  /// Internal/full distributed-row DOFs plus nonzero concentrated rows.
  int contribution_count = 0;
  int active_contribution_count =
      0;  ///< Denominator after active-mask selection.
  int maximum_node = 0;
  int maximum_label = 0;
  double maximum = 0;
  double spatial_mean = 0;
  double time_mean_used = .01;
  double maximum_flux = 0;
  double phi = 0;
  double increment_max = 0;
  double correction_max = 0;
  ThermalResidualBranch branch = ThermalResidualBranch::kNoDofs;
  bool passed = false;
};
/// Measures a successful assembly and stages one norm candidate. Borrows all
/// inputs, retaining none; controls must match this history's begun step.
/// Every failure discards the history candidate and preserves committed data.
/// Output may be invalid on failure. RateHistory is never changed.
StatusCode MeasureHeatFlux(const ThermalAssembly* assembly,
                           HeatFluxHistory* history,
                           const FieldControls* controls, HeatFluxResidual* out,
                           Diagnostic* d);
}  // namespace fesa
#endif  // FESA_ANALYSIS_HEAT_FLUX_RESIDUAL_H_
