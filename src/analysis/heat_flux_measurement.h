/// @file
/// @brief Thermal origin selection and norm staging, without equation
/// decisions.
#ifndef FESA_SRC_ANALYSIS_HEAT_FLUX_MEASUREMENT_H_
#define FESA_SRC_ANALYSIS_HEAT_FLUX_MEASUREMENT_H_
#include "fesa/analysis/heat_flux_residual.h"
#include "fesa/math/vector.h"
namespace fesa {
/// Borrows signed original element f_e and external Q_eff [E/T]. Null means an
/// empty list. Preserves element zeros/multiplicity, excludes external exact0.
/// Applies formulation7 active mask using the previous accepted qmax mean;
/// gathers/sums absolute values through the common Vector library operations.
/// Truly empty input gives qbar/qmax0; a nonempty but fully filtered list
/// fails. Success stages Preview exactly once and changes only origin counts,
/// qbar, qmax and time_mean_used in out. Failure invalidates any staged history
/// even before Preview, preserves committed history; out unchanged until
/// success. No retention; controls must match the history's active step.
StatusCode MeasureThermalFluxOrigins(const Vector* element_flux,
                                     const Vector* external_flux,
                                     HeatFluxHistory* history,
                                     const FieldControls* controls,
                                     HeatFluxResidual* out, Diagnostic* d);
}  // namespace fesa
#endif  // FESA_SRC_ANALYSIS_HEAT_FLUX_MEASUREMENT_H_
