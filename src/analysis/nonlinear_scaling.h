/// @file
/// @brief Simultaneous raw force/moment decision-scale conversion for sixDOF.
#ifndef FESA_ANALYSIS_NONLINEAR_SCALING_H_
#define FESA_ANALYSIS_NONLINEAR_SCALING_H_
#include "fesa/analysis/residual_report.h"
#include "fesa/io/control_settings.h"
namespace fesa {
/// SixDOF model with positive finite initial Lc. Computes force candidate
/// Cf_force*moment.raw/Lc and moment candidate Cf_moment*force.raw*Lc using
/// simultaneous RAW means. Applies only if own spatial<=raw and raw<candidate;
/// Explicit QU or Cf=0 disables the target conversion. Used scales start at
/// raw, flags false. An automatic other field can still use this raw QU.
/// No history mutation.
/// Borrows distinct fields and settings; raw means positive finite, spatial
/// nonnegative finite, Cf nonnegative finite. Owned field settings distinct.
/// Invalid input/any arithmetic failure preserves ALL field data. Success
/// changes only used_mean/converted. Skip irrelevant candidates;
/// exponent-normalized scalar arithmetic avoids intermediate overflow/underflow
/// of a representable decision scale. A used unrepresentable scale is fatal
/// kLibraryFailure, no clamp. Rounded equality with raw leaves converted=false.
/// d optional.
/// @see docs/MITC4/nonlinear_scaling_contract.md.
StatusCode ConvertNonlinearFieldMeans(const ControlSettings* controls,
                                      double characteristic_length,
                                      FieldResidual* force,
                                      FieldResidual* moment, Diagnostic* d);
}  // namespace fesa
#endif  // FESA_ANALYSIS_NONLINEAR_SCALING_H_
