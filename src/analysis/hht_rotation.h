/// @file
/// @brief Call-local R1/R2 implementation shared by HHT rotation entry points.
#ifndef FESA_ANALYSIS_HHT_ROTATION_H_
#define FESA_ANALYSIS_HHT_ROTATION_H_
#include "fesa/analysis/hht_integrator.h"
namespace fesa {
/// Selects predicted pose, caller-corrected pose update, or read-only R2
/// blocks.
enum class HhtRotationOperation { kPredict, kUpdate, kLinearize };
/// Borrows immutable parameters/state and caller-owned candidate buffers.
/// kPredict writes pose/omega/eta, kUpdate preserves pose and writes omega/eta,
/// kLinearize preserves all candidate data and writes required out n-by-9.
/// Other modes require null out. Does not set candidate.valid; public wrapper
/// owns validity. Required branch is cleared, chart failure sets true.
/// Validation preserves data; later failure invalidates requested outputs.
/// No retained scratch, body-spin change, convergence, recovery or commit.
/// @see docs/HHT/rotational_integrator.md R1/R2, units and chart guard.
StatusCode EvaluateHhtRotations(const HhtParameters* parameters,
                                const RotationalStateView* committed,
                                const RotationalCandidate* candidate, double dt,
                                HhtRotationOperation operation,
                                RotationalLinearization* out, bool* branch,
                                Diagnostic* d);
}  // namespace fesa
#endif  // FESA_ANALYSIS_HHT_ROTATION_H_
