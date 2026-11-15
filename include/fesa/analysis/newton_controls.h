/// @file
/// @brief Public nonlinear iteration control counts, separate from increments.
#ifndef FESA_ANALYSIS_NEWTON_CONTROLS_H_
#define FESA_ANALYSIS_NEWTON_CONTROLS_H_
namespace fesa {
/// Abaqus default/public counts plus the explicit FESA execution boundaries.
/// All counts positive; divergence_start>=3. Borrowed immutable during a solve.
/// No mass/time-integration or contact controls are implied.
/// @see docs/MITC4/abaqus_increment_control.md (F4 stop/counter rules).
/// @see docs/MITC4/newton_execution_contract.md.
struct NewtonControls {
  int divergence_start = 4;  ///< I0, first divergence check.
  int slow_start = 8;        ///< IR, first logarithmic convergence prediction.
  int nonquadratic_start = 9;  ///< IP, use Rp from this residual index.
  int prediction_limit = 16;   ///< IC, prediction threshold, not hard cap.
  int iteration_limit = 50;    ///< ISc, actual correction protection cap.
};
}  // namespace fesa
#endif  // FESA_ANALYSIS_NEWTON_CONTROLS_H_
