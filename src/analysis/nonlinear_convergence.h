/// @file
/// @brief Common field convergence and iteration-stop scalars, no state commit.
#ifndef FESA_ANALYSIS_NONLINEAR_CONVERGENCE_H_
#define FESA_ANALYSIS_NONLINEAR_CONVERGENCE_H_
#include "fesa/analysis/newton_controls.h"
#include "fesa/analysis/residual_report.h"
#include "fesa/io/control_settings.h"
namespace fesa {
/// Optional nonlinear branch data. Null in the unchanged legacy linear path.
struct FieldConvergenceOptions {
  bool correction_available = false;  ///< Actual correction measurement valid.
  bool estimate_available =
      false;                     ///< Optional general-branch exemption valid.
  double estimate_max = 0;       ///< Same length/radian units as correction.
  double residual_ratio = .005;  ///< Rn or Rp, nonnegative dimensionless.
};
/// Stop decision for one valid residual evaluation, evaluated after
/// convergence.
enum class NewtonDecision {
  kConverged,
  kContinue,
  kDivergence,
  kSlowConvergence,
  kIterationLimit
};
/// Earlier residual maxima, field0 force,1 moment. Caller owns scalar history.
/// maximum[field][lag]: lag0 previous,1 previous-previous; field flux units.
struct NewtonResidualHistory {
  int count = 0;  ///< Available preceding evaluations,0..2.
  double maximum[2][2] = {};
};
/// Caller-owned synchronous iteration result; no borrowed numerical buffers.
struct NewtonConvergenceResult {
  NewtonDecision decision = NewtonDecision::kContinue;
  bool estimate_available[2] = {};
  double estimate_max[2] = {};
  bool prediction_available[2] = {};
  double predicted_total[2] =
      {};  ///< Infinity means nonconvergent residual rate.
};
/// Uses measured/converted positive used_mean, physical
/// max/increment/correction from distinct force/moment fields. i>=0 counts
/// completed corrections, and correction_available must equal(i>0). Earlier
/// histories must cover min(i,2).
/// First-correction estimate uses the available initial maximum; later
/// estimates use the minimum of both preceding maxima. i0 has no estimate.
/// Valid decision (including nonconvergence) returns kOk. Mutates only branch,
/// phi, limits/pass and result. Preflight
/// preserves outputs; later failure invalidates decisions. No
/// state/DOF/model/library resource retained.
StatusCode EvaluateNewtonConvergence(
    const ControlSettings* controls, const NewtonControls* iteration_controls,
    double characteristic_length, int i, bool correction_available,
    const NewtonResidualHistory* history, FieldResidual* force,
    FieldResidual* moment, NewtonConvergenceResult* result, Diagnostic* d);
}  // namespace fesa
#endif  // FESA_ANALYSIS_NONLINEAR_CONVERGENCE_H_
