#include "analysis/field_decision.h"

#include <cmath>

#include "analysis/nonlinear_convergence.h"
#include "fesa/io/control_settings.h"
#include "model/model_storage.h"
namespace fesa {
namespace {
bool Nonnegative(double value) { return std::isfinite(value) && value >= 0; }
bool ValidControls(const FieldControls* controls) {
  return controls && Nonnegative(controls->rn) && Nonnegative(controls->cn) &&
         Nonnegative(controls->epsilon) && Nonnegative(controls->cepsilon) &&
         Nonnegative(controls->rl) && Nonnegative(controls->cf) &&
         Nonnegative(controls->epsilon_l) && Nonnegative(controls->epsilon_d) &&
         (controls->qu_explicit
              ? std::isfinite(controls->qu) && controls->qu > 0
              : std::isfinite(controls->q0) && controls->q0 > 0);
}
StatusCode NumericalFailure(Diagnostic* d) {
  return ReportFailure(StatusCode::kLibraryFailure,
                       DiagnosticStage::kComputation,
                       "Equilibrium threshold arithmetic", d);
}
}  // namespace
StatusCode EvaluateScalarFieldDecision(const FieldControls* controls,
                                       const FieldConvergenceOptions* options,
                                       bool displacement, double length,
                                       const ScalarFieldDecisionInput* input,
                                       ScalarFieldDecision* out,
                                       Diagnostic* d) {
  if (!ValidControls(controls) || !input || !out || !Nonnegative(length) ||
      input->dof_count < 0 || !Nonnegative(input->maximum) ||
      !Nonnegative(input->spatial_mean) || !Nonnegative(input->increment_max) ||
      !Nonnegative(input->correction_max))
    return model_storage::Invalid("Field criteria", d);
  if (options &&
      (!Nonnegative(options->residual_ratio) ||
       (options->estimate_available && !Nonnegative(options->estimate_max))))
    return model_storage::Invalid("Nonlinear field options", d);
  const bool actual_available = !options || options->correction_available;
  const double residual_ratio =
      options ? options->residual_ratio : controls->rn;
  out->branch = ScalarFieldBranch::kNoDofs;
  out->passed = true;
  out->residual_limit = out->correction_limit = 0;
  if (input->dof_count == 0) return StatusCode::kOk;
  if (!std::isfinite(input->used_mean) || input->used_mean <= 0 ||
      !std::isfinite(input->maximum) || !std::isfinite(input->spatial_mean))
    return model_storage::Invalid("Field criteria scale", d);
  out->phi = input->maximum / input->used_mean;
  if (!std::isfinite(out->phi))
    return ReportFailure(StatusCode::kLibraryFailure,
                         DiagnosticStage::kComputation,
                         "Equilibrium normalized residual", d);
  const double zero_limit = controls->epsilon * input->used_mean;
  if (!std::isfinite(zero_limit)) return NumericalFailure(d);
  if (input->spatial_mean <= zero_limit) {
    out->branch = ScalarFieldBranch::kZeroFlux;
    out->residual_limit = controls->epsilon;
    out->correction_limit = controls->cepsilon * input->increment_max;
    if (!std::isfinite(out->correction_limit)) return NumericalFailure(d);
    out->passed =
        out->phi <= out->residual_limit ||
        (actual_available && input->correction_max <= out->correction_limit);
  } else if (input->allow_linear && out->phi <= controls->rl) {
    out->branch = ScalarFieldBranch::kLinear;
    out->residual_limit = controls->rl;
  } else {
    out->branch = ScalarFieldBranch::kGeneral;
    out->residual_limit = residual_ratio;
    out->correction_limit = controls->cn * input->increment_max;
    if (!std::isfinite(out->correction_limit)) return NumericalFailure(d);
    const double small_limit = controls->epsilon_d * length;
    if (displacement && !std::isfinite(small_limit)) return NumericalFailure(d);
    // Before the first solve, an unchanged candidate has increment=0 but
    // supplies no evidence that the response to the new load is negligible.
    const bool skip =
        controls->cn == 1 || (actual_available && displacement && length > 0 &&
                              input->increment_max < small_limit);
    const bool estimate_passes = options && options->estimate_available &&
                                 options->estimate_max <= out->correction_limit;
    out->passed =
        out->phi <= residual_ratio &&
        (skip ||
         (actual_available && input->correction_max <= out->correction_limit) ||
         estimate_passes);
  }
  return out->passed ? StatusCode::kOk : StatusCode::kEquilibriumFailure;
}
}  // namespace fesa
