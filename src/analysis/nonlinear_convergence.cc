#include "analysis/nonlinear_convergence.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "analysis/equilibrium_check.h"
#include "model/model_storage.h"
namespace fesa {
namespace {
bool Nonnegative(double value) { return std::isfinite(value) && value >= 0; }
bool ValidControls(const FieldControls* c) {
  if (!c) return false;
  const double ratios[] = {c->rn, c->rp, c->cn,        c->epsilon,  c->cepsilon,
                           c->rl, c->cf, c->epsilon_l, c->epsilon_d};
  for (double v : ratios)
    if (!Nonnegative(v)) return false;
  return c->qu_explicit ? std::isfinite(c->qu) && c->qu > 0
                        : std::isfinite(c->q0) && c->q0 > 0;
}
bool ValidField(const FieldResidual* f) {
  return f && f->dof_count >= 0 && Nonnegative(f->maximum) &&
         Nonnegative(f->spatial_mean) && Nonnegative(f->increment_max) &&
         Nonnegative(f->correction_max) &&
         (f->dof_count == 0 ||
          (std::isfinite(f->used_mean) && f->used_mean > 0));
}
// Scalar product/ratio ordering avoids overflow of an otherwise finite
// estimate. A positive result below representable range is NA, never fake0.
bool Estimate(double a, double previous, double older, double correction,
              double* out) {
  const double denominator = std::min(previous, older);
  if (!(denominator > 0)) return false;
  if (a == 0 || correction == 0) {
    *out = 0;
    return true;
  }
  int ea = 0, ec = 0, ed = 0;
  const double ma = std::frexp(a, &ea), mc = std::frexp(correction, &ec),
               md = std::frexp(denominator, &ed);
  const double value = std::scalbn(ma * mc / md, ea + ec - ed);
  if (!std::isfinite(value) || !(value > 0)) return false;
  *out = value;
  return true;
}
double Predict(int i, const FieldResidual* f, double previous) {
  const double a = f->maximum;
  if (!(a > 0) || !(previous > 0) || a >= previous || f->residual_limit == 0)
    return std::numeric_limits<double>::infinity();
  // log1p preserves a rate near1; separate logs preserve tiny a/previous.
  const double log_rate = a > previous / 2
                              ? std::log1p((a - previous) / previous)
                              : std::log(a) - std::log(previous);
  if (!(log_rate < 0)) return std::numeric_limits<double>::infinity();
  const double log_target_ratio =
      std::log(f->residual_limit) + std::log(f->used_mean) - std::log(a);
  const double remaining = std::max(0.0, log_target_ratio / log_rate);
  return static_cast<double>(i) + std::ceil(remaining);
}
}  // namespace
StatusCode EvaluateNewtonConvergence(
    const ControlSettings* controls, const NewtonControls* c, double length,
    int i, bool correction_available, const NewtonResidualHistory* history,
    FieldResidual* force, FieldResidual* moment,
    NewtonConvergenceResult* result, Diagnostic* d) {
  if (!controls || !ValidControls(controls->displacement) ||
      !ValidControls(controls->rotation) || !c || c->divergence_start < 3 ||
      c->slow_start <= 0 || c->nonquadratic_start <= 0 ||
      c->prediction_limit <= 0 || c->iteration_limit <= 0 ||
      !Nonnegative(length) || i < 0 || correction_available != (i > 0) ||
      !history || history->count != std::min(i, 2) || !ValidField(force) ||
      !ValidField(moment) || force == moment || !result)
    return model_storage::Invalid("Newton convergence inputs", d);
  for (int field = 0; field < 2; ++field)
    for (int lag = 0; lag < history->count; ++lag)
      if (!Nonnegative(history->maximum[field][lag]))
        return model_storage::Invalid("Newton residual history", d);
  FieldResidual* fields[2] = {force, moment};
  const FieldControls* settings[2] = {controls->displacement,
                                      controls->rotation};
  FieldConvergenceOptions options[2];
  result->decision = NewtonDecision::kContinue;
  for (int field = 0; field < 2; ++field) {
    result->estimate_available[field] = false;
    result->estimate_max[field] = 0;
    result->prediction_available[field] = false;
    result->predicted_total[field] = 0;
    options[field].correction_available = correction_available;
    options[field].residual_ratio =
        i >= c->nonquadratic_start ? settings[field]->rp : settings[field]->rn;
    if (fields[field]->dof_count > 0 && i >= 1) {
      // At the first actual correction only a0 exists. Native first-increment
      // probes confirm a1/a0*c1; do not read unavailable lag1 or fabricate0.
      // Later corrections retain the published minimum of both prior maxima.
      // @see incremental_hexa8_start_estimate_fix.md.
      result->estimate_available[field] =
          Estimate(fields[field]->maximum, history->maximum[field][0],
                   history->maximum[field][i >= 2 ? 1 : 0],
                   fields[field]->correction_max,
                   &result->estimate_max[field]);
      options[field].estimate_available = result->estimate_available[field];
      options[field].estimate_max = result->estimate_max[field];
    }
  }
  Diagnostic local;
  const auto status = EvaluateJointFieldCriteriaWithOptions(
      controls, length, &options[0], &options[1], force, moment, &local);
  if (status != StatusCode::kOk && status != StatusCode::kEquilibriumFailure) {
    if (d)
      d->Record(local.code, local.stage, local.operation, local.vendor_code);
    return status;
  }
  if (force->passed && moment->passed) {
    result->decision = NewtonDecision::kConverged;
    return StatusCode::kOk;
  }
  bool diverging = false, slow = false;
  for (int field = 0; field < 2; ++field) {
    const FieldResidual* f = fields[field];
    if (f->dof_count == 0) continue;
    const double a = f->maximum;
    if (i >= c->divergence_start && history->count >= 2 &&
        std::min(a, history->maximum[field][0]) > history->maximum[field][1])
      diverging = true;
    if (i >= c->slow_start && history->count >= 1 &&
        f->phi > f->residual_limit) {
      result->prediction_available[field] = true;
      result->predicted_total[field] =
          Predict(i, f, history->maximum[field][0]);
      if (result->predicted_total[field] > c->prediction_limit) slow = true;
    }
  }
  if (diverging)
    result->decision = NewtonDecision::kDivergence;
  else if (slow)
    result->decision = NewtonDecision::kSlowConvergence;
  else if (i >= c->iteration_limit)
    result->decision = NewtonDecision::kIterationLimit;
  return StatusCode::kOk;
}
}  // namespace fesa
