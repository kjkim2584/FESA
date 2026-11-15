#include "fesa/analysis/hht_integrator.h"

#include <cmath>
#include <new>

#include "analysis/hht_rotation.h"

namespace fesa {
namespace {
StatusCode Invalid(Diagnostic* d) {
  return ReportFailure(StatusCode::kInvalidArgument,
                       DiagnosticStage::kValidation, "HHT kinematics input", d);
}
bool ValidParameters(const HhtParameters* parameters) {
  return parameters && std::isfinite(parameters->alpha) &&
         parameters->alpha >= -1. / 3 && parameters->alpha <= 0;
}
// Call-local scalars only. No vector or factorization reuse across attempts.
struct Coefficients {
  double predictor_acceleration = 0;
  double old_acceleration_velocity = 0;
  double new_acceleration_velocity = 0;
  double c0 = 0;
  double c1 = 0;
};
StatusCode ComputeCoefficients(const HhtParameters* parameters, double dt,
                               Coefficients* out, Diagnostic* d) {
  if (!ValidParameters(parameters) || !std::isfinite(dt) || dt <= 0)
    return Invalid(d);
  const double beta = (1 - parameters->alpha) * (1 - parameters->alpha) / 4;
  const double gamma = .5 - parameters->alpha;
  const double squared_time = dt * dt;
  out->c0 = (1 / beta / dt) / dt;
  out->c1 = gamma / beta / dt;
  out->predictor_acceleration = squared_time * (.5 - beta);
  out->old_acceleration_velocity = dt * (1 - gamma);
  out->new_acceleration_velocity = dt * gamma;
  if (!std::isfinite(squared_time) || squared_time <= 0 ||
      !std::isfinite(out->c0) || out->c0 <= 0 || !std::isfinite(out->c1) ||
      out->c1 <= 0)
    return ReportFailure(StatusCode::kEquilibriumFailure,
                         DiagnosticStage::kComputation,
                         "HHT coefficients not representable", d);
  return StatusCode::kOk;
}
bool ValidViews(const TranslationalStateView* state,
                const TranslationalCandidate* candidate, bool update) {
  if (!state || !candidate) return false;
  const Vector* vectors[] = {state->displacement, state->velocity,
                             state->acceleration, candidate->displacement,
                             candidate->velocity, candidate->acceleration};
  for (int i = 0; i < 6; ++i) {
    if (!vectors[i] || vectors[i]->Size() != vectors[0]->Size()) return false;
    for (int j = 0; j < i; ++j)
      if (vectors[i] == vectors[j]) return false;
    // Validation only; vector algebra is exclusively through common Vector.
    if (i < 3 || (update && i == 3))
      for (int entry = 0; entry < vectors[i]->Size(); ++entry)
        if (!std::isfinite(vectors[i]->Data()[entry])) return false;
  }
  return true;
}
class HhtIntegrator final : public TimeIntegrator {
 public:
  explicit HhtIntegrator(const HhtParameters* parameters)
      : parameters_(parameters) {}
  StatusCode PrepareKinematics(const TranslationalStateView* state, double dt,
                               TranslationalCandidate* candidate,
                               Diagnostic* d) const override {
    if (candidate) candidate->valid = false;
    if (!ValidViews(state, candidate, false)) return Invalid(d);
    Coefficients coefficients;
    auto status = ComputeCoefficients(parameters_, dt, &coefficients, d);
    if (status != StatusCode::kOk) return status;
    auto* u = candidate->displacement;
    status = u->CopyFrom(state->displacement, d);
    if (status == StatusCode::kOk)
      status = u->AddScaled(dt, state->velocity, d);
    if (status == StatusCode::kOk)
      status = u->AddScaled(coefficients.predictor_acceleration,
                            state->acceleration, d);
    auto* v = candidate->velocity;
    if (status == StatusCode::kOk) status = v->CopyFrom(state->velocity, d);
    if (status == StatusCode::kOk)
      status = v->AddScaled(coefficients.old_acceleration_velocity,
                            state->acceleration, d);
    auto* a = candidate->acceleration;
    if (status == StatusCode::kOk) status = a->CopyFrom(state->acceleration, d);
    if (status == StatusCode::kOk) status = a->Scale(0, d);
    candidate->valid = status == StatusCode::kOk;
    return status;
  }
  StatusCode UpdateKinematics(const TranslationalStateView* state, double dt,
                              TranslationalCandidate* candidate,
                              Diagnostic* d) const override {
    if (candidate) candidate->valid = false;
    if (!ValidViews(state, candidate, true)) return Invalid(d);
    Coefficients coefficients;
    auto status = ComputeCoefficients(parameters_, dt, &coefficients, d);
    if (status != StatusCode::kOk) return status;
    // H4 in displacement-difference form, without a retained affine vector.
    auto* a = candidate->acceleration;
    status = a->CopyFrom(candidate->displacement, d);
    if (status == StatusCode::kOk)
      status = a->AddScaled(-1, state->displacement, d);
    if (status == StatusCode::kOk)
      status = a->AddScaled(-dt, state->velocity, d);
    if (status == StatusCode::kOk)
      status = a->AddScaled(-coefficients.predictor_acceleration,
                            state->acceleration, d);
    if (status == StatusCode::kOk) status = a->Scale(coefficients.c0, d);
    auto* v = candidate->velocity;
    if (status == StatusCode::kOk) status = v->CopyFrom(state->velocity, d);
    if (status == StatusCode::kOk)
      status = v->AddScaled(coefficients.old_acceleration_velocity,
                            state->acceleration, d);
    if (status == StatusCode::kOk)
      status = v->AddScaled(coefficients.new_acceleration_velocity, a, d);
    candidate->valid = status == StatusCode::kOk;
    return status;
  }
  StatusCode EvaluateKinematicLinearization(double dt,
                                            TranslationalLinearization* out,
                                            Diagnostic* d) const override {
    if (!out) return Invalid(d);
    Coefficients coefficients;
    const auto status = ComputeCoefficients(parameters_, dt, &coefficients, d);
    if (status != StatusCode::kOk) return status;
    out->acceleration_per_displacement = coefficients.c0;
    out->velocity_per_displacement = coefficients.c1;
    return StatusCode::kOk;
  }
  StatusCode PrepareRotationalKinematics(const RotationalStateView* state,
                                         double dt,
                                         RotationalCandidate* candidate,
                                         bool* branch,
                                         Diagnostic* d) const override {
    if (candidate) candidate->valid = false;
    const auto status = EvaluateHhtRotations(parameters_, state, candidate, dt,
                                             HhtRotationOperation::kPredict,
                                             nullptr, branch, d);
    if (candidate) candidate->valid = status == StatusCode::kOk;
    return status;
  }
  StatusCode UpdateRotationalKinematics(const RotationalStateView* state,
                                        double dt,
                                        RotationalCandidate* candidate,
                                        bool* branch,
                                        Diagnostic* d) const override {
    if (candidate) candidate->valid = false;
    const auto status =
        EvaluateHhtRotations(parameters_, state, candidate, dt,
                             HhtRotationOperation::kUpdate, nullptr, branch, d);
    if (candidate) candidate->valid = status == StatusCode::kOk;
    return status;
  }
  StatusCode EvaluateRotationalLinearization(
      const RotationalStateView* state, const RotationalCandidate* candidate,
      double dt, RotationalLinearization* out, bool* branch,
      Diagnostic* d) const override {
    return EvaluateHhtRotations(parameters_, state, candidate, dt,
                                HhtRotationOperation::kLinearize, out, branch,
                                d);
  }

 private:
  const HhtParameters* parameters_;  // Borrowed, immutable for entire lifetime.
};
}  // namespace
StatusCode CreateHhtIntegrator(const HhtParameters* parameters,
                               TimeIntegrator** out, Diagnostic* d) {
  if (!out) return Invalid(d);
  *out = nullptr;
  if (!ValidParameters(parameters)) return Invalid(d);
  *out = new (std::nothrow) HhtIntegrator(parameters);
  if (!*out)
    return ReportFailure(StatusCode::kAllocationFailure,
                         DiagnosticStage::kAllocation, "CreateHhtIntegrator",
                         d);
  return StatusCode::kOk;
}
}  // namespace fesa
