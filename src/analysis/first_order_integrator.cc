#include "fesa/analysis/first_order_integrator.h"

#include <cmath>
#include <new>

namespace fesa {
StatusCode BackwardEulerIntegrator::Create(BackwardEulerIntegrator** out,
                                           Diagnostic* d) {
  if (!out)
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation,
                         "Null backward Euler output", d);
  *out = new (std::nothrow) BackwardEulerIntegrator;
  return *out ? StatusCode::kOk
              : ReportFailure(StatusCode::kAllocationFailure,
                              DiagnosticStage::kValidation,
                              "Backward Euler allocation", d);
}
StatusCode BackwardEulerIntegrator::CapacityCoefficient(double dt, double* out,
                                                        Diagnostic* d) const {
  if (!out || !std::isfinite(dt) || dt <= 0)
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation,
                         "Invalid backward Euler time increment", d);
  const double value = 1 / dt;
  if (!std::isfinite(value) || value <= 0)
    return ReportFailure(StatusCode::kLibraryFailure,
                         DiagnosticStage::kValidation,
                         "Unrepresentable backward Euler coefficient", d);
  *out = value;
  return StatusCode::kOk;
}
StatusCode BackwardEulerIntegrator::RecoverRate(const Vector* theta,
                                                const Vector* theta_old,
                                                double dt, Vector* rate,
                                                Diagnostic* d) const {
  double coefficient = 0;
  auto status = CapacityCoefficient(dt, &coefficient, d);
  if (status != StatusCode::kOk) return status;
  if (!theta || !theta_old || !rate || rate == theta || rate == theta_old ||
      theta->Size() != theta_old->Size() || theta->Size() != rate->Size())
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation,
                         "Invalid backward Euler rate vectors", d);
  status = rate->CopyFrom(theta, d);
  if (status == StatusCode::kOk) status = rate->AddScaled(-1, theta_old, d);
  if (status == StatusCode::kOk) status = rate->Scale(coefficient, d);
  return status;
}
}  // namespace fesa
