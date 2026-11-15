/// @file
/// @brief First-order thermal time integration, independent of mechanical
/// motion.
#ifndef FESA_ANALYSIS_FIRST_ORDER_INTEGRATOR_H_
#define FESA_ANALYSIS_FIRST_ORDER_INTEGRATOR_H_

#include "fesa/math/vector.h"

namespace fesa {
/// Pure numerical integration role; retains no state or borrowed arguments.
/// Time is in model T, temperature in theta, rate in theta/T.
class FirstOrderTimeIntegrator {
 public:
  virtual ~FirstOrderTimeIntegrator() = default;
  FirstOrderTimeIntegrator(const FirstOrderTimeIntegrator&) = delete;
  FirstOrderTimeIntegrator& operator=(const FirstOrderTimeIntegrator&) = delete;
  FirstOrderTimeIntegrator(FirstOrderTimeIntegrator&&) = delete;
  FirstOrderTimeIntegrator& operator=(FirstOrderTimeIntegrator&&) = delete;
  /// Writes capacity multiplier [1/T]. Nonpositive/nonfinite dt or null output
  /// is kInvalidArgument, preserving *out. Unrepresentable multiplier is
  /// kLibraryFailure. No state mutation, retention or allocation.
  virtual StatusCode CapacityCoefficient(double dt, double* out,
                                         Diagnostic* d) const = 0;
  /// Computes rate from borrowed same-size positive temperature vectors.
  /// rate must be distinct from both inputs; inputs may coincide. Invalid
  /// arguments preserve rate; a backend/numeric failure may invalidate it.
  /// Inputs/state are unchanged; no allocation or retained arguments.
  virtual StatusCode RecoverRate(const Vector* theta, const Vector* theta_old,
                                 double dt, Vector* rate,
                                 Diagnostic* d) const = 0;

 protected:
  FirstOrderTimeIntegrator() = default;
};
/// Backward Euler: C coefficient=1/dt, rate=(theta-theta_old)/dt.
/// @see docs/HeatTransfer/formulation.md section4.
class BackwardEulerIntegrator final : public FirstOrderTimeIntegrator {
 public:
  /// Caller owns success through this or the virtual base pointer. Failure
  /// nulls *out; out required and must not hold an existing object.
  static StatusCode Create(BackwardEulerIntegrator** out, Diagnostic* d);
  ~BackwardEulerIntegrator() override = default;
  BackwardEulerIntegrator(const BackwardEulerIntegrator&) = delete;
  BackwardEulerIntegrator& operator=(const BackwardEulerIntegrator&) = delete;
  BackwardEulerIntegrator(BackwardEulerIntegrator&&) = delete;
  BackwardEulerIntegrator& operator=(BackwardEulerIntegrator&&) = delete;
  StatusCode CapacityCoefficient(double dt, double* out,
                                 Diagnostic* d) const override;
  StatusCode RecoverRate(const Vector* theta, const Vector* theta_old,
                         double dt, Vector* rate, Diagnostic* d) const override;

 private:
  BackwardEulerIntegrator() = default;
};
}  // namespace fesa
#endif  // FESA_ANALYSIS_FIRST_ORDER_INTEGRATOR_H_
