/// @file
/// @brief Node-owned temperatures and reaction heat fluxes in thermal DOF
/// order.
#ifndef FESA_MODEL_NODAL_TEMPERATURE_RESULT_H_
#define FESA_MODEL_NODAL_TEMPERATURE_RESULT_H_

#include "fesa/math/vector.h"

namespace fesa {
/// Owns NT and RFL vectors. Index j corresponds to Abaqus DOF label 11+j.
/// NT uses the input temperature scale; RFL [E/T] is positive into the body.
/// The owning thermal Node controls lifetime and single-writer access.
class NodalTemperatureResult final {
 public:
  /// Creates zero results for 1..20 components; caller owns success.
  /// Failure nulls *out and releases partial allocations. d is borrowed.
  static StatusCode Create(int count, NodalTemperatureResult** out,
                           Diagnostic* d);
  ~NodalTemperatureResult();
  NodalTemperatureResult(const NodalTemperatureResult&) = delete;
  NodalTemperatureResult& operator=(const NodalTemperatureResult&) = delete;
  NodalTemperatureResult(NodalTemperatureResult&&) = delete;
  NodalTemperatureResult& operator=(NodalTemperatureResult&&) = delete;
  /// Fixed component count, shared by NT and RFL.
  int DofCount() const;
  /// Read-only borrowed values, valid until the owning Node is destroyed.
  const Vector* Temperatures() const;
  /// Read-only borrowed reaction fluxes, valid for the same lifetime.
  const Vector* ReactionFluxes() const;
  /// Copies borrowed, finite vectors of exactly DofCount() components.
  /// Rejects NT input==owned RFL or RFL input==owned NT before mutation;
  /// same-slot self-copy and a shared external input are allowed.
  /// Validation failure preserves both vectors; backend failure may leave
  /// partial results. Caller marks ThermalState unconfirmed before updating.
  /// No state commit or rollback occurs here (implementation_contract §3.1).
  StatusCode Update(const Vector* nt, const Vector* rfl, Diagnostic* d);

 private:
  NodalTemperatureResult() = default;
  Vector* temperatures_ = nullptr;     ///< Owned NT, initially zero.
  Vector* reaction_fluxes_ = nullptr;  ///< Owned RFL, initially zero.
};
}  // namespace fesa

#endif  // FESA_MODEL_NODAL_TEMPERATURE_RESULT_H_
