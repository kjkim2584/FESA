/// @file
/// @brief Thermal element response boundary independent of mechanical elements.
#ifndef FESA_MODEL_THERMAL_ELEMENT_H_
#define FESA_MODEL_THERMAL_ELEMENT_H_

#include <vector>

#include "fesa/math/matrix.h"
#include "fesa/model/node.h"
#include "fesa/model/thermal_condition.h"
#include "fesa/model/thermal_material.h"

namespace fesa {
class Domain;
/// Call-local distributed load, formulation8.4/8.5/8.6.1. All pointers are
/// borrowed only during evaluation. magnitude is inward power per volume/area
/// [E/T/L^3 or E/T/L^2]; face_label is BF or the normalized S1..S6/SNEG/SPOS
/// spelling.
struct ThermalLoadContext {
  ThermalConditionKind kind = ThermalConditionKind::kBodyFlux;
  const char* face_label = "BF";
  double magnitude = 0;
  double coefficient = 0;  ///< Effective scalar h, ignored for property table.
  const FilmPropertyTable* film_property =
      nullptr;                  ///< Borrowed immutable h(T).
  double sink_temperature = 0;  ///< End-time effective sink [theta].
  double sink_baseline = 0;  ///< Effective sink at current step start [theta].
  double ambient_temperature = 0;  ///< Radiation end-time ambient [theta].
  double ambient_baseline = 0;     ///< Radiation step-start ambient [theta].
  double emissivity = 0;           ///< Radiation emissivity, in [0,1].
  double sigma = 0;                ///< Positive Stefan-Boltzmann constant.
  double absolute_zero = 0;        ///< Temperature origin [theta].
  /// Borrowed local dof_count mask:1 only for free DOFs at an actual
  /// extrapolated iteration0 start (D13). Null means no negative-temperature
  /// exception. Prescribed/later-iterate/start-state checks remain strict.
  const int* extrapolated_free_dofs = nullptr;
  double step_scale = 1;     ///< A_step;1 also disables continuity correction.
  const Vector* step_start_temperature = nullptr;  ///< Borrowed local theta_s.
  /// Borrowed local first accepted T; null means h_c is current h(theta).
  const Vector* first_accepted_temperature = nullptr;
  const Vector* previous_outward_flux =
      nullptr;  ///< Borrowed local p[E/T],null0.
};
/// Overwrites caller-owned IDs with sorted unique thermal element node IDs.
/// Borrows Domain for this call; conditions never create participation.
/// Null arguments or allocation failure preserve IDs and all Domain state.
/// Used by input IC registration and thermal condition/DOF expansion.
StatusCode CollectThermalNodeIds(const Domain* domain, std::vector<int>* ids,
                                 Diagnostic* d);
/// Approved solid and shell thermal input identities.
enum class ThermalElementType { kDc3d4, kDc3d8, kDs4 };
/// Domain-owned immutable connectivity and response. No retained trial state.
/// Material and connected nodes are borrowed until element destruction.
class ThermalElement {
 public:
  virtual ~ThermalElement() = default;
  ThermalElement(const ThermalElement&) = delete;
  ThermalElement& operator=(const ThermalElement&) = delete;
  ThermalElement(ThermalElement&&) = delete;
  ThermalElement& operator=(ThermalElement&&) = delete;
  /// Original positive element identifier.
  int Id() const;
  /// Approved input identity and static keyword spelling.
  virtual ThermalElementType type() const = 0;
  virtual const char* type_name() const = 0;
  /// Connectivity size and thermal components per connected node.
  virtual int node_count() const = 0;
  virtual int dofs_per_node() const = 0;
  /// Total thermal components, ordered node then thickness point.
  int dof_count() const;
  /// Borrowed node at a zero-based local index; null outside connectivity.
  virtual const Node* GetNode(int local) const = 0;
  /// Borrowed immutable thermal material until Domain destruction.
  virtual const ThermalMaterial* material() const = 0;
  /// Shell-only section thickness [L]. Borrows out for this call; an override
  /// writes a positive finite value from its section without mutation or
  /// retention. Default kUnsupported preserves *out (including on solids).
  virtual StatusCode GetShellThickness(double* out, Diagnostic* d) const;
  /// Checks positive Jacobian/approved geometry without changing any state.
  virtual StatusCode ValidateGeometry(Diagnostic* d) const = 0;
  /// Overwrites dof_count-square symmetric K [E/T/theta], formulation §3/8.3.
  /// Null theta selects constant k; temperature-dependent conductivity
  /// requires a finite local temperature vector, borrowed only for this call.
  /// DC8 uses Gauss T, DC4 mean nodal k, DS4 natural-sign nodal k/Simpson.
  /// Omits dk/dtheta (modified Newton). Failure may invalidate out; no model
  /// mutation or retained input pointer.
  virtual StatusCode ComputeConduction(const Vector* element_theta, Matrix* out,
                                       Diagnostic* d) const = 0;
  /// Overwrites dof_count station volume weights [L^3], without rho*c.
  /// Failure may invalidate out. All output storage is caller-owned.
  virtual StatusCode ComputeCapacityWeights(Vector* out,
                                            Diagnostic* d) const = 0;
  /// Overwrites full element face-area weights[L²], zero on other DOFs.
  /// Normalized face labels S1..S6/SNEG/SPOS only; no input retention/model
  /// mutation. Default unsupported preserves out; supported failure may
  /// invalidate out. Active zero loads still validate the loaded geometry.
  virtual StatusCode ComputeSurfaceWeights(const char* face_label, Vector* out,
                                           Diagnostic* d) const;
  /// Distributed BF/face flux and nodal film/radiation hook.
  /// Borrows all inputs, retains none. Overwrites Q[E/T] and optional residual
  /// Jacobian d(-Q)/dtheta[E/T/theta], zero for flux kinds. Film uses a
  /// modified tangent that omits hprime from time correction (8.6.1).
  /// Film requires finite local theta and, when A<1, local theta_s. Null p
  /// means new0. Unsupported kind preserves outputs; supported numerical
  /// failure may invalidate outputs, never model/state data.
  /// Radiation uses nodal area and full 4*epsilon*sigma*(theta-Z)^3 tangent.
  /// Selected surface absolute temperatures theta-Z and theta_s-Z must be
  /// nonnegative; D13 permits a negative trial only where the call-local
  /// extrapolated_free_dofs mask is1. Start states remain nonnegative.
  /// Negative failures are never clamped. Zero is valid. Ambient>=Z.
  virtual StatusCode EvaluateLoadContribution(const ThermalLoadContext* load,
                                              const Vector* element_theta,
                                              Vector* external_flux,
                                              Matrix* jacobian,
                                              Diagnostic* d) const;

 protected:
  /// Sets immutable ID; derived factory validates definitions first.
  explicit ThermalElement(int id);

 private:
  int id_;
};
}  // namespace fesa

#endif  // FESA_MODEL_THERMAL_ELEMENT_H_
