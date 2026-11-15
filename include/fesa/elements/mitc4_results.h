/// @file
/// @brief One owned official MITC4 output payload, no numerical cache.
#ifndef FESA_ELEMENTS_MITC4_RESULTS_H_
#define FESA_ELEMENTS_MITC4_RESULTS_H_
#include "fesa/math/matrix.h"
#include "fesa/model/element_trial.h"
#include "fesa/model/node.h"
#include "fesa/model/shell_section.h"
namespace fesa {
/// Element owns this result set. Common AnalysisState owns validity: initial
/// zero buffers do not constitute confirmed results. No input/state retained.
class Mitc4Results final {
 public:
  /// Allocates fixed official output buffers. Caller deletes success; failure
  /// nulls*out and releases partial resources. d optional, all inputs borrowed.
  static StatusCode Create(Mitc4Results** out, Diagnostic* d);
  ~Mitc4Results();
  Mitc4Results(const Mitc4Results&) = delete;
  Mitc4Results& operator=(const Mitc4Results&) = delete;
  Mitc4Results(Mitc4Results&&) = delete;
  Mitc4Results& operator=(Mitc4Results&&) = delete;
  /// Recalculates direct force and linear S/LE at8IPs and8surface points,
  /// then16 element-node rows. No K/M or solver/commit/I/O. Caller marks
  /// results unconfirmed before this call. Preflight preserves payload;
  /// later failure may leave partial unconfirmed results, never rollback.
  /// u24 global node-major[U,UR], positive initial kd. Material may lackrho.
  StatusCode RecoverLinear(const Node* const* nodes, int count,
                           const ShellSection* section, const Vector* u,
                           double drilling, Diagnostic* d);
  /// Recomputes force only and finite Cauchy S/physical logarithmic LE at the
  /// same sample locations as linear recovery, then common layer extrapolation.
  /// Positive fixed initial kd, no density/K/M or state/history mutation.
  /// Inputs borrowed, failure required/cleared. Caller marks results
  /// unconfirmed first. Preflight preserves; later failure may leave partials.
  /// Current output basis and typed Gauss/tying/surface failures follow
  /// docs/MITC4/nonlinear_recovery_execution_contract.md.
  StatusCode RecoverNonlinear(const Node* const* nodes, int count,
                              const ShellSection* section,
                              const ElementTrialState* state, double drilling,
                              TrialEvaluationFailure* failure, Diagnostic* d);
  /// Borrowed global24 node-major force/moment, includes drilling springs.
  const Vector* internal_force() const;
  /// Borrowed3x3 common row basis[a,b,n], initial basis for linear recovery.
  const Matrix* output_axes() const;
  /// Borrowed8x6 S: [11,22,33,12,13,23], stress units, tensor shear.
  /// IP rows: plane1(-a,-a),2(a,-a),3(a,a),4(-a,a), each zeta[-a,+a].
  const Matrix* ip_stress() const;
  /// Borrowed8x6 LE: infinitesimal in linear, physical log in nonlinear;
  /// shear components engineering2epsilon. Uses physical gradient, not MITCB.
  const Matrix* ip_strain() const;
  /// Borrowed8x6 S/LE directly at zeta=-1 then+1, four planeGauss points each.
  const Matrix* surface_stress() const;
  const Matrix* surface_strain() const;
  /// Borrowed16x6 inplane extrapolated S/LE, no shared-node averaging.
  /// Layer order zeta=-a,+a,-1,+1, each four original connectivity corners.
  const Matrix* nodal_stress() const;
  const Matrix* nodal_strain() const;

 private:
  Mitc4Results() = default;
  Vector* internal_force_ = nullptr;  ///< Owned global resisting force.
  Matrix* axes_ = nullptr;            ///< Owned current official output axes.
  Matrix* ip_stress_ = nullptr;       ///< Owned8x6.
  Matrix* ip_strain_ = nullptr;       ///< Owned8x6.
  Matrix* surface_stress_ = nullptr;  ///< Owned8x6.
  Matrix* surface_strain_ = nullptr;  ///< Owned8x6.
  Matrix* nodal_stress_ = nullptr;    ///< Owned16x6.
  Matrix* nodal_strain_ = nullptr;    ///< Owned16x6.
};
}  // namespace fesa
#endif  // FESA_ELEMENTS_MITC4_RESULTS_H_
