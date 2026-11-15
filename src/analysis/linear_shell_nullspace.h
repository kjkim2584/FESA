/// @file
/// @brief Call-local MITC4 reference/current geometric mass kernel.
#ifndef FESA_ANALYSIS_LINEAR_SHELL_NULLSPACE_H_
#define FESA_ANALYSIS_LINEAR_SHELL_NULLSPACE_H_

#include "fesa/analysis/analysis_model.h"
#include "fesa/analysis/dof_manager.h"
#include "fesa/math/matrix.h"

namespace fesa {
/// Owns reference mass-null directions and optional finite-pose directions;
/// legacy name retained for existing linear callers. No model/state retained.
/// Construct/destroy within one initialization attempt, never cache between
/// increments. Field rows use original node-ID order; each direction is an
/// allowed global angular direction. No artificial mass or motion adjustment.
class LinearShellNullspace final {
 public:
  /// Requires a matching bound linear, undamped, positive-density S4 model and
  /// 6N mapping. Caller deletes success. Failure nulls *out and cleans scratch;
  /// model, conditions, mapping and results remain unchanged. Uses only the
  /// existing dimensionless dense-QR rank criterion for incident constraints.
  static StatusCode Create(const AnalysisModel* model, const DofManager* dofs,
                           LinearShellNullspace** out, Diagnostic* d);
  /// Builds a current kernel from the same reference incidence geometry and
  /// borrowed node-ID-ordered n-by-4 unit initial-to-current wxyz orientations.
  /// Requires existing nonlinear rotation subset (all3 free or all3 fixed);
  /// partial rotations are unsupported. Does not decide step execution support.
  /// Preserves reference directions for P17 body-spin; Directions is current.
  /// Same creation-only ownership/failure contract; no inputs/results mutated.
  static StatusCode CreateFinite(const AnalysisModel* model,
                                 const DofManager* dofs,
                                 const Matrix* orientations,
                                 LinearShellNullspace** out, Diagnostic* d);
  ~LinearShellNullspace();
  LinearShellNullspace(const LinearShellNullspace&) = delete;
  LinearShellNullspace& operator=(const LinearShellNullspace&) = delete;
  LinearShellNullspace(LinearShellNullspace&&) = delete;
  LinearShellNullspace& operator=(LinearShellNullspace&&) = delete;
  /// Number of independent allowed rotational null directions, possibly0.
  int Count() const;
  /// Original node ID, or-1 for invalid direction index.
  int NodeId(int index) const;
  /// Physical global angular DOF index for component0..2, or-1 if invalid.
  int GlobalDof(int index, int component) const;
  /// Corresponding free equation, -1 if constrained/invalid.
  int Equation(int index, int component) const;
  /// Borrowed owned node-capacity-by3 matrix; only first Count() rows valid.
  /// Each valid row is a unit direction in global angular component order.
  const Matrix* Directions() const;
  /// Borrowed initial unit directions, same row layout as Directions. Linear
  /// instances return Directions itself; finite instances own a separate copy.
  const Matrix* ReferenceDirections() const;
  /// Borrowed node-capacity vector; first Count() entries are positive sums
  /// of incident fixed drilling coefficients, with FL/rad units.
  const Vector* DrillingCoefficients() const;
  /// Projects only angular components of a borrowed global6N field onto each
  /// direction. Input/output disjoint; out sizeCount, or null when Count0.
  /// No retention or model changes. Invalid shape preserves out; after a
  /// backend failure discard out. Result units equal the input angular units.
  StatusCode ProjectRotations(const Vector* global, Vector* out,
                              Diagnostic* d) const;
  /// P14 dimensionless initial compatibility errors, one per null direction.
  /// Finite-pose instances return unsupported: caller uses P17/N6/N7.
  /// Borrowed 6N motion/moment mean UR[rad]/M[FL] or omega[rad/T]/M_dot[FL/T].
  /// Writes abs(k*n.motion-n.moment)/(k*norm(motion)+norm(moment)); zero
  /// denominator gives zero. Caller applies <=1e-10 to each entry. No state,
  /// model or input changes. Same size/alias/empty/output-failure contract as
  /// ProjectRotations; out must also differ from moment. Scratch is call-local.
  StatusCode EvaluateCompatibility(const Vector* motion, const Vector* moment,
                                   Vector* out, Diagnostic* d) const;
  /// Applies P14 <=1e-10 independently at each allowed null direction.
  /// velocity selects omega/M_dot units and failure label; otherwise UR/M.
  /// Borrows inputs, retains/changes no state. Failure reports first original
  /// node and kEquilibriumFailure for incompatibility; arithmetic/backend
  /// errors propagate. Empty kernel succeeds after input shape validation.
  StatusCode CheckCompatibility(const Vector* motion, const Vector* moment,
                                bool velocity, Diagnostic* d) const;

 private:
  LinearShellNullspace() = default;
  /// Populates this not-yet-published owner; caller cleans it on any failure.
  StatusCode Initialize(const AnalysisModel* model, const DofManager* dofs,
                        bool finite, Diagnostic* d);
  /// Call-local current direction construction; only on unpublished owner.
  StatusCode Rotate(const Matrix* orientations, Diagnostic* d);
  int count_ = 0;
  int global_count_ = 0;
  int* nodes_ = nullptr;      ///< Owned node-capacity original IDs.
  int* global_ = nullptr;     ///< Owned3*node-capacity physical DOF indices.
  int* equations_ = nullptr;  ///< Owned3*node-capacity free equation indices.
  Matrix* reference_directions_ = nullptr;  ///< Owned finite reference copy.
  Matrix* directions_ = nullptr;            ///< Owned node-capacity-by3.
  Vector* coefficients_ = nullptr;          ///< Owned node-capacity.
};
}  // namespace fesa
#endif  // FESA_ANALYSIS_LINEAR_SHELL_NULLSPACE_H_
