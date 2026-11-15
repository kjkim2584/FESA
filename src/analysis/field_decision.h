/// @file
/// @brief Shared scalar field criteria, with no mechanical or thermal state.
#ifndef FESA_SRC_ANALYSIS_FIELD_DECISION_H_
#define FESA_SRC_ANALYSIS_FIELD_DECISION_H_
#include "fesa/base/diagnostic.h"
namespace fesa {
struct FieldControls;
struct FieldConvergenceOptions;
/// Field-independent branch; consumers map it to their own report enum.
enum class ScalarFieldBranch { kNoDofs, kZeroFlux, kLinear, kGeneral };
/// Borrowed measurements in one consistent field's physical units.
struct ScalarFieldDecisionInput {
  int dof_count = 0;
  double maximum = 0;
  double spatial_mean = 0;
  double used_mean = 0;
  double increment_max = 0;
  double correction_max = 0;
  bool allow_linear = true;  ///< Mechanical all-active-field gate may disable.
};
/// Caller-owned scalar outputs; no numerical or physical state retained.
struct ScalarFieldDecision {
  double phi = 0;
  double residual_limit = 0;
  double correction_limit = 0;
  ScalarFieldBranch branch = ScalarFieldBranch::kNoDofs;
  bool passed = true;
};
/// Abaqus public field criteria, formulation7. Null options preserves legacy
/// actual-correction semantics. small_increment_exception applies only to
/// displacement with positive characteristic length, never temperature.
/// Borrows all inputs. Invalid preflight preserves out; later numerical error
/// invalidates the decision. Nonconvergence returns kEquilibriumFailure with
/// complete output. No allocation, model/state change, or library operation.
StatusCode EvaluateScalarFieldDecision(const FieldControls* controls,
                                       const FieldConvergenceOptions* options,
                                       bool small_increment_exception,
                                       double characteristic_length,
                                       const ScalarFieldDecisionInput* input,
                                       ScalarFieldDecision* out, Diagnostic* d);
}  // namespace fesa
#endif  // FESA_SRC_ANALYSIS_FIELD_DECISION_H_
