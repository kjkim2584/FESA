/// @file
/// @brief Non-installed verification seam for the approved P4b DL=.5 case.
#ifndef FESA_SRC_ANALYSIS_THERMAL_INCREMENT_POLICY_H_
#define FESA_SRC_ANALYSIS_THERMAL_INCREMENT_POLICY_H_
#include "fesa/analysis/thermal_increment_controller.h"
namespace fesa {
/// Call-local scalar policy. Only DL varies for verification; all other M1
/// coefficients and the IT3 ring stay fixed. No input/configuration path.
struct ThermalIncrementPolicy {
  double dl = .95;
};
/// Pure scalar endpoint arithmetic (formulation6.5). No state or ownership.
/// Invalid/non-progressing clock, or a nonterminal residual within the already
/// accumulated gamma bound, is kLibraryFailure. Failure preserves *dt/*end.
/// Distinct nonnull output pointers are required. Completion writes (0,period).
StatusCode PlanThermalIncrement(double period, double time, double proposal,
                                int accepted, double* dt, double* end,
                                Diagnostic* d);
/// Borrows parameters/policy only during creation, clones settings/copies DL.
/// Caller owns success; failure nulls output. Public Create uses default DL.
StatusCode CreateThermalIncrementControllerWithPolicy(
    const HeatTransferStepParameters* parameters,
    const ThermalIncrementPolicy* policy, ThermalIncrementController** out,
    Diagnostic* d);
}  // namespace fesa
#endif  // FESA_SRC_ANALYSIS_THERMAL_INCREMENT_POLICY_H_
