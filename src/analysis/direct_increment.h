/// @file
/// @brief Resolved nominal increment for fixed HHT time stepping.
#ifndef FESA_ANALYSIS_DIRECT_INCREMENT_H_
#define FESA_ANALYSIS_DIRECT_INCREMENT_H_

#include "analysis/increment_control.h"
#include "fesa/base/diagnostic.h"

namespace fesa {
/// Resolves DIRECT scalar settings, all in the model's time unit.
/// requested/period must be positive finite; minimum/maximum nonnegative
/// finite. Zero denotes unspecified. An explicit minimum cannot exceed the
/// requested increment or a nonzero maximum. The resolved nominal increment
/// is min(requested, period, nonzero maximum). Period clipping may put it
/// below minimum, as may a later terminal remainder; no cutback is implied.
/// Writes the required scalar only on success. No state, allocation or backend
/// work. Caller retains original settings/source and reports any reduction.
/// @see docs/HHT/implementation_contract.md
StatusCode ResolveDirectIncrement(double requested, double period,
                                  double minimum, double maximum,
                                  double* nominal, Diagnostic* d);
/// Initializes caller-owned scalar clock; initial frame consumes no increment.
/// All failures preserve out; no retained inputs, allocation or backend work.
StatusCode InitializeDirectIncrementClock(double requested, double period,
                                          double minimum, double maximum,
                                          IncrementClock* out, Diagnostic* d);
/// Proposes from the committed clock without changing it. Uses compensated
/// time addition and exact period at the terminal clip; no epsilon skip,
/// extrapolation or growth. A round-trip-compatible uniform h=T/N grid closes
/// at its Nth step with the actual remaining terminal size, including rounding.
/// INC counts successful advancing frames. See direct_endpoint_roundtrip.md.
/// Limit exhaustion is kEquilibriumFailure, unrepresentable advance is
/// kLibraryFailure. All failures preserve out. Inputs are call-scoped borrows.
StatusCode ProposeDirectIncrement(double period, int increment_limit,
                                  const IncrementClock* clock,
                                  IncrementTarget* out, Diagnostic* d);
/// Stages metadata only for kOk numerical outcome and the exact bound proposal.
/// Caller promotes this distinct output in the no-fail physical commit AFTER
/// recovery/output/flush. Discarding it preserves old time/counts. Failure
/// preserves out, propagates numerical status without retry. DIRECT has no
/// fabricated Newton iteration count or quick-success streak.
StatusCode AcceptDirectIncrement(double period, int increment_limit,
                                 const IncrementClock* clock,
                                 const IncrementTarget* target,
                                 StatusCode numerical_status,
                                 IncrementClock* out, Diagnostic* d);
}  // namespace fesa

#endif  // FESA_ANALYSIS_DIRECT_INCREMENT_H_
