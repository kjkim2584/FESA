/// @file
/// @brief Physical time-dependent CLOAD value/rate evaluation.
#ifndef FESA_ANALYSIS_DYNAMIC_LOAD_H_
#define FESA_ANALYSIS_DYNAMIC_LOAD_H_

#include "fesa/analysis/dof_manager.h"
#include "fesa/analysis/effective_conditions.h"
#include "fesa/math/vector.h"

namespace fesa {
/// Evaluates current resolved CLOAD origins at step/total time [T].
/// Requires a dynamic step, matching global 3N/6N mapping,
/// 0<=step_time<=period, total_time>=step_time, finite times; Beam line loads
/// are unsupported. Active values are reference force times relative amplitude;
/// no amplitude means constant force, removed origins are zero immediately.
/// STEP TIME holds the original definition step's end value after propagation;
/// TOTAL TIME uses total_time including prior static periods. No HHT weights or
/// static ramp. Creates global Cartesian external and unassembled origins in
/// the exact EffectiveConditions::Loads order, retaining fixed/cancelling/zero
/// entries. Components1..3 have F units;4..6 have FL units. Origins are decoded
/// with the corresponding load DOF, never reduced across force/moment units.
/// With no origins, *origins is null and *external is a zero vector.
/// Caller owns both successful vectors. Distinct creation-only outputs are
/// required, null on every failure. Inputs/state/history remain unchanged;
/// borrows inputs only for the call. Failure diagnostic may borrow Domain.
/// @see docs/HHT/io_contract.md
StatusCode EvaluateDynamicClLoads(const EffectiveConditions* conditions,
                                  const DofManager* dofs, double step_time,
                                  double total_time, Vector** external,
                                  Vector** origins, Diagnostic* d);
/// Same ownership, mapping, time and failure contract as the value evaluator,
/// returning the right-hand physical rate p_dot [F/T or FL/T]. Constant and
/// removed origins have zero rate; amplitude clocks/holds follow their existing
/// definition context. No HHT weighting, state mutation or velocity correction.
/// Both distinct creation-only outputs are null on every failure.
StatusCode EvaluateDynamicClLoadRates(const EffectiveConditions* conditions,
                                      const DofManager* dofs, double step_time,
                                      double total_time, Vector** external_rate,
                                      Vector** origin_rates, Diagnostic* d);
}  // namespace fesa
#endif  // FESA_ANALYSIS_DYNAMIC_LOAD_H_
