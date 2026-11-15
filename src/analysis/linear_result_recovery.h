/// @file
/// @brief Shared linear official result recovery after accepted computation.
#ifndef FESA_ANALYSIS_LINEAR_RESULT_RECOVERY_H_
#define FESA_ANALYSIS_LINEAR_RESULT_RECOVERY_H_

#include "analysis/analysis_access.h"
#include "fesa/analysis/analysis_model.h"
#include "fesa/analysis/analysis_state.h"
#include "fesa/analysis/dof_manager.h"
#include "fesa/elements/element_workspace.h"

namespace fesa {
/// Borrows model/mapping/state/workers for this call. Reads prepared physical
/// candidate displacement and reaction; never computes convergence/reactions.
/// Caller establishes node-ID/component correspondence; legacy 6N static
/// states without stored IDs remain supported. Physical sizes must match.
/// Marks official results unconfirmed with supplied step/increment before the
/// first node mutation, then updates node U/RF and element structural results.
/// Initial increment0 requires staged dynamic motion. No write/flush/commit,
/// extra geometry acceptance test or cutback. Failure may leave partial
/// unconfirmed results; caller discards candidate and terminates recovery.
/// Committed physical state is preserved. Workers remain alive through join;
/// optional nonthrowing observer has the existing node-update fault-seam role.
StatusCode RecoverLinearResults(const AnalysisModel* model,
                                const DofManager* dofs, AnalysisState* state,
                                int step, int increment,
                                ElementWorkspace** workers, int worker_count,
                                Diagnostic* d,
                                void (*observer)(AnalysisStage) = nullptr);
}  // namespace fesa
#endif  // FESA_ANALYSIS_LINEAR_RESULT_RECOVERY_H_
