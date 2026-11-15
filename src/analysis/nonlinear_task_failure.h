/// @file
/// @brief Shared post-join fatal/typed element-failure selection.
#ifndef FESA_ANALYSIS_NONLINEAR_TASK_FAILURE_H_
#define FESA_ANALYSIS_NONLINEAR_TASK_FAILURE_H_

#include "fesa/analysis/nonlinear_problem.h"
#include "fesa/model/element.h"
#include "fesa/model/element_trial.h"

namespace fesa {
/// Consumes the common executor's primary status/root after all workers join.
/// Borrows original-ID-sorted elements[count], per-item typed labels and
/// workers[worker_count] diagnostics. Elements are valid, IDs strictly ascend;
/// count0 permits null elements/labels, worker_count>0. selected is the
/// executor's separate required root diagnostic; required failure is cleared.
/// An observed fatal outranks retryable geometry/branch; lowest fatal ID wins.
/// All-retryable errors retain the executor's root/ID. No unexecuted item is
/// treated as observed failure. Returns selected status and typed cause.
/// No allocation/backend, callback, state/result mutation or retained pointers.
/// selected must not alias the borrowed worker diagnostics. Invalid metadata
/// reports validation failure there; worker arrays/labels remain unchanged.
StatusCode SelectNonlinearTaskFailure(StatusCode primary,
                                      const Element* const* elements, int count,
                                      const TrialEvaluationFailure* labels,
                                      const Diagnostic* workers,
                                      int worker_count, Diagnostic* selected,
                                      NonlinearFailure* failure);
}  // namespace fesa
#endif  // FESA_ANALYSIS_NONLINEAR_TASK_FAILURE_H_
