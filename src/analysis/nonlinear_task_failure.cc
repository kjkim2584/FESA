#include "analysis/nonlinear_task_failure.h"

#include <algorithm>

#include "base/diagnostic_copy.h"

namespace fesa {
StatusCode SelectNonlinearTaskFailure(StatusCode primary,
                                      const Element* const* elements, int count,
                                      const TrialEvaluationFailure* labels,
                                      const Diagnostic* workers,
                                      int worker_count, Diagnostic* selected,
                                      NonlinearFailure* failure) {
  if (failure) *failure = NonlinearFailure::kNone;
  if (!selected || !failure || count < 0 || worker_count <= 0 || !workers ||
      (count > 0 && (!elements || !labels)))
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation,
                         "Nonlinear task failure inputs", selected);
  if (primary != StatusCode::kEquilibriumFailure || count == 0) return primary;
  for (int i = 0; i < count; ++i)
    if (!elements[i] || (i > 0 && elements[i - 1]->Id() >= elements[i]->Id()))
      return ReportFailure(StatusCode::kInvalidArgument,
                           DiagnosticStage::kValidation,
                           "Nonlinear task failure element order", selected);
  const Diagnostic* fatal = nullptr;
  for (int worker = 0; worker < worker_count; ++worker) {
    const auto* candidate = workers + worker;
    if (candidate->code == StatusCode::kOk) continue;
    bool retryable = false;
    if (candidate->code == StatusCode::kEquilibriumFailure) {
      const auto found = std::lower_bound(
          elements, elements + count, candidate->element_id,
          [](const Element* element, int id) { return element->Id() < id; });
      if (found != elements + count &&
          (*found)->Id() == candidate->element_id) {
        const auto typed = labels[found - elements];
        retryable = typed == TrialEvaluationFailure::kGeometry ||
                    typed == TrialEvaluationFailure::kRotationBranch;
      }
    }
    if (!retryable && (!fatal || candidate->element_id < fatal->element_id))
      fatal = candidate;
  }
  if (fatal) {
    CopyDiagnosticMetadata(fatal, selected);
    return fatal->code;
  }
  for (int i = 0; i < count; ++i) {
    if (elements[i]->Id() != selected->element_id) continue;
    if (labels[i] == TrialEvaluationFailure::kGeometry)
      *failure = NonlinearFailure::kGeometry;
    else if (labels[i] == TrialEvaluationFailure::kRotationBranch)
      *failure = NonlinearFailure::kRotationBranch;
    break;
  }
  return primary;
}
}  // namespace fesa
