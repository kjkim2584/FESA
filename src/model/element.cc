#include "fesa/model/element.h"

namespace fesa {

Element::Element(int id) : id_(id) {}
int Element::Id() const { return id_; }
int Element::internal_dof_count() const { return 0; }
StatusCode Element::GetDofNodes(int* node_ids, int capacity,
                                Diagnostic* d) const {
  if (!node_ids || node_count() <= 0 || capacity < node_count())
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation, "Element::GetDofNodes",
                         d);
  for (int i = 0; i < node_count(); ++i)
    if (!GetNode(i))
      return ReportFailure(StatusCode::kInvalidArgument,
                           DiagnosticStage::kValidation,
                           "Element::GetDofNodes missing node", d);
  for (int i = 0; i < node_count(); ++i) node_ids[i] = GetNode(i)->Id();
  return StatusCode::kOk;
}
StatusCode Element::CharacteristicLength(double*, Diagnostic* d) const {
  return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kComputation,
                       "Element characteristic length unsupported", d);
}
bool Element::SupportsMass(MassMode) const { return false; }
StatusCode Element::ComputeMass(MassMode, ElementWorkspace*, Matrix*,
                                Diagnostic* d) const {
  return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kComputation,
                       "Element mass unsupported", d);
}
StatusCode Element::LumpFromConsistentMass(const Matrix*, Matrix*,
                                           Diagnostic* d) const {
  return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kComputation,
                       "Element mass lumping unsupported", d);
}
StatusCode Element::ComputeStiffness(ElementWorkspace*, Matrix*,
                                     Diagnostic* d) const {
  return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kComputation,
                       "Element::ComputeStiffness unsupported", d);
}
StatusCode Element::ComputeInternalForce(const Vector*, ElementWorkspace*,
                                         Vector*, Diagnostic* d) const {
  return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kComputation,
                       "Element::ComputeInternalForce unsupported", d);
}
bool Element::SupportsTrialEvaluation() const { return false; }
StatusCode Element::EvaluateTrial(const ElementTrialState*,
                                  const ElementWorkspace*,
                                  ElementEvaluationRequest, Vector*, Matrix*,
                                  TrialEvaluationFailure* failure,
                                  Diagnostic* d) const {
  if (failure) *failure = TrialEvaluationFailure::kNone;
  return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kComputation,
                       "Element finite-state trial unsupported", d);
}
StatusCode Element::UpdateTrialResults(const ElementTrialState*,
                                       const ElementWorkspace*,
                                       TrialEvaluationFailure* failure,
                                       Diagnostic* d) {
  if (failure) *failure = TrialEvaluationFailure::kNone;
  return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kComputation,
                       "Element finite-state recovery unsupported", d);
}
bool Element::SupportsTrialGeometryCheck() const { return false; }
StatusCode Element::CheckTrialGeometry(
    const ElementTrialState*, const ElementWorkspace*,
    TrialEvaluationFailure* failure, Diagnostic* d) const {
  if (failure) *failure = TrialEvaluationFailure::kNone;
  return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kComputation,
                       "Element trial geometry check unsupported", d);
}
bool Element::SupportsDynamicEvaluation(DynamicKinematics) const {
  return false;
}
StatusCode Element::EvaluateDynamic(
    const ElementMotionInput*, MassMode, const ElementWorkspace*,
    DynamicEvaluationRequest, ElementDynamicResult*,
    TrialEvaluationFailure* failure, Diagnostic* d) const {
  if (failure) *failure = TrialEvaluationFailure::kNone;
  return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kComputation,
                       "Element dynamic evaluation unsupported", d);
}
StatusCode Element::UpdateResults(const Vector*, ElementWorkspace*,
                                  Diagnostic* d) {
  return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kComputation,
                       "Element::UpdateResults unsupported", d);
}
StatusCode Element::ComputeLoadedInternalForce(const Vector* u,
                                               const Vector* load,
                                               ElementWorkspace* workspace,
                                               Vector* force,
                                               Diagnostic* d) const {
  if (!load) return ComputeInternalForce(u, workspace, force, d);
  return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kComputation,
                       "Element load-aware force unsupported", d);
}
StatusCode Element::UpdateLoadedResults(const Vector* u, const Vector* load,
                                        ElementWorkspace* workspace,
                                        Diagnostic* d) {
  if (!load) return UpdateResults(u, workspace, d);
  return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kComputation,
                       "Element load-aware results unsupported", d);
}

}  // namespace fesa
