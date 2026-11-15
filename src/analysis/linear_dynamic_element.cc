#include "analysis/linear_dynamic_element.h"

#include <cmath>
#include <memory>

#include "adapters/mkl/dense_ops.h"
#include "fesa/analysis/analysis_model.h"
#include "fesa/analysis/dof_manager.h"
#include "fesa/analysis/static_workspace.h"
#include "fesa/elements/beam3d.h"
#include "fesa/elements/mitc4_shell.h"
#include "fesa/elements/solid_element.h"
#include "model/model_storage.h"

namespace fesa {
StatusCode ResolveDynamicMassMode(const Element* element, MassMode requested,
                                  MassMode* out, Diagnostic* d) {
  if (!element || !out ||
      (requested != MassMode::kDefault && requested != MassMode::kConsistent &&
       requested != MassMode::kRowSum))
    return model_storage::Invalid("Dynamic mass request", d);
  const auto mode = requested == MassMode::kDefault
                        ? (element->SupportsMass(MassMode::kRowSum)
                               ? MassMode::kRowSum : MassMode::kConsistent)
                        : requested;
  if (!element->SupportsMass(mode)) {
    const auto status = ReportFailure(StatusCode::kUnsupported,
                                      DiagnosticStage::kValidation,
                                      "Unsupported dynamic mass mode", d);
    if (d) d->element_id = element->Id();
    return status;
  }
  *out = mode;
  return StatusCode::kOk;
}
bool HasMaterialRayleighDamping(const Element* element) {
  const auto* solid = dynamic_cast<const SolidElement*>(element);
  const auto* shell = dynamic_cast<const Mitc4Shell*>(element);
  const auto* beam = dynamic_cast<const Beam3D*>(element);
  const auto* material = solid ? solid->material()
                              : shell ? shell->material()
                              : beam ? beam->Section()->GetMaterial() : nullptr;
  return material && (material->rayleigh_alpha() != 0 ||
                      material->rayleigh_beta() != 0);
}
StatusCode ComputeLinearDynamicDamping(const Element* element,
                                       const ElementWorkspace* workspace,
                                       MassMode mode, Matrix* out,
                                       Diagnostic* d) {
  if (!element || !out)
    return model_storage::Invalid("Linear damping context", d);
  const auto resolved = ResolveDynamicMassMode(element, mode, &mode, d);
  if (resolved != StatusCode::kOk) return resolved;
  std::unique_ptr<Vector> vectors[3];
  for (auto& value : vectors) {
    Vector* raw = nullptr;
    const auto status = Vector::Create(element->dof_count(), &raw, d);
    value.reset(raw);
    if (status != StatusCode::kOk) return status;
  }
  std::unique_ptr<Matrix> matrices[4];
  for (auto& value : matrices) {
    Matrix* raw = nullptr;
    const auto status = Matrix::Create(element->dof_count(), element->dof_count(),
                                        &raw, d);
    value.reset(raw);
    if (status != StatusCode::kOk) return status;
  }
  ElementMotionInput input;
  input.linear_displacement = input.velocity = input.acceleration =
      vectors[0].get();
  ElementDynamicResult result;
  result.inertia_force = vectors[1].get();
  result.damping_force = vectors[2].get();
  result.inertia_displacement = matrices[0].get();
  result.inertia_velocity = matrices[1].get();
  result.inertia_acceleration = matrices[2].get();
  result.damping_displacement = matrices[3].get();
  result.damping_velocity = out;
  TrialEvaluationFailure failure = TrialEvaluationFailure::kNone;
  return element->EvaluateDynamic(&input, mode, workspace,
      DynamicEvaluationRequest::kForcesAndLinearization, &result, &failure, d);
}
StatusCode ValidateLinearDynamicModel(const AnalysisModel* model, MassMode mode,
                                      bool* shell, Diagnostic* d) {
  if (!model || !shell || !model->Conditions() || model->Elements()->empty() ||
      (mode != MassMode::kConsistent && mode != MassMode::kRowSum &&
       mode != MassMode::kDefault))
    return model_storage::Invalid("Linear dynamic model inputs", d);
  const auto family = model->Elements()->begin()->second->family();
  const bool is_shell = family == ElementFamily::kShell;
  const StepDefinition* step = nullptr;
  auto status = model->Conditions()->SourceDomain()->FindStepDefinition(
      model->Conditions()->StepNumber(), &step, d);
  if (status != StatusCode::kOk) return status;
  if (step->HasFiniteGeometry())
    return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                         "Linear dynamic kinematics only", d);
  for (const auto& entry : *model->Elements()) {
    const auto* element = entry.second;
    const auto* solid = dynamic_cast<const SolidElement*>(element);
    const auto* plate = dynamic_cast<const Mitc4Shell*>(element);
    if (element->family() != family ||
        element->dofs_per_node() != (is_shell ? 6 : 3) ||
        !(is_shell ? plate != nullptr : solid != nullptr) ||
        !element->SupportsDynamicEvaluation(DynamicKinematics::kLinear)) {
      status = ReportFailure(
          StatusCode::kUnsupported, DiagnosticStage::kValidation,
          "Linear dynamics requires homogeneous supported Solid or MITC4", d);
      if (d) d->element_id = entry.first;
      return status;
    }
    MassMode physical;
    status = ResolveDynamicMassMode(element, mode, &physical, d);
    if (status != StatusCode::kOk) return status;
    double density = 0;
    status = (is_shell ? plate->material() : solid->material())
                 ->Density(&density, d);
    if (status != StatusCode::kOk) {
      if (d) d->element_id = entry.first;
      return status;
    }
  }
  for (const auto& entry : *model->Nodes())
    if (entry.second->dof_count() != (is_shell ? 6 : 3))
      return model_storage::Invalid("Linear dynamic physical layout", d);
  *shell = is_shell;
  return StatusCode::kOk;
}
StatusCode ValidateLinearDynamicWorkers(bool shell, ElementWorkspace** workers,
                                        int count, Diagnostic* d) {
  if (count <= 0 || (!shell && !workers))
    return model_storage::Invalid("Linear dynamic workers", d);
  for (int i = 0; workers && i < count; ++i) {
    if (!workers[i] && !shell)
      return model_storage::Invalid("Linear dynamic worker", d);
    for (int j = 0; workers[i] && j < i; ++j)
      if (workers[i] == workers[j])
        return model_storage::Invalid("Aliased linear dynamic workers", d);
  }
  return StatusCode::kOk;
}
StatusCode ValidateHhtSolidAssembly(const AnalysisModel* model,
                                    const DofManager* dofs,
                                    const EffectiveConditions* conditions,
                                    StaticWorkspace* workspace,
                                    const HhtSolidAssemblyInput* dynamic,
                                    ElementWorkspace** workers,
                                    int worker_count, Diagnostic* d) {
  if (!model || !dofs || !conditions || !workspace || !dynamic)
    return model_storage::Invalid("HHT assembly context", d);
  if (worker_count <= 0 || !std::isfinite(dynamic->alpha) ||
      dynamic->alpha < -1. / 3 || dynamic->alpha > 0 ||
      !std::isfinite(dynamic->c0) || dynamic->c0 <= 0 ||
      (dynamic->mode != MassMode::kConsistent &&
       dynamic->mode != MassMode::kRowSum &&
       dynamic->mode != MassMode::kDefault) ||
      conditions->LineLoadState())
    return model_storage::Invalid("HHT reduced assembly inputs", d);
  bool shell = false;
  auto status = ValidateLinearDynamicModel(model, dynamic->mode, &shell, d);
  if (status == StatusCode::kOk)
    status = ValidateLinearDynamicWorkers(shell, workers, worker_count, d);
  if (status != StatusCode::kOk) return status;
  for (const auto* input :
       {dynamic->committed, dynamic->predictor, dynamic->external}) {
    if (!input || input->Size() != dofs->GlobalCount() ||
        input == workspace->Increment() ||
        input == workspace->ExternalForce() ||
        input == workspace->InternalForce() ||
        input == workspace->ReducedRhs() ||
        input == workspace->ReducedSolution())
      return model_storage::Invalid("HHT assembly motion/force shape/alias", d);
    for (int i = 0; i < input->Size(); ++i)
      if (!std::isfinite(input->Data()[i]))
        return model_storage::Invalid("HHT assembly finite input", d);
  }
  for (const auto& entry : *model->Elements()) {
    if (!HasMaterialRayleighDamping(entry.second)) continue;
    if (!std::isfinite(dynamic->c1) || dynamic->c1 <= 0)
      return model_storage::Invalid("HHT damping velocity coefficient", d);
    for (const auto* input : {dynamic->previous_velocity,
                             dynamic->predicted_velocity})
      if (!input || input->Size() != dofs->GlobalCount() ||
          !dense_ops::AllFinite(input->Size(), input->Data()) ||
          input == workspace->Increment() ||
          input == workspace->ExternalForce() ||
          input == workspace->InternalForce() ||
          input == workspace->ReducedRhs() ||
          input == workspace->ReducedSolution())
        return model_storage::Invalid("HHT damping velocity input", d);
  }
  for (const auto& entry : *model->Nodes())
    if (dofs->NodeDofCount(entry.first) != (shell ? 6 : 3))
      return model_storage::Invalid("HHT assembly requires 3N mapping", d);
  return StatusCode::kOk;
}
StatusCode EvaluateHhtSolidContribution(
    const Element* element, ElementWorkspace* workspace, MassMode mode,
    double alpha, double c0, const Vector* committed, const Vector* predictor,
    const Vector* prescribed, Matrix* effective, Vector* rhs, Diagnostic* d,
    const Vector* previous_velocity, const Vector* predicted_velocity,
    double c1) {
  const int size = element ? element->dof_count() : 0;
  if (!element || (!workspace && element->family() != ElementFamily::kShell) ||
      !effective || !rhs || effective->Rows() != size ||
      effective->Columns() != size || rhs->Size() != size ||
      !std::isfinite(alpha) || alpha < -1. / 3 || alpha > 0 ||
      !std::isfinite(c0) || c0 <= 0 ||
      (mode != MassMode::kConsistent && mode != MassMode::kRowSum &&
       mode != MassMode::kDefault))
    return model_storage::Invalid("HHT element contribution inputs", d);
  for (const auto* input : {committed, predictor, prescribed}) {
    if (!input || input == rhs || input->Size() != size ||
        !dense_ops::AllFinite(size, input->Data()))
      return model_storage::Invalid("HHT element contribution displacement", d);
  }
  const bool damped = HasMaterialRayleighDamping(element);
  if (damped) {
    if (!std::isfinite(c1) || c1 <= 0)
      return model_storage::Invalid("HHT damping coefficient", d);
    for (const auto* input : {previous_velocity, predicted_velocity})
      if (!input || input == rhs || input->Size() != size ||
          !dense_ops::AllFinite(size, input->Data()))
        return model_storage::Invalid("HHT damping history", d);
  }
  if ((element->family() != ElementFamily::kSolid &&
       element->family() != ElementFamily::kShell) ||
      !element->SupportsDynamicEvaluation(DynamicKinematics::kLinear))
    return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                         "HHT contribution requires supported linear element", d);
  auto status = ResolveDynamicMassMode(element, mode, &mode, d);
  if (status != StatusCode::kOk) return status;
  Matrix* raw_mass = nullptr;
  status = Matrix::Create(size, size, &raw_mass, d);
  std::unique_ptr<Matrix> mass(raw_mass);
  if (status != StatusCode::kOk) return status;
  Vector* raw_force = nullptr;
  status = Vector::Create(size, &raw_force, d);
  std::unique_ptr<Vector> force(raw_force);
  if (status != StatusCode::kOk) return status;
  status = element->ComputeStiffness(workspace, effective, d);
  if (status != StatusCode::kOk) return status;
  status = element->ComputeMass(mode, workspace, mass.get(), d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::MultiplyVector(effective, false, committed, rhs, d);
  if (status != StatusCode::kOk) return status;
  status = rhs->Scale(alpha, d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::MultiplyVector(mass.get(), false, predictor, force.get(), d);
  if (status != StatusCode::kOk) return status;
  status = rhs->AddScaled(c0, force.get(), d);
  if (status != StatusCode::kOk) return status;
  status = effective->Scale(1 + alpha, d);
  if (status != StatusCode::kOk) return status;
  status = effective->AddScaled(c0, mass.get(), d);
  if (status != StatusCode::kOk) return status;
  if (damped) {
    Matrix* raw_c = nullptr;
    status = Matrix::Create(size, size, &raw_c, d);
    std::unique_ptr<Matrix> damping(raw_c);
    if (status != StatusCode::kOk) return status;
    status = ComputeLinearDynamicDamping(element, workspace, mode, raw_c, d);
    Vector* raw_v = nullptr;
    if (status == StatusCode::kOk) status = Vector::Create(size, &raw_v, d);
    std::unique_ptr<Vector> drive(raw_v);
    if (status == StatusCode::kOk) status = drive->CopyFrom(predictor, d);
    if (status == StatusCode::kOk) status = drive->Scale(c1, d);
    if (status == StatusCode::kOk)
      status = drive->AddScaled(-1, predicted_velocity, d);
    if (status == StatusCode::kOk) status = drive->Scale(1 + alpha, d);
    if (status == StatusCode::kOk)
      status = drive->AddScaled(alpha, previous_velocity, d);
    if (status == StatusCode::kOk)
      status = Matrix::MultiplyVector(raw_c, false, raw_v, force.get(), d);
    if (status == StatusCode::kOk) status = rhs->AddScaled(1, force.get(), d);
    if (status == StatusCode::kOk)
      status = effective->AddScaled((1 + alpha) * c1, raw_c, d);
    if (status != StatusCode::kOk) return status;
  }
  status = Matrix::MultiplyVector(effective, false, prescribed, force.get(), d);
  if (status != StatusCode::kOk) return status;
  return rhs->AddScaled(-1, force.get(), d);
}
}  // namespace fesa
