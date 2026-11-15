#include "fesa/elements/incremental_solid_element.h"

#include <memory>
#include <new>

#include "elements/incremental_solid_geometry.h"
#include "elements/incremental_hexa8_ops.h"
#include "elements/incremental_hexa8_workspace.h"
#include "elements/incremental_solid_ops.h"
#include "elements/incremental_solid_recovery.h"
#include "elements/solid_element_read_view.h"
#include "elements/solid_trial_geometry.h"
#include "fesa/model/incremental_solid_trial_context.h"
#include "model/model_storage.h"

namespace fesa {
IncrementalSolidElement::IncrementalSolidElement(
    int id, const SolidSection* section, SolidFormulation formulation)
    : Element(id), section_(section), formulation_(formulation),
      node_count_(formulation == SolidFormulation::kTetra4Full ? 4 :
          formulation == SolidFormulation::kTetra10Full ? 10 :
          formulation == SolidFormulation::kWedge6 ? 6 :
          formulation == SolidFormulation::kPyramid5Full ? 5 :
          formulation == SolidFormulation::kHexa8Full ? 8 : 0),
      point_count_(formulation == SolidFormulation::kTetra4Full ? 1 :
          formulation == SolidFormulation::kTetra10Full ? 4 :
          formulation == SolidFormulation::kWedge6 ? 2 :
          formulation == SolidFormulation::kPyramid5Full ? 5 :
          formulation == SolidFormulation::kHexa8Full ? 8 : 0) {}
IncrementalSolidElement::~IncrementalSolidElement() {
  delete[] nodes_;
  if (gradients_)
    for (int g = 0; g < point_count_; ++g) delete gradients_[g];
  delete[] gradients_;
  delete weight_;
  delete force_;
  delete stress_;
  delete strain_;
  delete nodal_stress_;
  delete nodal_strain_;
}
ElementFamily IncrementalSolidElement::family() const {
  return ElementFamily::kSolid;
}
const char* IncrementalSolidElement::type_name() const {
  switch (formulation_) {
    case SolidFormulation::kTetra4Full: return "C3D4";
    case SolidFormulation::kTetra10Full: return "C3D10";
    case SolidFormulation::kWedge6: return "C3D6";
    case SolidFormulation::kPyramid5Full: return "C3D5";
    case SolidFormulation::kHexa8Full: return "C3D8";
    default: return "unsupported";
  }
}
SolidFormulation IncrementalSolidElement::formulation() const {
  return formulation_;
}
int IncrementalSolidElement::node_count() const { return node_count_; }
int IncrementalSolidElement::dofs_per_node() const { return 3; }
int IncrementalSolidElement::dof_count() const { return 3 * node_count_; }
int IncrementalSolidElement::integration_point_count() const {
  return point_count_;
}
const Node* IncrementalSolidElement::GetNode(int local) const {
  return nodes_ && local >= 0 && local < node_count_ ? nodes_[local] : nullptr;
}
const Vector* IncrementalSolidElement::internal_force() const { return force_; }
const Matrix* IncrementalSolidElement::ip_stress() const { return stress_; }
const Matrix* IncrementalSolidElement::ip_strain() const { return strain_; }
const Matrix* IncrementalSolidElement::nodal_stress() const {
  return nodal_stress_;
}
const Matrix* IncrementalSolidElement::nodal_strain() const {
  return nodal_strain_;
}
const SolidSection* IncrementalSolidElement::section() const {
  return section_;
}
const ElasticMaterial* IncrementalSolidElement::material() const {
  return section_->material();
}
StatusCode IncrementalSolidElement::Initialize(const Node* const* nodes,
                                               Diagnostic* d) {
  if (!nodes || node_count_ <= 0 || nodes_)
    return model_storage::Invalid("Incremental solid nodes", d);
  for (int i = 0; i < node_count_; ++i) {
    if (!nodes[i] || nodes[i]->dof_count() != 3)
      return model_storage::Invalid("Incremental solid node layout", d);
    for (int j = 0; j < i; ++j)
      if (nodes[i]->Id() == nodes[j]->Id())
        return model_storage::Invalid("Incremental solid repeated node", d);
  }
  nodes_ = new (std::nothrow) const Node*[node_count_];
  gradients_ = new (std::nothrow) Matrix*[point_count_]();
  if (!nodes_ || !gradients_)
    return model_storage::AllocationFailure("Incremental connectivity/H", d);
  for (int i = 0; i < node_count_; ++i) nodes_[i] = nodes[i];
  StatusCode status = StatusCode::kOk;
  for (int g = 0; g < point_count_ && status == StatusCode::kOk; ++g)
    status = Matrix::Create(3, node_count_, &gradients_[g], d);
  if (status == StatusCode::kOk)
    status = Vector::Create(point_count_, &weight_, d);
  if (status == StatusCode::kOk)
    status = Vector::Create(dof_count(), &force_, d);
  if (status == StatusCode::kOk)
    status = Matrix::Create(point_count_, 6, &stress_, d);
  if (status == StatusCode::kOk)
    status = Matrix::Create(point_count_, 6, &strain_, d);
  if (status == StatusCode::kOk)
    status = Matrix::Create(node_count_, 6, &nodal_stress_, d);
  if (status == StatusCode::kOk)
    status = Matrix::Create(node_count_, 6, &nodal_strain_, d);
  return status;
}
StatusCode IncrementalSolidElement::ValidateGeometry(Diagnostic* d) const {
  SolidReferenceData* raw = nullptr;
  auto status = SolidReferenceData::Create(formulation_,
      IntegrationPurpose::kStiffness, &raw, d);
  std::unique_ptr<SolidReferenceData> data(raw);
  if (status != StatusCode::kOk) return status;
  if (data->node_count() != node_count_ || data->point_count() != point_count_)
    return model_storage::Invalid("Incremental immutable layout", d);
  return BuildIncrementalInitialGeometry(this, data.get(), gradients_,
                                         weight_, d);
}
StatusCode IncrementalSolidElement::CharacteristicLength(double* out,
                                                         Diagnostic* d) const {
  return ComputeIncrementalCharacteristicLength(this, out, d);
}
bool IncrementalSolidElement::SupportsTrialEvaluation() const { return true; }
bool IncrementalSolidElement::SupportsTrialGeometryCheck() const {
  return true;
}
StatusCode IncrementalSolidElement::EvaluateTrial(
    const ElementTrialState* state, const ElementWorkspace* workspace,
    ElementEvaluationRequest request, Vector* force, Matrix* tangent,
    TrialEvaluationFailure* failure, Diagnostic* d) const {
  if (failure) *failure = TrialEvaluationFailure::kNone;
  const auto* context = state ? state->incremental_solid : nullptr;
  const auto* reference =
      workspace
          ? workspace->reference(formulation(), IntegrationPurpose::kStiffness)
          : nullptr;
  if (!failure || !context || !context->workspace || !reference ||
      reference->formulation() != formulation() ||
      reference->purpose() != IntegrationPurpose::kStiffness ||
      reference->node_count() != node_count_ ||
      reference->point_count() != point_count_ ||
      context->workspace->node_count() != node_count_ ||
      context->workspace->point_count() != point_count_ ||
      state->orientations || state->committed_orientations || state->body_spin)
    return model_storage::Invalid("Incremental trial context", d);
  int failed_point = -1;
  const auto evaluate = formulation() == SolidFormulation::kHexa8Full ?
      EvaluateIncrementalHexa8Points : EvaluateIncrementalSolidPoints;
  auto status = evaluate(
      gradients_, weight_, context->committed_translations, state->translations,
      context->committed_stress, material(), context->workspace, request, force,
      nullptr, tangent, context->candidate_stress, failure, d,
      context->use_committed_geometric_stress, &failed_point);
  if (status != StatusCode::kOk && failed_point >= 0 && d) {
    const double* x = reference->natural_coordinates()->Data() +
                      3 * failed_point;
    d->geometry_point_available = true;
    d->natural_xi = x[0]; d->natural_eta = x[1]; d->natural_zeta = x[2];
  }
  // The common dispatcher retries only equilibrium failures with typed labels.
  // Committed singularities and backend faults retain their fatal status.
  if (status == StatusCode::kSingular &&
      *failure != TrialEvaluationFailure::kNone) {
    status = StatusCode::kEquilibriumFailure;
    if (d) d->code = status;
  }
  return status;
}
StatusCode IncrementalSolidElement::CheckTrialGeometry(
    const ElementTrialState* state, const ElementWorkspace* workspace,
    TrialEvaluationFailure* failure, Diagnostic* d) const {
  SolidElementReadView view;
  const auto status = SolidElementReadView::Bind(this, &view, d);
  if (status != StatusCode::kOk) return status;
  return CheckSolidTrialGeometry(&view, state, workspace, failure, d);
}
StatusCode IncrementalSolidElement::UpdateTrialResults(
    const ElementTrialState* state, const ElementWorkspace* workspace,
    TrialEvaluationFailure* failure, Diagnostic* d) {
  auto status =
      EvaluateTrial(state, workspace, ElementEvaluationRequest::kInternalForce,
                    force_, nullptr, failure, d);
  if (status != StatusCode::kOk) return status;
  const auto* context = state->incremental_solid;
  const Matrix* const* recovery_f = context->workspace->candidate_f_data();
  if (formulation() == SolidFormulation::kHexa8Full)
    recovery_f = context->workspace->hexa8_workspace()->modified_f_data();
  return RecoverIncrementalSolidPoints(
      workspace->reference(formulation(), IntegrationPurpose::kStiffness),
      recovery_f, context->candidate_stress,
      context->workspace, stress_, strain_, nodal_stress_, nodal_strain_,
      failure, d);
}
}  // namespace fesa
