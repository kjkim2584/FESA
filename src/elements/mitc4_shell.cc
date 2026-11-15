#include "fesa/elements/mitc4_shell.h"

#include <cmath>
#include <memory>
#include <new>

#include "elements/mitc4_initial.h"
#include "elements/mitc4_trial.h"
#include "model/model_storage.h"

namespace fesa {

Mitc4Shell::Mitc4Shell(int id, const ShellSection* section)
    : Element(id), section_(section), material_(section->material()) {}

Mitc4Shell::~Mitc4Shell() {
  delete results_;
  delete initial_axes_;
}

StatusCode Mitc4Shell::Create(int id, const Node* const* nodes, int count,
                              const ShellSection* section, Mitc4Shell** out,
                              Diagnostic* d) {
  if (!out) return model_storage::Invalid("Mitc4Shell::Create output", d);
  *out = nullptr;
  if (!model_storage::ValidId(id) || !nodes || count != 4 || !section)
    return model_storage::Invalid("Mitc4Shell::Create definition", d);
  std::unique_ptr<Mitc4Shell> shell(new (std::nothrow) Mitc4Shell(id, section));
  if (!shell) return model_storage::AllocationFailure("Mitc4Shell::Create", d);
  auto status = Matrix::Create(3, 3, &shell->initial_axes_, d);
  if (status != StatusCode::kOk) return status;
  // The initial derivation exports model-output axes and kd together; K is
  // discarded inside it. Numerical evaluations still reconstruct their axes.
  status = ComputeMitc4InitialProperties(
      nodes, count, section, shell->initial_axes_, &shell->drilling_, d);
  if (status != StatusCode::kOk) return status;
  status = Mitc4Results::Create(&shell->results_, d);
  if (status != StatusCode::kOk) return status;
  for (int i = 0; i < 4; ++i) shell->nodes_[i] = nodes[i];
  *out = shell.release();
  return StatusCode::kOk;
}

ElementFamily Mitc4Shell::family() const { return ElementFamily::kShell; }
const char* Mitc4Shell::type_name() const { return "S4"; }
int Mitc4Shell::node_count() const { return 4; }
int Mitc4Shell::dofs_per_node() const { return 6; }
int Mitc4Shell::dof_count() const { return 24; }

const Node* Mitc4Shell::GetNode(int local) const {
  return local >= 0 && local < 4 ? nodes_[local] : nullptr;
}

StatusCode Mitc4Shell::CharacteristicLength(double* out, Diagnostic* d) const {
  if (!out)
    return model_storage::Invalid("Mitc4Shell characteristic length output", d);
  Vector* raw = nullptr;
  auto status = Vector::Create(3, &raw, d);
  std::unique_ptr<Vector> first(raw);
  if (status != StatusCode::kOk) return status;
  raw = nullptr;
  status = Vector::Create(3, &raw, d);
  std::unique_ptr<Vector> second(raw);
  if (status != StatusCode::kOk) return status;
  double mean = 0.0;
  // Lc is the mean of the four initial edges, not the geometry guard diameter.
  // Factory-validated immutable nodes need no full geometry or K evaluation.
  for (int i = 0; i < 4; ++i) {
    for (int axis = 0; axis < 3; ++axis) {
      double coordinate = 0.0;
      status = nodes_[i]->Coordinate(axis, &coordinate, d);
      if (status != StatusCode::kOk) return status;
      status = first->Set(axis, coordinate, d);
      if (status != StatusCode::kOk) return status;
      status = nodes_[(i + 1) % 4]->Coordinate(axis, &coordinate, d);
      if (status != StatusCode::kOk) return status;
      status = second->Set(axis, coordinate, d);
      if (status != StatusCode::kOk) return status;
    }
    status = first->AddScaled(-1.0, second.get(), d);
    if (status != StatusCode::kOk) return status;
    double length = 0.0;
    status = first->Norm2(&length, d);
    if (status != StatusCode::kOk) return status;
    mean += length / 4.0;
  }
  if (!std::isfinite(mean) || mean <= 0.0)
    return ReportFailure(StatusCode::kLibraryFailure,
                         DiagnosticStage::kComputation,
                         "Mitc4Shell characteristic length range", d);
  *out = mean;
  return StatusCode::kOk;
}

StatusCode Mitc4Shell::ComputeStiffness(ElementWorkspace*, Matrix* out,
                                        Diagnostic* d) const {
  double drilling = 0.0;
  return ComputeMitc4InitialStiffness(nodes_, 4, section_, out, &drilling, d);
}

StatusCode Mitc4Shell::ComputeInternalForce(const Vector* u, ElementWorkspace*,
                                            Vector* out, Diagnostic* d) const {
  return ComputeMitc4LinearInternalForce(nodes_, 4, section_, u, drilling_, out,
                                         d);
}

bool Mitc4Shell::SupportsTrialEvaluation() const { return true; }
StatusCode Mitc4Shell::EvaluateTrial(const ElementTrialState* state,
                                     const ElementWorkspace*,
                                     ElementEvaluationRequest request,
                                     Vector* force, Matrix* tangent,
                                     TrialEvaluationFailure* failure,
                                     Diagnostic* d) const {
  return EvaluateMitc4Trial(nodes_, 4, section_, state, request, drilling_,
                            force, tangent, failure, d);
}

bool Mitc4Shell::SupportsMass(MassMode mode) const {
  return mode == MassMode::kConsistent || mode == MassMode::kRowSum;
}

StatusCode Mitc4Shell::ComputeMass(MassMode mode, ElementWorkspace*,
                                   Matrix* out, Diagnostic* d) const {
  if (!SupportsMass(mode))
    return ReportFailure(StatusCode::kUnsupported,
                         DiagnosticStage::kComputation,
                         "Mitc4Shell mass mode unsupported", d);
  if (mode == MassMode::kConsistent)
    return ComputeMitc4InitialMass(nodes_, 4, section_, out, d);
  if (!out || out->Rows() != 24 || out->Columns() != 24)
    return model_storage::Invalid("Mitc4Shell mass output/shape", d);
  Matrix* raw = nullptr;
  auto status = Matrix::Create(24, 24, &raw, d);
  std::unique_ptr<Matrix> consistent(raw);
  if (status != StatusCode::kOk) return status;
  status = ComputeMitc4InitialMass(nodes_, 4, section_, consistent.get(), d);
  if (status != StatusCode::kOk) return status;
  return LumpMitc4InitialMass(nodes_, 4, consistent.get(), out, d);
}

StatusCode Mitc4Shell::LumpFromConsistentMass(const Matrix* consistent,
                                              Matrix* lumped,
                                              Diagnostic* d) const {
  return LumpMitc4InitialMass(nodes_, 4, consistent, lumped, d);
}

StatusCode Mitc4Shell::UpdateResults(const Vector* u, ElementWorkspace*,
                                     Diagnostic* d) {
  return results_->RecoverLinear(nodes_, 4, section_, u, drilling_, d);
}
StatusCode Mitc4Shell::UpdateTrialResults(const ElementTrialState* state,
                                          const ElementWorkspace*,
                                          TrialEvaluationFailure* failure,
                                          Diagnostic* d) {
  const auto status = results_->RecoverNonlinear(nodes_, 4, section_, state,
                                                 drilling_, failure, d);
  if (status != StatusCode::kOk && d) d->element_id = Id();
  return status;
}
const Matrix* Mitc4Shell::initial_output_axes() const { return initial_axes_; }
const Mitc4Results* Mitc4Shell::results() const { return results_; }

const ShellSection* Mitc4Shell::section() const { return section_; }
const ElasticMaterial* Mitc4Shell::material() const { return material_; }
double Mitc4Shell::drilling_stiffness() const { return drilling_; }

}  // namespace fesa
