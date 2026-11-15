#include "analysis/linear_shell_nullspace.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <memory>
#include <new>
#include <vector>

#include "base/diagnostic_copy.h"
#include "fesa/elements/mitc4_shell.h"
#include "fesa/math/dense_factorization.h"
#include "fesa/math/rotation.h"
#include "model/model_storage.h"

namespace fesa {
namespace {
StatusCode Frame(const Mitc4Shell* shell, Matrix* frame, Diagnostic* d) {
  Diagnostic probe;
  CopyDiagnosticMetadata(d, &probe);
  probe.element_id = shell->Id();
  const auto status = shell->ComputeInitialFrame(frame, &probe);
  if (status != StatusCode::kOk) CopyDiagnosticMetadata(&probe, d);
  return status;
}
StatusCode Independent(const Matrix* candidate, bool* independent,
                       Diagnostic* d) {
  *independent = false;
  Diagnostic probe;
  CopyDiagnosticMetadata(d, &probe);
  DenseQrFactorization* raw = nullptr;
  const auto status = DenseQrFactorization::Create(candidate, &raw, &probe);
  std::unique_ptr<DenseQrFactorization> owner(raw);
  if (status == StatusCode::kSingular) return StatusCode::kOk;
  if (status != StatusCode::kOk) {
    CopyDiagnosticMetadata(&probe, d);
    return status;
  }
  *independent = true;
  return StatusCode::kOk;
}
}  // namespace
LinearShellNullspace::~LinearShellNullspace() {
  delete[] nodes_;
  delete[] global_;
  delete[] equations_;
  delete directions_;
  delete reference_directions_;
  delete coefficients_;
}
StatusCode LinearShellNullspace::Create(const AnalysisModel* model,
                                        const DofManager* dofs,
                                        LinearShellNullspace** out,
                                        Diagnostic* d) {
  if (!out) return model_storage::Invalid("Shell nullspace output", d);
  *out = nullptr;
  try {
    std::unique_ptr<LinearShellNullspace> result(new LinearShellNullspace);
    const auto status = result->Initialize(model, dofs, false, d);
    if (status != StatusCode::kOk) return status;
    *out = result.release();
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return model_storage::AllocationFailure("Linear shell nullspace", d);
  }
}

StatusCode LinearShellNullspace::CreateFinite(const AnalysisModel* model,
                                              const DofManager* dofs,
                                              const Matrix* orientations,
                                              LinearShellNullspace** out,
                                              Diagnostic* d) {
  if (out) *out = nullptr;
  if (!out || !orientations)
    return model_storage::Invalid("Finite shell nullspace arguments", d);
  try {
    std::unique_ptr<LinearShellNullspace> result(new LinearShellNullspace);
    auto status = result->Initialize(model, dofs, true, d);
    if (status == StatusCode::kOk) status = result->Rotate(orientations, d);
    if (status != StatusCode::kOk) return status;
    *out = result.release();
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return model_storage::AllocationFailure("Finite shell nullspace", d);
  }
}
StatusCode LinearShellNullspace::Rotate(const Matrix* orientations,
                                        Diagnostic* d) {
  const int nodes = global_count_ / 6;
  if (orientations->Rows() != nodes || orientations->Columns() != 4)
    return model_storage::Invalid("Finite shell nullspace pose shape", d);
  auto status = Matrix::Create(nodes, 3, &reference_directions_, d);
  if (status == StatusCode::kOk)
    status = reference_directions_->CopyFrom(directions_, d);
  if (status != StatusCode::kOk) return status;
  Matrix* raw_rotation = nullptr;
  status = Matrix::Create(3, 3, &raw_rotation, d);
  std::unique_ptr<Matrix> rotation(raw_rotation);
  if (status != StatusCode::kOk) return status;
  Vector *raw_q = nullptr, *raw_initial = nullptr, *raw_current = nullptr;
  status = Vector::Create(4, &raw_q, d);
  std::unique_ptr<Vector> q(raw_q);
  if (status == StatusCode::kOk) status = Vector::Create(3, &raw_initial, d);
  std::unique_ptr<Vector> initial(raw_initial);
  if (status == StatusCode::kOk) status = Vector::Create(3, &raw_current, d);
  std::unique_ptr<Vector> current(raw_current);
  if (status != StatusCode::kOk) return status;
  // Validate the complete input pose, including fully constrained/unused rows.
  for (int node = 0; node < nodes; ++node) {
    const int indices[] = {4 * node, 4 * node + 1, 4 * node + 2, 4 * node + 3};
    status = orientations->GatherTo(indices, 4, q.get(), 0, d);
    if (status == StatusCode::kOk)
      status = RotationMatrix(q.get(), rotation.get(), d);
    if (status != StatusCode::kOk) return status;
  }
  for (int i = 0; i < count_; ++i) {
    const int node = global_[3 * i] / 6;
    const int q_indices[] = {4 * node, 4 * node + 1, 4 * node + 2,
                             4 * node + 3};
    const int indices[] = {3 * i, 3 * i + 1, 3 * i + 2};
    status = orientations->GatherTo(q_indices, 4, q.get(), 0, d);
    if (status == StatusCode::kOk)
      status = RotationMatrix(q.get(), rotation.get(), d);
    if (status == StatusCode::kOk)
      status = reference_directions_->GatherTo(indices, 3, initial.get(), 0, d);
    if (status == StatusCode::kOk)
      status = Matrix::MultiplyVector(rotation.get(), false, initial.get(),
                                      current.get(), d);
    double norm = 0;
    if (status == StatusCode::kOk) status = current->Norm2(&norm, d);
    if (status == StatusCode::kOk) status = current->Scale(1 / norm, d);
    if (status == StatusCode::kOk)
      status = directions_->SetRow(i, current.get(), d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}
StatusCode LinearShellNullspace::Initialize(const AnalysisModel* model,
                                            const DofManager* dofs, bool finite,
                                            Diagnostic* d) {
  if (!model || !dofs || !model->Conditions() || model->Nodes()->empty() ||
      model->Elements()->empty() || model->Nodes()->size() > INT_MAX / 6 ||
      model->Nodes()->size() != static_cast<size_t>(dofs->NodeCount()) ||
      dofs->GlobalCount() != 6 * dofs->NodeCount())
    return model_storage::Invalid("Linear shell nullspace model/mapping", d);
  const StepDefinition* step = nullptr;
  auto status = model->Conditions()->SourceDomain()->FindStepDefinition(
      model->Conditions()->StepNumber(), &step, d);
  if (status != StatusCode::kOk) return status;
  if (!finite && step->HasFiniteGeometry())
    return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                         "Initial linear shell nullspace only", d);
  for (const auto& entry : *model->Nodes())
    if (dofs->NodeDofCount(entry.first) != 6)
      return model_storage::Invalid("Shell nullspace node mapping", d);
  if (finite) {
    for (int row = 0; row < dofs->NodeCount(); ++row) {
      int free = 0;
      for (int c = 4; c <= 6; ++c)
        if (dofs->Equation(dofs->GlobalDof(dofs->NodeId(row), c)) >= 0) ++free;
      if (free != 0 && free != 3)
        return ReportFailure(
            StatusCode::kUnsupported, DiagnosticStage::kValidation,
            "Finite shell rotations require all3 components", d);
    }
  }

  const int nodes = dofs->NodeCount();
  global_count_ = dofs->GlobalCount();
  nodes_ = new int[nodes];
  global_ = new int[3 * nodes];
  equations_ = new int[3 * nodes];
  status = Matrix::Create(nodes, 3, &directions_, d);
  if (status == StatusCode::kOk)
    status = Vector::Create(nodes, &coefficients_, d);
  if (status != StatusCode::kOk) return status;
  // Explicit size constructors avoid throwing debug-proxy allocations inside
  // MSVC's noexcept default vector constructor. All adjacency is call-local.
  using Incidents = std::vector<const Mitc4Shell*>;
  std::vector<std::unique_ptr<Incidents>> incidence(nodes);
  for (const auto& entry : *model->Elements()) {
    const auto* shell = dynamic_cast<const Mitc4Shell*>(entry.second);
    if (!shell || !shell->SupportsDynamicEvaluation(
                      finite ? DynamicKinematics::kFiniteRotation
                             : DynamicKinematics::kLinear)) {
      if (d) d->element_id = entry.first;
      return ReportFailure(StatusCode::kUnsupported,
                           DiagnosticStage::kValidation,
                           "Nullspace requires supported S4 kinematics", d);
    }
    double density = 0;
    status = shell->material()->Density(&density, d);
    if (status != StatusCode::kOk) return status;
    for (int i = 0; i < 4; ++i) {
      const int global = dofs->GlobalDof(shell->GetNode(i)->Id(), 1);
      if (global < 0 || !dofs->Participates(shell->GetNode(i)->Id()))
        return model_storage::Invalid("Shell incidence mapping", d);
      auto* owner = &incidence[global / 6];
      if (!*owner) owner->reset(new Incidents(0));
      (*owner)->push_back(shell);
    }
  }
  std::unique_ptr<Matrix> matrices[3];
  for (auto* owner : {&matrices[0], &matrices[1], &matrices[2]}) {
    Matrix* raw = nullptr;
    status = Matrix::Create(3, 3, &raw, d);
    owner->reset(raw);
    if (status != StatusCode::kOk) return status;
  }
  auto* first = matrices[0].get();
  auto* other = matrices[1].get();
  auto* candidate = matrices[2].get();
  Vector* raw = nullptr;
  status = Vector::Create(3, &raw, d);
  std::unique_ptr<Vector> direction(raw);
  if (status != StatusCode::kOk) return status;
  for (int row = 0; row < nodes; ++row) {
    const auto* incident = incidence[row].get();
    if (!incident) continue;
    const int node = dofs->NodeId(row);
    status = Frame(incident->front(), first, d);
    if (status == StatusCode::kOk)
      status = candidate->CopyBlockFrom(first, 0, 0, 2, 3, 0, 0, d);
    if (status != StatusCode::kOk) return status;
    bool full = false;
    for (size_t i = 1; i < incident->size() && !full; ++i) {
      status = Frame((*incident)[i], other, d);
      if (status != StatusCode::kOk) return status;
      for (int tangent = 0; tangent < 2 && !full; ++tangent) {
        status = candidate->CopyBlockFrom(other, tangent, 0, 1, 3, 2, 0, d);
        if (status == StatusCode::kOk)
          status = Independent(candidate, &full, d);
        if (status != StatusCode::kOk) return status;
      }
    }
    for (int component = 0; component < 3 && !full; ++component) {
      if (dofs->Equation(dofs->GlobalDof(node, component + 4)) >= 0) continue;
      status = direction->Scale(0, d);
      if (status == StatusCode::kOk) status = direction->Set(component, 1, d);
      if (status == StatusCode::kOk)
        status = candidate->SetRow(2, direction.get(), d);
      if (status == StatusCode::kOk) status = Independent(candidate, &full, d);
      if (status != StatusCode::kOk) return status;
    }
    if (full) continue;
    const int normal_indices[] = {6, 7, 8};
    status = first->GatherTo(normal_indices, 3, direction.get(), 0, d);
    if (status != StatusCode::kOk) return status;
    for (int component = 0; component < 3; ++component)
      if (dofs->Equation(dofs->GlobalDof(node, component + 4)) < 0) {
        status = direction->Set(component, 0, d);
        if (status != StatusCode::kOk) return status;
      }
    double norm = 0;
    status = direction->Norm2(&norm, d);
    if (status != StatusCode::kOk) return status;
    if (norm == 0)
      return model_storage::Invalid("Empty shell null direction", d);
    status = direction->Scale(1 / norm, d);
    if (status == StatusCode::kOk)
      status = directions_->SetRow(count_, direction.get(), d);
    if (status != StatusCode::kOk) return status;
    Vector* raw_coefficients = nullptr;
    status = Vector::Create(static_cast<int>(incident->size()),
                            &raw_coefficients, d);
    std::unique_ptr<Vector> coefficients(raw_coefficients);
    if (status != StatusCode::kOk) return status;
    for (size_t i = 0; i < incident->size(); ++i) {
      status = coefficients->Set(static_cast<int>(i),
                                 (*incident)[i]->drilling_stiffness(), d);
      if (status != StatusCode::kOk) return status;
    }
    double sum = 0;
    status = coefficients->AbsoluteSum(&sum, d);
    if (status == StatusCode::kOk) status = coefficients_->Set(count_, sum, d);
    if (status != StatusCode::kOk) return status;
    nodes_[count_] = node;
    for (int component = 0; component < 3; ++component) {
      global_[3 * count_ + component] = dofs->GlobalDof(node, component + 4);
      equations_[3 * count_ + component] =
          dofs->Equation(global_[3 * count_ + component]);
    }
    ++count_;
  }
  return StatusCode::kOk;
}
int LinearShellNullspace::Count() const { return count_; }
int LinearShellNullspace::NodeId(int index) const {
  return index >= 0 && index < count_ ? nodes_[index] : -1;
}
int LinearShellNullspace::GlobalDof(int index, int component) const {
  return index >= 0 && index < count_ && component >= 0 && component < 3
             ? global_[3 * index + component]
             : -1;
}
int LinearShellNullspace::Equation(int index, int component) const {
  return index >= 0 && index < count_ && component >= 0 && component < 3
             ? equations_[3 * index + component]
             : -1;
}
const Matrix* LinearShellNullspace::Directions() const { return directions_; }
const Vector* LinearShellNullspace::DrillingCoefficients() const {
  return coefficients_;
}
const Matrix* LinearShellNullspace::ReferenceDirections() const {
  return reference_directions_ ? reference_directions_ : directions_;
}
StatusCode LinearShellNullspace::ProjectRotations(const Vector* global,
                                                  Vector* out,
                                                  Diagnostic* d) const {
  if (!global || global->Size() != global_count_ || global == out ||
      (count_ == 0 ? out != nullptr : !out || out->Size() != count_))
    return model_storage::Invalid("Shell nullspace projection shape", d);
  if (count_ == 0) return StatusCode::kOk;
  Vector *raw_local = nullptr, *raw_direction = nullptr;
  auto status = Vector::Create(3, &raw_local, d);
  std::unique_ptr<Vector> local(raw_local);
  if (status == StatusCode::kOk) status = Vector::Create(3, &raw_direction, d);
  std::unique_ptr<Vector> direction(raw_direction);
  if (status != StatusCode::kOk) return status;
  for (int i = 0; i < count_; ++i) {
    const int indices[] = {3 * i, 3 * i + 1, 3 * i + 2};
    status = local->GatherFrom(global, global_ + 3 * i, 3, 0, d);
    if (status == StatusCode::kOk)
      status = directions_->GatherTo(indices, 3, direction.get(), 0, d);
    double value = 0;
    if (status == StatusCode::kOk)
      status = local->Dot(direction.get(), &value, d);
    if (status == StatusCode::kOk) status = out->Set(i, value, d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}
StatusCode LinearShellNullspace::EvaluateCompatibility(const Vector* motion,
                                                       const Vector* moment,
                                                       Vector* out,
                                                       Diagnostic* d) const {
  if (reference_directions_)
    return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                         "Finite shell compatibility requires P17", d);
  if (!motion || !moment || motion->Size() != global_count_ ||
      moment->Size() != global_count_ || motion == out || moment == out ||
      (count_ == 0 ? out != nullptr : !out || out->Size() != count_))
    return model_storage::Invalid("Shell initial compatibility shape", d);
  if (count_ == 0) return StatusCode::kOk;
  Vector *raw_x = nullptr, *raw_p = nullptr, *raw_n = nullptr;
  auto status = Vector::Create(3, &raw_x, d);
  std::unique_ptr<Vector> x(raw_x);
  if (status == StatusCode::kOk) status = Vector::Create(3, &raw_p, d);
  std::unique_ptr<Vector> p(raw_p);
  if (status == StatusCode::kOk) status = Vector::Create(3, &raw_n, d);
  std::unique_ptr<Vector> n(raw_n);
  if (status != StatusCode::kOk) return status;
  for (int i = 0; i < count_; ++i) {
    status = x->GatherFrom(motion, global_ + 3 * i, 3, 0, d);
    if (status == StatusCode::kOk)
      status = p->GatherFrom(moment, global_ + 3 * i, 3, 0, d);
    const int indices[] = {3 * i, 3 * i + 1, 3 * i + 2};
    if (status == StatusCode::kOk)
      status = directions_->GatherTo(indices, 3, n.get(), 0, d);
    double maxima[2] = {}, norms[2] = {}, projections[2] = {};
    int exponents[2] = {};
    double mantissas[2] = {};
    Vector* fields[] = {x.get(), p.get()};
    for (int field = 0; field < 2 && status == StatusCode::kOk; ++field) {
      int index = 0;
      status = fields[field]->MaximumAbs(&maxima[field], &index, d);
      if (status != StatusCode::kOk) return status;
      if (maxima[field] == 0) continue;
      mantissas[field] = std::frexp(maxima[field], &exponents[field]);
      // Normalize by max without forming its possibly overflowing reciprocal.
      int remaining = -exponents[field];
      while (remaining != 0 && status == StatusCode::kOk) {
        const int shift = std::max(-512, std::min(512, remaining));
        status = fields[field]->Scale(std::scalbn(1.0, shift), d);
        remaining -= shift;
      }
      if (status == StatusCode::kOk)
        status = fields[field]->Scale(1 / mantissas[field], d);
      if (status == StatusCode::kOk)
        status = fields[field]->Norm2(&norms[field], d);
      if (status == StatusCode::kOk)
        status = fields[field]->Dot(n.get(), &projections[field], d);
    }
    if (status != StatusCode::kOk) return status;
    int spring_exponent = 0;
    const double spring =
        std::frexp(coefficients_->Data()[i], &spring_exponent);
    exponents[0] += spring_exponent;
    mantissas[0] *= spring;
    const int common = maxima[0] == 0   ? exponents[1]
                       : maxima[1] == 0 ? exponents[0]
                                        : std::max(exponents[0], exponents[1]);
    double weights[2] = {};
    for (int field = 0; field < 2; ++field) {
      if (maxima[field] == 0) continue;
      weights[field] = std::scalbn(mantissas[field], exponents[field] - common);
      if (weights[field] == 0)
        return model_storage::Invalid("Shell compatibility scale underflow", d);
    }
    const double denominator = weights[0] * norms[0] + weights[1] * norms[1];
    const double error = denominator == 0
                             ? 0
                             : std::abs(weights[0] * projections[0] -
                                        weights[1] * projections[1]) /
                                   denominator;
    status = out->Set(i, error, d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}
StatusCode LinearShellNullspace::CheckCompatibility(const Vector* motion,
                                                    const Vector* moment,
                                                    bool velocity,
                                                    Diagnostic* d) const {
  Vector* raw = nullptr;
  auto status = count_ == 0 ? StatusCode::kOk : Vector::Create(count_, &raw, d);
  std::unique_ptr<Vector> errors(raw);
  if (status == StatusCode::kOk)
    status = EvaluateCompatibility(motion, moment, raw, d);
  if (status != StatusCode::kOk) return status;
  for (int i = 0; i < count_; ++i) {
    if (raw->Data()[i] <= 1e-10) continue;
    status = ReportFailure(StatusCode::kEquilibriumFailure,
                           DiagnosticStage::kComputation,
                           velocity ? "Shell initial velocity compatibility"
                                    : "Shell initial rotation compatibility",
                           d);
    if (d) d->node_id = nodes_[i];
    return status;
  }
  return StatusCode::kOk;
}
}  // namespace fesa
