#include "elements/solid_trial_geometry.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>

#include "fesa/math/dense_factorization.h"
#include "elements/solid_element_read_view.h"

namespace fesa {
namespace {
StatusCode NewMatrix(int rows, int columns, std::unique_ptr<Matrix>* owner,
                     Diagnostic* d) {
  Matrix* raw = nullptr;
  const auto status = Matrix::Create(rows, columns, &raw, d);
  owner->reset(raw);
  return status;
}
StatusCode NewVector(std::unique_ptr<Vector>* owner, Diagnostic* d) {
  Vector* raw = nullptr;
  const auto status = Vector::Create(3, &raw, d);
  owner->reset(raw);
  return status;
}
StatusCode Invalid(Diagnostic* d) {
  return ReportFailure(StatusCode::kInvalidArgument,
                       DiagnosticStage::kValidation, "Solid trial geometry", d);
}
StatusCode GeometryFailure(const SolidReferenceData* data, int point,
    bool measured, double value, TrialEvaluationFailure* failure,
    Diagnostic* d) {
  *failure = TrialEvaluationFailure::kGeometry;
  const auto status = ReportFailure(StatusCode::kEquilibriumFailure,
      DiagnosticStage::kComputation, "Solid current mapping guard", d);
  if (d) {
    const double* p = data->natural_coordinates()->Data() + 3 * point;
    d->geometry_point_available = true;
    d->natural_xi = p[0];
    d->natural_eta = p[1];
    d->natural_zeta = p[2];
    d->geometry_value_available = measured;
    d->geometry_quantity = GeometryQuantity::kNormalizedCurrentJacobian;
    d->geometry_value = value;
    d->geometry_limit = 256 * std::numeric_limits<double>::epsilon();
  }
  return status;
}
}  // namespace
StatusCode CheckSolidTrialGeometry(const SolidElement* element,
    const ElementTrialState* state, const ElementWorkspace* workspace,
    TrialEvaluationFailure* failure, Diagnostic* d) {
  if (failure) *failure = TrialEvaluationFailure::kNone;
  SolidElementReadView view;
  const auto status = SolidElementReadView::Bind(element, &view, d);
  if (status != StatusCode::kOk) return status;
  return CheckSolidTrialGeometry(&view, state, workspace, failure, d);
}
StatusCode CheckSolidTrialGeometry(const SolidElementReadView* element,
    const ElementTrialState* state, const ElementWorkspace* workspace,
    TrialEvaluationFailure* failure, Diagnostic* d) {
  if (failure) *failure = TrialEvaluationFailure::kNone;
  if (!element || !failure || !state || !state->translations || !workspace ||
      state->orientations || state->committed_orientations || state->body_spin ||
      (!element->element() ||
       state->translations->Size() != element->element()->dof_count())) return Invalid(d);
  const auto* data = workspace->reference(element->formulation(),
                                           IntegrationPurpose::kTrialGeometry);
  if (!data || data->formulation() != element->formulation() ||
      data->purpose() != IntegrationPurpose::kTrialGeometry ||
      data->node_count() != element->node_count()) return Invalid(d);
  for (int i = 0; i < state->translations->Size(); ++i)
    if (!std::isfinite(state->translations->Data()[i])) return Invalid(d);
  const int n = element->node_count();
  std::unique_ptr<Matrix> current, displacement, ones, anchor, shift, jacobian;
  std::unique_ptr<Vector> first, second;
  auto status = NewMatrix(n, 3, &current, d);
  if (status != StatusCode::kOk) return status;
  status = NewMatrix(n, 3, &displacement, d);
  if (status != StatusCode::kOk) return status;
  status = NewVector(&first, d);
  if (status != StatusCode::kOk) return status;
  status = NewVector(&second, d);
  if (status != StatusCode::kOk) return status;
  for (int node = 0; node < n; ++node) {
    for (int axis = 0; axis < 3; ++axis) {
      double x = 0;
      status = element->GetNode(node)->Coordinate(axis, &x, d);
      if (status != StatusCode::kOk) return status;
      status = current->Set(node, axis, x, d);
      if (status != StatusCode::kOk) return status;
    }
    const int indices[] = {3 * node, 3 * node + 1, 3 * node + 2};
    status = first->GatherFrom(state->translations, indices, 3, 0, d);
    if (status != StatusCode::kOk) return status;
    status = displacement->SetRow(node, first.get(), d);
    if (status != StatusCode::kOk) return status;
  }
  status = current->AddScaled(1, displacement.get(), d);
  if (status != StatusCode::kOk) return status;
  double h = 0;
  for (int i = 0; i < n; ++i) {
    const int a[] = {3 * i, 3 * i + 1, 3 * i + 2};
    for (int j = i + 1; j < n; ++j) {
      const int b[] = {3 * j, 3 * j + 1, 3 * j + 2};
      status = current->GatherTo(a, 3, first.get(), 0, d);
      if (status != StatusCode::kOk) return status;
      status = current->GatherTo(b, 3, second.get(), 0, d);
      if (status != StatusCode::kOk) return status;
      status = first->AddScaled(-1, second.get(), d);
      if (status != StatusCode::kOk) return status;
      double length = 0;
      status = first->Norm2(&length, d);
      if (status != StatusCode::kOk) return status;
      h = std::max(h, length);
    }
  }
  if (h == 0) return GeometryFailure(data, 0, false, 0, failure, d);
  if (!std::isfinite(h) || !std::isfinite(1 / h))
    return ReportFailure(StatusCode::kLibraryFailure,
        DiagnosticStage::kComputation, "Solid current geometry scale range", d);
  status = NewMatrix(n, 1, &ones, d);
  if (status != StatusCode::kOk) return status;
  status = NewMatrix(1, 3, &anchor, d);
  if (status != StatusCode::kOk) return status;
  status = NewMatrix(n, 3, &shift, d);
  if (status != StatusCode::kOk) return status;
  status = NewMatrix(3, 3, &jacobian, d);
  if (status != StatusCode::kOk) return status;
  for (int i = 0; i < n; ++i) {
    status = ones->Set(i, 0, 1, d);
    if (status != StatusCode::kOk) return status;
  }
  status = anchor->CopyBlockFrom(current.get(), 0, 0, 1, 3, 0, 0, d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Multiply(ones.get(), false, anchor.get(), false,
                             shift.get(), d);
  if (status != StatusCode::kOk) return status;
  status = current->AddScaled(-1, shift.get(), d);
  if (status != StatusCode::kOk) return status;
  status = current->Scale(1 / h, d);
  if (status != StatusCode::kOk) return status;
  for (int point = 0; point < data->point_count(); ++point) {
    status = Matrix::Multiply(current.get(), true,
        data->natural_gradients(point), true, jacobian.get(), d);
    if (status != StatusCode::kOk) return status;
    DenseLuFactorization* raw = nullptr;
    status = DenseLuFactorization::Create(jacobian.get(), &raw, d);
    std::unique_ptr<DenseLuFactorization> factor(raw);
    if (status == StatusCode::kSingular)
      return GeometryFailure(data, point, false, 0, failure, d);
    if (status != StatusCode::kOk) return status;
    double determinant = 0;
    status = factor->Determinant(&determinant, d);
    if (status != StatusCode::kOk) return status;
    if (!(determinant > 256 * std::numeric_limits<double>::epsilon()))
      return GeometryFailure(data, point, true, determinant, failure, d);
  }
  return StatusCode::kOk;
}
}  // namespace fesa
