#include "elements/incremental_solid_geometry.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>

#include "fesa/math/dense_factorization.h"
namespace fesa {
namespace {
StatusCode CoordinateVector(const Node* node, Vector* out, Diagnostic* d) {
  for (int axis = 0; axis < 3; ++axis) {
    double value = 0;
    auto status = node->Coordinate(axis, &value, d);
    if (status == StatusCode::kOk) status = out->Set(axis, value, d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}
StatusCode NewVector(std::unique_ptr<Vector>* out, Diagnostic* d) {
  Vector* raw = nullptr;
  const auto status = Vector::Create(3, &raw, d);
  out->reset(raw);
  return status;
}
StatusCode NewMatrix(int rows, int columns, std::unique_ptr<Matrix>* out,
                     Diagnostic* d) {
  Matrix* raw = nullptr;
  const auto status = Matrix::Create(rows, columns, &raw, d);
  out->reset(raw);
  return status;
}
StatusCode InvalidGeometry(Diagnostic* d) {
  return ReportFailure(StatusCode::kInputError, DiagnosticStage::kValidation,
                       "Incremental normalized initial Jacobian", d);
}
StatusCode CenterJacobianLength(const IncrementalSolidElement* element, double* out,
                                Diagnostic* d) {
  const bool wedge = element->formulation() == SolidFormulation::kWedge6;
  const int nodes = wedge ? 6 : 8;
  Matrix *coordinates = nullptr, *gradient = nullptr, *jacobian = nullptr;
  Vector *anchor = nullptr, *row = nullptr, *lengths = nullptr,
         *weights = nullptr;
  auto status = Matrix::Create(nodes, 3, &coordinates, d);
  std::unique_ptr<Matrix> co(coordinates);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Create(3, nodes, &gradient, d);
  std::unique_ptr<Matrix> go(gradient);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Create(3, 3, &jacobian, d);
  std::unique_ptr<Matrix> jo(jacobian);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(3, &anchor, d);
  std::unique_ptr<Vector> ao(anchor);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(3, &row, d);
  std::unique_ptr<Vector> ro(row);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(3, &lengths, d);
  std::unique_ptr<Vector> lo(lengths);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(3, &weights, d);
  std::unique_ptr<Vector> wo(weights);
  if (status != StatusCode::kOk) return status;
  const int signs[8][3] = {{-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1},
                           {-1, -1, 1},  {1, -1, 1},  {1, 1, 1},  {-1, 1, 1}};
  const int triangle_derivative[2][3] = {{-1, 1, 0}, {-1, 0, 1}};
  for (int axis = 0; axis < 3; ++axis) {
    double value = 0;
    status = element->GetNode(0)->Coordinate(axis, &value, d);
    if (status != StatusCode::kOk) return status;
    status = anchor->Set(axis, value, d);
    if (status != StatusCode::kOk) return status;
    status = weights->Set(axis, wedge && axis < 2 ? 1. / 3 : 2. / 3, d);
    if (status != StatusCode::kOk) return status;
  }
  for (int node = 0; node < nodes; ++node) {
    for (int axis = 0; axis < 3; ++axis) {
      double value = 0;
      status = element->GetNode(node)->Coordinate(axis, &value, d);
      if (status != StatusCode::kOk) return status;
      status = row->Set(axis, value, d);
      if (status != StatusCode::kOk) return status;
      // F4 at (1/3,1/3,0) or F3 at the Hexa center; only dN is needed.
      const double natural =
          wedge ? (axis < 2 ? triangle_derivative[axis][node % 3] / 2.
                            : (node < 3 ? -1. : 1.) / 6)
                : signs[node][axis] / 8.;
      status = gradient->Set(axis, node, natural, d);
      if (status != StatusCode::kOk) return status;
    }
    status = row->AddScaled(-1, anchor, d);
    if (status != StatusCode::kOk) return status;
    status = coordinates->SetRow(node, row, d);
    if (status != StatusCode::kOk) return status;
  }
  status = Matrix::Multiply(coordinates, true, gradient, true, jacobian, d);
  if (status != StatusCode::kOk) return status;
  for (int axis = 0; axis < 3; ++axis) {
    const int indices[] = {axis, 3 + axis, 6 + axis};
    status = jacobian->GatherTo(indices, 3, row, 0, d);
    if (status != StatusCode::kOk) return status;
    double norm = 0;
    status = row->Norm2(&norm, d);
    if (status != StatusCode::kOk) return status;
    status = lengths->Set(axis, norm, d);
    if (status != StatusCode::kOk) return status;
  }
  double length = 0;
  status = lengths->Dot(weights, &length, d);
  if (status != StatusCode::kOk) return status;
  if (!(length > 0) || !std::isfinite(length))
    return ReportFailure(StatusCode::kLibraryFailure,
                         DiagnosticStage::kComputation,
                         "Solid center-J characteristic length range", d);
  *out = length;
  return StatusCode::kOk;
}
// User-adopted empirical Pyramid length: four base samples at a,b=+/-0.51.
// Only bilinear base tangents and apex distances are needed; no apex gradient.
StatusCode PyramidLength(const IncrementalSolidElement* element, double* out,
                         Diagnostic* d) {
  Matrix *coordinates = nullptr, *coefficients = nullptr, *values = nullptr;
  Vector *anchor = nullptr, *apex = nullptr, *row = nullptr, *lengths = nullptr,
         *weights = nullptr;
  auto status = Matrix::Create(4, 3, &coordinates, d);
  std::unique_ptr<Matrix> co(coordinates);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Create(12, 4, &coefficients, d);
  std::unique_ptr<Matrix> cfo(coefficients);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Create(12, 3, &values, d);
  std::unique_ptr<Matrix> vo(values);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(3, &anchor, d);
  std::unique_ptr<Vector> ao(anchor);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(3, &apex, d);
  std::unique_ptr<Vector> apo(apex);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(3, &row, d);
  std::unique_ptr<Vector> ro(row);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(12, &lengths, d);
  std::unique_ptr<Vector> lo(lengths);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(12, &weights, d);
  std::unique_ptr<Vector> wo(weights);
  if (status != StatusCode::kOk) return status;
  for (int axis = 0; axis < 3; ++axis) {
    double value = 0;
    status = element->GetNode(0)->Coordinate(axis, &value, d);
    if (status != StatusCode::kOk) return status;
    status = anchor->Set(axis, value, d);
    if (status != StatusCode::kOk) return status;
    status = element->GetNode(4)->Coordinate(axis, &value, d);
    if (status != StatusCode::kOk) return status;
    status = apex->Set(axis, value, d);
    if (status != StatusCode::kOk) return status;
  }
  status = apex->AddScaled(-1, anchor, d);
  if (status != StatusCode::kOk) return status;
  // Base node zero is the translation anchor, hence its relative row is zero.
  for (int node = 1; node < 4; ++node) {
    for (int axis = 0; axis < 3; ++axis) {
      double value = 0;
      status = element->GetNode(node)->Coordinate(axis, &value, d);
      if (status != StatusCode::kOk) return status;
      status = row->Set(axis, value, d);
      if (status != StatusCode::kOk) return status;
    }
    status = row->AddScaled(-1, anchor, d);
    if (status != StatusCode::kOk) return status;
    status = coordinates->SetRow(node, row, d);
    if (status != StatusCode::kOk) return status;
  }
  const int signs[4][2] = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}};
  for (int point = 0; point < 4; ++point) {
    const double a = .51 * signs[point][0];
    const double b = .51 * signs[point][1];
    for (int node = 0; node < 4; ++node) {
      const double coefficients_at_point[] = {
          signs[node][0] * (1 + signs[node][1] * b) / 4,
          signs[node][1] * (1 + signs[node][0] * a) / 4,
          (1 + signs[node][0] * a) * (1 + signs[node][1] * b) / 4};
      for (int component = 0; component < 3; ++component) {
        status = coefficients->Set(3 * point + component, node,
                                   coefficients_at_point[component], d);
        if (status != StatusCode::kOk) return status;
      }
    }
  }
  status = Matrix::Multiply(coefficients, false, coordinates, false, values, d);
  if (status != StatusCode::kOk) return status;
  for (int sample = 0; sample < 12; ++sample) {
    const int indices[] = {3 * sample, 3 * sample + 1, 3 * sample + 2};
    status = values->GatherTo(indices, 3, row, 0, d);
    if (status != StatusCode::kOk) return status;
    if (sample % 3 == 2) {
      status = row->AddScaled(-1, apex, d);
      if (status != StatusCode::kOk) return status;
    }
    double norm = 0;
    status = row->Norm2(&norm, d);
    if (status != StatusCode::kOk) return status;
    status = lengths->Set(sample, norm, d);
    if (status != StatusCode::kOk) return status;
    status = weights->Set(sample, sample % 3 == 2 ? 1. / 12 : .273 / 4, d);
    if (status != StatusCode::kOk) return status;
  }
  double length = 0;
  status = lengths->Dot(weights, &length, d);
  if (status != StatusCode::kOk) return status;
  if (!(length > 0) || !std::isfinite(length))
    return ReportFailure(StatusCode::kLibraryFailure,
                         DiagnosticStage::kComputation,
                         "Pyramid characteristic length range", d);
  *out = length;
  return StatusCode::kOk;
}
StatusCode TetraLength(const IncrementalSolidElement* element, double* out,
                       Diagnostic* d) {
  if (!out) return ReportFailure(StatusCode::kInvalidArgument,
      DiagnosticStage::kValidation, "Incremental solid length", d);
  std::unique_ptr<Vector> anchor, edge, lengths, weights;
  auto status = NewVector(&anchor, d);
  if (status == StatusCode::kOk) status = NewVector(&edge, d);
  if (status == StatusCode::kOk) status = NewVector(&lengths, d);
  if (status == StatusCode::kOk) status = NewVector(&weights, d);
  if (status == StatusCode::kOk)
    status = CoordinateVector(element->GetNode(0), anchor.get(), d);
  if (status != StatusCode::kOk) return status;
  for (int node = 1; node < 4; ++node) {
    status = CoordinateVector(element->GetNode(node), edge.get(), d);
    if (status == StatusCode::kOk)
      status = edge->AddScaled(-1, anchor.get(), d);
    double length = 0;
    if (status == StatusCode::kOk) status = edge->Norm2(&length, d);
    if (status == StatusCode::kOk) status = lengths->Set(node - 1, length, d);
    if (status == StatusCode::kOk)
      status = weights->Set(node - 1, node == 3 ? 2. / 3 : 1. / 3, d);
    if (status != StatusCode::kOk) return status;
  }
  double value = 0;
  status = lengths->Dot(weights.get(), &value, d);
  if (status != StatusCode::kOk) return status;
  if (!(value > 0) || !std::isfinite(value))
    return ReportFailure(StatusCode::kLibraryFailure,
                         DiagnosticStage::kComputation,
                         "Incremental solid length range", d);
  *out = value;
  return StatusCode::kOk;
}
}  // namespace
// Relative coordinates and normalized J share one call lifetime. Matrix math
// stays in the common library; only scalar quadrature/guard coefficients are
// formed here. No constitutive state or legacy SolidElement is instantiated.
StatusCode BuildIncrementalInitialGeometry(const Element* element,
    const SolidReferenceData* stiffness, Matrix* const* gradients,
    Vector* weights, Diagnostic* d) {
  auto invalid = [d]() {
    return ReportFailure(StatusCode::kInvalidArgument,
        DiagnosticStage::kValidation, "Incremental initial geometry shape", d);
  };
  if (!element || !stiffness || !gradients || !weights ||
      stiffness->purpose() != IntegrationPurpose::kStiffness ||
      element->node_count() != stiffness->node_count() ||
      weights->Size() != stiffness->point_count()) return invalid();
  const int n = stiffness->node_count(), ng = stiffness->point_count();
  const auto type = stiffness->formulation();
  const bool tetra = type == SolidFormulation::kTetra4Full ||
                     type == SolidFormulation::kTetra10Full;
  for (int g = 0; g < ng; ++g) {
    if (!gradients[g] || gradients[g]->Rows() != 3 ||
        gradients[g]->Columns() != n) return invalid();
    for (int j = 0; j < g; ++j)
      if (gradients[g] == gradients[j]) return invalid();
  }
  for (int i = 0; i < n; ++i)
    if (!element->GetNode(i)) return invalid();
  std::unique_ptr<Vector> anchor, edge;
  auto status = NewVector(&anchor, d);
  if (status == StatusCode::kOk) status = NewVector(&edge, d);
  if (status != StatusCode::kOk) return status;
  double h = 0;
  for (int i = 0; i < n; ++i) {
    status = CoordinateVector(element->GetNode(i), anchor.get(), d);
    if (status != StatusCode::kOk) return status;
    for (int j = i + 1; j < n; ++j) {
      status = CoordinateVector(element->GetNode(j), edge.get(), d);
      if (status == StatusCode::kOk)
        status = edge->AddScaled(-1, anchor.get(), d);
      double length = 0;
      if (status == StatusCode::kOk) status = edge->Norm2(&length, d);
      if (status != StatusCode::kOk) return status;
      h = std::max(h, length);
    }
  }
  if (!(h > 0)) return InvalidGeometry(d);
  if (!std::isfinite(1 / h))
    return ReportFailure(StatusCode::kLibraryFailure,
        DiagnosticStage::kComputation, "Incremental initial scale range", d);
  std::unique_ptr<Matrix> coordinates, jacobian, corners, natural;
  status = NewMatrix(n, 3, &coordinates, d);
  if (status == StatusCode::kOk) status = NewMatrix(3, 3, &jacobian, d);
  if (status == StatusCode::kOk)
    status = CoordinateVector(element->GetNode(0), anchor.get(), d);
  if (status != StatusCode::kOk) return status;
  for (int node = 0; node < n; ++node) {
    status = CoordinateVector(element->GetNode(node), edge.get(), d);
    if (status == StatusCode::kOk)
      status = edge->AddScaled(-1, anchor.get(), d);
    if (status == StatusCode::kOk)
      status = coordinates->SetRow(node, edge.get(), d);
    if (status != StatusCode::kOk) return status;
  }
  status = coordinates->Scale(1 / h, d);
  if (status != StatusCode::kOk) return status;
  if (tetra) {
    status = NewMatrix(4, 3, &corners, d);
    if (status == StatusCode::kOk) status = NewMatrix(3, 4, &natural, d);
    if (status == StatusCode::kOk)
      status = corners->CopyBlockFrom(coordinates.get(), 0, 0, 4, 3, 0, 0, d);
    for (int axis = 0; axis < 3 && status == StatusCode::kOk; ++axis) {
      status = natural->Set(axis, 0, -1, d);
      if (status == StatusCode::kOk) status = natural->Set(axis, axis+1, 1, d);
    }
    if (status != StatusCode::kOk) return status;
  }
  // Initial Tetra10 midsides are checked without projecting the input. Mapping
  // uses the four corner affine J, while all ten displacement shapes remain.
  if (type == SolidFormulation::kTetra10Full) {
    const int edges[6][2] = {{0,1},{1,2},{2,0},{0,3},{1,3},{2,3}};
    std::unique_ptr<Vector> midpoint;
    status = NewVector(&midpoint, d);
    for (int k = 0; k < 6 && status == StatusCode::kOk; ++k) {
      const int a = 3 * edges[k][0], b = 3 * edges[k][1], m = 3 * (4+k);
      const int ai[] = {a,a+1,a+2}, bi[] = {b,b+1,b+2};
      const int mi[] = {m,m+1,m+2};
      status = coordinates->GatherTo(ai, 3, anchor.get(), 0, d);
      if (status == StatusCode::kOk)
        status = coordinates->GatherTo(bi, 3, edge.get(), 0, d);
      if (status == StatusCode::kOk)
        status = coordinates->GatherTo(mi, 3, midpoint.get(), 0, d);
      if (status == StatusCode::kOk) status = anchor->Scale(.5, d);
      if (status == StatusCode::kOk)
        status = anchor->AddScaled(.5, edge.get(), d);
      if (status == StatusCode::kOk)
        status = midpoint->AddScaled(-1, anchor.get(), d);
      double deviation = 0;
      if (status == StatusCode::kOk) status = midpoint->Norm2(&deviation, d);
      if (status != StatusCode::kOk) return status;
      if (deviation > 1e-10) {
        if (d) d->node_id = element->GetNode(4+k)->Id();
        return ReportFailure(StatusCode::kInputError,
            DiagnosticStage::kValidation, "Incremental curved Tetra10", d);
      }
    }
    if (status != StatusCode::kOk) return status;
  }
  auto point = [&](const Matrix* derivative, Matrix* gradient,
                   double* determinant) {
    auto result = Matrix::Multiply(tetra ? corners.get() : coordinates.get(),
        true, tetra ? natural.get() : derivative, true, jacobian.get(), d);
    DenseLuFactorization* raw = nullptr;
    if (result == StatusCode::kOk)
      result = DenseLuFactorization::Create(jacobian.get(), &raw, d);
    std::unique_ptr<DenseLuFactorization> factor(raw);
    if (result == StatusCode::kSingular) return InvalidGeometry(d);
    if (result == StatusCode::kOk)
      result = factor->Determinant(determinant, d);
    if (result != StatusCode::kOk) return result;
    if (!(*determinant > 256 * std::numeric_limits<double>::epsilon()))
      return InvalidGeometry(d);
    if (gradient) {
      result = factor->Solve(derivative, true, gradient, d);
      if (result == StatusCode::kOk) result = gradient->Scale(1 / h, d);
    }
    return result;
  };
  // Initial guards use the same immutable ordered non-apex sample set as the
  // current geometry contract; initial failures retain input-error semantics.
  SolidReferenceData* raw = nullptr;
  status = SolidReferenceData::Create(type, IntegrationPurpose::kTrialGeometry,
      &raw, d);
  std::unique_ptr<SolidReferenceData> samples(raw);
  if (status != StatusCode::kOk) return status;
  auto location = [d](const SolidReferenceData* data, int g) {
    if (!d) return;
    const double* x = data->natural_coordinates()->Data() + 3*g;
    d->geometry_point_available = true;
    d->natural_xi = x[0]; d->natural_eta = x[1]; d->natural_zeta = x[2];
  };
  if (!tetra)
    for (int g = 0; g < samples->point_count(); ++g) {
      double determinant = 0;
      status = point(samples->natural_gradients(g), nullptr, &determinant);
      if (status != StatusCode::kOk) { location(samples.get(), g); return status; }
    }
  for (int g = 0; g < ng; ++g) {
    double determinant = 0;
    status = point(stiffness->natural_gradients(g), gradients[g], &determinant);
    if (status != StatusCode::kOk) { location(stiffness, g); return status; }
    const double q = type == SolidFormulation::kTetra4Full
        ? determinant * h * h * h / 6
        : determinant * h * h * h * stiffness->weights()->Data()[g];
    if (!(q > 0) || !std::isfinite(q))
      return ReportFailure(StatusCode::kLibraryFailure,
          DiagnosticStage::kComputation, "Incremental initial volume range", d);
    status = weights->Set(g, q, d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}
StatusCode ComputeIncrementalCharacteristicLength(
    const IncrementalSolidElement* element, double* out, Diagnostic* d) {
  if (!element || !out)
    return ReportFailure(StatusCode::kInvalidArgument,
        DiagnosticStage::kValidation, "Incremental length input", d);
  const auto type = element->formulation();
  if (type == SolidFormulation::kTetra4Full ||
      type == SolidFormulation::kTetra10Full) return TetraLength(element, out, d);
  if (type == SolidFormulation::kPyramid5Full)
    return PyramidLength(element, out, d);
  return CenterJacobianLength(element, out, d);
}
}  // namespace fesa
