#include "elements/thermal_face_geometry.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>

#include "fesa/math/dense_factorization.h"
#include "fesa/model/thermal_element.h"
#include "model/model_storage.h"

namespace fesa {
namespace thermal_face_geometry {
namespace {
StatusCode GeometryError(Diagnostic* d) {
  return ReportFailure(StatusCode::kInputError, DiagnosticStage::kValidation,
                       "Thermal face degenerate or nonpositive Gauss geometry",
                       d);
}
StatusCode Row(const Matrix* source, int row, Vector* out, Diagnostic* d) {
  const int indices[3] = {3 * row, 3 * row + 1, 3 * row + 2};
  return source->GatherTo(indices, 3, out, 0, d);
}
StatusCode Cross(const Vector* a, const Vector* b, Vector* out, Diagnostic* d) {
  Matrix* raw = nullptr;
  auto status = Matrix::Create(3, 3, &raw, d);
  std::unique_ptr<Matrix> skew(raw);
  if (status != StatusCode::kOk) return status;
  const double* v = a->Data();
  const double c[9] = {0, -v[2], v[1], v[2], 0, -v[0], -v[1], v[0], 0};
  for (int i = 0; i < 9; ++i) {
    status = raw->Set(i / 3, i % 3, c[i], d);
    if (status != StatusCode::kOk) return status;
  }
  return Matrix::MultiplyVector(raw, false, b, out, d);
}
}  // namespace
// Caller owns the returned normalized 4x2 local coordinates. All intermediate
// geometry is call-local. Roundoff guards are formulation 3.3, not tolerances
// for convergence or acceptance.
StatusCode CreateCoordinates(const Node* const* nodes, Matrix** out,
                             double* scale, Diagnostic* d) {
  *out = nullptr;
  std::unique_ptr<Matrix> matrices[6];
  const int rows[6] = {4, 4, 1, 4, 2, 4};
  const int cols[6] = {3, 1, 3, 3, 3, 2};
  for (int i = 0; i < 6; ++i) {
    Matrix* raw = nullptr;
    const auto status = Matrix::Create(rows[i], cols[i], &raw, d);
    matrices[i].reset(raw);
    if (status != StatusCode::kOk) return status;
  }
  std::unique_ptr<Vector> vectors[4];
  for (auto& owner : vectors) {
    Vector* raw = nullptr;
    const auto status = Vector::Create(3, &raw, d);
    owner.reset(raw);
    if (status != StatusCode::kOk) return status;
  }
  Matrix* xyz = matrices[0].get();
  Matrix* ones = matrices[1].get();
  Matrix* anchor = matrices[2].get();
  Matrix* translation = matrices[3].get();
  Matrix* basis = matrices[4].get();
  Matrix* local = matrices[5].get();
  Vector* first = vectors[0].get();
  Vector* second = vectors[1].get();
  Vector* normal = vectors[2].get();
  Vector* point = vectors[3].get();
  auto status = StatusCode::kOk;
  for (int a = 0; a < 4; ++a) {
    for (int i = 0; i < 3; ++i) {
      double coordinate = 0;
      status = nodes[a]->Coordinate(i, &coordinate, d);
      if (status == StatusCode::kOk) status = xyz->Set(a, i, coordinate, d);
      if (status != StatusCode::kOk) return status;
    }
    status = ones->Set(a, 0, 1, d);
    if (status != StatusCode::kOk) return status;
  }
  status = anchor->CopyBlockFrom(xyz, 0, 0, 1, 3, 0, 0, d);
  if (status == StatusCode::kOk)
    status = Matrix::Multiply(ones, false, anchor, false, translation, d);
  if (status == StatusCode::kOk) status = xyz->AddScaled(-1, translation, d);
  if (status != StatusCode::kOk) return status;
  double diameter = 0;
  for (int a = 0; a < 4; ++a)
    for (int b = a + 1; b < 4; ++b) {
      status = Row(xyz, a, first, d);
      if (status == StatusCode::kOk) status = Row(xyz, b, second, d);
      if (status == StatusCode::kOk) status = first->AddScaled(-1, second, d);
      double length = 0;
      if (status == StatusCode::kOk) status = first->Norm2(&length, d);
      if (status != StatusCode::kOk) return status;
      diameter = std::max(diameter, length);
    }
  if (!(diameter > 0)) return GeometryError(d);
  if (!std::isfinite(1 / diameter))
    return ReportFailure(StatusCode::kLibraryFailure,
                         DiagnosticStage::kComputation,
                         "Thermal face geometry scale range", d);
  status = xyz->Scale(1 / diameter, d);
  // Diagonal normal: (x3-x1) cross (x4-x2), not a mechanical edge normal.
  if (status == StatusCode::kOk) status = Row(xyz, 2, first, d);
  if (status == StatusCode::kOk) status = Row(xyz, 3, second, d);
  if (status == StatusCode::kOk) status = Row(xyz, 1, point, d);
  if (status == StatusCode::kOk) status = second->AddScaled(-1, point, d);
  if (status == StatusCode::kOk) status = Cross(first, second, normal, d);
  double norm = 0;
  if (status == StatusCode::kOk) status = normal->Norm2(&norm, d);
  if (status != StatusCode::kOk) return status;
  const double epsilon = std::numeric_limits<double>::epsilon();
  if (!(norm > 256 * epsilon)) return GeometryError(d);
  status = normal->Scale(1 / norm, d);
  if (status == StatusCode::kOk) status = Row(xyz, 1, first, d);
  double projection = 0;
  if (status == StatusCode::kOk) status = first->Dot(normal, &projection, d);
  if (status == StatusCode::kOk)
    status = first->AddScaled(-projection, normal, d);
  if (status == StatusCode::kOk) status = first->Norm2(&norm, d);
  if (status != StatusCode::kOk) return status;
  if (!(norm > 256 * epsilon)) return GeometryError(d);
  status = first->Scale(1 / norm, d);
  if (status == StatusCode::kOk) status = Cross(normal, first, second, d);
  if (status == StatusCode::kOk) status = basis->SetRow(0, first, d);
  if (status == StatusCode::kOk) status = basis->SetRow(1, second, d);
  if (status != StatusCode::kOk) return status;
  for (int a = 0; a < 4; ++a) {
    status = Row(xyz, a, point, d);
    double height = 0;
    if (status == StatusCode::kOk) status = normal->Dot(point, &height, d);
    if (status != StatusCode::kOk) return status;
    if (std::abs(height) > 64 * epsilon)
      return ReportFailure(StatusCode::kUnsupported,
                           DiagnosticStage::kValidation,
                           "Non-planar thermal face geometry", d);
  }
  status = Matrix::Multiply(xyz, false, basis, true, local, d);
  if (status == StatusCode::kOk) {
    *scale = diameter;
    *out = matrices[5].release();
  }
  return status;
}
constexpr int kSigns[4][2] = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}};
// Normalized J[alpha,i]=dxhat_i/dxi_alpha; J grad_xhat(N)=grad_xi(N).
// Physical length factors are restored only in the scalar Gram/volume terms.
// shape and gradient are borrowed outputs. Null gradient skips the solve.
StatusCode MapPoint(const Matrix* xy, int point, Matrix* gradient,
                    double* shape, double* determinant, Diagnostic* d) {
  const double xi = kSigns[point][0] / std::sqrt(3.0);
  const double eta = kSigns[point][1] / std::sqrt(3.0);
  Matrix* raw = nullptr;
  auto status = Matrix::Create(2, 4, &raw, d);
  std::unique_ptr<Matrix> natural(raw);
  if (status != StatusCode::kOk) return status;
  for (int a = 0; a < 4; ++a) {
    const int sx = kSigns[a][0], sy = kSigns[a][1];
    if (shape) shape[a] = (1 + sx * xi) * (1 + sy * eta) / 4;
    status = natural->Set(0, a, sx * (1 + sy * eta) / 4, d);
    if (status == StatusCode::kOk)
      status = natural->Set(1, a, sy * (1 + sx * xi) / 4, d);
    if (status != StatusCode::kOk) return status;
  }
  status = Matrix::Create(2, 2, &raw, d);
  std::unique_ptr<Matrix> jacobian(raw);
  if (status == StatusCode::kOk)
    status = Matrix::Multiply(natural.get(), false, xy, false, raw, d);
  if (status != StatusCode::kOk) return status;
  DenseLuFactorization* factor = nullptr;
  status = DenseLuFactorization::Create(raw, &factor, d);
  std::unique_ptr<DenseLuFactorization> owner(factor);
  if (status == StatusCode::kSingular) return GeometryError(d);
  if (status != StatusCode::kOk) return status;
  double det = 0;
  status = factor->Determinant(&det, d);
  if (status != StatusCode::kOk) return status;
  if (!(det > 0)) return GeometryError(d);
  if (gradient) status = factor->Solve(natural.get(), false, gradient, d);
  if (status == StatusCode::kOk) *determinant = det;
  return status;
}

StatusCode AssignWeights(const ThermalElement* element, const int* local_nodes,
                         int count, int section_point, Vector* out,
                         Diagnostic* d) {
  if (!element || !local_nodes || !out || (count != 3 && count != 4) ||
      out->Size() != element->dof_count() || section_point < 0 ||
      section_point >= element->dofs_per_node())
    return model_storage::Invalid("Thermal surface weight arguments", d);
  const Node* nodes[4] = {};
  int indices[4] = {};
  for (int i = 0; i < count; ++i) {
    nodes[i] = element->GetNode(local_nodes[i]);
    if (!nodes[i]) return model_storage::Invalid("Thermal face node", d);
    indices[i] = local_nodes[i] * element->dofs_per_node() + section_point;
  }
  Vector* raw = nullptr;
  auto status = Vector::Create(count, &raw, d);
  std::unique_ptr<Vector> area(raw);
  if (status != StatusCode::kOk) return status;
  if (count == 4) {
    Matrix* raw_xy = nullptr;
    double length = 0;
    status = CreateCoordinates(nodes, &raw_xy, &length, d);
    std::unique_ptr<Matrix> xy(raw_xy);
    if (status != StatusCode::kOk) return status;
    status = Vector::Create(4, &raw, d);
    std::unique_ptr<Vector> station(raw);
    if (status != StatusCode::kOk) return status;
    for (int g = 0; g < 4; ++g) {
      double shape[4] = {}, determinant = 0;
      status = MapPoint(xy.get(), g, nullptr, shape, &determinant, d);
      for (int i = 0; i < 4 && status == StatusCode::kOk; ++i) {
        const double weight = determinant * length * shape[i] * length;
        if (!std::isfinite(weight) || weight <= 0)
          return ReportFailure(StatusCode::kLibraryFailure,
                               DiagnosticStage::kComputation,
                               "Thermal face area range", d);
        status = station->Set(i, weight, d);
      }
      if (status == StatusCode::kOk)
        status = area->AddScaled(1, station.get(), d);
      if (status != StatusCode::kOk) return status;
    }
  } else {
    std::unique_ptr<Vector> vectors[4];
    for (auto& owner : vectors) {
      status = Vector::Create(3, &raw, d);
      owner.reset(raw);
      if (status != StatusCode::kOk) return status;
    }
    for (int i = 0; i < 3; ++i)
      for (int axis = 0; axis < 3; ++axis) {
        double coordinate = 0;
        status = nodes[i]->Coordinate(axis, &coordinate, d);
        if (status == StatusCode::kOk)
          status = vectors[i]->Set(axis, coordinate, d);
        if (status != StatusCode::kOk) return status;
      }
    status = vectors[1]->AddScaled(-1, vectors[0].get(), d);
    if (status == StatusCode::kOk)
      status = vectors[2]->AddScaled(-1, vectors[0].get(), d);
    if (status == StatusCode::kOk)
      status = vectors[3]->CopyFrom(vectors[1].get(), d);
    if (status == StatusCode::kOk)
      status = vectors[3]->AddScaled(-1, vectors[2].get(), d);
    double length = 0;
    for (int i = 1; i < 4 && status == StatusCode::kOk; ++i) {
      double value = 0;
      status = vectors[i]->Norm2(&value, d);
      length = std::max(length, value);
    }
    if (status != StatusCode::kOk) return status;
    if (!(length > 0) || !std::isfinite(1 / length)) return GeometryError(d);
    status = vectors[1]->Scale(1 / length, d);
    if (status == StatusCode::kOk) status = vectors[2]->Scale(1 / length, d);
    if (status == StatusCode::kOk)
      status = Cross(vectors[1].get(), vectors[2].get(), vectors[3].get(), d);
    double norm = 0;
    if (status == StatusCode::kOk) status = vectors[3]->Norm2(&norm, d);
    if (status != StatusCode::kOk) return status;
    const double weight = (norm / 6) * length * length;
    if (!(weight > 0) || !std::isfinite(weight))
      return ReportFailure(StatusCode::kLibraryFailure,
                           DiagnosticStage::kComputation,
                           "Thermal triangular face area range", d);
    for (int i = 0; i < 3 && status == StatusCode::kOk; ++i)
      status = area->Set(i, weight, d);
  }
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(out->Size(), &raw, d);
  std::unique_ptr<Vector> zero(raw);
  if (status == StatusCode::kOk) status = out->CopyFrom(zero.get(), d);
  if (status == StatusCode::kOk)
    status = out->ScatterAdd(area.get(), indices, count, d);
  return status;
}
}  // namespace thermal_face_geometry
}  // namespace fesa
