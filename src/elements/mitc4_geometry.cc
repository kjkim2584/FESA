#include "elements/mitc4_geometry.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>

#include "model/model_storage.h"
namespace fesa {
namespace {
StatusCode GeometryError(Diagnostic* d) {
  return ReportFailure(StatusCode::kInputError, DiagnosticStage::kValidation,
                       "MITC4 warped, folded or degenerate initial geometry",
                       d);
}
StatusCode Cross(const Vector* a, const Vector* b, Vector* out, Diagnostic* d) {
  Matrix* raw = nullptr;
  auto status = Matrix::Create(3, 3, &raw, d);
  std::unique_ptr<Matrix> skew(raw);
  if (status != StatusCode::kOk) return status;
  const double* v = a->Data();
  const double coefficients[9] = {0,     -v[2], v[1], v[2], 0,
                                  -v[0], -v[1], v[0], 0};
  for (int i = 0; i < 9; ++i) {
    status = skew->Set(i / 3, i % 3, coefficients[i], d);
    if (status != StatusCode::kOk) return status;
  }
  return Matrix::MultiplyVector(skew.get(), false, b, out, d);
}
StatusCode Row(const Matrix* source, int row, Vector* out, Diagnostic* d) {
  const int indices[3] = {3 * row, 3 * row + 1, 3 * row + 2};
  return source->GatherTo(indices, 3, out, 0, d);
}
}  // namespace
StatusCode ValidateMitc4Geometry(const Node* const* nodes, int count,
                                 Matrix* local_coordinates, Matrix* frame,
                                 double* diameter, double* length,
                                 Diagnostic* d) {
  if (!nodes || count != 4 || !local_coordinates || !frame || !diameter ||
      !length || diameter == length || local_coordinates == frame ||
      local_coordinates->Rows() != 4 || local_coordinates->Columns() != 3 ||
      frame->Rows() != 3 || frame->Columns() != 3)
    return model_storage::Invalid("MITC4 geometry output/shape", d);
  for (int i = 0; i < 4; ++i) {
    if (!nodes[i] || nodes[i]->dof_count() != 6)
      return model_storage::Invalid("MITC4 geometry connectivity", d);
    for (int j = 0; j < i; ++j)
      if (nodes[i] == nodes[j] || nodes[i]->Id() == nodes[j]->Id())
        return model_storage::Invalid("MITC4 duplicate node", d);
  }
  std::unique_ptr<Matrix> matrices[8];
  const int rows[8] = {4, 4, 1, 4, 3, 2, 2, 4};
  const int cols[8] = {3, 1, 3, 3, 3, 4, 3, 3};
  for (int i = 0; i < 8; ++i) {
    Matrix* raw = nullptr;
    const auto code = Matrix::Create(rows[i], cols[i], &raw, d);
    matrices[i].reset(raw);
    if (code != StatusCode::kOk) return code;
  }
  std::unique_ptr<Vector> vectors[4];
  for (int i = 0; i < 4; ++i) {
    Vector* raw = nullptr;
    const auto code = Vector::Create(3, &raw, d);
    vectors[i].reset(raw);
    if (code != StatusCode::kOk) return code;
  }
  Matrix* coordinates = matrices[0].get();
  // Borrowed aliases of the call-local owners above; no cross-call storage.
  Matrix* ones = matrices[1].get();
  Matrix* anchor = matrices[2].get();
  Matrix* translation = matrices[3].get();
  Matrix* basis = matrices[4].get();
  Matrix* derivatives = matrices[5].get();
  Matrix* tangents = matrices[6].get();
  Matrix* local = matrices[7].get();
  Vector* first = vectors[0].get();
  Vector* second = vectors[1].get();
  Vector* normal = vectors[2].get();
  Vector* point = vectors[3].get();
  StatusCode status = StatusCode::kOk;
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 3; ++j) {
      double value = 0;
      status = nodes[i]->Coordinate(j, &value, d);
      if (status != StatusCode::kOk) return status;
      status = coordinates->Set(i, j, value, d);
      if (status != StatusCode::kOk) return status;
    }
  double scale = 0;
  double perimeter = 0;
  for (int i = 0; i < 4; ++i)
    for (int j = i + 1; j < 4; ++j) {
      status = Row(coordinates, i, first, d);
      if (status != StatusCode::kOk) return status;
      status = Row(coordinates, j, second, d);
      if (status != StatusCode::kOk) return status;
      status = first->AddScaled(-1, second, d);
      if (status != StatusCode::kOk) return status;
      double distance = 0;
      status = first->Norm2(&distance, d);
      if (status != StatusCode::kOk) return status;
      scale = std::max(scale, distance);
      if (j == i + 1 || (i == 0 && j == 3)) perimeter += distance / 4;
    }
  if (!(scale > 0)) return GeometryError(d);
  if (!std::isfinite(1 / scale) || !std::isfinite(perimeter))
    return ReportFailure(StatusCode::kLibraryFailure,
                         DiagnosticStage::kComputation,
                         "MITC4 geometry scale range", d);
  for (int i = 0; i < 4; ++i) {
    status = ones->Set(i, 0, 1, d);
    if (status != StatusCode::kOk) return status;
  }
  status = anchor->CopyBlockFrom(coordinates, 0, 0, 1, 3, 0, 0, d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Multiply(ones, false, anchor, false, translation, d);
  if (status != StatusCode::kOk) return status;
  status = coordinates->AddScaled(-1, translation, d);
  if (status != StatusCode::kOk) return status;
  status = coordinates->Scale(1 / scale, d);
  if (status != StatusCode::kOk) return status;
  status = Row(coordinates, 1, first, d);
  if (status != StatusCode::kOk) return status;
  status = Row(coordinates, 3, second, d);
  if (status != StatusCode::kOk) return status;
  status = Cross(first, second, normal, d);
  if (status != StatusCode::kOk) return status;
  double norm = 0;
  status = normal->Norm2(&norm, d);
  if (status != StatusCode::kOk) return status;
  const double epsilon = std::numeric_limits<double>::epsilon();
  if (!(norm > 256 * epsilon)) return GeometryError(d);
  status = normal->Scale(1 / norm, d);
  if (status != StatusCode::kOk) return status;
  status = first->Norm2(&norm, d);
  if (status != StatusCode::kOk) return status;
  status = first->Scale(1 / norm, d);
  if (status != StatusCode::kOk) return status;
  status = Cross(normal, first, second, d);
  if (status != StatusCode::kOk) return status;
  for (int i = 0; i < 3; ++i) {
    status = basis->SetRow(i, vectors[i].get(), d);
    if (status != StatusCode::kOk) return status;
  }
  for (int i = 0; i < 4; ++i) {
    status = Row(coordinates, i, point, d);
    if (status != StatusCode::kOk) return status;
    double height = 0;
    status = normal->Dot(point, &height, d);
    if (status != StatusCode::kOk) return status;
    if (std::abs(height) > 64 * epsilon) return GeometryError(d);
  }
  const int signs[4][2] = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}};
  for (int corner = 0; corner < 4; ++corner) {
    const double r = signs[corner][0];
    const double s = signs[corner][1];
    for (int i = 0; i < 4; ++i) {
      status =
          derivatives->Set(0, i, signs[i][0] * (1 + signs[i][1] * s) / 4.0, d);
      if (status != StatusCode::kOk) return status;
      status =
          derivatives->Set(1, i, signs[i][1] * (1 + signs[i][0] * r) / 4.0, d);
      if (status != StatusCode::kOk) return status;
    }
    status =
        Matrix::Multiply(derivatives, false, coordinates, false, tangents, d);
    if (status != StatusCode::kOk) return status;
    status = Row(tangents, 0, first, d);
    if (status != StatusCode::kOk) return status;
    status = Row(tangents, 1, second, d);
    if (status != StatusCode::kOk) return status;
    status = Cross(first, second, point, d);
    if (status != StatusCode::kOk) return status;
    double determinant = 0;
    status = normal->Dot(point, &determinant, d);
    if (status != StatusCode::kOk) return status;
    if (!(determinant > 256 * epsilon)) return GeometryError(d);
  }
  status = Matrix::Multiply(coordinates, false, basis, true, local, d);
  if (status != StatusCode::kOk) return status;
  status = local->Scale(scale, d);
  if (status != StatusCode::kOk) return status;
  status = frame->CopyFrom(basis, d);
  if (status != StatusCode::kOk) return status;
  status = local_coordinates->CopyFrom(local, d);
  if (status != StatusCode::kOk) return status;
  *diameter = scale;
  *length = perimeter;
  return StatusCode::kOk;
}
}  // namespace fesa
