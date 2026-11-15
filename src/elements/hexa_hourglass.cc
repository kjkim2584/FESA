#include "elements/hexa_hourglass.h"

#include <cmath>
#include <memory>
#include <new>

#include "fesa/model/elastic_material.h"
#include "model/model_storage.h"
namespace fesa {
HexaHourglass::~HexaHourglass() {
  delete h_;
  delete c_;
  delete weighted_h_;
  delete q_;
  delete resistance_;
  delete force_;
}
StatusCode HexaHourglass::Create(const Matrix* coordinates, double scale,
                                 const Matrix* gradient, double volume,
                                 const ElasticMaterial* material,
                                 bool stiffness, HexaHourglass** out,
                                 Diagnostic* d) {
  if (!out) return model_storage::Invalid("Hourglass output", d);
  *out = nullptr;
  if (!coordinates || coordinates->Rows() != 8 || coordinates->Columns() != 3 ||
      !gradient || gradient->Rows() != 3 || gradient->Columns() != 8 ||
      !material || !(scale > 0) || !std::isfinite(scale) || !(volume > 0) ||
      !std::isfinite(volume))
    return model_storage::Invalid("Hourglass geometry/material arguments", d);
  std::unique_ptr<HexaHourglass> result(new (std::nothrow) HexaHourglass);
  if (!result) return model_storage::AllocationFailure("Hourglass object", d);
  auto status = result->Initialize(coordinates, scale, gradient, volume,
                                   material, stiffness, d);
  if (status == StatusCode::kOk) *out = result.release();
  return status;
}
StatusCode HexaHourglass::Initialize(const Matrix* coordinates, double scale,
                                     const Matrix* gradient, double volume,
                                     const ElasticMaterial* material,
                                     bool stiffness, Diagnostic* d) {
  Matrix *modes = nullptr, *gamma_x = nullptr, *projection = nullptr,
         *gram = nullptr, *padded = nullptr;
  Vector *row = nullptr, *diagonal = nullptr;
  auto status = Matrix::Create(4, 8, &modes, d);
  std::unique_ptr<Matrix> modes_owner(modes);
  if (status != StatusCode::kOk) return status;
  const int gamma[4][8] = {{1, 1, -1, -1, -1, -1, 1, 1},
                           {1, -1, -1, 1, -1, 1, 1, -1},
                           {1, -1, 1, -1, 1, -1, 1, -1},
                           {-1, 1, -1, 1, 1, -1, 1, -1}};
  for (int mode = 0; mode < 4; ++mode)
    for (int node = 0; node < 8; ++node) {
      status = modes->Set(mode, node, gamma[mode][node], d);
      if (status != StatusCode::kOk) return status;
    }
  status = Matrix::Create(4, 3, &gamma_x, d);
  std::unique_ptr<Matrix> gamma_x_owner(gamma_x);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Create(4, 8, &projection, d);
  std::unique_ptr<Matrix> projection_owner(projection);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Multiply(modes, false, coordinates, false, gamma_x, d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Multiply(gamma_x, false, gradient, false, projection, d);
  if (status != StatusCode::kOk) return status;
  status = projection->Scale(scale, d);
  if (status != StatusCode::kOk) return status;
  status = modes->AddScaled(-1, projection, d);
  if (status != StatusCode::kOk) return status;
  // Lift H to H3=H tensor I3 with a zero mode row and external gather/copy.
  status = Matrix::Create(5, 8, &padded, d);
  std::unique_ptr<Matrix> padded_owner(padded);
  if (status != StatusCode::kOk) return status;
  status = padded->CopyBlockFrom(modes, 0, 0, 4, 8, 0, 0, d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Create(12, 24, &h_, d);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(24, &row, d);
  std::unique_ptr<Vector> row_owner(row);
  if (status != StatusCode::kOk) return status;
  for (int mode = 0; mode < 4; ++mode)
    for (int axis = 0; axis < 3; ++axis) {
      int indices[24] = {};
      for (int node = 0; node < 8; ++node)
        for (int component = 0; component < 3; ++component)
          indices[3 * node + component] =
              8 * (component == axis ? mode : 4) + node;
      status = padded->GatherTo(indices, 24, row, 0, d);
      if (status != StatusCode::kOk) return status;
      status = h_->SetRow(3 * mode + axis, row, d);
      if (status != StatusCode::kOk) return status;
    }
  status = Matrix::Create(3, 3, &gram, d);
  std::unique_ptr<Matrix> gram_owner(gram);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Multiply(gradient, false, gradient, true, gram, d);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(12, &diagonal, d);
  std::unique_ptr<Vector> diagonal_owner(diagonal);
  if (status != StatusCode::kOk) return status;
  const double mu = (material->young() / 2) / (1 + material->poisson());
  for (int i = 0; i < 12; ++i) {
    const int axis = i % 3;
    const double coefficient =
        .2357 * .005 * mu * volume * gram->Data()[4 * axis];
    if (!(coefficient > 0) || !std::isfinite(coefficient))
      return ReportFailure(StatusCode::kLibraryFailure,
                           DiagnosticStage::kComputation,
                           "Hourglass coefficient range", d);
    status = diagonal->Set(i, coefficient, d);
    if (status != StatusCode::kOk) return status;
  }
  status = Matrix::Create(12, 12, &c_, d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Diagonal(diagonal, c_, d);
  if (status != StatusCode::kOk) return status;
  if (stiffness) return Matrix::Create(12, 24, &weighted_h_, d);
  status = Vector::Create(12, &q_, d);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(12, &resistance_, d);
  if (status != StatusCode::kOk) return status;
  return Vector::Create(24, &force_, d);
}
StatusCode HexaHourglass::AddStiffness(Matrix* out, Diagnostic* d) {
  if (!weighted_h_ || !out || out->Rows() != 24 || out->Columns() != 24)
    return model_storage::Invalid("Hourglass stiffness output/purpose", d);
  auto status = Matrix::Multiply(c_, false, h_, false, weighted_h_, d);
  if (status != StatusCode::kOk) return status;
  return Matrix::MultiplyAdd(h_, true, weighted_h_, false, 1, 1, out, d);
}
StatusCode HexaHourglass::AddForce(const Vector* u, Vector* out,
                                   Diagnostic* d) {
  if (!q_ || !u || !out || u == out || u->Size() != 24 || out->Size() != 24)
    return model_storage::Invalid("Hourglass force input/output/purpose", d);
  auto status = Matrix::MultiplyVector(h_, false, u, q_, d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::MultiplyVector(c_, false, q_, resistance_, d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::MultiplyVector(h_, true, resistance_, force_, d);
  if (status != StatusCode::kOk) return status;
  return out->AddScaled(1, force_, d);
}
}  // namespace fesa
