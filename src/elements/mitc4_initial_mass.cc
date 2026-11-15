#include "elements/mitc4_initial.h"

#include <cmath>
#include <memory>

#include "elements/mitc4_geometry.h"
#include "fesa/math/dense_factorization.h"
#include "model/model_storage.h"

namespace fesa {
namespace {

bool MassShape(const Matrix* matrix) {
  return matrix && matrix->Rows() == 24 && matrix->Columns() == 24;
}

StatusCode CreateMatrix(int rows, int columns, std::unique_ptr<Matrix>* owner,
                        Diagnostic* d) {
  Matrix* raw = nullptr;
  const auto status = Matrix::Create(rows, columns, &raw, d);
  owner->reset(raw);
  return status;
}

StatusCode FillTransform(const Matrix* frame, Matrix* transform,
                         Diagnostic* d) {
  // q_local=T*q_global: identical row basis for translations and rotations.
  for (int block = 0; block < 8; ++block) {
    const auto status = transform->CopyBlockFrom(
        frame, 0, 0, 3, 3, 3 * block, 3 * block, d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}

StatusCode SetDerivatives(double r, double s, Matrix* derivatives,
                          Diagnostic* d) {
  const int signs[4][2] = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}};
  for (int i = 0; i < 4; ++i) {
    auto status = derivatives->Set(
        0, i, signs[i][0] * (1.0 + signs[i][1] * s) / 4.0, d);
    if (status != StatusCode::kOk) return status;
    status = derivatives->Set(
        1, i, signs[i][1] * (1.0 + signs[i][0] * r) / 4.0, d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}

StatusCode SetVelocityInterpolation(double r, double s, double z, Matrix* h,
                                     Diagnostic* d) {
  const int signs[4][2] = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}};
  for (int i = 0; i < 4; ++i) {
    const double shape = (1.0 + signs[i][0] * r) *
                         (1.0 + signs[i][1] * s) / 4.0;
    for (int component = 0; component < 3; ++component) {
      const auto status = h->Set(component, 6 * i + component, shape, d);
      if (status != StatusCode::kOk) return status;
    }
    // Ni[I,-z[n]x], n=(0,0,1): vx+=Ni*z*omega_y, vy-=Ni*z*omega_x.
    auto status = h->Set(0, 6 * i + 4, shape * z, d);
    if (status != StatusCode::kOk) return status;
    status = h->Set(1, 6 * i + 3, -shape * z, d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}

StatusCode IntegrateMass(const Matrix* coordinates, double thickness,
                         double density, Matrix* local_mass, Diagnostic* d) {
  std::unique_ptr<Matrix> xy;
  std::unique_ptr<Matrix> derivatives;
  std::unique_ptr<Matrix> jacobian;
  std::unique_ptr<Matrix> interpolation;
  auto status = CreateMatrix(4, 2, &xy, d);
  if (status != StatusCode::kOk) return status;
  status = CreateMatrix(2, 4, &derivatives, d);
  if (status != StatusCode::kOk) return status;
  status = CreateMatrix(2, 2, &jacobian, d);
  if (status != StatusCode::kOk) return status;
  status = CreateMatrix(3, 24, &interpolation, d);
  if (status != StatusCode::kOk) return status;
  status = xy->CopyBlockFrom(coordinates, 0, 0, 4, 2, 0, 0, d);
  if (status != StatusCode::kOk) return status;
  const double gauss = 1.0 / std::sqrt(3.0);
  for (double r : {-gauss, gauss}) {
    for (double s : {-gauss, gauss}) {
      status = SetDerivatives(r, s, derivatives.get(), d);
      if (status != StatusCode::kOk) return status;
      status = Matrix::Multiply(derivatives.get(), false, xy.get(), false,
                                 jacobian.get(), d);
      if (status != StatusCode::kOk) return status;
      DenseLuFactorization* raw_factor = nullptr;
      status = DenseLuFactorization::Create(jacobian.get(), &raw_factor, d);
      std::unique_ptr<DenseLuFactorization> factor(raw_factor);
      if (status != StatusCode::kOk) return status;
      double determinant = 0.0;
      status = factor->Determinant(&determinant, d);
      if (status != StatusCode::kOk) return status;
      if (!(determinant > 0.0)) {
        return ReportFailure(StatusCode::kInputError,
                             DiagnosticStage::kValidation,
                             "MITC4 initial mass nonpositive Jacobian", d);
      }
      const double weight = density * determinant * (thickness / 2.0);
      if (!std::isfinite(weight) || !(weight > 0.0)) {
        return ReportFailure(StatusCode::kLibraryFailure,
                             DiagnosticStage::kComputation,
                             "MITC4 initial mass integration weight range", d);
      }
      for (double zeta : {-gauss, gauss}) {
        const double z = (thickness / 2.0) * zeta;
        status = SetVelocityInterpolation(r, s, z, interpolation.get(), d);
        if (status != StatusCode::kOk) return status;
        // The eight full-volume Gauss weights are one. All matrix products
        // and integration accumulation use the common MKL-backed operations.
        status = Matrix::MultiplyAdd(interpolation.get(), true,
                                      interpolation.get(), false, weight, 1.0,
                                      local_mass, d);
        if (status != StatusCode::kOk) return status;
      }
    }
  }
  return StatusCode::kOk;
}

}  // namespace

StatusCode ComputeMitc4InitialMass(const Node* const* nodes, int count,
                                  const ShellSection* section, Matrix* mass,
                                  Diagnostic* d) {
  if (!nodes || count != 4 || !section || !MassShape(mass))
    return model_storage::Invalid("MITC4 initial mass inputs/shape", d);
  double density = 0.0;
  auto status = section->material()->Density(&density, d);
  if (status != StatusCode::kOk) return status;
  const double thickness = section->thickness();
  if (!std::isfinite(density) || !(density > 0.0) ||
      !std::isfinite(thickness) || !(thickness > 0.0)) {
    return ReportFailure(StatusCode::kInputError, DiagnosticStage::kValidation,
                         "MITC4 initial mass density/thickness", d);
  }
  std::unique_ptr<Matrix> coordinates;
  std::unique_ptr<Matrix> frame;
  std::unique_ptr<Matrix> transform;
  std::unique_ptr<Matrix> local_mass;
  std::unique_ptr<Matrix> product;
  status = CreateMatrix(4, 3, &coordinates, d);
  if (status != StatusCode::kOk) return status;
  status = CreateMatrix(3, 3, &frame, d);
  if (status != StatusCode::kOk) return status;
  double diameter = 0.0;
  double length = 0.0;
  status = ValidateMitc4Geometry(nodes, count, coordinates.get(), frame.get(),
                                 &diameter, &length, d);
  if (status != StatusCode::kOk) return status;
  status = CreateMatrix(24, 24, &transform, d);
  if (status != StatusCode::kOk) return status;
  status = CreateMatrix(24, 24, &local_mass, d);
  if (status != StatusCode::kOk) return status;
  status = CreateMatrix(24, 24, &product, d);
  if (status != StatusCode::kOk) return status;
  status = FillTransform(frame.get(), transform.get(), d);
  if (status != StatusCode::kOk) return status;
  status = IntegrateMass(coordinates.get(), thickness, density,
                          local_mass.get(), d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Multiply(local_mass.get(), false, transform.get(), false,
                             product.get(), d);
  if (status != StatusCode::kOk) return status;
  return Matrix::Multiply(transform.get(), true, product.get(), false, mass, d);
}

StatusCode LumpMitc4InitialMass(const Node* const* nodes, int count,
                               const Matrix* consistent, Matrix* lumped,
                               Diagnostic* d) {
  if (!nodes || count != 4 || !MassShape(consistent) || !MassShape(lumped) ||
      consistent == lumped) {
    return model_storage::Invalid("MITC4 initial lumping inputs/shape/alias",
                                  d);
  }
  for (int i = 0; i < 24 * 24; ++i) {
    if (!std::isfinite(consistent->Data()[i]))
      return model_storage::Invalid("MITC4 initial lumping nonfinite input", d);
  }
  std::unique_ptr<Matrix> matrices[8];
  const int rows[8] = {4, 3, 24, 24, 24, 24, 24, 24};
  const int columns[8] = {3, 3, 24, 24, 24, 4, 24, 24};
  for (int i = 0; i < 8; ++i) {
    const auto status = CreateMatrix(rows[i], columns[i], &matrices[i], d);
    if (status != StatusCode::kOk) return status;
  }
  Matrix* coordinates = matrices[0].get();
  Matrix* frame = matrices[1].get();
  Matrix* transform = matrices[2].get();
  Matrix* product = matrices[3].get();
  Matrix* local = matrices[4].get();
  Matrix* selected = matrices[5].get();
  Matrix* diagonal = matrices[6].get();
  Matrix* transformed = matrices[7].get();
  double diameter = 0.0;
  double length = 0.0;
  auto status = ValidateMitc4Geometry(nodes, count, coordinates, frame,
                                      &diameter, &length, d);
  if (status != StatusCode::kOk) return status;
  status = FillTransform(frame, transform, d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Multiply(transform, false, consistent, false, product, d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Multiply(product, false, transform, true, local, d);
  if (status != StatusCode::kOk) return status;
  Vector* raw_row = nullptr;
  status = Vector::Create(4, &raw_row, d);
  std::unique_ptr<Vector> row_values(raw_row);
  if (status != StatusCode::kOk) return status;
  Vector* raw_sums = nullptr;
  status = Vector::Create(24, &raw_sums, d);
  std::unique_ptr<Vector> sums(raw_sums);
  if (status != StatusCode::kOk) return status;
  // mass_lumping_contract.md: sum only the four nodal columns of the same
  // local component. The drilling rows stay exactly zero by definition.
  for (int row = 0; row < 24; ++row) {
    const int component = row % 6;
    if (component == 5) continue;
    int indices[4] = {};
    for (int node = 0; node < 4; ++node)
      indices[node] = 24 * row + 6 * node + component;
    status = local->GatherTo(indices, 4, row_values.get(), 0, d);
    if (status != StatusCode::kOk) return status;
    status = selected->SetRow(row, row_values.get(), d);
    if (status != StatusCode::kOk) return status;
  }
  status = selected->RowSums(sums.get(), d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Diagonal(sums.get(), diagonal, d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Multiply(diagonal, false, transform, false, transformed, d);
  if (status != StatusCode::kOk) return status;
  return Matrix::Multiply(transform, true, transformed, false, lumped, d);
}

}  // namespace fesa
