#include "elements/solid_stvk.h"

#include <cmath>
#include <memory>

#include "elements/solid_access.h"
#include "fesa/math/dense_factorization.h"

namespace fesa {
namespace {
StatusCode MatrixBuffer(int size, std::unique_ptr<Matrix>* out, Diagnostic* d) {
  Matrix* raw = nullptr;
  const auto status = Matrix::Create(size, size, &raw, d);
  out->reset(raw);
  return status;
}
StatusCode VectorBuffer(int size, std::unique_ptr<Vector>* out, Diagnostic* d) {
  Vector* raw = nullptr;
  const auto status = Vector::Create(size, &raw, d);
  out->reset(raw);
  return status;
}
bool Shape(const Matrix* matrix, int size) {
  return matrix && matrix->Rows() == size && matrix->Columns() == size;
}
StatusCode Geometry(bool measured, double determinant, Diagnostic* d) {
  const auto status = ReportFailure(StatusCode::kEquilibriumFailure,
                       DiagnosticStage::kComputation,
                       "StVK nonpositive determinant", d);
  if (d) {
    d->geometry_quantity = GeometryQuantity::kDeformationJacobian;
    d->geometry_value_available = measured;
    d->geometry_value = determinant;
    d->geometry_limit = 0;
  }
  return status;
}
}  // namespace

StatusCode EvaluateStvkPoint(const Matrix* f, const ElasticMaterial* material,
                            double* energy, Matrix* second_piola,
                            Matrix* first_piola, Matrix* tangent, Diagnostic* d) {
  const Matrix* outputs[] = {second_piola, first_piola, tangent};
  bool valid = Shape(f, 3) && material &&
               (energy || second_piola || first_piola || tangent);
  for (int i = 0; i < 3; ++i) {
    if (!outputs[i]) continue;
    valid = valid && Shape(outputs[i], i == 2 ? 9 : 3) && outputs[i] != f;
    for (int j = 0; j < i; ++j) valid = valid && outputs[i] != outputs[j];
  }
  if (valid) {
    for (int i = 0; i < 9; ++i) valid = valid && std::isfinite(f->Data()[i]);
  }
  if (!valid) {
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation, "EvaluateStvkPoint", d);
  }
  DenseLuFactorization* raw_factor = nullptr;
  auto status = DenseLuFactorization::Create(f, &raw_factor, d);
  std::unique_ptr<DenseLuFactorization> factor(raw_factor);
  if (status == StatusCode::kSingular) return Geometry(false, 0, d);
  if (status != StatusCode::kOk) return status;
  double determinant = 0;
  status = factor->Determinant(&determinant, d);
  if (status != StatusCode::kOk) return status;
  if (!(determinant > 0)) return Geometry(true, determinant, d);

  ObserveSolidEvaluation(SolidEvaluationEvent::kConstitutive);
  const double nu = material->poisson();
  const double mu = material->young() / (2 * (1 + nu));
  const double lambda = material->young() * nu / ((1 + nu) * (1 - 2 * nu));
  if (!std::isfinite(mu) || !std::isfinite(lambda)) {
    return ReportFailure(StatusCode::kLibraryFailure,
                         DiagnosticStage::kComputation, "StVK Lame overflow", d);
  }
  std::unique_ptr<Matrix> identity, strain, stress;
  for (auto* owner : {&identity, &strain, &stress}) {
    status = MatrixBuffer(3, owner, d);
    if (status != StatusCode::kOk) return status;
  }
  for (int i = 0; i < 3; ++i) {
    status = identity->Set(i, i, 1, d);
    if (status != StatusCode::kOk) return status;
  }
  status = Matrix::Multiply(f, true, f, false, strain.get(), d);
  if (status != StatusCode::kOk) return status;
  status = strain->AddScaled(-1, identity.get(), d);
  if (status != StatusCode::kOk) return status;
  status = strain->Scale(.5, d);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<Vector> diagonal, ones;
  status = VectorBuffer(3, &diagonal, d);
  if (status != StatusCode::kOk) return status;
  status = VectorBuffer(3, &ones, d);
  if (status != StatusCode::kOk) return status;
  for (int i = 0; i < 3; ++i) {
    status = ones->Set(i, 1, d);
    if (status != StatusCode::kOk) return status;
  }
  const int diagonal_indices[] = {0, 4, 8};
  status = strain->GatherTo(diagonal_indices, 3, diagonal.get(), 0, d);
  if (status != StatusCode::kOk) return status;
  double trace = 0;
  status = diagonal->Dot(ones.get(), &trace, d);
  if (status != StatusCode::kOk) return status;
  const double twice_mu = 2 * mu;
  const double volumetric_stress = lambda * trace;
  if (!std::isfinite(twice_mu) || !std::isfinite(volumetric_stress)) {
    return ReportFailure(StatusCode::kLibraryFailure,
                         DiagnosticStage::kComputation,
                         "StVK stress coefficient overflow", d);
  }
  status = stress->CopyFrom(strain.get(), d);
  if (status != StatusCode::kOk) return status;
  status = stress->Scale(twice_mu, d);
  if (status != StatusCode::kOk) return status;
  status = stress->AddScaled(volumetric_stress, identity.get(), d);
  if (status != StatusCode::kOk) return status;

  double density = 0;
  if (energy) {
    std::unique_ptr<Vector> entries;
    status = VectorBuffer(9, &entries, d);
    if (status != StatusCode::kOk) return status;
    const int indices[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
    status = strain->GatherTo(indices, 9, entries.get(), 0, d);
    if (status != StatusCode::kOk) return status;
    double square = 0;
    status = entries->Dot(entries.get(), &square, d);
    if (status != StatusCode::kOk) return status;
    density = .5 * lambda * trace * trace + mu * square;
    if (!std::isfinite(density)) {
      return ReportFailure(StatusCode::kLibraryFailure,
                           DiagnosticStage::kComputation, "StVK energy", d);
    }
  }
  if (second_piola) {
    status = second_piola->CopyFrom(stress.get(), d);
    if (status != StatusCode::kOk) return status;
  }
  if (first_piola) {
    status = Matrix::Multiply(f, false, stress.get(), false, first_piola, d);
    if (status != StatusCode::kOk) return status;
  }
  if (tangent) {
    ObserveSolidEvaluation(SolidEvaluationEvent::kMaterialTangent);
    std::unique_ptr<Matrix> left;
    status = MatrixBuffer(3, &left, d);
    if (status != StatusCode::kOk) return status;
    status = Matrix::Multiply(f, false, f, true, left.get(), d);
    if (status != StatusCode::kOk) return status;
    // M3 defines each constitutive coefficient; no matrix product is expanded.
    for (int row = 0; row < 9; ++row) {
      const int i = row / 3, first = row % 3;
      for (int column = 0; column < 9; ++column) {
        const int j = column / 3, second = column % 3;
        const double coefficient =
            (i == j ? stress->Data()[3 * second + first] : 0) +
            lambda * f->Data()[row] * f->Data()[column] +
            mu * f->Data()[3 * i + second] * f->Data()[3 * j + first] +
            (first == second ? mu * left->Data()[3 * i + j] : 0);
        if (!std::isfinite(coefficient)) {
          return ReportFailure(StatusCode::kLibraryFailure,
                               DiagnosticStage::kComputation,
                               "StVK tangent overflow", d);
        }
        status = tangent->Set(row, column, coefficient, d);
        if (status != StatusCode::kOk) return status;
      }
    }
  }
  if (energy) *energy = density;
  return StatusCode::kOk;
}
}  // namespace fesa
