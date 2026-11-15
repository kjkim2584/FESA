#include "elements/solid_nonlinear_recovery.h"

#include <cmath>
#include <memory>

#include "elements/solid_recovery.h"
#include "elements/solid_stvk.h"
#include "fesa/math/dense_factorization.h"

namespace fesa {
namespace {
bool PointShape(const Matrix* matrix) {
  return matrix && matrix->Rows() == 3 && matrix->Columns() == 3;
}
StatusCode InvalidOutput(Diagnostic* d) {
  return ReportFailure(StatusCode::kInvalidArgument,
                       DiagnosticStage::kValidation, "Solid point output", d);
}
StatusCode PointMatrix(std::unique_ptr<Matrix>* owner, Diagnostic* d) {
  Matrix* raw = nullptr;
  const auto status = Matrix::Create(3, 3, &raw, d);
  owner->reset(raw);
  return status;
}
}  // namespace

StatusCode ComputeSolidPointOutput(const Matrix* f,
    const ElasticMaterial* material, Matrix* cauchy, Matrix* logarithmic_strain,
    TrialEvaluationFailure* failure, Diagnostic* d) {
  if (failure) *failure = TrialEvaluationFailure::kNone;
  if (!failure || !material || !PointShape(f) || !PointShape(cauchy) ||
      !PointShape(logarithmic_strain) || cauchy == logarithmic_strain ||
      cauchy == f || logarithmic_strain == f) return InvalidOutput(d);
  for (int i = 0; i < 9; ++i)
    if (!std::isfinite(f->Data()[i])) return InvalidOutput(d);
  std::unique_ptr<Matrix> first_piola, left;
  auto status = PointMatrix(&first_piola, d);
  if (status != StatusCode::kOk) return status;
  status = EvaluateStvkPoint(f, material, nullptr, nullptr, first_piola.get(),
                             nullptr, d);
  if (status == StatusCode::kEquilibriumFailure)
    *failure = TrialEvaluationFailure::kGeometry;
  if (status != StatusCode::kOk) return status;
  DenseLuFactorization* raw = nullptr;
  status = DenseLuFactorization::Create(f, &raw, d);
  std::unique_ptr<DenseLuFactorization> factor(raw);
  if (status != StatusCode::kOk) return status;
  double determinant = 0;
  status = factor->Determinant(&determinant, d);
  if (status != StatusCode::kOk) return status;
  const double inverse_j = 1 / determinant;
  if (!std::isfinite(inverse_j)) {
    return ReportFailure(StatusCode::kLibraryFailure,
                         DiagnosticStage::kComputation,
                         "Solid output inverse determinant overflow", d);
  }
  status = Matrix::Multiply(first_piola.get(), false, f, true, cauchy, d);
  if (status != StatusCode::kOk) return status;
  status = cauchy->Scale(inverse_j, d);
  if (status != StatusCode::kOk) return status;
  status = PointMatrix(&left, d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Multiply(f, false, f, true, left.get(), d);
  if (status != StatusCode::kOk) return status;
  bool not_positive = false;
  status = Matrix::SymmetricPositiveDefiniteLog(left.get(), logarithmic_strain,
                                                &not_positive, d);
  if (not_positive) *failure = TrialEvaluationFailure::kGeometry;
  if (status != StatusCode::kOk) return status;
  return logarithmic_strain->Scale(.5, d);
}

StatusCode RecoverSolidPointOutputs(const SolidReferenceData* data,
    const Matrix* const* material_gradients, const ElasticMaterial* material,
    Matrix* ip_stress, Matrix* ip_strain, Matrix* nodal_stress,
    Matrix* nodal_strain, TrialEvaluationFailure* failure, Diagnostic* d) {
  if (failure) *failure = TrialEvaluationFailure::kNone;
  if (!data || !material || !material_gradients || !failure)
    return InvalidOutput(d);
  switch (data->formulation()) {
    case SolidFormulation::kTetra4Full:
    case SolidFormulation::kTetra10Full:
    case SolidFormulation::kHexa8Full:
    case SolidFormulation::kWedge6:
    case SolidFormulation::kPyramid5Full:
      break;
    default:
      return ReportFailure(StatusCode::kUnsupported,
                           DiagnosticStage::kValidation,
                           "Nonlinear solid recovery formulation", d);
  }
  auto status = ValidateSolidTensorRecovery(data, ip_stress, ip_strain,
                                             nodal_stress, nodal_strain, d);
  if (status != StatusCode::kOk) return status;
  for (int g = 0; g < data->point_count(); ++g) {
    const auto* f = material_gradients[g];
    if (!PointShape(f) || f == ip_stress || f == ip_strain ||
        f == nodal_stress || f == nodal_strain) return InvalidOutput(d);
    for (int i = 0; i < 9; ++i)
      if (!std::isfinite(f->Data()[i])) return InvalidOutput(d);
  }
  std::unique_ptr<Matrix> cauchy, logarithmic_strain;
  status = PointMatrix(&cauchy, d);
  if (status != StatusCode::kOk) return status;
  status = PointMatrix(&logarithmic_strain, d);
  if (status != StatusCode::kOk) return status;
  Vector* row_raw = nullptr;
  status = Vector::Create(6, &row_raw, d);
  std::unique_ptr<Vector> row(row_raw);
  if (status != StatusCode::kOk) return status;
  const int order[] = {0, 4, 8, 1, 2, 5};
  for (int g = 0; g < data->point_count(); ++g) {
    status = ComputeSolidPointOutput(material_gradients[g], material,
        cauchy.get(), logarithmic_strain.get(), failure, d);
    if (status != StatusCode::kOk) return status;
    status = cauchy->GatherTo(order, 6, row.get(), 0, d);
    if (status != StatusCode::kOk) return status;
    status = ip_stress->SetRow(g, row.get(), d);
    if (status != StatusCode::kOk) return status;
    status = logarithmic_strain->GatherTo(order, 6, row.get(), 0, d);
    if (status != StatusCode::kOk) return status;
    status = ip_strain->SetRow(g, row.get(), d);
    if (status != StatusCode::kOk) return status;
  }
  return ExtrapolateSolidTensors(data, ip_stress, ip_strain, nodal_stress,
                                  nodal_strain, d);
}
}  // namespace fesa
