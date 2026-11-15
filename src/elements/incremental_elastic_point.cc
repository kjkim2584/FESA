#include "elements/incremental_elastic_point.h"

#include <cmath>
#include <memory>

#include "fesa/math/dense_factorization.h"
#include "model/model_storage.h"

namespace fesa {
namespace {
// IE2-IE4 buffers; identity is immutable. No operand/result shares storage.
enum Buffer {
  kI, kDifference, kSum, kInverseSum, kG, kGt, kStrain, kSpin, kB, kC,
  kRotation, kRotatedLeft, kStress, kInverseF, kPk1,
  kTwoIMinusG, kRotationPlusI, kDirection, kDg, kDgt, kDd, kDw,
  kDrRhs, kDr, kDs, kTerm, kOtherTerm, kProduct, kOtherProduct, kDp, kIso
};
static_assert(kIso + 1 == IncrementalElasticWorkspace::kMatrixCount);
bool Shape(const Matrix* matrix, int size) {
  return matrix && matrix->Rows() == size && matrix->Columns() == size;
}
StatusCode Range(const char* operation, Diagnostic* d) {
  return ReportFailure(StatusCode::kLibraryFailure,
                       DiagnosticStage::kComputation, operation, d);
}
StatusCode Geometry(const char* operation, TrialEvaluationFailure kind,
    TrialEvaluationFailure* failure, Diagnostic* d, bool measured = false,
    double determinant = 0, int vendor = 0) {
  *failure = kind;
  if (d) {
    d->Record(StatusCode::kSingular, DiagnosticStage::kComputation,
              operation, vendor);
    if (kind == TrialEvaluationFailure::kGeometry && measured) {
      d->geometry_quantity = GeometryQuantity::kDeformationJacobian;
      d->geometry_value_available = true;
      d->geometry_value = determinant;
      d->geometry_limit = 0;
    }
  }
  return StatusCode::kSingular;
}
StatusCode Factor(const Matrix* matrix,
    std::unique_ptr<DenseLuFactorization>* out, Diagnostic* d) {
  DenseLuFactorization* raw = nullptr;
  const auto status = DenseLuFactorization::Create(matrix, &raw, d);
  out->reset(raw);
  return status;
}
StatusCode Combine(const Matrix* a, double alpha, const Matrix* b, double beta,
                    Matrix* out, Diagnostic* d) {
  auto status = out->CopyFrom(a, d);
  if (status == StatusCode::kOk) status = out->Scale(alpha, d);
  if (status == StatusCode::kOk) status = out->AddScaled(beta, b, d);
  return status;
}
StatusCode Trace(const Matrix* input, IncrementalElasticWorkspace* workspace,
                  double* out, Diagnostic* d) {
  const int indices[] = {0, 4, 8};
  auto status = input->GatherTo(indices, 3, workspace->diagonal(), 0, d);
  if (status == StatusCode::kOk)
    status = workspace->diagonal()->Dot(workspace->ones(), out, d);
  return status;
}
StatusCode StoreSymmetricTensor(Matrix* tensor,
    IncrementalElasticWorkspace* workspace, Diagnostic* d) {
  // Explicit packed-six history roundtrip. Select the independent upper
  // entries; no averaging, clipping or constitutive correction is applied.
  const int tensor_indices[] = {0, 1, 2, 1, 4, 5, 2, 5, 8};
  auto status = tensor->GatherTo(tensor_indices, 9, workspace->flat(), 0, d);
  if (status != StatusCode::kOk) return status;
  for (int row = 0; row < 3; ++row) {
    const int indices[] = {3 * row, 3 * row + 1, 3 * row + 2};
    status = workspace->diagonal()->GatherFrom(workspace->flat(), indices,
                                               3, 0, d);
    if (status == StatusCode::kOk)
      status = tensor->SetRow(row, workspace->diagonal(), d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}
StatusCode AddElastic(const Matrix* strain, double lambda, double mu,
    IncrementalElasticWorkspace* workspace, Matrix* stress, Diagnostic* d) {
  double trace = 0;
  auto status = Trace(strain, workspace, &trace, d);
  if (status != StatusCode::kOk) return status;
  if (!std::isfinite(lambda * trace))
    return Range("Incremental isotropic trace range", d);
  auto* iso = workspace->matrix(kIso);
  status = iso->CopyFrom(workspace->matrix(kI), d);
  if (status == StatusCode::kOk) status = iso->Scale(lambda * trace, d);
  if (status == StatusCode::kOk) status = stress->AddScaled(2 * mu, strain, d);
  if (status == StatusCode::kOk) status = stress->AddScaled(1, iso, d);
  if (status == StatusCode::kOk) status = StoreSymmetricTensor(stress, workspace, d);
  return status;
}
StatusCode Derivative(const Matrix* sn, double lambda, double mu, double j,
    const DenseLuFactorization* rotation_factor,
    IncrementalElasticWorkspace* w, Matrix* out, Diagnostic* d) {
  const auto* identity = w->matrix(kI);
  auto status = Combine(identity, 2, w->matrix(kG), -1,
                         w->matrix(kTwoIMinusG), d);
  if (status == StatusCode::kOk)
    status = Combine(w->matrix(kRotation), 1, identity, 1,
                      w->matrix(kRotationPlusI), d);
  if (status != StatusCode::kOk) return status;
  const int flat_indices[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  for (int column = 0; column < 9; ++column) {
    // Unit direction dF=e_k e_B^T, IE4; scalar basis coefficients only.
    status = w->matrix(kDirection)->Scale(0, d);
    if (status == StatusCode::kOk)
      status = w->matrix(kDirection)->Set(column / 3, column % 3, 1, d);
    if (status == StatusCode::kOk)
      status = Matrix::Multiply(w->matrix(kTwoIMinusG), false,
          w->matrix(kDirection), false, w->matrix(kProduct), d);
    if (status == StatusCode::kOk)
      status = Matrix::Multiply(w->matrix(kProduct), false,
          w->matrix(kInverseSum), false, w->matrix(kDg), d);
    if (status == StatusCode::kOk)
      status = Matrix::Multiply(w->matrix(kDg), true, identity, false,
                                w->matrix(kDgt), d);
    if (status == StatusCode::kOk)
      status = Combine(w->matrix(kDg), .5, w->matrix(kDgt), .5,
                        w->matrix(kDd), d);
    if (status == StatusCode::kOk)
      status = Combine(w->matrix(kDg), .5, w->matrix(kDgt), -.5,
                        w->matrix(kDw), d);
    if (status == StatusCode::kOk)
      status = Matrix::Multiply(w->matrix(kDw), false,
          w->matrix(kRotationPlusI), false, w->matrix(kDrRhs), d);
    if (status == StatusCode::kOk) status = w->matrix(kDrRhs)->Scale(.5, d);
    if (status == StatusCode::kOk)
      status = rotation_factor->Solve(w->matrix(kDrRhs), false,
                                      w->matrix(kDr), d);
    // dσ=(dR)σn R^T + (Rσn)(dR)^T + C:dΔd, including both rotations.
    if (status == StatusCode::kOk)
      status = Matrix::Multiply(w->matrix(kDr), false, sn, false,
                                w->matrix(kProduct), d);
    if (status == StatusCode::kOk)
      status = Matrix::Multiply(w->matrix(kProduct), false,
          w->matrix(kRotation), true, w->matrix(kDs), d);
    if (status == StatusCode::kOk)
      status = Matrix::Multiply(w->matrix(kRotatedLeft), false,
          w->matrix(kDr), true, w->matrix(kOtherTerm), d);
    if (status == StatusCode::kOk)
      status = w->matrix(kDs)->AddScaled(1, w->matrix(kOtherTerm), d);
    if (status == StatusCode::kOk)
      status = AddElastic(w->matrix(kDd), lambda, mu, w, w->matrix(kDs), d);
    if (status != StatusCode::kOk) return status;
    // dP=tr(F^-1 H)P + J(dσ)F^-T - P H^T F^-T, IE4.
    status = Matrix::Multiply(w->matrix(kInverseF), false,
        w->matrix(kDirection), false, w->matrix(kProduct), d);
    double volume_trace = 0;
    if (status == StatusCode::kOk)
      status = Trace(w->matrix(kProduct), w, &volume_trace, d);
    if (status == StatusCode::kOk)
      status = w->matrix(kDp)->CopyFrom(w->matrix(kPk1), d);
    if (status == StatusCode::kOk)
      status = w->matrix(kDp)->Scale(volume_trace, d);
    if (status == StatusCode::kOk)
      status = Matrix::Multiply(w->matrix(kDs), false,
          w->matrix(kInverseF), true, w->matrix(kTerm), d);
    if (status == StatusCode::kOk)
      status = w->matrix(kDp)->AddScaled(j, w->matrix(kTerm), d);
    if (status == StatusCode::kOk)
      status = Matrix::Multiply(w->matrix(kPk1), false,
          w->matrix(kDirection), true, w->matrix(kProduct), d);
    if (status == StatusCode::kOk)
      status = Matrix::Multiply(w->matrix(kProduct), false,
          w->matrix(kInverseF), true, w->matrix(kOtherProduct), d);
    if (status == StatusCode::kOk)
      status = w->matrix(kDp)->AddScaled(-1, w->matrix(kOtherProduct), d);
    if (status == StatusCode::kOk)
      status = w->matrix(kDp)->GatherTo(flat_indices, 9, w->flat(), 0, d);
    if (status == StatusCode::kOk)
      status = w->derivative_rows()->SetRow(column, w->flat(), d);
    if (status != StatusCode::kOk) return status;
  }
  return Matrix::Multiply(w->derivative_rows(), true, w->identity9(), false,
                           out, d);
}
}  // namespace
StatusCode EvaluateIncrementalElasticPoint(const Matrix* fn, const Matrix* f,
    const Matrix* sn, const ElasticMaterial* material,
    IncrementalElasticWorkspace* w, Matrix* stress, Matrix* pk1, Matrix* dp,
    TrialEvaluationFailure* failure, Diagnostic* d) {
  if (failure) *failure = TrialEvaluationFailure::kNone;
  if (!failure || !w || !material || !Shape(fn, 3) || !Shape(f, 3) ||
      !Shape(sn, 3) || !Shape(stress, 3) || !Shape(pk1, 3) ||
      (dp && !Shape(dp, 9)))
    return model_storage::Invalid("Incremental point dimensions", d);
  const Matrix* inputs[] = {fn, f, sn};
  const Matrix* outputs[] = {stress, pk1, dp};
  for (const auto* input : inputs) {
    for (int i = 0; i < 9; ++i)
      if (!std::isfinite(input->Data()[i]))
        return model_storage::Invalid("Incremental point finite input", d);
    for (int i = 0; i < IncrementalElasticWorkspace::kMatrixCount; ++i)
      if (input == w->matrix(i))
        return model_storage::Invalid("Incremental point scratch alias", d);
  }
  for (int i = 0; i < 3; ++i) {
    if (!outputs[i]) continue;
    for (const auto* input : inputs)
      if (outputs[i] == input)
        return model_storage::Invalid("Incremental point input/output alias", d);
    for (int j = 0; j < i; ++j)
      if (outputs[i] == outputs[j])
        return model_storage::Invalid("Incremental point output alias", d);
    for (int j = 0; j < IncrementalElasticWorkspace::kMatrixCount; ++j)
      if (outputs[i] == w->matrix(j))
        return model_storage::Invalid("Incremental point scratch alias", d);
    if (outputs[i] == w->derivative_rows() || outputs[i] == w->identity9())
      return model_storage::Invalid("Incremental point derivative alias", d);
  }
  for (int i = 0; i < 3; ++i)
    for (int j = i + 1; j < 3; ++j)
      if (sn->Data()[3 * i + j] != sn->Data()[3 * j + i])
        return model_storage::Invalid("Incremental committed tensor symmetry", d);
  const double mu = material->young() / (2 * (1 + material->poisson()));
  const double lambda = material->young() * material->poisson() /
      ((1 + material->poisson()) * (1 - 2 * material->poisson()));
  if (!(mu > 0) || !std::isfinite(2 * mu) || !std::isfinite(lambda))
    return Range("Incremental elastic coefficient range", d);
  std::unique_ptr<DenseLuFactorization> committed_factor, candidate_factor,
      sum_factor, rotation_factor;
  auto status = Factor(fn, &committed_factor, d);
  if (status != StatusCode::kOk) return status;  // Committed faults stay fatal.
  double jn = 0, j = 0;
  status = committed_factor->Determinant(&jn, d);
  if (status != StatusCode::kOk) return status;
  if (!(jn > 0))
    return model_storage::Invalid("Incremental invalid committed determinant", d);
  status = Factor(f, &candidate_factor, d);
  if (status == StatusCode::kSingular)
    return Geometry("Incremental singular candidate F",
        TrialEvaluationFailure::kGeometry, failure, d, false, 0,
        d ? d->vendor_code : 0);
  if (status != StatusCode::kOk) return status;
  status = candidate_factor->Determinant(&j, d);
  if (status != StatusCode::kOk) return status;
  if (!(j > 0))
    return Geometry("Incremental nonpositive candidate J",
                    TrialEvaluationFailure::kGeometry, failure, d, true, j);
  const auto* identity = w->matrix(kI);
  status = Combine(f, 1, fn, -1, w->matrix(kDifference), d);
  if (status == StatusCode::kOk)
    status = Combine(f, 1, fn, 1, w->matrix(kSum), d);
  if (status == StatusCode::kOk)
    status = Factor(w->matrix(kSum), &sum_factor, d);
  if (status == StatusCode::kSingular)
    return Geometry("Incremental singular midpoint F+Fn",
        TrialEvaluationFailure::kGeometry, failure, d, false, 0,
        d ? d->vendor_code : 0);
  if (status != StatusCode::kOk) return status;
  status = sum_factor->Solve(identity, false, w->matrix(kInverseSum), d);
  if (status == StatusCode::kOk)
    status = Matrix::Multiply(w->matrix(kDifference), false,
        w->matrix(kInverseSum), false, w->matrix(kG), d);
  if (status == StatusCode::kOk) status = w->matrix(kG)->Scale(2, d);
  if (status == StatusCode::kOk)
    status = Matrix::Multiply(w->matrix(kG), true, identity, false,
                              w->matrix(kGt), d);
  if (status == StatusCode::kOk)
    status = Combine(w->matrix(kG), .5, w->matrix(kGt), .5,
                      w->matrix(kStrain), d);
  if (status == StatusCode::kOk)
    status = Combine(w->matrix(kG), .5, w->matrix(kGt), -.5,
                      w->matrix(kSpin), d);
  if (status == StatusCode::kOk)
    status = Combine(identity, 1, w->matrix(kSpin), -.5, w->matrix(kB), d);
  if (status == StatusCode::kOk)
    status = Combine(identity, 1, w->matrix(kSpin), .5, w->matrix(kC), d);
  if (status == StatusCode::kOk)
    status = Factor(w->matrix(kB), &rotation_factor, d);
  if (status == StatusCode::kSingular)
    return Geometry("Incremental singular rotation solve",
        TrialEvaluationFailure::kRotationBranch, failure, d, false, 0,
        d ? d->vendor_code : 0);
  if (status != StatusCode::kOk) return status;
  status = rotation_factor->Solve(w->matrix(kC), false,
                                   w->matrix(kRotation), d);
  if (status == StatusCode::kOk)
    status = Matrix::Multiply(w->matrix(kRotation), false, sn, false,
                              w->matrix(kRotatedLeft), d);
  if (status == StatusCode::kOk)
    status = Matrix::Multiply(w->matrix(kRotatedLeft), false,
        w->matrix(kRotation), true, w->matrix(kStress), d);
  if (status == StatusCode::kOk)
    status = AddElastic(w->matrix(kStrain), lambda, mu, w, w->matrix(kStress), d);
  if (status == StatusCode::kOk)
    status = candidate_factor->Solve(identity, false, w->matrix(kInverseF), d);
  if (status == StatusCode::kOk)
    status = Matrix::Multiply(w->matrix(kStress), false,
        w->matrix(kInverseF), true, w->matrix(kPk1), d);
  if (status == StatusCode::kOk) status = w->matrix(kPk1)->Scale(j, d);
  if (status == StatusCode::kOk && dp)
    status = Derivative(sn, lambda, mu, j, rotation_factor.get(), w, dp, d);
  if (status == StatusCode::kOk) status = stress->CopyFrom(w->matrix(kStress), d);
  if (status == StatusCode::kOk) status = pk1->CopyFrom(w->matrix(kPk1), d);
  return status;
}
}  // namespace fesa
