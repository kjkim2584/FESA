#include "elements/mitc4_nonlinear_ops.h"

#include <cmath>
#include <initializer_list>
#include <limits>
#include <memory>

#include "fesa/math/dense_factorization.h"
#include "model/model_storage.h"

namespace fesa {
namespace {
StatusCode CreateMatrix(int rows, int columns, std::unique_ptr<Matrix>* owner,
                        Diagnostic* d) {
  Matrix* raw = nullptr;
  const auto status = Matrix::Create(rows, columns, &raw, d);
  owner->reset(raw);
  return status;
}

StatusCode GeometryFailure(TrialEvaluationFailure* failure, Diagnostic* d) {
  *failure = TrialEvaluationFailure::kGeometry;
  return ReportFailure(StatusCode::kEquilibriumFailure,
                       DiagnosticStage::kComputation,
                       "MITC4 trial collapsed or inverted volume basis", d);
}

StatusCode NumericalFailure(const char* operation, Diagnostic* d) {
  return ReportFailure(StatusCode::kLibraryFailure,
                       DiagnosticStage::kComputation, operation, d);
}

StatusCode Skew(const Vector* vector, Matrix* out, Diagnostic* d) {
  const double* v = vector->Data();
  const double coefficients[9] = {0.0,   -v[2], v[1], v[2], 0.0,
                                  -v[0], -v[1], v[0], 0.0};
  for (int i = 0; i < 9; ++i) {
    const auto status = out->Set(i / 3, i % 3, coefficients[i], d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}

StatusCode CrossWithScratch(const Vector* a, const Vector* b, Vector* out,
                            Matrix* skew, Diagnostic* d) {
  const auto status = Skew(a, skew, d);
  if (status != StatusCode::kOk) return status;
  return Matrix::MultiplyVector(skew, false, b, out, d);
}

StatusCode ValidateVolumeArguments(const Vector* const* initial,
                                   const Vector* const* current,
                                   double* initial_determinant,
                                   double* jacobian,
                                   TrialEvaluationFailure* failure,
                                   Diagnostic* d) {
  if (!failure) return model_storage::Invalid("MITC4 volume failure flag", d);
  *failure = TrialEvaluationFailure::kNone;
  if (!initial || !current || !initial_determinant || !jacobian ||
      initial_determinant == jacobian)
    return model_storage::Invalid("MITC4 volume basis inputs/outputs", d);
  for (int axis = 0; axis < 3; ++axis)
    for (const Vector* vector : {initial[axis], current[axis]}) {
      if (!vector || vector->Size() != 3)
        return model_storage::Invalid("MITC4 volume basis vector shape", d);
      for (int i = 0; i < 3; ++i)
        if (!std::isfinite(vector->Data()[i]))
          return model_storage::Invalid("MITC4 volume basis finite input", d);
    }

  return StatusCode::kOk;
}

}  // namespace
StatusCode CrossMitc4Vectors(const Vector* a, const Vector* b, Vector* out,
                             Diagnostic* d) {
  if (!a || !b || !out || a->Size() != 3 || b->Size() != 3 ||
      out->Size() != 3 || out == a || out == b)
    return model_storage::Invalid("MITC4 cross input/shape/alias", d);
  for (int i = 0; i < 3; ++i)
    if (!std::isfinite(a->Data()[i]) || !std::isfinite(b->Data()[i]))
      return model_storage::Invalid("MITC4 cross finite input", d);
  std::unique_ptr<Matrix> skew;
  auto status = CreateMatrix(3, 3, &skew, d);
  if (status != StatusCode::kOk) return status;
  return CrossWithScratch(a, b, out, skew.get(), d);
}
void Mitc4VolumeWorkspace::Clear() {
  delete cross_;
  cross_ = nullptr;
  delete skew_;
  skew_ = nullptr;
  for (int axis = 0; axis < 3; ++axis) {
    delete normalized_[axis];
    normalized_[axis] = nullptr;
  }
}

Mitc4VolumeWorkspace::~Mitc4VolumeWorkspace() { Clear(); }

StatusCode Mitc4VolumeWorkspace::Initialize(Diagnostic* d) {
  if (skew_) return StatusCode::kOk;
  auto status = Vector::Create(3, &cross_, d);
  for (int axis = 0; axis < 3 && status == StatusCode::kOk; ++axis)
    status = Vector::Create(3, &normalized_[axis], d);
  if (status == StatusCode::kOk) status = Matrix::Create(3, 3, &skew_, d);
  if (status != StatusCode::kOk) Clear();
  return status;
}

StatusCode Mitc4VolumeWorkspace::Validate(
    const Vector* const* initial, const Vector* const* current,
    double* initial_determinant, double* jacobian,
    TrialEvaluationFailure* failure, Diagnostic* d) {
  auto status = ValidateVolumeArguments(initial, current, initial_determinant,
                                        jacobian, failure, d);
  if (status != StatusCode::kOk) return status;
  if (!skew_)
    return model_storage::Invalid("MITC4 volume workspace not initialized", d);
  // Check the reference before any current-geometry classification. Reducing
  // a load increment cannot repair an unrepresentable reference volume.
  status = CrossWithScratch(initial[0], initial[1], cross_, skew_, d);
  if (status != StatusCode::kOk) return status;
  double det_initial = 0.0;
  status = cross_->Dot(initial[2], &det_initial, d);
  if (status != StatusCode::kOk) return status;
  if (!(det_initial > 0.0))
    return ReportFailure(StatusCode::kInputError, DiagnosticStage::kValidation,
                         "MITC4 initial volume basis range", d);
  for (int axis = 0; axis < 3; ++axis) {
    double norm = 0.0;
    status = current[axis]->Norm2(&norm, d);
    if (status != StatusCode::kOk) return status;
    if (!(norm > 0.0) || !std::isfinite(1.0 / norm))
      return GeometryFailure(failure, d);
    status = normalized_[axis]->CopyFrom(current[axis], d);
    if (status != StatusCode::kOk) return status;
    status = normalized_[axis]->Scale(1.0 / norm, d);
    if (status != StatusCode::kOk) return status;
  }
  status = CrossWithScratch(normalized_[0], normalized_[1], cross_, skew_, d);
  if (status != StatusCode::kOk) return status;
  double triple = 0.0;
  status = cross_->Dot(normalized_[2], &triple, d);
  if (status != StatusCode::kOk) return status;
  if (!(triple > 256.0 * std::numeric_limits<double>::epsilon()))
    return GeometryFailure(failure, d);
  status = CrossWithScratch(current[0], current[1], cross_, skew_, d);
  if (status != StatusCode::kOk) return status;
  double det_current = 0.0;
  status = cross_->Dot(current[2], &det_current, d);
  if (status != StatusCode::kOk) return status;
  const double j = det_current / det_initial;
  if (!(j > 0.0) || !std::isfinite(j)) return GeometryFailure(failure, d);
  *initial_determinant = det_initial;
  *jacobian = j;
  return StatusCode::kOk;
}

StatusCode ValidateMitc4VolumeBasis(const Vector* const* initial,
                                    const Vector* const* current,
                                    double* initial_determinant,
                                    double* jacobian,
                                    TrialEvaluationFailure* failure,
                                    Diagnostic* d) {
  auto status = ValidateVolumeArguments(initial, current, initial_determinant,
                                        jacobian, failure, d);
  if (status != StatusCode::kOk) return status;
  Mitc4VolumeWorkspace workspace;
  status = workspace.Initialize(d);
  if (status != StatusCode::kOk) return status;
  return workspace.Validate(initial, current, initial_determinant, jacobian,
                             failure, d);
}

StatusCode BuildMitc4NonlinearCartesianTransform(const Matrix* initial_basis,
                                                 const Matrix* frame_transpose,
                                                 Matrix* p, Diagnostic* d) {
  if (!initial_basis || !frame_transpose || !p || initial_basis->Rows() != 3 ||
      initial_basis->Columns() != 3 || frame_transpose->Rows() != 3 ||
      frame_transpose->Columns() != 3 || p->Rows() != 5 || p->Columns() != 5 ||
      p == initial_basis || p == frame_transpose)
    return model_storage::Invalid("MITC4 Cartesian transform shape/alias", d);
  for (const Matrix* matrix : {initial_basis, frame_transpose})
    for (int i = 0; i < 9; ++i)
      if (!std::isfinite(matrix->Data()[i]))
        return model_storage::Invalid("MITC4 Cartesian transform finite input",
                                      d);

  // initial_basis stores G_r,G_s,G_zeta as rows. Solve its transpose:
  // H=G_columns^-1*Q0^T. Keep all components of the validated initial data.
  DenseLuFactorization* raw = nullptr;
  auto status = DenseLuFactorization::Create(initial_basis, &raw, d);
  std::unique_ptr<DenseLuFactorization> factor(raw);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<Matrix> h;
  status = CreateMatrix(3, 3, &h, d);
  if (status != StatusCode::kOk) return status;
  status = factor->Solve(frame_transpose, true, h.get(), d);
  if (status != StatusCode::kOk) return status;
  return BuildMitc4NonlinearTransformFromH(h.get(), p, d);
}
StatusCode BuildMitc4NonlinearTransformFromH(const Matrix* h, Matrix* p,
                                             Diagnostic* d) {
  if (!h || !p || h->Rows() != 3 || h->Columns() != 3 || p->Rows() != 5 ||
      p->Columns() != 5 || h == p)
    return model_storage::Invalid("MITC4 H transform shape/alias", d);
  for (int i = 0; i < 9; ++i)
    if (!std::isfinite(h->Data()[i]))
      return model_storage::Invalid("MITC4 H transform finite input", d);
  constexpr int kA[5] = {0, 1, 0, 1, 2};
  constexpr int kB[5] = {0, 1, 1, 2, 0};
  const double* values = h->Data();
  StatusCode status = StatusCode::kOk;
  for (int row = 0; row < 5; ++row) {
    const int a = kA[row];
    const int b = kB[row];
    const double scale = a == b ? 1.0 : 2.0;
    // E_local=H^T E_cov H; E_zeta,zeta is excluded from shell energy.
    const double coefficients[5] = {
        scale * values[a] * values[b], scale * values[3 + a] * values[3 + b],
        scale * (values[a] * values[3 + b] + values[3 + a] * values[b]),
        scale * (values[a] * values[6 + b] + values[6 + a] * values[b]),
        scale *
            (values[3 + a] * values[6 + b] + values[6 + a] * values[3 + b])};
    for (int column = 0; column < 5; ++column) {
      if (!std::isfinite(coefficients[column]))
        return NumericalFailure("MITC4 trial Cartesian transform range", d);
      status = p->Set(row, column, coefficients[column], d);
      if (status != StatusCode::kOk) return status;
    }
  }
  return StatusCode::kOk;
}

StatusCode BuildMitc4NonlinearTyingSelector(double r, double s, Matrix* out,
                                            Diagnostic* d) {
  if (!out || out->Rows() != 5 || out->Columns() != 7 || !std::isfinite(r) ||
      !std::isfinite(s) || std::abs(r) > 1 || std::abs(s) > 1)
    return model_storage::Invalid("MITC4 tying selector input", d);
  const double coefficients[35] = {
      1,           0,           0, 0, 0, 0, 0, 0, 1, 0,           0,          0,
      0,           0,           0, 0, 1, 0, 0, 0, 0, 0,           0,          0,
      (1 - s) / 2, (1 + s) / 2, 0, 0, 0, 0, 0, 0, 0, (1 - r) / 2, (1 + r) / 2};
  for (int i = 0; i < 35; ++i) {
    const auto status = out->Set(i / 7, i % 7, coefficients[i], d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}

}  // namespace fesa
