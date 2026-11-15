#include "fesa/math/matrix.h"

#include <cmath>
#include <limits>
#include <memory>

#include "adapters/mkl/dense_ops.h"

namespace fesa {
namespace {
StatusCode Invalid(Diagnostic* d) {
  return ReportFailure(StatusCode::kInvalidArgument,
                       DiagnosticStage::kValidation,
                       "Symmetric matrix logarithm input/shape/alias", d);
}
StatusCode Backend(const char* operation, int native, Diagnostic* d) {
  const bool allocation = dense_ops::LapackAllocationFailure(native);
  const auto code =
      allocation ? StatusCode::kAllocationFailure : StatusCode::kLibraryFailure;
  if (d)
    d->Record(code,
              allocation ? DiagnosticStage::kAllocation
                         : DiagnosticStage::kComputation,
              operation, native);
  return code;
}
StatusCode CreateSquare(int size, std::unique_ptr<Matrix>* out, Diagnostic* d) {
  Matrix* raw = nullptr;
  const auto status = Matrix::Create(size, size, &raw, d);
  out->reset(raw);
  return status;
}
StatusCode CheckSymmetric(const Matrix* input, Diagnostic* d) {
  const int size = input->Rows();
  const double* data = input->Data();
  if (!dense_ops::BeginOperation())
    return Backend("Symmetric logarithm norm preflight", 0, d);
  const int maximum = dense_ops::MaximumAbsIndex(size * size, data);
  const double scale = std::abs(data[maximum]);
  const double tolerance = 64 * std::numeric_limits<double>::epsilon() * scale;
  for (int i = 0; i < size; ++i)
    for (int j = i + 1; j < size; ++j)
      if (std::abs(data[i * size + j] - data[j * size + i]) > tolerance)
        return Invalid(d);
  return StatusCode::kOk;
}
}  // namespace
StatusCode Matrix::SymmetricPositiveDefiniteLog(const Matrix* input,
                                                Matrix* out, bool* not_positive,
                                                Diagnostic* d) {
  if (!not_positive) return Invalid(d);
  *not_positive = false;
  if (!input || !out || input == out || input->rows_ != input->columns_ ||
      out->rows_ != input->rows_ || out->columns_ != input->columns_ ||
      !dense_ops::AllFinite(input->rows_ * input->columns_, input->data_))
    return Invalid(d);
  const int size = input->rows_;
  auto status = CheckSymmetric(input, d);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<Matrix> q;
  status = CreateSquare(size, &q, d);
  if (status != StatusCode::kOk) return status;
  status = q->CopyFrom(input, d);
  if (status != StatusCode::kOk) return status;
  Vector* raw = nullptr;
  status = Vector::Create(size, &raw, d);
  std::unique_ptr<Vector> eigenvalues(raw);
  if (status != StatusCode::kOk) return status;
  if (!dense_ops::BeginOperation())
    return Backend("Symmetric eigensolver preflight", 0, d);
  q->finite_verified_ = false;
  eigenvalues->finite_verified_ = false;
  const int native =
      dense_ops::SymmetricEigen(size, q->data_, eigenvalues->data_);
  if (native != 0) return Backend("Symmetric eigensolver", native, d);
  if (!dense_ops::AllFinite(size * size, q->data_) ||
      !dense_ops::AllFinite(size, eigenvalues->data_))
    return Backend("Nonfinite symmetric eigensystem", 0, d);
  for (int i = 0; i < size; ++i) {
    if (eigenvalues->Data()[i] <= 0) {
      *not_positive = true;
      return ReportFailure(StatusCode::kEquilibriumFailure,
                           DiagnosticStage::kComputation,
                           "Matrix logarithm nonpositive eigenvalue", d);
    }
  }
  if (!dense_ops::BeginOperation())
    return Backend("Positive spectral logarithm preflight", 0, d);
  eigenvalues->finite_verified_ = false;
  dense_ops::NaturalLog(size, eigenvalues->data_, eigenvalues->data_);
  if (!dense_ops::AllFinite(size, eigenvalues->data_))
    return Backend("Nonfinite spectral logarithm", 0, d);
  std::unique_ptr<Matrix> diagonal, product, result;
  status = CreateSquare(size, &diagonal, d);
  if (status == StatusCode::kOk) status = CreateSquare(size, &product, d);
  if (status == StatusCode::kOk) status = CreateSquare(size, &result, d);
  if (status == StatusCode::kOk)
    status = Diagonal(eigenvalues.get(), diagonal.get(), d);
  if (status == StatusCode::kOk)
    status = Multiply(q.get(), false, diagonal.get(), false, product.get(), d);
  if (status == StatusCode::kOk)
    status = Multiply(product.get(), false, q.get(), true, result.get(), d);
  return status == StatusCode::kOk ? out->CopyFrom(result.get(), d) : status;
}
}  // namespace fesa
