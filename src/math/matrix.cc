#include "fesa/math/matrix.h"

#include <cmath>
#include <memory>
#include <new>

#include "adapters/mkl/dense_ops.h"

namespace fesa {
namespace {

StatusCode Invalid(const char* operation, Diagnostic* d) {
  return ReportFailure(StatusCode::kInvalidArgument,
                       DiagnosticStage::kValidation, operation, d);
}

StatusCode BackendFailure(const char* operation, Diagnostic* d) {
  return ReportFailure(StatusCode::kLibraryFailure,
                       DiagnosticStage::kComputation, operation, d);
}

}  // namespace

StatusCode Matrix::GatherTo(const int* indices, int count, Vector* out,
                            int offset, Diagnostic* d) const {
  if (!out || count < 0 || offset < 0 || offset > out->size_ ||
      count > out->size_ - offset || (count > 0 && !indices))
    return Invalid("Matrix::GatherTo", d);
  for (int i = 0; i < count; ++i)
    if (indices[i] < 0 || indices[i] >= rows_ * columns_)
      return Invalid("Matrix::GatherTo", d);
  if (count == 0) return StatusCode::kOk;
  if (!dense_ops::BeginOperation())
    return BackendFailure("Matrix::GatherTo", d);
  const bool verified = out->finite_verified_ && finite_verified_;
  out->finite_verified_ = false;
  dense_ops::Gather(count, data_, indices, out->data_ + offset);
  out->finite_verified_ = verified;
  return StatusCode::kOk;
}

Matrix::Matrix(int rows, int columns, double* data)
    : rows_(rows), columns_(columns), data_(data) {}
Matrix::~Matrix() { dense_ops::Release(data_); }
bool Matrix::IsFinite() const {
  return finite_verified_ || dense_ops::AllFinite(rows_ * columns_, data_);
}
bool Matrix::VerifyFinite() {
  finite_verified_ = dense_ops::AllFinite(rows_ * columns_, data_);
  return finite_verified_;
}

StatusCode Matrix::Create(int rows, int columns, Matrix** out, Diagnostic* d) {
  if (out == nullptr) return Invalid("Matrix::Create", d);
  *out = nullptr;
  if (!dense_ops::ValidDimensions(rows, columns)) {
    return Invalid("Matrix::Create", d);
  }
  double* data = dense_ops::Allocate(static_cast<std::size_t>(rows) * columns);
  if (data != nullptr) {
    if (dense_ops::AllowObjectAllocation()) {
      *out = new (std::nothrow) Matrix(rows, columns, data);
    }
    if (*out != nullptr) return StatusCode::kOk;
    dense_ops::Release(data);
  }
  return ReportFailure(StatusCode::kAllocationFailure,
                       DiagnosticStage::kAllocation, "Matrix::Create", d);
}

StatusCode Matrix::Get(int row, int column, double* out, Diagnostic* d) const {
  if (row < 0 || row >= rows_ || column < 0 || column >= columns_
      || out == nullptr) return Invalid("Matrix::Get", d);
  *out = data_[row * columns_ + column];
  return StatusCode::kOk;
}

StatusCode Matrix::Set(int row, int column, double value, Diagnostic* d) {
  if (row < 0 || row >= rows_ || column < 0 || column >= columns_
      || !std::isfinite(value)) return Invalid("Matrix::Set", d);
  data_[row * columns_ + column] = value;
  return StatusCode::kOk;
}

StatusCode Matrix::CopyFrom(const Matrix* input, Diagnostic* d) {
  if (input == nullptr || input->rows_ != rows_
      || input->columns_ != columns_) {
    return Invalid("Matrix::CopyFrom", d);
  }
  if (input == this) return StatusCode::kOk;
  if (!dense_ops::BeginOperation()) {
    return BackendFailure("Matrix::CopyFrom", d);
  }
  finite_verified_ = false;
  dense_ops::Copy(rows_ * columns_, input->data_, data_);
  finite_verified_ = input->finite_verified_;
  return StatusCode::kOk;
}

StatusCode Matrix::Multiply(const Matrix* a, bool transpose_a,
                            const Matrix* b, bool transpose_b,
                            Matrix* out, Diagnostic* d) {
  if (a == nullptr || b == nullptr || out == nullptr || out == a || out == b) {
    return Invalid("Matrix::Multiply", d);
  }
  const int rows = transpose_a ? a->columns_ : a->rows_;
  const int inner = transpose_a ? a->rows_ : a->columns_;
  const int b_rows = transpose_b ? b->columns_ : b->rows_;
  const int columns = transpose_b ? b->rows_ : b->columns_;
  if (inner != b_rows || out->rows_ != rows || out->columns_ != columns) {
    return Invalid("Matrix::Multiply", d);
  }
  if (!dense_ops::BeginOperation()) {
    return BackendFailure("Matrix::Multiply", d);
  }
  out->finite_verified_ = false;
  dense_ops::Multiply(rows, columns, inner, transpose_a, transpose_b,
                      a->data_, a->columns_, b->data_, b->columns_, out->data_);
  return out->VerifyFinite() ? StatusCode::kOk
                             : BackendFailure("Matrix::Multiply", d);
}

StatusCode Matrix::MultiplyVector(const Matrix* a, bool transpose_a,
                                  const Vector* x, Vector* out, Diagnostic* d) {
  if (a == nullptr || x == nullptr || out == nullptr || x == out) {
    return Invalid("Matrix::MultiplyVector", d);
  }
  const int rows = transpose_a ? a->columns_ : a->rows_;
  const int inner = transpose_a ? a->rows_ : a->columns_;
  if (x->size_ != inner || out->size_ != rows) {
    return Invalid("Matrix::MultiplyVector", d);
  }
  if (!dense_ops::BeginOperation()) {
    return BackendFailure("Matrix::MultiplyVector", d);
  }
  out->finite_verified_ = false;
  dense_ops::MultiplyVector(a->rows_, a->columns_, transpose_a,
                            a->data_, x->data_, out->data_);
  return out->VerifyFinite() ? StatusCode::kOk
                             : BackendFailure("Matrix::MultiplyVector", d);
}

const double* Matrix::Data() const { return data_; }
int Matrix::Rows() const { return rows_; }
int Matrix::Columns() const { return columns_; }

StatusCode Matrix::AddScaled(double alpha, const Matrix* input, Diagnostic* d) {
  if (!input || input == this || input->rows_ != rows_ ||
      input->columns_ != columns_ || !std::isfinite(alpha) || !IsFinite() ||
      !input->IsFinite())
    return Invalid("Matrix::AddScaled", d);
  if (!dense_ops::BeginOperation())
    return BackendFailure("Matrix::AddScaled", d);
  finite_verified_ = false;
  dense_ops::AddScaled(rows_ * columns_, alpha, input->data_, data_);
  return VerifyFinite() ? StatusCode::kOk
                        : BackendFailure("Matrix::AddScaled", d);
}
StatusCode Matrix::Scale(double alpha, Diagnostic* d) {
  if (!std::isfinite(alpha) || !IsFinite()) return Invalid("Matrix::Scale", d);
  if (!dense_ops::BeginOperation()) return BackendFailure("Matrix::Scale", d);
  finite_verified_ = false;
  dense_ops::Scale(rows_ * columns_, alpha, data_);
  return VerifyFinite() ? StatusCode::kOk : BackendFailure("Matrix::Scale", d);
}
StatusCode Matrix::RowSums(Vector* out, Diagnostic* d) const {
  if (!out || out->Size() != rows_ || !IsFinite())
    return Invalid("Matrix::RowSums", d);
  Vector* raw = nullptr;
  auto status = Vector::Create(columns_, &raw, d);
  std::unique_ptr<Vector> ones(raw);
  if (status != StatusCode::kOk) return status;
  for (int i = 0; i < columns_; ++i) {
    status = ones->Set(i, 1, d);  // Constant selector coefficient definition.
    if (status != StatusCode::kOk) return status;
  }
  return MultiplyVector(this, false, ones.get(), out, d);
}
StatusCode Matrix::SetRow(int row, const Vector* values, Diagnostic* d) {
  if (row < 0 || row >= rows_ || !values || values->Size() != columns_ ||
      !values->IsFinite())
    return Invalid("Matrix::SetRow", d);
  if (!dense_ops::BeginOperation()) return BackendFailure("Matrix::SetRow", d);
  const bool verified = finite_verified_;
  finite_verified_ = false;
  dense_ops::Copy(columns_, values->Data(), data_ + row * columns_);
  finite_verified_ = verified;
  return StatusCode::kOk;
}
StatusCode Matrix::CopyBlockFrom(const Matrix* source, int source_row,
                                 int source_column, int row_count,
                                 int column_count, int target_row,
                                 int target_column, Diagnostic* d) {
  if (!source || source == this || source_row < 0 || source_column < 0 ||
      target_row < 0 || target_column < 0 || row_count <= 0 ||
      column_count <= 0 || row_count > source->rows_ ||
      column_count > source->columns_ || row_count > rows_ ||
      column_count > columns_ || source_row > source->rows_ - row_count ||
      source_column > source->columns_ - column_count ||
      target_row > rows_ - row_count ||
      target_column > columns_ - column_count || !source->IsFinite())
    return Invalid("Matrix::CopyBlockFrom", d);
  if (!dense_ops::BeginOperation())
    return BackendFailure("Matrix::CopyBlockFrom", d);
  const bool verified = finite_verified_;
  finite_verified_ = false;
  for (int row = 0; row < row_count; ++row)
    dense_ops::Copy(
        column_count,
        source->data_ + (source_row + row) * source->columns_ + source_column,
        data_ + (target_row + row) * columns_ + target_column);
  finite_verified_ = verified;
  return StatusCode::kOk;
}
StatusCode Matrix::MultiplyAdd(const Matrix* a, bool transpose_a,
                               const Matrix* b, bool transpose_b, double alpha,
                               double beta, Matrix* out, Diagnostic* d) {
  if (!a || !b || !out || out == a || out == b || !std::isfinite(alpha) ||
      !std::isfinite(beta))
    return Invalid("Matrix::MultiplyAdd", d);
  const int rows = transpose_a ? a->columns_ : a->rows_;
  const int inner = transpose_a ? a->rows_ : a->columns_;
  const int b_rows = transpose_b ? b->columns_ : b->rows_;
  const int columns = transpose_b ? b->rows_ : b->columns_;
  if (inner != b_rows || out->rows_ != rows || out->columns_ != columns ||
      !a->IsFinite() || !b->IsFinite() || (beta != 0 && !out->IsFinite()))
    return Invalid("Matrix::MultiplyAdd", d);
  if (!dense_ops::BeginOperation())
    return BackendFailure("Matrix::MultiplyAdd", d);
  out->finite_verified_ = false;
  dense_ops::MultiplyAdd(rows, columns, inner, transpose_a, transpose_b, alpha,
                         beta, a->data_, a->columns_, b->data_, b->columns_,
                         out->data_);
  return out->VerifyFinite() ? StatusCode::kOk
                             : BackendFailure("Matrix::MultiplyAdd", d);
}
StatusCode Matrix::Diagonal(const Vector* diagonal, Matrix* out,
                            Diagnostic* d) {
  if (!diagonal || !out || out->rows_ != diagonal->Size() ||
      out->columns_ != diagonal->Size() || !diagonal->IsFinite())
    return Invalid("Matrix::Diagonal", d);
  if (!dense_ops::BeginOperation())
    return BackendFailure("Matrix::Diagonal", d);
  out->finite_verified_ = false;
  const int native = dense_ops::SetZero(out->rows_, out->columns_, out->data_);
  if (native != 0) {
    if (d)
      d->Record(StatusCode::kLibraryFailure, DiagnosticStage::kComputation,
                "Matrix::Diagonal zero initialization", native);
    return StatusCode::kLibraryFailure;
  }
  dense_ops::CopyStrided(diagonal->Size(), diagonal->Data(), 1, out->data_,
                         out->columns_ + 1);
  out->finite_verified_ = true;
  return StatusCode::kOk;
}

}  // namespace fesa
