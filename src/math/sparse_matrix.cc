#include "fesa/math/sparse_matrix.h"

#include <cmath>
#include <memory>
#include <new>

#include "adapters/mkl/dense_ops.h"
#include "adapters/mkl/sparse_ops.h"
#include "model/model_storage.h"
namespace fesa {
namespace {
StatusCode Backend(int status, Diagnostic* d) {
  const StatusCode code = status == 2 ? StatusCode::kAllocationFailure
                                      : StatusCode::kLibraryFailure;
  if (d)
    d->Record(code,
              status == 2 ? DiagnosticStage::kAllocation
                          : DiagnosticStage::kComputation,
              "SparseMatrix", status);
  return code;
}
}  // namespace
SparseMatrix::~SparseMatrix() {
  sparse_ops::Release(handle_);
  delete[] offsets_;
}
StatusCode SparseMatrix::Export(Diagnostic* d) {
  int size = 0;
  int* starts = nullptr;
  int* ends = nullptr;
  int* columns = nullptr;
  double* values = nullptr;
  const int status =
      sparse_ops::Export(handle_, &size, &starts, &ends, &columns, &values);
  if (status != 0) return Backend(status, d);
  if (size != size_ || !starts || !ends) return Backend(3, d);
  offsets_ = new (std::nothrow) int[static_cast<std::size_t>(size_) + 1];
  if (!offsets_) return model_storage::AllocationFailure("SparseMatrix", d);
  int previous_end = 0;
  for (int row = 0; row < size_; ++row) {
    if (starts[row] != previous_end || ends[row] < starts[row])
      return Backend(3, d);
    offsets_[row] = starts[row];
    int previous_column =
        storage_ == SparseStorage::kSymmetricUpper ? row - 1 : -1;
    for (int j = starts[row]; j < ends[row]; ++j) {
      if (!columns || !values || columns[j] <= previous_column ||
          columns[j] >= size_ || !std::isfinite(values[j]))
        return Backend(3, d);
      previous_column = columns[j];
    }
    previous_end = ends[row];
  }
  offsets_[size_] = previous_end;
  nonzeros_ = previous_end;
  columns_ = columns;
  values_ = values;
  return StatusCode::kOk;
}
StatusCode SparseMatrix::Create(int size, int count, const int* rows,
                                const int* columns, const Vector* values,
                                SparseMatrix** out, Diagnostic* d) {
  return CreateWithStorage(size, count, rows, columns, values,
                           SparseStorage::kSymmetricUpper, out, d);
}
StatusCode SparseMatrix::CreateGeneral(int size, int count, const int* rows,
                                       const int* columns, const Vector* values,
                                       SparseMatrix** out, Diagnostic* d) {
  return CreateWithStorage(size, count, rows, columns, values,
                           SparseStorage::kGeneral, out, d);
}
StatusCode SparseMatrix::CreateWithStorage(int size, int count, const int* rows,
                                           const int* columns,
                                           const Vector* values,
                                           SparseStorage storage,
                                           SparseMatrix** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("SparseMatrix::Create", d);
  *out = nullptr;
  if (size <= 0 || count < 0 || !sparse_ops::ValidSizes(size, count) ||
      (count > 0 && (!rows || !columns || !values || count > values->Size())))
    return model_storage::Invalid("SparseMatrix::Create", d);
  for (int i = 0; i < count; ++i)
    if (rows[i] < 0 || rows[i] >= size || columns[i] < 0 ||
        (storage == SparseStorage::kSymmetricUpper && columns[i] < rows[i]) ||
        columns[i] >= size)
      return model_storage::Invalid("SparseMatrix::Create", d);
  std::unique_ptr<SparseMatrix> result(new (std::nothrow) SparseMatrix);
  if (!result || !dense_ops::AllowObjectAllocation())
    return model_storage::AllocationFailure("SparseMatrix::Create", d);
  result->size_ = size;
  result->storage_ = storage;
  std::unique_ptr<int[]> diagonal;
  std::unique_ptr<Vector> zero;
  if (count == 0) {
    diagonal.reset(new (std::nothrow) int[size]);
    Vector* buffer = nullptr;
    const StatusCode status = Vector::Create(size, &buffer, d);
    zero.reset(buffer);
    if (status != StatusCode::kOk) return status;
    if (!diagonal) return model_storage::AllocationFailure("SparseMatrix", d);
    for (int i = 0; i < size; ++i) diagonal[i] = i;
    rows = diagonal.get();
    columns = diagonal.get();
    values = zero.get();
    count = size;
  }
  if (!dense_ops::BeginOperation()) return Backend(4, d);
  const int native = sparse_ops::Convert(size, count, rows, columns,
                                         values->Data(), &result->handle_);
  if (native != 0) return Backend(native, d);
  const StatusCode status = result->Export(d);
  if (status != StatusCode::kOk) return status;
  *out = result.release();
  return StatusCode::kOk;
}
StatusCode SparseMatrix::Add(const SparseMatrix* a, const SparseMatrix* b,
                             SparseMatrix** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("SparseMatrix::Add", d);
  *out = nullptr;
  if (!a || !b || a->size_ != b->size_ || a->storage_ != b->storage_ ||
      !sparse_ops::ValidAdditionSizes(a->size_, a->nonzeros_, b->nonzeros_))
    return model_storage::Invalid("SparseMatrix::Add", d);
  std::unique_ptr<SparseMatrix> result(new (std::nothrow) SparseMatrix);
  if (!result || !dense_ops::AllowObjectAllocation())
    return model_storage::AllocationFailure("SparseMatrix::Add", d);
  result->size_ = a->size_;
  result->storage_ = a->storage_;
  if (!dense_ops::BeginOperation()) return Backend(4, d);
  const int native = sparse_ops::Add(a->handle_, b->handle_, &result->handle_);
  if (native != 0) return Backend(native, d);
  const StatusCode status = result->Export(d);
  if (status != StatusCode::kOk) return status;
  *out = result.release();
  return StatusCode::kOk;
}
int SparseMatrix::Size() const { return size_; }
int SparseMatrix::Nonzeros() const { return nonzeros_; }
SparseStorage SparseMatrix::Storage() const { return storage_; }
const int* SparseMatrix::RowOffsets() const { return offsets_; }
const int* SparseMatrix::Columns() const { return columns_; }
const double* SparseMatrix::Values() const { return values_; }
}  // namespace fesa
