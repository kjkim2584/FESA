#include "adapters/mkl/sparse_ops.h"

#include <limits>

#include "mkl.h"
namespace fesa {
namespace sparse_ops {
static_assert(sizeof(MKL_INT) == sizeof(int), "LP64 required");
bool ValidSizes(std::size_t size, std::size_t count) {
  const auto lp64 = static_cast<std::size_t>(std::numeric_limits<int>::max());
  const auto bytes = std::numeric_limits<std::size_t>::max();
  return size > 0 && size <= lp64 && count <= lp64 &&
         size < bytes / sizeof(int) && count <= bytes / 16;
}
bool ValidCsrStorageSize(std::size_t size, std::size_t count) {
  constexpr std::size_t kEntryBytes = sizeof(double) + sizeof(int);
  constexpr std::size_t kRowBytes = 3 * sizeof(int);
  const auto bytes = std::numeric_limits<std::size_t>::max();
  if (size == 0 || size > (bytes - sizeof(int)) / kRowBytes) return false;
  const std::size_t row_bytes = kRowBytes * size + sizeof(int);
  return count <= (bytes - row_bytes) / kEntryBytes;
}
bool ValidAdditionSizes(std::size_t size, std::size_t count_a,
                        std::size_t count_b) {
  if (!ValidSizes(size, count_a) || !ValidSizes(size, count_b)) return false;
  const auto lp64 = static_cast<std::size_t>(std::numeric_limits<int>::max());
  if (count_a > lp64 - count_b) return false;
  return ValidCsrStorageSize(size, count_a + count_b);
}
int Convert(int size, int count, const int* rows, const int* columns,
            const double* values, void** out) {
  *out = nullptr;
  sparse_matrix_t coo = nullptr;
  sparse_status_t status = mkl_sparse_d_create_coo(
      &coo, SPARSE_INDEX_BASE_ZERO, size, size, count, const_cast<int*>(rows),
      const_cast<int*>(columns), const_cast<double*>(values));
  if (status != SPARSE_STATUS_SUCCESS) return static_cast<int>(status);
  sparse_matrix_t csr = nullptr;
  status = mkl_sparse_convert_csr(coo, SPARSE_OPERATION_NON_TRANSPOSE, &csr);
  const sparse_status_t released = mkl_sparse_destroy(coo);
  if (status != SPARSE_STATUS_SUCCESS) return static_cast<int>(status);
  if (released != SPARSE_STATUS_SUCCESS) {
    mkl_sparse_destroy(csr);
    return static_cast<int>(released);
  }
  status = mkl_sparse_order(csr);
  if (status != SPARSE_STATUS_SUCCESS) {
    mkl_sparse_destroy(csr);
    return static_cast<int>(status);
  }
  *out = csr;
  return 0;
}
int Add(void* a, void* b, void** out) {
  *out = nullptr;
  sparse_matrix_t result = nullptr;
  sparse_status_t status = mkl_sparse_d_add(
      SPARSE_OPERATION_NON_TRANSPOSE, static_cast<sparse_matrix_t>(a), 1,
      static_cast<sparse_matrix_t>(b), &result);
  if (status != SPARSE_STATUS_SUCCESS) return static_cast<int>(status);
  status = mkl_sparse_order(result);
  if (status != SPARSE_STATUS_SUCCESS) {
    mkl_sparse_destroy(result);
    return static_cast<int>(status);
  }
  *out = result;
  return 0;
}
int Export(void* handle, int* size, int** starts, int** ends, int** columns,
           double** values) {
  sparse_index_base_t base = SPARSE_INDEX_BASE_ONE;
  int cols = 0;
  const int status = static_cast<int>(
      mkl_sparse_d_export_csr(static_cast<sparse_matrix_t>(handle), &base, size,
                              &cols, starts, ends, columns, values));
  if (status != 0) return status;
  return base == SPARSE_INDEX_BASE_ZERO && *size == cols ? 0 : 3;
}
void Release(void* handle) {
  if (handle) mkl_sparse_destroy(static_cast<sparse_matrix_t>(handle));
}
}  // namespace sparse_ops
}  // namespace fesa
