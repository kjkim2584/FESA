#include "adapters/mkl/dense_ops.h"

#include <cmath>
#include <cstring>
#include <limits>
#include <memory>
#include <new>

// Keep facade validation; see docs/Performance/mkl_direct_call.md.
// Only this adapter opts into Intel kernels; threaded fallback is unchanged.
#define MKL_DIRECT_CALL
#include "mkl.h"

namespace fesa {
namespace dense_ops {
namespace {

thread_local bool fail_allocation = false;
thread_local bool fail_operation = false;
thread_local bool fail_object_allocation = false;
thread_local int copy_calls = 0;
thread_local int live_allocations = 0;
thread_local std::uint64_t allocation_calls = 0;
static_assert(sizeof(MKL_INT) == 4, "Only oneMKL LP64 is supported");

}  // namespace

bool ValidDimensions(int rows, int columns) {
  if (rows <= 0 || columns <= 0) return false;
  const auto r = static_cast<std::size_t>(rows);
  const auto c = static_cast<std::size_t>(columns);
  return r <= static_cast<std::size_t>(std::numeric_limits<MKL_INT>::max()) / c
      && r <= std::numeric_limits<std::size_t>::max() / sizeof(double) / c;
}

double* Allocate(std::size_t count) {
  if (allocation_calls != std::numeric_limits<std::uint64_t>::max())
    ++allocation_calls;
  if (fail_allocation) {
    fail_allocation = false;
    return nullptr;
  }
  double* data = static_cast<double*>(mkl_calloc(count, sizeof(double), 64));
  if (data != nullptr) ++live_allocations;
  return data;
}

void Release(double* data) {
  if (data != nullptr) --live_allocations;
  mkl_free(data);
}

bool BeginOperation() {
  const bool failed = fail_operation;
  fail_operation = false;
  return !failed;
}

bool AllFinite(int count, const double* data) {
  static_assert(sizeof(double) == sizeof(std::uint64_t) &&
                    std::numeric_limits<double>::is_iec559,
                "Finite classification requires IEEE binary64");
  for (int i = 0; i < count; ++i) {
    // Read-only representation classification, not a numerical array copy.
    std::uint64_t bits = 0;
    std::memcpy(&bits, data + i, sizeof(bits));
    if ((bits & UINT64_C(0x7ff0000000000000)) == UINT64_C(0x7ff0000000000000))
      return false;
  }
  return true;
}

void Copy(int count, const double* input, double* out) {
  if (copy_calls != std::numeric_limits<int>::max()) ++copy_calls;
  cblas_dcopy(count, input, 1, out, 1);
}

void Multiply(int rows, int columns, int inner, bool transpose_a,
              bool transpose_b, const double* a, int stride_a,
              const double* b, int stride_b, double* out) {
  cblas_dgemm(CblasRowMajor, transpose_a ? CblasTrans : CblasNoTrans,
              transpose_b ? CblasTrans : CblasNoTrans, rows, columns, inner,
              1.0, a, stride_a, b, stride_b, 0.0, out, columns);
}

void MultiplyVector(int rows, int columns, bool transpose,
                    const double* a, const double* x, double* out) {
  cblas_dgemv(CblasRowMajor, transpose ? CblasTrans : CblasNoTrans,
              rows, columns, 1.0, a, columns, x, 1, 0.0, out, 1);
}

void AddScaled(int count, double alpha, const double* x, double* y) {
  cblas_daxpy(count, alpha, x, 1, y, 1);
}

void Scale(int count, double alpha, double* x) {
  cblas_dscal(count, alpha, x, 1);
}

double Dot(int count, const double* a, const double* b) {
  // The direct-call macro takes the address of its stride arguments.
  constexpr int kStride = 1;
  return cblas_ddot(count, a, kStride, b, kStride);
}

double Norm2(int count, const double* a) { return cblas_dnrm2(count, a, 1); }

void Gather(int count, const double* source, const int* indices, double* out) {
  cblas_dgthr(count, source, out, indices);
}
double AbsoluteSum(int count, const double* source) {
  return cblas_dasum(count, source, 1);
}
int MaximumAbsIndex(int count, const double* source) {
  return static_cast<int>(cblas_idamax(count, source, 1));
}
int ScatterAdd(int size, int count, const int* indices, const double* source,
               double* out) {
  const std::size_t rows = static_cast<std::size_t>(count) + 1;
  if (rows > std::numeric_limits<std::size_t>::max() / sizeof(int)) return -1;
  std::unique_ptr<int[]> offsets(new (std::nothrow) int[rows]);
  double* ones = Allocate(count);
  if (!offsets || !ones) {
    Release(ones);
    return -1;
  }
  for (int i = 0; i < count; ++i) {
    offsets[i] = i;
    ones[i] = 1;  // Selector coefficient definition, not numerical summation.
  }
  offsets[count] = count;
  sparse_matrix_t selector = nullptr;
  const sparse_status_t created = mkl_sparse_d_create_csr(
      &selector, SPARSE_INDEX_BASE_ZERO, count, size, offsets.get(),
      offsets.get() + 1, const_cast<int*>(indices), ones);
  int status = static_cast<int>(created);
  if (created == SPARSE_STATUS_SUCCESS) {
    matrix_descr descriptor = {};
    descriptor.type = SPARSE_MATRIX_TYPE_GENERAL;
    // Keep duplicate-index reduction independent of the caller's thread count.
    // See docs/Linear static/parallel_reproducibility.md.
    const int previous_threads = mkl_set_num_threads_local(1);
    status = static_cast<int>(mkl_sparse_d_mv(
        SPARSE_OPERATION_TRANSPOSE, 1, selector, descriptor, source, 1, out));
    mkl_set_num_threads_local(previous_threads);
    const int released = static_cast<int>(mkl_sparse_destroy(selector));
    if (status == 0) status = released;
  }
  Release(ones);
  return status;
}

void CopyStrided(int count, const double* input, int input_stride, double* out,
                 int output_stride) {
  if (copy_calls != std::numeric_limits<int>::max()) ++copy_calls;
  cblas_dcopy(count, input, input_stride, out, output_stride);
}
void MultiplyAdd(int rows, int columns, int inner, bool transpose_a,
                 bool transpose_b, double alpha, double beta, const double* a,
                 int stride_a, const double* b, int stride_b, double* out) {
  cblas_dgemm(CblasRowMajor, transpose_a ? CblasTrans : CblasNoTrans,
              transpose_b ? CblasTrans : CblasNoTrans, rows, columns, inner,
              alpha, a, stride_a, b, stride_b, beta, out, columns);
}
int SetZero(int rows, int columns, double* out) {
  return LAPACKE_dlaset(LAPACK_ROW_MAJOR, 'A', rows, columns, 0, 0, out,
                        columns);
}
int LuFactor(int size, double* factors, int* pivots) {
  return LAPACKE_dgetrf(LAPACK_ROW_MAJOR, size, size, factors, size, pivots);
}
int LuSolve(int size, int rhs_count, bool transpose, const double* factors,
            const int* pivots, double* rhs) {
  return LAPACKE_dgetrs(LAPACK_ROW_MAJOR, transpose ? 'T' : 'N', size,
                        rhs_count, factors, size, pivots, rhs, rhs_count);
}
int QrFactor(int rows, int columns, double* factors, double* tau) {
  return LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, rows, columns, factors, columns, tau);
}
int ApplyQt(int rows, int columns, int rhs_count, const double* factors,
            const double* tau, double* rhs) {
  return LAPACKE_dormqr(LAPACK_ROW_MAJOR, 'L', 'T', rows, rhs_count, columns,
                        factors, columns, tau, rhs, rhs_count);
}
int SolveUpper(int size, int rhs_count, const double* factors, int stride,
               double* rhs) {
  return LAPACKE_dtrtrs(LAPACK_ROW_MAJOR, 'U', 'N', 'N', size, rhs_count,
                        factors, stride, rhs, rhs_count);
}
bool LapackAllocationFailure(int info) {
  return info == LAPACK_WORK_MEMORY_ERROR ||
         info == LAPACK_TRANSPOSE_MEMORY_ERROR;
}
void LogAbs(int count, const double* input, double* out) {
  vdAbs(count, input, out);
  vdLn(count, out, out);  // VM explicitly permits exact in-place arrays.
}

void FailNextAllocation() { fail_allocation = true; }
void FailNextObjectAllocation() { fail_object_allocation = true; }
bool AllowObjectAllocation() {
  const bool failed = fail_object_allocation;
  fail_object_allocation = false;
  return !failed;
}
void FailNextOperation() { fail_operation = true; }
int CopyCallCount() { return copy_calls; }
void SetCopyCallCountForTest(int count) { copy_calls = count; }
int LiveAllocationCount() { return live_allocations; }
std::uint64_t AllocationCallCount() { return allocation_calls; }

int SymmetricEigen(int size, double* eigenvectors, double* eigenvalues) {
  return LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', size, eigenvectors, size,
                       eigenvalues);
}
void NaturalLog(int count, const double* input, double* out) {
  vdLn(count, input, out);
}
}  // namespace dense_ops
}  // namespace fesa
