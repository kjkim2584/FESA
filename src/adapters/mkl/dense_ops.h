/// @file
/// @brief Internal dense BLAS boundary and thread-local test seams.
#ifndef FESA_ADAPTERS_MKL_DENSE_OPS_H_
#define FESA_ADAPTERS_MKL_DENSE_OPS_H_

#include <cstddef>
#include <cstdint>

namespace fesa {
namespace dense_ops {

/// Checks positive LP64 element count and byte size before allocation.
bool ValidDimensions(int rows, int columns);
/// Returns owned 64-byte aligned zero memory; Release must free it.
double* Allocate(std::size_t count);
/// Releases only adapter-owned memory, including null.
void Release(double* data);
/// Preflight for operations: consumes the thread-local failure seam.
bool BeginOperation();
/// Read-only binary64 classification, not a numerical vector operation.
/// Returns false at the first NaN/infinity. Does not change FP flags/masks;
/// a nonpositive count does not access data. Positive count requires data.
bool AllFinite(int count, const double* data);
/// BLAS operations borrow their buffers for the call. Sizes are prevalidated.
void Copy(int count, const double* input, double* out);
/// Computes op(a)*op(b), row-major, into distinct out.
void Multiply(int rows, int columns, int inner, bool transpose_a,
              bool transpose_b, const double* a, int stride_a,
              const double* b, int stride_b, double* out);
/// Computes op(a)*x, row-major, into distinct out.
void MultiplyVector(int rows, int columns, bool transpose,
                    const double* a, const double* x, double* out);
/// Computes y <- y + alpha*x on distinct buffers.
void AddScaled(int count, double alpha, const double* x, double* y);
/// Computes x <- alpha*x in place.
void Scale(int count, double alpha, double* x);
/// Returns the dot product of borrowed buffers.
double Dot(int count, const double* a, const double* b);
/// Returns the Euclidean norm of a borrowed buffer.
double Norm2(int count, const double* a);
/// Indexed read, indices validated by common math; disjoint buffers.
void Gather(int count, const double* source, const int* indices, double* out);
/// Absolute sum and maximum index via BLAS, positive count only.
double AbsoluteSum(int count, const double* source);
int MaximumAbsIndex(int count, const double* source);
/// General selector transpose accumulation. Returns 0 on success, -1 for
/// adapter allocation failure, otherwise the Sparse BLAS status. Only
/// successfully created handles are released; arrays live through release.
/// The accumulation uses one local MKL thread for reproducible duplicate-index
/// reduction, restoring the caller's previous local setting on every return.
int ScatterAdd(int size, int count, const int* indices, const double* source,
               double* out);

/// Strided numerical copy, disjoint buffers with positive validated strides.
void CopyStrided(int count, const double* input, int input_stride, double* out,
                 int output_stride);
/// General row-major gemm accumulation; output disjoint from input arrays.
void MultiplyAdd(int rows, int columns, int inner, bool transpose_a,
                 bool transpose_b, double alpha, double beta, const double* a,
                 int stride_a, const double* b, int stride_b, double* out);
/// Overwrites all entries with zero without reading previous values.
int SetZero(int rows, int columns, double* out);
/// Row-major LAPACK factor/solve routines; return native info. Dimensions and
/// disjoint array lengths are validated by common math before entry.
int LuFactor(int size, double* factors, int* pivots);
/// Solves with existing LU into separate RHS work storage.
int LuSolve(int size, int rhs_count, bool transpose, const double* factors,
            const int* pivots, double* rhs);
/// Stores QR and implicit reflectors into separate owned arrays.
int QrFactor(int rows, int columns, double* factors, double* tau);
/// Applies implicit Q^T to m by rhs_count work storage.
int ApplyQt(int rows, int columns, int rhs_count, const double* factors,
            const double* tau, double* rhs);
/// Solves upper triangular factors, with the supplied row stride.
int SolveUpper(int size, int rhs_count, const double* factors, int stride,
               double* rhs);
/// Classifies native LAPACKE work/transpose allocation failures internally.
bool LapackAllocationFailure(int info);
/// Row-major symmetric eigenvectors in columns, upper triangle input, ascending
/// eigenvalues. Overwrites input; native LAPACKE info, no vendor public types.
int SymmetricEigen(int size, double* eigenvectors, double* eigenvalues);
/// VM natural logarithm of positive input. Supports in-place input/output.
void NaturalLog(int count, const double* input, double* out);
/// VM log(abs(input)); disjoint input/output, VM-supported in-place log stage.
void LogAbs(int count, const double* input, double* out);

/// Internal tests only: fail next allocation/operation on this thread.
void FailNextAllocation();
/// Internal tests only: fail object creation after numerical buffer allocation.
void FailNextObjectAllocation();
/// Object factories consume this seam before nothrow new.
bool AllowObjectAllocation();
/// Internal tests only: fail before the next numerical operation mutates data.
void FailNextOperation();
/// Internal tests only: backend copy calls, saturating at INT_MAX (a lower
/// bound thereafter); saturation never skips or changes the actual BLAS copy.
int CopyCallCount();
/// Internal tests only: sets count in [0, INT_MAX] for boundary tests.
void SetCopyCallCountForTest(int count);
/// Internal tests only: all numerical allocation attempts on this thread,
/// including failed and subsequently freed allocations; saturates at UINT64_MAX.
std::uint64_t AllocationCallCount();
/// Internal tests only: current number of unreleased adapter allocations.
int LiveAllocationCount();

}  // namespace dense_ops
}  // namespace fesa

#endif  // FESA_ADAPTERS_MKL_DENSE_OPS_H_
