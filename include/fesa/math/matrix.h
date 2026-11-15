/// @file
/// @brief Fixed row-major MKL-backed dense matrix with explicit copying.
#ifndef FESA_MATH_MATRIX_H_
#define FESA_MATH_MATRIX_H_

#include "fesa/base/diagnostic.h"
#include "fesa/math/vector.h"

namespace fesa {

/// Owns a 64-byte aligned buffer; no copy/move. Single-writer object.
/// Coefficients are unitless here; consumers define their physical meaning.
class Matrix {
 public:
  /// Creates a zero matrix; caller owns *out and deletes it. Failure sets null.
  /// out is required and must not hold an existing allocation. d is optional.
  static StatusCode Create(int rows, int columns, Matrix** out, Diagnostic* d);
  ~Matrix();
  Matrix(const Matrix&) = delete;
  Matrix& operator=(const Matrix&) = delete;
  Matrix(Matrix&&) = delete;
  Matrix& operator=(Matrix&&) = delete;

  /// Borrows this/out/d; invalid index/null out preserves scalar output.
  StatusCode Get(int row, int column, double* out, Diagnostic* d) const;
  /// Sets one finite coefficient; invalid arguments preserve this.
  StatusCode Set(int row, int column, double value, Diagnostic* d);
  /// Copies borrowed input with matching dimensions; validation failure
  /// preserves this. Self-copy skips the backend.
  StatusCode CopyFrom(const Matrix* input, Diagnostic* d);
  /// Computes op(a)*op(b) into caller-owned out, which must be distinct.
  /// All pointers are borrowed; invalid arguments/backend preflight preserve
  /// out. Nonfinite BLAS output reports failure after mutation; discard out.
  static StatusCode Multiply(const Matrix* a, bool transpose_a,
                             const Matrix* b, bool transpose_b,
                             Matrix* out, Diagnostic* d);
  /// Computes op(a)*x into distinct caller-owned out; same failure contract.
  static StatusCode MultiplyVector(const Matrix* a, bool transpose_a,
                                   const Vector* x, Vector* out, Diagnostic* d);
  /// Computes log(A)=Q diag(log(mu)) Q^T for matching disjoint square matrices.
  /// Borrows all objects. Requires finite symmetric input within64eps*maxAbs;
  /// the upper triangle defines the matrix inside this roundoff tolerance.
  /// Required not_positive is cleared on entry; nonpositive eigenvalue sets
  /// it and returns kEquilibriumFailure, with no abs/clamp. Invalid inputs
  /// return kInvalidArgument; native eigensolver errors retain vendor code.
  /// All temporary buffers are per-call. Input is never changed. Failures
  /// preserve out except a backend failure during the final explicit copy.
  /// @see docs/Math/symmetric_log_contract.md.
  static StatusCode SymmetricPositiveDefiniteLog(const Matrix* input,
                                                 Matrix* out,
                                                 bool* not_positive,
                                                 Diagnostic* d);
  /// Borrowed read-only row-major data valid until destruction; never free it.
  const double* Data() const;
  /// Fixed positive dimensions.
  int Rows() const;
  /// Fixed positive column count.
  int Columns() const;
  /// Gathers zero-based row-major flattened indices into out[offset..].
  /// Repeated indices allowed. Preflight preserves out; backend failure
  /// invalidates it. count==0 skips backend. All pointers borrowed for call.
  StatusCode GatherTo(const int* indices, int count, Vector* out, int offset,
                      Diagnostic* d) const;

  /// Adds alpha*input in place; separate matching inputs required.
  /// Validation preserves this; backend/nonfinite-output failure discards it.
  StatusCode AddScaled(double alpha, const Matrix* input, Diagnostic* d);
  /// Scales in place; finite scalar/input required, same failure contract.
  StatusCode Scale(double alpha, Diagnostic* d);
  /// Writes row sums to matching out without changing this. Validation
  /// preserves output; backend failure invalidates output. Pointers borrowed
  /// for call.
  StatusCode RowSums(Vector* out, Diagnostic* d) const;
  /// Copies a matching finite row vector using external numerical copying.
  /// Invalid row/shape preserves this; backend failure discards this.
  StatusCode SetRow(int row, const Vector* values, Diagnostic* d);
  /// Copies a nonempty rectangular block of a separate source. Invalid range,
  /// alias or nonfinite source preserves this; backend failure discards this.
  StatusCode CopyBlockFrom(const Matrix* source, int source_row,
                           int source_column, int row_count, int column_count,
                           int target_row, int target_column, Diagnostic* d);
  /// Computes alpha*op(a)*op(b)+beta*out on disjoint matching buffers.
  /// Finite scalars/inputs required; beta=0 ignores previous output contents.
  /// Validation preserves output; backend failure invalidates output.
  static StatusCode MultiplyAdd(const Matrix* a, bool transpose_a,
                                const Matrix* b, bool transpose_b, double alpha,
                                double beta, Matrix* out, Diagnostic* d);
  /// Overwrites square out with diagonal and zero off-diagonal entries.
  /// Shapes must match. Validation preserves out; backend failure discards it.
  static StatusCode Diagonal(const Vector* diagonal, Matrix* out,
                             Diagnostic* d);

 private:
  friend class FiniteValidationPeer;
  /// Reads the proof or scans unknown data; never changes shared const state.
  bool IsFinite() const;
  /// Always scans the entire owned output and updates its proof after mutation.
  bool VerifyFinite();
  /// Whole-buffer finite proof; zero creation proves it, false means unknown.
  bool finite_verified_ = true;
  friend class DenseLuFactorization;
  friend class DenseQrFactorization;
  Matrix(int rows, int columns, double* data);
  int rows_;
  int columns_;
  double* data_;  ///< Owned aligned allocation; released by destructor.
};

}  // namespace fesa

#endif  // FESA_MATH_MATRIX_H_
