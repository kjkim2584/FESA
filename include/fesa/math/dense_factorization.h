/// @file
/// @brief Owned dense factorizations with vendor-independent solve contracts.
#ifndef FESA_MATH_DENSE_FACTORIZATION_H_
#define FESA_MATH_DENSE_FACTORIZATION_H_

#include "fesa/math/matrix.h"

namespace fesa {

/// Owns an explicitly copied square LU and pivot metadata. No cross-call cache.
class DenseLuFactorization {
 public:
  /// Finite square input required, borrowed only for creation and unchanged.
  /// Caller owns successful factor; failure nulls *out and frees resources.
  static StatusCode Create(const Matrix* input, DenseLuFactorization** out,
                           Diagnostic* d);
  ~DenseLuFactorization();
  DenseLuFactorization(const DenseLuFactorization&) = delete;
  DenseLuFactorization& operator=(const DenseLuFactorization&) = delete;
  DenseLuFactorization(DenseLuFactorization&&) = delete;
  DenseLuFactorization& operator=(DenseLuFactorization&&) = delete;
  /// Solves A*x=rhs or A^T*x=rhs for all RHS columns. Inputs/out must be
  /// separate and matching. Validation preserves out; backend failure discards
  /// out. Factors and RHS stay unchanged; pointers borrowed for this call.
  StatusCode Solve(const Matrix* rhs, bool transpose, Matrix* out,
                   Diagnostic* d) const;
  /// Writes signed determinant. Null/overflow/underflow/backend error preserves
  /// the scalar output. Uses LU and external diagonal/log reduction operations.
  StatusCode Determinant(double* out, Diagnostic* d) const;

 private:
  /// Adopts already validated factor storage and pivot metadata.
  DenseLuFactorization(Matrix* factors, int* pivots);
  Matrix* factors_;  ///< Owned copied LU, unchanged by solves.
  int* pivots_;      ///< Owned one-based pivot metadata, backend-neutral int.
};

/// Owns copied m>=n QR and implicit Q reflectors; full numerical column rank.
/// Rank policy and units follow the dense math contract in
/// docs/SolidElement/architecture.md.
class DenseQrFactorization {
 public:
  /// Input must be finite, m>=n and pass the documented R-diagonal threshold.
  /// Input is unchanged. Caller owns success; failure nulls output and cleans.
  static StatusCode Create(const Matrix* input, DenseQrFactorization** out,
                           Diagnostic* d);
  ~DenseQrFactorization();
  DenseQrFactorization(const DenseQrFactorization&) = delete;
  DenseQrFactorization& operator=(const DenseQrFactorization&) = delete;
  DenseQrFactorization(DenseQrFactorization&&) = delete;
  DenseQrFactorization& operator=(DenseQrFactorization&&) = delete;
  /// Minimizes ||A*x-rhs||2 for multiple RHS columns using existing factors.
  /// RHS is m by r and out n by r, disjoint. Validation preserves out;
  /// backend failure discards out. Factor/RHS inputs remain unchanged.
  StatusCode SolveLeastSquares(const Matrix* rhs, Matrix* out,
                               Diagnostic* d) const;

 private:
  /// Adopts already validated factor storage and reflector coefficients.
  DenseQrFactorization(Matrix* factors, Vector* tau);
  Matrix* factors_;  ///< Owned QR, unchanged by solves.
  Vector* tau_;      ///< Owned implicit reflector coefficients.
};

}  // namespace fesa

#endif  // FESA_MATH_DENSE_FACTORIZATION_H_
