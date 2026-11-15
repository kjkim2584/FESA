/// @file
/// @brief Fixed-size MKL-backed dense vector with explicit copying.
#ifndef FESA_MATH_VECTOR_H_
#define FESA_MATH_VECTOR_H_

#include "fesa/base/diagnostic.h"

namespace fesa {

/// Owns its buffer; copy/move are forbidden. Single-writer object.
/// Dimensions/components are unitless; callers define physical conventions.
class Vector {
 public:
  /// Creates a zero vector. Caller owns *out and deletes it. Failure sets null.
  /// out is required and must not hold an existing allocation. d is optional.
  static StatusCode Create(int size, Vector** out, Diagnostic* d);
  ~Vector();
  Vector(const Vector&) = delete;
  Vector& operator=(const Vector&) = delete;
  Vector(Vector&&) = delete;
  Vector& operator=(Vector&&) = delete;

  /// Borrows this/out/d. Invalid index/null out preserves *out.
  StatusCode Get(int index, double* out, Diagnostic* d) const;
  /// Sets one finite coefficient; failed validation preserves the vector.
  StatusCode Set(int index, double value, Diagnostic* d);
  /// Explicit copy from borrowed input; size mismatch preserves this.
  /// Self-copy succeeds without a backend call.
  StatusCode CopyFrom(const Vector* input, Diagnostic* d);
  /// this <- this + alpha*x. Borrowed x must be distinct and same size.
  /// Invalid arguments/backend preflight preserve this; numeric overflow
  /// reports failure after mutation, so the caller must discard this result.
  StatusCode AddScaled(double alpha, const Vector* x, Diagnostic* d);
  /// Explicit in-place scaling; same failure contract as AddScaled.
  StatusCode Scale(double alpha, Diagnostic* d);
  /// Writes dot product on success; any failure preserves scalar output.
  StatusCode Dot(const Vector* x, double* out, Diagnostic* d) const;
  /// Writes Euclidean norm on success; any failure preserves scalar output.
  StatusCode Norm2(double* out, Diagnostic* d) const;
  /// Gathers count indexed values into this[offset..]. Duplicates allowed.
  /// All indices are zero-based, checked before mutation. Inputs must be
  /// distinct. count==0 is a no-op; no backend call. Failure invalidates output
  /// only after a numerical backend operation has begun.
  StatusCode GatherFrom(const Vector* source, const int* indices, int count,
                        int offset, Diagnostic* d);
  /// Adds source[0..count) to selected destinations, summing duplicate indices
  /// through a library selector. Inputs distinct; same validation contract.
  StatusCode ScatterAdd(const Vector* source, const int* indices, int count,
                        Diagnostic* d);
  /// Sum of absolute components; scalar output preserved on failure.
  StatusCode AbsoluteSum(double* out, Diagnostic* d) const;
  /// Largest absolute component and first corresponding index via BLAS.
  StatusCode MaximumAbs(double* out, int* index, Diagnostic* d) const;
  /// Read-only borrowed contiguous storage, valid until this is destroyed.
  /// Caller must not free or cast away const from this pointer.
  const double* Data() const;
  /// Fixed positive component count.
  int Size() const;

 private:
  friend class FiniteValidationPeer;
  /// Reads the proof or scans unknown data; never changes shared const state.
  bool IsFinite() const;
  /// Always scans the entire owned output and updates its proof after mutation.
  bool VerifyFinite();
  /// Whole-buffer finite proof; zero creation proves it, false means unknown.
  bool finite_verified_ = true;
  // Same math module; permits checked matrix-vector output.
  friend class Matrix;
  friend class DenseLuFactorization;
  friend class DenseQrFactorization;
  friend class PardisoSolver;
  Vector(int size, double* data);
  int size_;
  double* data_;  ///< Owned aligned allocation, released by destructor.
};

}  // namespace fesa

#endif  // FESA_MATH_VECTOR_H_
