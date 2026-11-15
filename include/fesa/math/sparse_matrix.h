/// @file
/// @brief Immutable zero-based sparse matrix with explicit storage semantics.
#ifndef FESA_MATH_SPARSE_MATRIX_H_
#define FESA_MATH_SPARSE_MATRIX_H_
#include "fesa/math/vector.h"
namespace fesa {
/// Numerical interpretation of the stored entries; never inferred from values.
enum class SparseStorage { kSymmetricUpper, kGeneral };
/// Owns canonical CSR: sorted unique columns, explicit zeros retained.
/// Values follow the consumer's units. No copy/move or mutable storage access.
class SparseMatrix {
 public:
  /// Converts borrowed COO upper entries using library duplicate summation.
  /// count may be zero (zero diagonal matrix); size must be positive.
  /// values has at least count entries. Caller deletes successful *out;
  /// failure nulls it, preserving inputs. Inputs need live only for this call.
  static StatusCode Create(int size, int count, const int* rows,
                           const int* columns, const Vector* values,
                           SparseMatrix** out, Diagnostic* d);
  /// Converts full general COO; lower and upper entries remain independent.
  /// Same size, count, ownership and failure contract as Create. No symmetry
  /// expansion or diagonal regularization. Empty input retains a zero diagonal.
  static StatusCode CreateGeneral(int size, int count, const int* rows,
                                  const int* columns, const Vector* values,
                                  SparseMatrix** out, Diagnostic* d);
  /// Library sum of equally sized, equally stored borrowed matrices.
  /// Same ownership contract; mixed storage is kInvalidArgument.
  /// Before allocation/backend work, rejects kInvalidArgument if the sum of
  /// stored input counts exceeds LP64 or the managed CSR byte bound overflows
  /// size_t. This conservative bound can reject very large overlapping inputs
  /// whose exact union would fit. Rejection nulls *out and preserves inputs.
  static StatusCode Add(const SparseMatrix* a, const SparseMatrix* b,
                        SparseMatrix** out, Diagnostic* d);
  ~SparseMatrix();
  SparseMatrix(const SparseMatrix&) = delete;
  SparseMatrix& operator=(const SparseMatrix&) = delete;
  SparseMatrix(SparseMatrix&&) = delete;
  SparseMatrix& operator=(SparseMatrix&&) = delete;
  /// Positive order and stored count, respectively.
  int Size() const;
  int Nonzeros() const;
  /// Interpretation fixed at construction, preserved by Add.
  SparseStorage Storage() const;
  /// Borrowed zero-based n+1 row offsets; lifetime is this object.
  const int* RowOffsets() const;
  /// Borrowed sorted columns and numerical values; never mutate/free them.
  const int* Columns() const;
  const double* Values() const;

 private:
  SparseMatrix() = default;
  /// Shared factory; failure nulls *out and preserves all borrowed inputs.
  static StatusCode CreateWithStorage(int size, int count, const int* rows,
                                      const int* columns, const Vector* values,
                                      SparseStorage storage, SparseMatrix** out,
                                      Diagnostic* d);
  /// Validates/exports owned backend CSR; failure leaves this destructible.
  StatusCode Export(Diagnostic* d);
  int size_ = 0;
  int nonzeros_ = 0;
  SparseStorage storage_ = SparseStorage::kSymmetricUpper;
  void* handle_ = nullptr;          ///< Owned opaque adapter allocation.
  int* offsets_ = nullptr;          ///< Owned three-array CSR index projection.
  const int* columns_ = nullptr;    ///< Borrowed from handle_.
  const double* values_ = nullptr;  ///< Borrowed from handle_.
};
}  // namespace fesa
#endif  // FESA_MATH_SPARSE_MATRIX_H_
