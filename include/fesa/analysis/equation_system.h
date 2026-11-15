/// @file
/// @brief Immutable reduced linear equation system.
#ifndef FESA_ANALYSIS_EQUATION_SYSTEM_H_
#define FESA_ANALYSIS_EQUATION_SYSTEM_H_
#include "fesa/math/sparse_matrix.h"
namespace fesa {
/// Mathematical operator type, explicitly selected by the analysis caller.
/// Symmetric types require upper CSR; kGeneral requires full general CSR.
enum class MatrixType {
  kSymmetricPositiveDefinite,
  kSymmetricIndefinite,
  kGeneral
};
/// Owns the reduced operator and RHS. The caller defines whether the solution
/// is a total displacement (Kff, pf-Kfc*uc) or a Newton correction (tangent,
/// residual). Storage mode is explicit in SparseMatrix; no full K is stored.
class EquationSystem {
 public:
  /// Transfers matrix/rhs ownership ONLY on success. Nulls *out on failure.
  /// Both must be positive-dimensional and have equal size. Rejects invalid
  /// type/storage combinations before ownership transfer; inputs are preserved.
  static StatusCode Create(SparseMatrix* matrix, Vector* rhs,
                           MatrixType type,
                           EquationSystem** out, Diagnostic* d);
  ~EquationSystem();
  EquationSystem(const EquationSystem&) = delete;
  EquationSystem& operator=(const EquationSystem&) = delete;
  EquationSystem(EquationSystem&&) = delete;
  EquationSystem& operator=(EquationSystem&&) = delete;
  /// Borrowed read-only reduced storage, valid until destruction.
  const SparseMatrix* Stiffness() const;
  const Vector* Rhs() const;
  /// Caller-selected mathematical type, never inferred from CSR or entries.
  MatrixType Type() const;

 private:
  EquationSystem() = default;
  SparseMatrix* matrix_ = nullptr;  ///< Owned.
  Vector* rhs_ = nullptr;           ///< Owned.
  MatrixType type_ = MatrixType::kGeneral;  ///< Set by the successful factory.
};
}  // namespace fesa
#endif  // FESA_ANALYSIS_EQUATION_SYSTEM_H_
