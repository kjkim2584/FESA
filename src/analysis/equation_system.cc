#include "fesa/analysis/equation_system.h"

#include <new>

#include "model/model_storage.h"
namespace fesa {
StatusCode EquationSystem::Create(SparseMatrix* matrix, Vector* rhs,
                                  MatrixType type,
                                  EquationSystem** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("EquationSystem::Create", d);
  *out = nullptr;
  if (!matrix || !rhs || matrix->Size() != rhs->Size())
    return model_storage::Invalid("EquationSystem::Create", d);
  const bool symmetric = type == MatrixType::kSymmetricPositiveDefinite ||
                         type == MatrixType::kSymmetricIndefinite;
  if ((!symmetric && type != MatrixType::kGeneral) ||
      (symmetric && matrix->Storage() != SparseStorage::kSymmetricUpper) ||
      (!symmetric && matrix->Storage() != SparseStorage::kGeneral)) {
    return model_storage::Invalid("EquationSystem::Create.type", d);
  }
  auto* result = new (std::nothrow) EquationSystem;
  if (!result) return model_storage::AllocationFailure("EquationSystem", d);
  result->matrix_ = matrix;
  result->rhs_ = rhs;
  result->type_ = type;
  *out = result;
  return StatusCode::kOk;
}
EquationSystem::~EquationSystem() {
  delete matrix_;
  delete rhs_;
}
const SparseMatrix* EquationSystem::Stiffness() const { return matrix_; }
const Vector* EquationSystem::Rhs() const { return rhs_; }
MatrixType EquationSystem::Type() const { return type_; }
}  // namespace fesa
