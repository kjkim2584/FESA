#include "elements/incremental_elastic_workspace.h"
#include <new>
#include <memory>
#include "model/model_storage.h"
#include "elements/incremental_hexa8_workspace.h"
namespace fesa {
StatusCode IncrementalElasticWorkspace::Create(IncrementalElasticWorkspace** out,
                                               Diagnostic* d) {
  return CreateForSolid(4, 1, out, d);
}
StatusCode IncrementalElasticWorkspace::CreateForSolid(int node_count,
    int point_count, IncrementalElasticWorkspace** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("Incremental workspace output", d);
  *out = nullptr;
  if (!((node_count == 4 && point_count == 1) ||
        (node_count == 10 && point_count == 4) ||
        (node_count == 6 && point_count == 2) ||
        (node_count == 5 && point_count == 5) ||
        (node_count == 8 && point_count == 8)))
    return model_storage::Invalid("Incremental workspace node/IP pair", d);
  std::unique_ptr<IncrementalElasticWorkspace> workspace(
      new (std::nothrow) IncrementalElasticWorkspace);
  if (!workspace)
    return model_storage::AllocationFailure("Incremental workspace", d);
  workspace->node_count_ = node_count;
  workspace->point_count_ = point_count;
  for (auto*& matrix : workspace->matrices_) {
    const auto status = Matrix::Create(3, 3, &matrix, d);
    if (status != StatusCode::kOk) return status;
  }
  auto status = Matrix::Create(9, 9, &workspace->derivative_rows_, d);
  if (status == StatusCode::kOk)
    status = Matrix::Create(9, 9, &workspace->identity9_, d);
  if (status == StatusCode::kOk)
    status = Vector::Create(3, &workspace->diagonal_, d);
  if (status == StatusCode::kOk)
    status = Vector::Create(3, &workspace->ones_, d);
  if (status == StatusCode::kOk)
    status = Vector::Create(9, &workspace->flat_, d);
  if (status == StatusCode::kOk)
    status = Vector::Create(6, &workspace->tensor_, d);
  if (status == StatusCode::kOk)
    status = Vector::Create(3 * node_count, &workspace->point_force_, d);
  if (status != StatusCode::kOk) return status;
  const int nd = 3 * node_count;
  const int rows[] = {9, 9, 9, nd, 3, 3, 3, 3, 3, 3};
  const int columns[] = {nd, 9, nd, nd, node_count, 3, 3, 3, 3, 3};
  for (int i = 0; i < kElementMatrixCount; ++i) {
    status = Matrix::Create(rows[i], columns[i],
                            &workspace->element_matrices_[i], d);
    if (status != StatusCode::kOk) return status;
  }
  workspace->candidate_f_ = new (std::nothrow) Matrix*[point_count]();
  if (!workspace->candidate_f_)
    return model_storage::AllocationFailure("Incremental candidate F array", d);
  for (int g = 0; g < point_count; ++g) {
    status = Matrix::Create(3, 3, &workspace->candidate_f_[g], d);
    if (status != StatusCode::kOk) return status;
  }
  for (int i = 0; i < 9; ++i) {
    status = workspace->identity9_->Set(i, i, 1, d);
    if (status != StatusCode::kOk) return status;
    if (i >= 3) continue;
    status = workspace->matrices_[0]->Set(i, i, 1, d);
    if (status == StatusCode::kOk) status = workspace->ones_->Set(i, 1, d);
    if (status != StatusCode::kOk) return status;
  }
  if (node_count == 8) {
    status = IncrementalHexa8Workspace::Create(&workspace->hexa8_, d);
    if (status != StatusCode::kOk) return status;
  }
  *out = workspace.release();
  return StatusCode::kOk;
}
IncrementalElasticWorkspace::~IncrementalElasticWorkspace() {
  delete hexa8_;
  for (const auto* matrix : matrices_) delete matrix;
  delete derivative_rows_;
  delete identity9_;
  delete diagonal_;
  delete ones_;
  delete flat_;
  for (const auto* matrix : element_matrices_) delete matrix;
  delete tensor_;
  if (candidate_f_)
    for (int g = 0; g < point_count_; ++g) delete candidate_f_[g];
  delete[] candidate_f_;
  delete point_force_;
}
Matrix* IncrementalElasticWorkspace::matrix(int index) const {
  return index >= 0 && index < kMatrixCount ? matrices_[index] : nullptr;
}
Matrix* IncrementalElasticWorkspace::derivative_rows() const {
  return derivative_rows_;
}
Matrix* IncrementalElasticWorkspace::identity9() const { return identity9_; }
Vector* IncrementalElasticWorkspace::diagonal() const { return diagonal_; }
Vector* IncrementalElasticWorkspace::ones() const { return ones_; }
Vector* IncrementalElasticWorkspace::flat() const { return flat_; }
Matrix* IncrementalElasticWorkspace::element_matrix(int index) const {
  return index >= 0 && index < kElementMatrixCount
             ? element_matrices_[index] : nullptr;
}
Vector* IncrementalElasticWorkspace::tensor() const { return tensor_; }
Matrix* IncrementalElasticWorkspace::candidate_f(int point) const {
  return point >= 0 && point < point_count_ ? candidate_f_[point] : nullptr;
}
}  // namespace fesa
