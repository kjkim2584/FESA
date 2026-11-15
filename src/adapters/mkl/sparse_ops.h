/// @file
/// @brief Internal opaque Sparse BLAS boundary, successful handles only.
#ifndef FESA_ADAPTERS_MKL_SPARSE_OPS_H_
#define FESA_ADAPTERS_MKL_SPARSE_OPS_H_
#include <cstddef>
namespace fesa {
namespace sparse_ops {
/// Count/byte preflight, independently testable without allocating huge
/// buffers.
bool ValidSizes(std::size_t size, std::size_t count);
/// Checks managed CSR bytes: 12*count+12*size+4, without imposing LP64.
/// Includes values/columns, two vendor row arrays and owned n+1 offsets;
/// excludes opaque backend scratch. No allocation or numerical operation.
bool ValidCsrStorageSize(std::size_t size, std::size_t count);
/// Conservative addition preflight: each input and count_a+count_b must fit
/// LP64, and the result's managed CSR bytes must fit size_t. May reject an
/// overlapping union that would fit exactly. No allocation/backend call.
bool ValidAdditionSizes(std::size_t size, std::size_t count_a,
                        std::size_t count_b);
/// Inputs borrowed for call; successful output owns independent CSR storage.
/// count>0; zero-based indices and storage mode are validated by SparseMatrix.
int Convert(int size, int count, const int* rows, const int* columns,
            const double* values, void** out);
/// New owned sorted CSR sum; inputs remain valid and unchanged.
/// Caller must first validate result count/storage with ValidAdditionSizes.
int Add(void* a, void* b, void** out);
/// Exports borrowed arrays, expiring when the handle is destroyed.
int Export(void* handle, int* size, int** starts, int** ends, int** columns,
           double** values);
/// Releases a successfully created handle (null allowed).
void Release(void* handle);
}  // namespace sparse_ops
}  // namespace fesa
#endif  // FESA_ADAPTERS_MKL_SPARSE_OPS_H_
