/// @file
/// @brief Internal allocation helpers for owned textual model metadata.
#ifndef FESA_MODEL_MODEL_STORAGE_H_
#define FESA_MODEL_MODEL_STORAGE_H_

#include "fesa/base/diagnostic.h"

namespace fesa {
namespace model_storage {

/// Copies text to a new caller-owned char[] or returns null on allocation
/// error.
/// Null input is rejected. Input is borrowed only for this call.
char* CopyText(const char* text);
/// Original IDs use the approved Abaqus subset's positive nine-digit range.
bool ValidId(int id);
/// Records a model argument error in optional borrowed d.
StatusCode Invalid(const char* operation, Diagnostic* d);
/// Records allocation failure in optional borrowed d.
StatusCode AllocationFailure(const char* operation, Diagnostic* d);

}  // namespace model_storage
}  // namespace fesa

#endif  // FESA_MODEL_MODEL_STORAGE_H_
