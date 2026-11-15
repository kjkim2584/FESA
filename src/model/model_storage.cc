#include "model/model_storage.h"

#include <cstring>
#include <new>

namespace fesa {
namespace model_storage {

char* CopyText(const char* text) {
  if (text == nullptr) return nullptr;
  const std::size_t length = std::strlen(text) + 1;
  char* copy = new (std::nothrow) char[length];
  if (copy != nullptr) std::memcpy(copy, text, length);
  return copy;
}

bool ValidId(int id) { return id > 0 && id <= 999999999; }

StatusCode Invalid(const char* operation, Diagnostic* d) {
  return ReportFailure(StatusCode::kInvalidArgument,
                       DiagnosticStage::kValidation, operation, d);
}

StatusCode AllocationFailure(const char* operation, Diagnostic* d) {
  return ReportFailure(StatusCode::kAllocationFailure,
                       DiagnosticStage::kAllocation, operation, d);
}

}  // namespace model_storage
}  // namespace fesa
