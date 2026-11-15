/// @file
/// @brief Vendor-independent result codes for fallible solver operations.
#ifndef FESA_BASE_STATUS_CODE_H_
#define FESA_BASE_STATUS_CODE_H_

namespace fesa {

/// Failure categories; kOk alone denotes success.
enum class StatusCode {
  kOk,
  kInvalidArgument,
  kInputError,
  kUnsupported,
  kAllocationFailure,
  kLibraryFailure,
  kSingular,
  kEquilibriumFailure,
  kIoFailure
};

}  // namespace fesa

#endif  // FESA_BASE_STATUS_CODE_H_
