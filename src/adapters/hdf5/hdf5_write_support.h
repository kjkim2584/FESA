/// @file
/// @brief Internal shared native HDF5 ownership, fault and error boundaries.
#ifndef FESA_ADAPTERS_HDF5_HDF5_WRITE_SUPPORT_H_
#define FESA_ADAPTERS_HDF5_HDF5_WRITE_SUPPORT_H_
#include "adapters/hdf5/hdf5_result_writer.h"
#include "fesa/base/diagnostic.h"
#include "hdf5.h"
namespace fesa {
namespace hdf5 {
namespace internal {
/// Borrowed scalar metadata for one numeric dataset, valid for the call.
struct DatasetMetadata {
  const char* unit_dimension;
  const char* evaluation;
  const char* coordinate_system;
  const char* const* components;
  int component_count;
};
/// Checks native status and records the shared native error stack on failure.
bool Check(hid_t result, const char* operation, Diagnostic* d);
/// Writes a borrowed UTF-8 scalar attribute. Failure may leave partial output.
bool ScalarText(hid_t object, const char* name, const char* value,
                Diagnostic* d);
/// Writes borrowed UTF-8 strings. Failure may leave partial output.
bool StringDataset(hid_t group, const char* name, int count,
                   const char* const* values, Diagnostic* d);
/// Writes borrowed numeric storage through the shared fault/close boundary.
/// Optional metadata overrides legacy name-based metadata for this call only.
/// Failure may leave partial native output; caller aborts the writer.
bool Dataset(hid_t group, const char* name, hid_t file_type, hid_t memory_type,
             int rank, const hsize_t* dims, const void* data, Diagnostic* d,
             const DatasetMetadata* metadata = nullptr);
/// Consumes the same thread-local one-shot fault used by the writer.
bool Inject(FaultPoint point);
/// Captures native error depth before cleanup. Borrows stable operation/d for
/// the call, changes only optional diagnostic and thread-local error depth.
StatusCode Failure(const char* operation, Diagnostic* d,
                   StatusCode code = StatusCode::kIoFailure,
                   bool native = false);
/// Owns one native handle. No HDF5 type escapes this adapter-only header.
class Handle {
 public:
  /// Owns id until checked Close or destruction; close must match its kind.
  Handle(hid_t id, herr_t (*close)(hid_t)) : id_(id), close_(close) {}
  /// Best-effort cleanup only; success paths must call Close explicitly.
  ~Handle() {
    if (id_ >= 0) close_(id_);
  }
  Handle(const Handle&) = delete;
  Handle& operator=(const Handle&) = delete;
  Handle(Handle&&) = delete;
  Handle& operator=(Handle&&) = delete;
  /// Borrowed native ID, invalid after successful Close/Release.
  hid_t get() const { return id_; }
  /// Transfers the primitive native ID; caller immediately adopts it.
  hid_t Release() {
    const hid_t id = id_;
    id_ = -1;
    return id;
  }
  /// Native close then the shared one-shot child-close fault. Failure retains
  /// an ID only if native close failed and it is still valid. No double close.
  bool Close(Diagnostic* d) {
    if (id_ < 0) return true;
    const hid_t id = id_;
    if (close_(id) < 0) {
      Failure("HDF5 child close", d, StatusCode::kIoFailure, true);
      if (H5Iis_valid(id) <= 0) id_ = -1;
      return false;
    }
    id_ = -1;
    if (Inject(FaultPoint::kChildClose)) {
      Failure("HDF5 child close", d);
      return false;
    }
    return true;
  }

 private:
  hid_t id_;
  herr_t (*close_)(hid_t);
};
}  // namespace internal
}  // namespace hdf5
}  // namespace fesa
#endif  // FESA_ADAPTERS_HDF5_HDF5_WRITE_SUPPORT_H_
