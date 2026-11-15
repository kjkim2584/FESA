/// @file
/// @brief Internal sequential HDF5 adapter fault boundaries for verification.
#ifndef FESA_ADAPTERS_HDF5_HDF5_RESULT_WRITER_H_
#define FESA_ADAPTERS_HDF5_HDF5_RESULT_WRITER_H_
namespace fesa {
namespace hdf5 {
/// One-shot faults before the native boundary, except child/file close which
/// perform native close first then report failure (no double close).
enum class FaultPoint {
  kNone,
  kDataset,
  kPayloadFlush,
  kComplete,
  kHighWater,
  kCompletionFlush,
  kFinishState,
  kFinishFlush,
  kChildClose,
  kFileClose,
  kLateReopen,
  kLateState
};
/// Internal tests only; thread-local and consumed once. No public solver flags.
void FailNext(FaultPoint point);
/// Number of captured native error stack entries on this thread, latest error.
int NativeErrorDepth();
}  // namespace hdf5
}  // namespace fesa
#endif  // FESA_ADAPTERS_HDF5_HDF5_RESULT_WRITER_H_
