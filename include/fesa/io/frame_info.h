/// @file
/// @brief Scalar identity and physical times for one output frame.
#ifndef FESA_IO_FRAME_INFO_H_
#define FESA_IO_FRAME_INFO_H_
#include <cstdint>
namespace fesa {
/// Physical solution procedure, independent of element family.
enum class FrameProcedure {
  kStatic,
  kDynamic,
  kHeatTransferSteady,
  kHeatTransferTransient
};
/// Kinematic formulation used to prepare this frame's physical results.
enum class FrameGeometry { kLinear, kNonlinear };
/// Initialization output is distinct from a successful advancing increment.
enum class FrameKind { kInitial, kIncrementEnd };
/// Borrowed read-only for WriteFrame. step_name is valid UTF-8 (empty allowed)
/// and must outlive the call. Times are physical model times, not wall time.
/// Frame numbers start at 1 across steps. Static frames have positive
/// increments; dynamic initial frames use increment/time 0 and never complete a
/// step.
struct FrameInfo {
  std::int64_t frame_number = 1;
  std::int64_t step_number = 1;
  std::int64_t increment_number = 1;
  const char* step_name = "";
  double step_time = 1;
  double total_time = 1;
  FrameProcedure procedure = FrameProcedure::kStatic;
  FrameGeometry geometry = FrameGeometry::kLinear;
  FrameKind kind = FrameKind::kIncrementEnd;
  /// Whole-step completion, distinct from increment/frame completion.
  /// Incremental callers set false until the successful step endpoint.
  bool completes_step = true;
};
}  // namespace fesa
#endif  // FESA_IO_FRAME_INFO_H_
