/// @file
/// @brief Internal analysis-log stdio fault boundaries; never public flags.
#ifndef FESA_IO_LOG_STDIO_H_
#define FESA_IO_LOG_STDIO_H_
namespace fesa {
namespace io {
/// Simulated native error at checked fputs/fflush/fclose boundary. Close always
/// performs native close first; stream ownership is relinquished exactly once.
enum class LogFault { kNone, kWrite, kFlush, kClose };
/// One-shot per-thread injection after skip matching calls. Default skip=0.
void FailLogAfter(LogFault fault, int skip = 0);
/// Read-only count of checked Write boundaries on this control thread.
/// Allows measured ordinal fault replay without a new fault option.
unsigned long long LogWriteCallCount();
}  // namespace io
}  // namespace fesa
#endif  // FESA_IO_LOG_STDIO_H_
