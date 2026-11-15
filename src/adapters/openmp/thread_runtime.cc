#include "adapters/openmp/thread_runtime.h"

#define NOMINMAX
#include <Windows.h>

#include <algorithm>
#include <atomic>
#include <chrono>

#include "base/diagnostic_copy.h"
#include "mkl.h"
#include "omp.h"
namespace fesa {
namespace {
thread_local int configured_threads = 1;
// Zero means this coordinator has not captured native affinity yet.
thread_local int pardiso_cnr_thread_ceiling = 0;
thread_local bool element_execution = false;
thread_local int observed_team_maximum = 0;
thread_local double dispatch_seconds = 0;
thread_local void (*failure_observer)(ElementFailureEvent, int, int, void*) =
    nullptr;
thread_local void* failure_observer_context = nullptr;

}  // namespace
StatusCode QueryProcessAvailableCpus(int* count, Diagnostic* d) {
  if (!count)
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation,
                         "QueryProcessAvailableCpus", d);
  DWORD_PTR process_mask = 0, system_mask = 0;
  if (!GetProcessAffinityMask(GetCurrentProcess(), &process_mask,
                              &system_mask)) {
    const DWORD error = GetLastError();
    if (d)
      d->Record(StatusCode::kLibraryFailure, DiagnosticStage::kValidation,
                "QueryProcessAvailableCpus", static_cast<int>(error));
    return StatusCode::kLibraryFailure;
  }
  int available = 0;
  for (; process_mask; process_mask >>= 1)
    if (process_mask & 1) ++available;
  if (available == 0)
    return ReportFailure(StatusCode::kLibraryFailure,
                         DiagnosticStage::kValidation,
                         "QueryProcessAvailableCpus: empty affinity mask", d);
  *count = available;
  return StatusCode::kOk;
}
StatusCode ConfigureThreadRuntime(int requested_threads, Diagnostic* d) {
  if (requested_threads <= 0 || InElementExecution())
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation, "ConfigureThreadRuntime",
                         d);
  int available = 0;
  const auto status = QueryProcessAvailableCpus(&available, d);
  if (status != StatusCode::kOk) return status;
  omp_set_dynamic(0);
  omp_set_num_threads(requested_threads);
  mkl_set_dynamic(0);
  mkl_set_num_threads_local(0);
  mkl_set_num_threads(requested_threads);
  configured_threads = requested_threads;
  pardiso_cnr_thread_ceiling = available;
  observed_team_maximum = 0;
  dispatch_seconds = 0;
  return StatusCode::kOk;
}
StatusCode PardisoCnrThreadCeiling(int* ceiling, Diagnostic* d) {
  if (!ceiling || InElementExecution())
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation,
                         "PardisoCnrThreadCeiling", d);
  if (pardiso_cnr_thread_ceiling == 0) {
    int available = 0;
    const auto status = QueryProcessAvailableCpus(&available, d);
    if (status != StatusCode::kOk) return status;
    pardiso_cnr_thread_ceiling = available;
  }
  *ceiling = pardiso_cnr_thread_ceiling;
  return StatusCode::kOk;
}
int ConfiguredThreadCount() { return configured_threads; }
int ObservedElementTeamMaximum() { return observed_team_maximum; }
double ElementDispatchSeconds() { return dispatch_seconds; }
bool InElementExecution() { return element_execution || omp_in_parallel(); }
void SetElementFailureObserver(
    void (*callback)(ElementFailureEvent, int, int, void*), void* context) {
  failure_observer = callback;
  failure_observer_context = context;
}
StatusCode ExecuteElementTasks(int count, int worker_count, ElementTask task,
                               void* context, Diagnostic* diagnostics,
                               int* observed, Diagnostic* d) {
  if (observed) *observed = 0;
  if (count < 0 || worker_count <= 0 || !task || !diagnostics ||
      InElementExecution())
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation, "ExecuteElementTasks",
                         d);
  for (int i = 0; i < worker_count; ++i) {
    CopyDiagnosticMetadata(d, &diagnostics[i]);
    diagnostics[i].Record(StatusCode::kOk, DiagnosticStage::kNone, nullptr);
    diagnostics[i].element_id = diagnostics[i].node_id = 0;
  }
  if (count == 0) return StatusCode::kOk;
  const auto observer = failure_observer;
  void* const observer_context = failure_observer_context;
  std::atomic<bool> stop(false);
  int team = 0;
  const auto start = std::chrono::steady_clock::now();
#pragma omp parallel num_threads(worker_count)
  {
    const int worker = omp_get_thread_num();
    const int previous = mkl_set_num_threads_local(1);
    element_execution = true;
#pragma omp single
    {
      team = omp_get_num_threads();
    }
#pragma omp for schedule(static)
    for (int item = 0; item < count; ++item) {
      if (stop.load(std::memory_order_acquire)) {
        if (observer)
          observer(ElementFailureEvent::kStopObserved, item, worker,
                   observer_context);
        continue;
      }
      const StatusCode status =
          task(item, worker, context, &diagnostics[worker]);
      if (status != StatusCode::kOk) {
        if (diagnostics[worker].code == StatusCode::kOk)
          diagnostics[worker].Record(status, DiagnosticStage::kComputation,
                                     "Element task");
        stop.store(true, std::memory_order_release);
        if (observer)
          observer(ElementFailureEvent::kStopPublished, item, worker,
                   observer_context);
      }
    }
    mkl_set_num_threads_local(previous);
    element_execution = false;
  }
  dispatch_seconds += std::chrono::duration<double>(
      std::chrono::steady_clock::now() - start).count();
  observed_team_maximum = std::max(observed_team_maximum, team);
  if (observed) *observed = team;
  const Diagnostic* failure = nullptr;
  for (int i = 0; i < worker_count; ++i)
    if (diagnostics[i].code != StatusCode::kOk &&
        (!failure || diagnostics[i].element_id < failure->element_id))
      failure = &diagnostics[i];
  if (!failure) return StatusCode::kOk;
  CopyDiagnosticMetadata(failure, d);
  return failure->code;
}
}  // namespace fesa
