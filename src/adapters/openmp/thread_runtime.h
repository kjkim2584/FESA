/// @file
/// @brief Process-local Intel OpenMP/MKL configuration and joined dispatch.
#ifndef FESA_ADAPTERS_OPENMP_THREAD_RUNTIME_H_
#define FESA_ADAPTERS_OPENMP_THREAD_RUNTIME_H_
#include "fesa/base/diagnostic.h"
namespace fesa {
/// Returns the positive CPU count in the native current-process affinity mask.
/// Borrows count/d for the call; null count is invalid. Native failure reports
/// kLibraryFailure and vendor_code; every failure preserves *count and runtime.
StatusCode QueryProcessAvailableCpus(int* count, Diagnostic* d);
/// Sets OpenMP/MKL dynamic=false and requested positive thread counts.
/// Call before numerical work, outside parallel regions; never changes the
/// persistent environment. After native affinity discovery succeeds, replaces
/// the coordinator's requested count and independent CNR ceiling snapshot and
/// resets observations/time. Any failure preserves all prior runtime/TLS state.
/// Positive requests above affinity count remain valid for direct API callers.
StatusCode ConfigureThreadRuntime(int requested_threads, Diagnostic* d);
/// Returns this coordinator's last successful positive affinity snapshot.
/// Before configure, lazily queries/caches affinity without changing requested
/// count (default one), MKL/OpenMP settings or observations/time. Borrows
/// output and optional d for the call. Null output/element execution is
/// invalid; any failure preserves output/cache. This is a CNR setting, not
/// observed workers. Assumes stable affinity and sequential PARDISO analysis
/// without RNG interference. See
/// docs/Linear static/parallel_reproducibility.md.
StatusCode PardisoCnrThreadCeiling(int* ceiling, Diagnostic* d);
/// Requested count last configured on this coordinator; defaults to one.
int ConfiguredThreadCount();
/// Coordinator observations since successful ConfigureThreadRuntime: maximum
/// actual outer team (zero before work), and summed dispatch wall seconds.
/// Includes fork/join and packing kernels, not a sum of worker CPU durations.
/// Nested within assembly/equilibrium/recovery timings; never add both totals.
int ObservedElementTeamMaximum();
double ElementDispatchSeconds();
/// True inside FESA element dispatch (including a serialized team) or any
/// active outer OpenMP region. Used to guard coordinator-only numerical work.
bool InElementExecution();
/// Callback borrows context and exclusive worker diagnostic until return.
/// It must not throw or log; writes stay within its assigned disjoint output
/// region/object. Checked call-local scratch may be allocated and released
/// before return; no cross-call workspace/reuse policy is implied.
using ElementTask = StatusCode (*)(int item, int worker, void* context,
                                   Diagnostic* d);
/// Failure-only test observation after publishing stop or observing its load.
enum class ElementFailureEvent { kStopPublished, kStopObserved };
/// Borrows nonthrowing callback/context on this coordinator until reset. The
/// dispatch snapshots them before entry; callbacks can run concurrently and
/// must not log, throw, exit, mutate numerical output, or outlive the join.
/// No callback/branch is added to successful task execution. Null clears it.
void SetElementFailureObserver(
    void (*callback)(ElementFailureEvent, int item, int worker, void* context),
    void* context);
/// Executes disjoint tasks with worker_count preallocated diagnostics. Each
/// worker saves its MKL local setting, uses local=1, restores it and joins.
/// On failure workers skip new work after observing stop; in-flight work ends.
/// Returns the smallest recorded failing element ID and its cause after join.
/// Caller discards ALL task output on failure. observed receives actual team
/// size (zero for no work), not a claim about internal MKL/PARDISO teams.
/// Borrowed context/diagnostics remain valid through the join; d may be null.
StatusCode ExecuteElementTasks(int count, int worker_count, ElementTask task,
                               void* context, Diagnostic* diagnostics,
                               int* observed, Diagnostic* d);
}  // namespace fesa
#endif  // FESA_ADAPTERS_OPENMP_THREAD_RUNTIME_H_
