/// @file
/// @brief Shared bounded scratch and checked index mapping for E1/E2.
#ifndef FESA_ANALYSIS_NUMERICAL_SUPPORT_H_
#define FESA_ANALYSIS_NUMERICAL_SUPPORT_H_
#include "analysis/static_assembler.h"
#include "fesa/elements/element_workspace.h"
#include "fesa/loads/beam_line_load.h"
namespace fesa {
namespace numerical_support {
/// Small exclusive caller-side numerical scratch; element is borrowed from
/// StaticAnalysis or the standalone Buffers owner. No global arrays per worker.
class WorkerBuffers {
 public:
  /// Creates exact local-shape scratch; caller deletes success. Borrows
  /// optional element context through destruction; failure nulls and cleans.
  /// Numerical calls overwrite scratch, including on failure.
  static StatusCode Create(ElementWorkspace* element, WorkerBuffers** out,
                           Diagnostic* d, int element_dofs = 12);
  ~WorkerBuffers();
  WorkerBuffers(const WorkerBuffers&) = delete;
  WorkerBuffers& operator=(const WorkerBuffers&) = delete;
  WorkerBuffers(WorkerBuffers&&) = delete;
  WorkerBuffers& operator=(WorkerBuffers&&) = delete;
  ElementWorkspace* element = nullptr;  ///< Borrowed through joined dispatch.
  Matrix* stiffness = nullptr;  ///< Owned local K, or M in initial mass assembly.
  Vector* u = nullptr;                  ///< Owned exact local candidate.
  Vector* force = nullptr;  ///< Owned exact local generalized force.
  BeamLineLoadWorkspace* line_workspace = nullptr;  ///< Owned when enabled.
  Vector* line_density = nullptr;  ///< Owned global F/L, three components.
  Vector* line_global = nullptr;   ///< Owned equivalent global [F,M], 12.
  Vector* line_local = nullptr;    ///< Owned equivalent local [F,M], 12.
 private:
  WorkerBuffers() = default;
};
/// Bounded contribution storage; no full dense stiffness or per-element global
/// vector. Local owner deletes this after the numerical procedure returns.
class Buffers {
 public:
  static constexpr int kElements = 16384;
  static constexpr int kMaxElementDofs = 30;
  static constexpr int kMaxUpperEntries = 465;
  /// Caller deletes successful *out; failure cleans partial allocations.
  /// element_dofs is a capacity/stride, not a cache. Legacy Beam fallback
  /// creates its required contexts; other families receive supplied contexts
  /// or null and decide whether a context is required for their operation.
  static StatusCode Create(bool matrix_storage, int capacity, Buffers** out,
                           Diagnostic* d, ElementWorkspace** workers = nullptr,
                           int worker_count = 1, int element_dofs = 12,
                           bool beam_workspace = true);
  ~Buffers();
  Buffers(const Buffers&) = delete;
  Buffers& operator=(const Buffers&) = delete;
  Buffers(Buffers&&) = delete;
  Buffers& operator=(Buffers&&) = delete;
  WorkerBuffers** workers = nullptr;  ///< Owned exclusive scratch per worker.
  Diagnostic* diagnostics = nullptr;  ///< Owned per-worker records.
  int element_dofs = 12;  ///< Physical stride/capacity, never a result cache.
  int worker_count = 0;  ///< Prepared count, not actual observed team size.
  Element** elements = nullptr;     ///< Owned array; borrowed element entries.
  int* global_indices =
      nullptr;  ///< Owned element_dofs indices per chunk element.
  int* matrix_offsets = nullptr;    ///< Owned disjoint packed prefix offsets.
  int* vector_offsets = nullptr;    ///< Owned disjoint packed prefix offsets.
  Vector* vector_values = nullptr;  ///< Owned bounded contribution values.
  Vector* matrix_values = nullptr;  ///< Owned bounded upper COO values.
  int* vector_indices = nullptr;    ///< Owned destination equations/global.
  int* rows = nullptr;              ///< Owned upper COO rows, assembly only.
  int* columns = nullptr;           ///< Owned upper COO columns, assembly only.
 private:
  Buffers() = default;
  ElementWorkspace** owned_elements_ = nullptr;  ///< Standalone fallback only.
};
/// Checks positive model/global workspace shapes and matching bound conditions.
StatusCode Validate(const AnalysisModel* model, const DofManager* dofs,
                    const EffectiveConditions* conditions,
                    StaticWorkspace* workspace, Diagnostic* d);
/// Writes element's actual physical global indices (at most 30); fails for
/// invalid/missing node mapping.
StatusCode ElementIndices(const Element* element, const DofManager* dofs,
                          int* indices, Diagnostic* d);
/// Overwrites full external force from already library-summed definitions.
StatusCode ExternalForce(const DofManager* dofs,
                         const EffectiveConditions* conditions, Vector* out,
                         Diagnostic* d);
/// Prepares optional scratch on fresh buffers before joined dispatch. null
/// state is a no-op and preserves the old model's allocations. On failure the
/// caller discards Buffers; its destructor releases partial preparation.
StatusCode PrepareLineLoads(const BeamLineLoadState* state, Buffers* buffers,
                            Diagnostic* d);
/// L1 evaluation into this worker's scratch; caller provides a prepared worker
/// and nonnull state. Model/State remain unchanged. Failure discards scratch.
StatusCode EvaluateLineLoad(const BeamLineLoadState* state,
                            const Element* element, WorkerBuffers* worker,
                            Diagnostic* d);
}  // namespace numerical_support
}  // namespace fesa
#endif  // FESA_ANALYSIS_NUMERICAL_SUPPORT_H_
