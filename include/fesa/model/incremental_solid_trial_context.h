/// @file
/// @brief Borrowed per-call incremental-solid baseline, outputs and scratch.
#ifndef FESA_MODEL_INCREMENTAL_SOLID_TRIAL_CONTEXT_H_
#define FESA_MODEL_INCREMENTAL_SOLID_TRIAL_CONTEXT_H_
#include "fesa/math/matrix.h"
namespace fesa {
class IncrementalElasticWorkspace;
/// Pose/history/analysis own these buffers and outlive this call. The element
/// never retains pointers. Committed U3n/global L and stress ng×6/global F/L²
/// use [11,22,33,12,13,23], tensor shear. Candidate stress is a separate
/// ng×6 numerical output; failure invalidates it, never committed values.
struct IncrementalSolidTrialContext {
  const Vector* committed_translations = nullptr;  ///< Borrowed pose copy.
  const Matrix* committed_stress = nullptr;  ///< Borrowed history baseline.
  Matrix* candidate_stress = nullptr;  ///< Borrowed pose-owned output.
  IncrementalElasticWorkspace* workspace = nullptr;  ///< Borrowed worker scratch.
  /// Call-local IE5-R geometric stress source; never material/history state.
  /// false uses candidate stress, true uses the fixed committed baseline.
  bool use_committed_geometric_stress = false;
};
}  // namespace fesa
#endif  // FESA_MODEL_INCREMENTAL_SOLID_TRIAL_CONTEXT_H_
