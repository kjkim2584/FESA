/// @file
/// @brief Worker-owned numerical scratch, never material history/results.
#ifndef FESA_ELEMENTS_INCREMENTAL_ELASTIC_WORKSPACE_H_
#define FESA_ELEMENTS_INCREMENTAL_ELASTIC_WORKSPACE_H_
#include "fesa/math/matrix.h"
namespace fesa {
class IncrementalHexa8Workspace;
/// Single-worker buffers for IE2-IE4. No cross-call input pointers retained.
class IncrementalElasticWorkspace {
 public:
  static constexpr int kMatrixCount = 31;  ///< Required3×3 IE scratch slots.
  static constexpr int kElementMatrixCount = 10;  ///< Point/operator buffers.
  /// Preallocates scratch; caller deletes success. Failure cleans/nulls out.
  static StatusCode Create(IncrementalElasticWorkspace** out, Diagnostic* d);
  /// Exact full-type storage for4/1,10/4,6/2,5/5,8/8 node/IP pairs. Hexa8
  /// storage does not enable ordinary pointwise evaluation. Caller owns success;
  /// failure nulls out and destroys partial buffers. Borrows out/d for the call.
  static StatusCode CreateForSolid(int node_count, int point_count,
      IncrementalElasticWorkspace** out, Diagnostic* d);
  ~IncrementalElasticWorkspace();
  IncrementalElasticWorkspace(const IncrementalElasticWorkspace&) = delete;
  IncrementalElasticWorkspace& operator=(const IncrementalElasticWorkspace&) = delete;
  IncrementalElasticWorkspace(IncrementalElasticWorkspace&&) = delete;
  IncrementalElasticWorkspace& operator=(IncrementalElasticWorkspace&&) = delete;
  /// Borrowed internal buffers, valid until destruction; never share workers.
  Matrix* matrix(int index) const;
  Matrix* derivative_rows() const;
  Matrix* identity9() const;
  Vector* diagonal() const;
  Vector* ones() const;
  Vector* flat() const;
  /// Buffers: map9×3n, dP/A9×9, product9×3n, point K3n×3n,
  /// current spatial gradient3×n, then
  /// committed F/candidate F/committed σ/candidate σ/P, each3×3.
  Matrix* element_matrix(int index) const;
  /// Borrowed6-component pack/unpack transport, never an owned IP history.
  Vector* tensor() const;
  /// Immutable physical sizes and borrowed scratch, never committed history.
  int node_count() const { return node_count_; }
  int point_count() const { return point_count_; }
  Matrix* candidate_f(int point) const;
  const Matrix* const* candidate_f_data() const { return candidate_f_; }
  Vector* point_force() const { return point_force_; }
  /// Borrowed owned coupled scratch only for8/8; other types return null.
  IncrementalHexa8Workspace* hexa8_workspace() const { return hexa8_; }
 private:
  IncrementalElasticWorkspace() = default;
  Matrix* matrices_[kMatrixCount] = {};  ///< Owned; index0 immutable I.
  Matrix* derivative_rows_ = nullptr;  ///< Owned9×9 transposed dP/dF.
  Matrix* identity9_ = nullptr;  ///< Owned immutable9×9 identity.
  Vector* diagonal_ = nullptr;  ///< Owned3-entry diagonal gather.
  Vector* ones_ = nullptr;  ///< Owned immutable3-entry trace coefficients.
  Vector* flat_ = nullptr;  ///< Owned9-entry derivative column transport.
  /// Owned exact-sized element scratch.
  Matrix* element_matrices_[kElementMatrixCount] = {};
  Vector* tensor_ = nullptr;  ///< Owned6-entry tensor transport.
  int node_count_ = 0, point_count_ = 0;
  /// Owned array and ng owned3×3 candidate buffers. Invalid until recomputed;
  /// callers borrow only until the next worker evaluation or destruction.
  Matrix** candidate_f_ = nullptr;
  Vector* point_force_ = nullptr;  ///< Owned3n point contribution scratch.
  IncrementalHexa8Workspace* hexa8_ = nullptr;  ///< Owned only for8/8.
};
}  // namespace fesa
#endif  // FESA_ELEMENTS_INCREMENTAL_ELASTIC_WORKSPACE_H_
