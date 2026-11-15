/// @file
/// @brief Worker-only storage for coupled incremental HEX8 virtual work.
#ifndef FESA_ELEMENTS_INCREMENTAL_HEXA8_WORKSPACE_H_
#define FESA_ELEMENTS_INCREMENTAL_HEXA8_WORKSPACE_H_

#include "fesa/math/matrix.h"

namespace fesa {
/// Owns every buffer below; one worker, never history. No input pointers kept.
/// Candidate quantities are invalid after any failed evaluation and rebuilt
/// on every call. Borrowed getters last until next evaluation/destruction.
class IncrementalHexa8Workspace {
 public:
  /// Eight independent point buffers. All3x3 except Jacobian9x9/gradient3x8.
  enum Point {
    kFn,
    kInverse,
    kOrdinaryStress,
    kOrdinaryPk1,
    kOrdinaryJacobian,
    kG,
    kMidInverse,
    kModifiedF,
    kSpatialGradient,
    kPointCount
  };
  /// Fixed operation storage: first7 are3x3; then explicit operator/row sizes.
  enum Buffer {
    kDelta,
    kMid,
    kTemp,
    kTemp2,
    kStressTemp,
    kWeightedStress,
    kEffectivePk1,
    kLocal,
    kMap,
    kProduct,
    kThetaRow,
    kMidRow,
    kDivRow,
    kGlobalTheta,
    kGlobalDiv,
    kRankRow,
    kBufferCount
  };
  /// Eight-entry scalar coefficient vectors; ones is immutable after creation.
  enum Coefficient { kJ, kM, kW, kTheta, kPn, kOnes, kCoefficientCount };
  /// Caller owns success. Failure destroys partial buffers and nulls out.
  static StatusCode Create(IncrementalHexa8Workspace** out, Diagnostic* d);
  ~IncrementalHexa8Workspace();
  IncrementalHexa8Workspace(const IncrementalHexa8Workspace&) = delete;
  IncrementalHexa8Workspace& operator=(const IncrementalHexa8Workspace&) =
      delete;
  IncrementalHexa8Workspace(IncrementalHexa8Workspace&&) = delete;
  IncrementalHexa8Workspace& operator=(IncrementalHexa8Workspace&&) = delete;
  /// Borrowed mutable scratch; callers must respect enum dimensions/semantics.
  Matrix* point(int g, Point slot) const;
  Matrix* buffer(Buffer slot) const;
  Vector* coefficient(Coefficient slot) const;
  /// Borrowed array of modified F buffers for LE recovery after success only.
  const Matrix* const* modified_f_data() const { return modified_f_; }

 private:
  IncrementalHexa8Workspace() = default;
  Matrix* points_[8][kPointCount] = {};  ///< Owned per-IP quantities.
  Matrix* buffers_[kBufferCount] = {};   ///< Owned fixed operator scratch.
  Vector* coefficients_[kCoefficientCount] = {};  ///< Owned scalar arrays.
  const Matrix* modified_f_[8] = {};  ///< Non-owning aliases into points_.
};
}  // namespace fesa
#endif  // FESA_ELEMENTS_INCREMENTAL_HEXA8_WORKSPACE_H_
