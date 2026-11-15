#include "elements/incremental_hexa8_workspace.h"

#include <memory>
#include <new>

#include "model/model_storage.h"

namespace fesa {
StatusCode IncrementalHexa8Workspace::Create(IncrementalHexa8Workspace** out,
                                             Diagnostic* d) {
  if (!out) return model_storage::Invalid("HEX8 workspace output", d);
  *out = nullptr;
  std::unique_ptr<IncrementalHexa8Workspace> w(new (std::nothrow)
                                                   IncrementalHexa8Workspace);
  if (!w) return model_storage::AllocationFailure("HEX8 workspace", d);
  for (int g = 0; g < 8; ++g) {
    for (int i = 0; i < kPointCount; ++i) {
      const int r = i == kOrdinaryJacobian ? 9 : 3;
      const int c = i == kOrdinaryJacobian ? 9 : i == kSpatialGradient ? 8 : 3;
      auto status = Matrix::Create(r, c, &w->points_[g][i], d);
      if (status != StatusCode::kOk) return status;
    }
    w->modified_f_[g] = w->points_[g][kModifiedF];
  }
  const int rows[] = {3, 3, 3, 3, 3, 3, 3, 9, 9, 9, 1, 1, 1, 1, 1, 1};
  const int cols[] = {3, 3, 3, 3, 3, 3, 3, 9, 24, 24, 9, 9, 9, 24, 24, 9};
  static_assert(sizeof(rows) / sizeof(rows[0]) == kBufferCount);
  for (int i = 0; i < kBufferCount; ++i) {
    auto status = Matrix::Create(rows[i], cols[i], &w->buffers_[i], d);
    if (status != StatusCode::kOk) return status;
  }
  for (int i = 0; i < kCoefficientCount; ++i) {
    auto status = Vector::Create(8, &w->coefficients_[i], d);
    if (status != StatusCode::kOk) return status;
  }
  for (int g = 0; g < 8; ++g) {
    auto status = w->coefficients_[kOnes]->Set(g, 1, d);
    if (status != StatusCode::kOk) return status;
  }
  *out = w.release();
  return StatusCode::kOk;
}
IncrementalHexa8Workspace::~IncrementalHexa8Workspace() {
  for (auto& point : points_)
    for (const auto* x : point) delete x;
  for (const auto* x : buffers_) delete x;
  for (const auto* x : coefficients_) delete x;
}
Matrix* IncrementalHexa8Workspace::point(int g, Point slot) const {
  return g >= 0 && g < 8 && slot >= 0 && slot < kPointCount ? points_[g][slot]
                                                            : nullptr;
}
Matrix* IncrementalHexa8Workspace::buffer(Buffer slot) const {
  return slot >= 0 && slot < kBufferCount ? buffers_[slot] : nullptr;
}
Vector* IncrementalHexa8Workspace::coefficient(Coefficient slot) const {
  return slot >= 0 && slot < kCoefficientCount ? coefficients_[slot] : nullptr;
}
}  // namespace fesa
