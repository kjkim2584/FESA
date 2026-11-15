#include "analysis/finite_hht_tangent.h"

#include <cmath>
#include <memory>

#include "adapters/mkl/dense_ops.h"
#include "model/model_storage.h"

namespace fesa {
namespace {
StatusCode Compose(int internal_count, double alpha,
                   const TranslationalLinearization* translation,
                   const RotationalLinearization* rotation,
                   const Matrix* structural,
                   const ElementDynamicResult* physical, Matrix* out,
                   Diagnostic* d) {
  if (!translation || !rotation || !structural || !physical || !out ||
      !std::isfinite(alpha) || alpha < -1. / 3 || alpha > 0 ||
      !std::isfinite(translation->velocity_per_displacement) ||
      !std::isfinite(translation->acceleration_per_displacement) ||
      translation->velocity_per_displacement <= 0 ||
      translation->acceleration_per_displacement <= 0 ||
      (internal_count == 0 ? structural->Rows() % 6 != 0
                           : structural->Rows() != 14) ||
      structural->Rows() != structural->Columns())
    return model_storage::Invalid("Finite HHT tangent inputs", d);
  const int size = structural->Rows(), nodes = (size - internal_count) / 6;
  const Matrix* physical_blocks[] = {structural,
                                     physical->inertia_displacement,
                                     physical->inertia_velocity,
                                     physical->inertia_acceleration,
                                     physical->damping_displacement,
                                     physical->damping_velocity};
  for (const auto* matrix : physical_blocks)
    if (!matrix || matrix == out || matrix->Rows() != size ||
        matrix->Columns() != size ||
        !dense_ops::AllFinite(size * size, matrix->Data()))
      return model_storage::Invalid("Finite HHT tangent physical blocks", d);
  if (out->Rows() != size || out->Columns() != size)
    return model_storage::Invalid("Finite HHT tangent output shape", d);
  for (const auto* matrix :
       {rotation->velocity_per_rotation, rotation->acceleration_per_rotation})
    if (!matrix || matrix == out || matrix->Rows() != nodes ||
        matrix->Columns() != 9 ||
        !dense_ops::AllFinite(9 * nodes, matrix->Data()))
      return model_storage::Invalid("Finite HHT tangent R2 blocks", d);
  std::unique_ptr<Matrix> buffers[3];
  for (auto& buffer : buffers) {
    Matrix* raw = nullptr;
    const auto status = Matrix::Create(size, size, &raw, d);
    buffer.reset(raw);
    if (status != StatusCode::kOk) return status;
  }
  auto* velocity = buffers[0].get();
  auto* acceleration = buffers[1].get();
  auto* result = buffers[2].get();
  for (int node = 0; node < nodes; ++node) {
    for (int c = 0; c < 3; ++c) {
      auto status = velocity->Set(6 * node + c, 6 * node + c,
                                  translation->velocity_per_displacement, d);
      if (status == StatusCode::kOk)
        status =
            acceleration->Set(6 * node + c, 6 * node + c,
                              translation->acceleration_per_displacement, d);
      if (status == StatusCode::kOk)
        status = velocity->CopyBlockFrom(rotation->velocity_per_rotation, node,
                                         3 * c, 1, 3, 6 * node + 3 + c,
                                         6 * node + 3, d);
      if (status == StatusCode::kOk)
        status = acceleration->CopyBlockFrom(
            rotation->acceleration_per_rotation, node, 3 * c, 1, 3,
            6 * node + 3 + c, 6 * node + 3, d);
      if (status != StatusCode::kOk) return status;
    }
  }
  for (int i = size - internal_count; i < size; ++i) {
    auto status =
        velocity->Set(i, i, translation->velocity_per_displacement, d);
    if (status == StatusCode::kOk)
      status = acceleration->Set(i, i,
                                 translation->acceleration_per_displacement, d);
    if (status != StatusCode::kOk) return status;
  }
  auto status = result->CopyFrom(structural, d);
  if (status == StatusCode::kOk) status = result->Scale(1 + alpha, d);
  if (status == StatusCode::kOk)
    status = result->AddScaled(1 + alpha, physical->damping_displacement, d);
  if (status == StatusCode::kOk)
    status = result->AddScaled(1, physical->inertia_displacement, d);
  // Right multiplication is essential: columns are spatial corrections.
  if (status == StatusCode::kOk)
    status = Matrix::MultiplyAdd(physical->inertia_velocity, false, velocity,
                                 false, 1, 1, result, d);
  if (status == StatusCode::kOk)
    status = Matrix::MultiplyAdd(physical->inertia_acceleration, false,
                                 acceleration, false, 1, 1, result, d);
  if (status == StatusCode::kOk)
    status = Matrix::MultiplyAdd(physical->damping_velocity, false, velocity,
                                 false, 1 + alpha, 1, result, d);
  return status == StatusCode::kOk ? out->CopyFrom(result, d) : status;
}
}  // namespace
StatusCode ComposeFiniteHhtTangent(
    double alpha, const TranslationalLinearization* translation,
    const RotationalLinearization* rotation, const Matrix* structural,
    const ElementDynamicResult* physical, Matrix* out, Diagnostic* d) {
  return Compose(0, alpha, translation, rotation, structural, physical, out, d);
}
StatusCode ComposeBeamFiniteHhtTangent(
    double alpha, const TranslationalLinearization* translation,
    const RotationalLinearization* rotation, const Matrix* structural,
    const ElementDynamicResult* physical, Matrix* out, Diagnostic* d) {
  return Compose(2, alpha, translation, rotation, structural, physical, out, d);
}
}  // namespace fesa
