#include "fesa/analysis/rotation_state.h"

#include <cassert>
#include <cmath>
#include <limits>
#include <memory>
#include <new>
#include <utility>

#include "fesa/math/rotation.h"
#include "model/model_storage.h"

namespace fesa {
namespace {

StatusCode GatherOrientation(const Matrix* orientations, int node, Vector* out,
                             Diagnostic* d) {
  const int indices[4] = {4 * node, 4 * node + 1, 4 * node + 2, 4 * node + 3};
  return orientations->GatherTo(indices, 4, out, 0, d);
}

// Always derive trial a from the accepted increment's reference. Newton
// corrections are spatial and must not be accumulated as body-spin history.
StatusCode UpdateBodySpin(const Matrix* committed_q, const Vector* committed_a,
                          int node, const Vector* orientation,
                          Vector* candidate_a, bool* branch, Diagnostic* d) {
  std::unique_ptr<Vector> owners[3];
  const int sizes[3] = {4, 3, 3};
  for (int i = 0; i < 3; ++i) {
    Vector* raw = nullptr;
    const auto status = Vector::Create(sizes[i], &raw, d);
    owners[i].reset(raw);
    if (status != StatusCode::kOk) return status;
  }
  auto status = GatherOrientation(committed_q, node, owners[0].get(), d);
  if (status != StatusCode::kOk) return status;
  status = BodyRotationIncrement(owners[0].get(), orientation, owners[1].get(),
                                 branch, d);
  if (status != StatusCode::kOk) return status;
  const int node_indices[3] = {3 * node, 3 * node + 1, 3 * node + 2};
  status = owners[2]->GatherFrom(committed_a, node_indices, 3, 0, d);
  if (status != StatusCode::kOk) return status;
  status = owners[2]->AddScaled(1.0, owners[1].get(), d);
  if (status != StatusCode::kOk) return status;
  const int components[3] = {0, 1, 2};
  return candidate_a->GatherFrom(owners[2].get(), components, 3, 3 * node, d);
}

}  // namespace

StatusCode RotationState::Create(int node_count, RotationState** out,
                                 Diagnostic* d) {
  if (!out) return model_storage::Invalid("RotationState::Create output", d);
  *out = nullptr;
  if (node_count <= 0 || node_count > std::numeric_limits<int>::max() / 4)
    return model_storage::Invalid("RotationState::Create dimensions", d);
  std::unique_ptr<RotationState> state(new (std::nothrow) RotationState);
  if (!state)
    return model_storage::AllocationFailure("RotationState::Create", d);
  state->node_count_ = node_count;
  auto status = Matrix::Create(node_count, 4, &state->committed_q_, d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Create(node_count, 4, &state->candidate_q_, d);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(3 * node_count, &state->committed_a_, d);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(3 * node_count, &state->candidate_a_, d);
  if (status != StatusCode::kOk) return status;
  for (int node = 0; node < node_count; ++node) {
    status = state->committed_q_->Set(node, 0, 1.0, d);
    if (status != StatusCode::kOk) return status;
  }
  *out = state.release();
  return StatusCode::kOk;
}

RotationState::~RotationState() {
  delete candidate_a_;
  delete committed_a_;
  delete candidate_q_;
  delete committed_q_;
}

StatusCode RotationState::PrepareCandidate(Diagnostic* d) {
  candidate_ready_ = false;
  auto status = candidate_q_->CopyFrom(committed_q_, d);
  if (status != StatusCode::kOk) return status;
  status = candidate_a_->CopyFrom(committed_a_, d);
  if (status != StatusCode::kOk) return status;
  candidate_ready_ = true;
  return StatusCode::kOk;
}

StatusCode RotationState::ApplyCorrection(int node, const Vector* spatial_delta,
                                          bool* branch_exceeded,
                                          Diagnostic* d) {
  if (branch_exceeded) *branch_exceeded = false;
  const bool ready = candidate_ready_;
  candidate_ready_ = false;
  if (!ready || node < 0 || node >= node_count_ || !branch_exceeded ||
      !spatial_delta || spatial_delta->Size() != 3 ||
      spatial_delta == committed_a_ || spatial_delta == candidate_a_)
    return model_storage::Invalid("RotationState correction inputs", d);
  double angle = 0.0;
  auto status = spatial_delta->Norm2(&angle, d);
  if (status != StatusCode::kOk) return status;
  constexpr double kPi = 3.141592653589793238462643383279502884;
  if (angle >= kPi ||
      std::cos(angle / 2.0) <= 256 * std::numeric_limits<double>::epsilon()) {
    *branch_exceeded = true;
    return ReportFailure(StatusCode::kEquilibriumFailure,
                         DiagnosticStage::kComputation,
                         "RotationState applied correction crosses pi", d);
  }
  std::unique_ptr<Vector> owners[3];
  for (auto& owner : owners) {
    Vector* raw = nullptr;
    status = Vector::Create(4, &raw, d);
    owner.reset(raw);
    if (status != StatusCode::kOk) return status;
  }
  status = GatherOrientation(candidate_q_, node, owners[0].get(), d);
  if (status != StatusCode::kOk) return status;
  status = RotationExp(spatial_delta, owners[1].get(), d);
  if (status != StatusCode::kOk) return status;
  status =
      ComposeRotations(owners[1].get(), owners[0].get(), owners[2].get(), d);
  if (status != StatusCode::kOk) return status;
  status = UpdateBodySpin(committed_q_, committed_a_, node, owners[2].get(),
                          candidate_a_, branch_exceeded, d);
  if (status != StatusCode::kOk) return status;
  status = candidate_q_->SetRow(node, owners[2].get(), d);
  if (status != StatusCode::kOk) return status;
  candidate_ready_ = true;
  return StatusCode::kOk;
}

StatusCode RotationState::SetCandidateOrientation(int node, const Vector* q,
                                                  bool* branch_exceeded,
                                                  Diagnostic* d) {
  if (branch_exceeded) *branch_exceeded = false;
  const bool ready = candidate_ready_;
  candidate_ready_ = false;
  if (!ready || node < 0 || node >= node_count_ || !branch_exceeded || !q ||
      q->Size() != 4 || q == committed_a_ || q == candidate_a_)
    return model_storage::Invalid("RotationState prescribed orientation", d);
  Vector* raw = nullptr;
  auto status = Vector::Create(4, &raw, d);
  std::unique_ptr<Vector> orientation(raw);
  if (status != StatusCode::kOk) return status;
  status = NormalizeQuaternion(q, orientation.get(), d);
  if (status != StatusCode::kOk) return status;
  status = UpdateBodySpin(committed_q_, committed_a_, node, orientation.get(),
                          candidate_a_, branch_exceeded, d);
  if (status != StatusCode::kOk) return status;
  status = candidate_q_->SetRow(node, orientation.get(), d);
  if (status != StatusCode::kOk) return status;
  candidate_ready_ = true;
  return StatusCode::kOk;
}

StatusCode RotationState::PrincipalOutput(Vector* out, Diagnostic* d) const {
  if (!candidate_ready_ || !out || out->Size() != 3 * node_count_ ||
      out == committed_a_ || out == candidate_a_)
    return model_storage::Invalid("RotationState principal output", d);
  Vector* raw = nullptr;
  auto status = Vector::Create(4, &raw, d);
  std::unique_ptr<Vector> orientation(raw);
  if (status != StatusCode::kOk) return status;
  raw = nullptr;
  status = Vector::Create(3, &raw, d);
  std::unique_ptr<Vector> principal(raw);
  if (status != StatusCode::kOk) return status;
  const int components[3] = {0, 1, 2};
  for (int node = 0; node < node_count_; ++node) {
    status = GatherOrientation(candidate_q_, node, orientation.get(), d);
    if (status != StatusCode::kOk) return status;
    status = PrincipalRotationVector(orientation.get(), principal.get(), d);
    if (status != StatusCode::kOk) return status;
    status = out->GatherFrom(principal.get(), components, 3, 3 * node, d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}

void RotationState::Commit() noexcept {
  assert(candidate_ready_);
  std::swap(committed_q_, candidate_q_);
  std::swap(committed_a_, candidate_a_);
  candidate_ready_ = false;
}

void RotationState::DiscardCandidate() noexcept { candidate_ready_ = false; }
int RotationState::NodeCount() const { return node_count_; }
const Matrix* RotationState::CommittedOrientations() const {
  return committed_q_;
}
const Vector* RotationState::CommittedBodySpin() const { return committed_a_; }
const Matrix* RotationState::CandidateOrientations() const {
  return candidate_ready_ ? candidate_q_ : nullptr;
}
const Vector* RotationState::CandidateBodySpin() const {
  return candidate_ready_ ? candidate_a_ : nullptr;
}

}  // namespace fesa
