#include "fesa/model/node.h"

#include <cmath>
#include <memory>
#include <new>

#include "model/model_storage.h"

namespace fesa {

Node::Node(int id, double x, double y, double z, Vector* results)
    : id_(id), x_(x), y_(y), z_(z), results_(results) {}
Node::~Node() {
  delete results_;
  delete temperature_results_;
}

StatusCode Node::Create(int id, double x, double y, double z,
                        Node** out, Diagnostic* d) {
  return Create(id, x, y, z, 6, out, d);
}
StatusCode Node::Create(int id, double x, double y, double z, int dof_count,
                        Node** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("Node::Create", d);
  *out = nullptr;
  if (!model_storage::ValidId(id) || !std::isfinite(x) || !std::isfinite(y) ||
      !std::isfinite(z) || (dof_count != 3 && dof_count != 6))
    return model_storage::Invalid("Node::Create", d);
  Vector* raw = nullptr;
  const auto status = Vector::Create(2 * dof_count, &raw, d);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<Vector> results(raw);
  *out = new (std::nothrow) Node(id, x, y, z, results.get());
  if (!*out) return model_storage::AllocationFailure("Node::Create", d);
  results.release();
  return StatusCode::kOk;
}
int Node::dof_count() const { return results_ ? results_->Size() / 2 : 0; }
int Node::Id() const { return id_; }
StatusCode Node::Coordinate(int axis, double* out, Diagnostic* d) const {
  if (axis < 0 || axis > 2 || out == nullptr) {
    return model_storage::Invalid("Node::Coordinate", d);
  }
  *out = axis == 0 ? x_ : (axis == 1 ? y_ : z_);
  return StatusCode::kOk;
}
StatusCode Node::UpdateResults(const Vector* values, Diagnostic* d) {
  if (!results_)
    return model_storage::Invalid("Thermal node has no mechanical results", d);
  return results_->CopyFrom(values, d);
}
const Vector* Node::Results() const { return results_; }
StatusCode Node::CreateThermal(int id, double x, double y, double z, int count,
                               Node** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("Node::CreateThermal output", d);
  *out = nullptr;
  if (!model_storage::ValidId(id) || !std::isfinite(x) || !std::isfinite(y) ||
      !std::isfinite(z))
    return model_storage::Invalid("Node::CreateThermal ID/coordinates", d);
  NodalTemperatureResult* raw = nullptr;
  const auto status = NodalTemperatureResult::Create(count, &raw, d);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<NodalTemperatureResult> results(raw);
  *out = new (std::nothrow) Node(id, x, y, z, nullptr);
  if (!*out) return model_storage::AllocationFailure("Node::CreateThermal", d);
  (*out)->temperature_results_ = results.release();
  return StatusCode::kOk;
}
int Node::temperature_dof_count() const {
  return temperature_results_ ? temperature_results_->DofCount() : 0;
}
const NodalTemperatureResult* Node::TemperatureResults() const {
  return temperature_results_;
}
StatusCode Node::UpdateTemperatureResults(const Vector* nt, const Vector* rfl,
                                         Diagnostic* d) {
  if (!temperature_results_)
    return model_storage::Invalid("Mechanical node has no thermal results", d);
  return temperature_results_->Update(nt, rfl, d);
}

}  // namespace fesa
