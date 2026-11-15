#include "fesa/model/load_definition.h"

#include <cmath>
#include <memory>
#include <new>

#include "model/model_storage.h"

namespace fesa {

LoadDefinition::LoadDefinition(int id, char* scope, ConditionTarget target,
    int node_id, char* name, int dof, double value, char* source_file,
    int source_line, int step, ConditionOperation operation)
    : id_(id), scope_(scope), target_(target), node_id_(node_id), name_(name),
      dof_(dof), value_(value), source_file_(source_file),
      source_line_(source_line), step_(step), operation_(operation) {}
LoadDefinition::~LoadDefinition() {
  delete[] original_target_;
  delete[] scope_;
  delete[] name_;
  delete[] source_file_;
}
StatusCode LoadDefinition::Create(int id, const char* scope,
    ConditionTarget target, int node_id, const char* name, int dof,
    double value, const char* source_file, int source_line, int step,
    ConditionOperation operation, LoadDefinition** out, Diagnostic* d) {
  if (out == nullptr) {
    return model_storage::Invalid("LoadDefinition::Create", d);
  }
  *out = nullptr;
  if (!model_storage::ValidId(id) || scope == nullptr || name == nullptr
      || source_file == nullptr || source_line <= 0 || step < 0
      || dof < 1 || dof > 6 || !std::isfinite(value)
      || (target != ConditionTarget::kNode
          && target != ConditionTarget::kNodeSet)
      || (operation != ConditionOperation::kMod
          && operation != ConditionOperation::kNew)
      || (target == ConditionTarget::kNode && !model_storage::ValidId(node_id))
      || (target == ConditionTarget::kNodeSet && (node_id != 0 || !name[0]))) {
    return model_storage::Invalid("LoadDefinition::Create", d);
  }
  std::unique_ptr<char[]> scope_copy(model_storage::CopyText(scope));
  std::unique_ptr<char[]> name_copy(model_storage::CopyText(name));
  std::unique_ptr<char[]> file_copy(model_storage::CopyText(source_file));
  if (scope_copy && name_copy && file_copy) {
    *out = new (std::nothrow) LoadDefinition(id, scope_copy.get(), target,
        node_id, name_copy.get(), dof, value, file_copy.get(), source_line,
        step, operation);
  }
  if (*out == nullptr) {
    return model_storage::AllocationFailure("LoadDefinition::Create", d);
  }
  scope_copy.release();
  name_copy.release();
  file_copy.release();
  return StatusCode::kOk;
}
int LoadDefinition::Id() const { return id_; }
const char* LoadDefinition::Scope() const { return scope_; }
ConditionTarget LoadDefinition::Target() const { return target_; }
int LoadDefinition::NodeId() const { return node_id_; }
const char* LoadDefinition::Name() const { return name_; }
int LoadDefinition::Dof() const { return dof_; }
double LoadDefinition::Value() const { return value_; }
const char* LoadDefinition::SourceFile() const { return source_file_; }
int LoadDefinition::SourceLine() const { return source_line_; }
int LoadDefinition::Step() const { return step_; }
ConditionOperation LoadDefinition::Operation() const { return operation_; }

StatusCode LoadDefinition::SetOriginalTarget(const char* text, Diagnostic* d) {
  if (!text)
    return model_storage::Invalid("LoadDefinition::SetOriginalTarget", d);
  std::unique_ptr<char[]> copy(model_storage::CopyText(text));
  if (!copy)
    return model_storage::AllocationFailure("LoadDefinition::SetOriginalTarget",
                                            d);
  delete[] original_target_;
  original_target_ = copy.release();
  return StatusCode::kOk;
}
const char* LoadDefinition::OriginalTarget() const {
  return original_target_ ? original_target_ : "";
}

StatusCode LoadDefinition::BindDynamicAmplitude(const Amplitude* amplitude,
                                                Diagnostic* d) {
  if (!amplitude || amplitude_)
    return model_storage::Invalid("LoadDefinition amplitude binding", d);
  amplitude_ = amplitude;
  return StatusCode::kOk;
}
const Amplitude* LoadDefinition::DynamicAmplitude() const { return amplitude_; }

}  // namespace fesa
