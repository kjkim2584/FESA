#include "fesa/model/boundary_definition.h"

#include <cmath>
#include <memory>
#include <new>

#include "model/model_storage.h"

namespace fesa {

BoundaryDefinition::BoundaryDefinition(int id, char* scope,
    ConditionTarget target,
    int node_id, char* name, int dof, double value, char* source_file,
    int source_line, int step, ConditionOperation operation)
    : id_(id), scope_(scope), target_(target), node_id_(node_id), name_(name),
      dof_(dof), value_(value), source_file_(source_file),
      source_line_(source_line), step_(step), operation_(operation) {}
BoundaryDefinition::~BoundaryDefinition() {
  delete[] original_target_;
  delete[] scope_;
  delete[] name_;
  delete[] source_file_;
}
StatusCode BoundaryDefinition::Create(int id, const char* scope,
    ConditionTarget target, int node_id, const char* name, int dof,
    double value, const char* source_file, int source_line, int step,
    ConditionOperation operation, BoundaryDefinition** out, Diagnostic* d) {
  if (out == nullptr) {
    return model_storage::Invalid("BoundaryDefinition::Create", d);
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
    return model_storage::Invalid("BoundaryDefinition::Create", d);
  }
  std::unique_ptr<char[]> scope_copy(model_storage::CopyText(scope));
  std::unique_ptr<char[]> name_copy(model_storage::CopyText(name));
  std::unique_ptr<char[]> file_copy(model_storage::CopyText(source_file));
  if (scope_copy && name_copy && file_copy) {
    *out = new (std::nothrow) BoundaryDefinition(id, scope_copy.get(), target,
        node_id, name_copy.get(), dof, value, file_copy.get(), source_line,
        step, operation);
  }
  if (*out == nullptr) {
    return model_storage::AllocationFailure("BoundaryDefinition::Create", d);
  }
  scope_copy.release();
  name_copy.release();
  file_copy.release();
  return StatusCode::kOk;
}
int BoundaryDefinition::Id() const { return id_; }
const char* BoundaryDefinition::Scope() const { return scope_; }
ConditionTarget BoundaryDefinition::Target() const { return target_; }
int BoundaryDefinition::NodeId() const { return node_id_; }
const char* BoundaryDefinition::Name() const { return name_; }
int BoundaryDefinition::Dof() const { return dof_; }
double BoundaryDefinition::Value() const { return value_; }
const char* BoundaryDefinition::SourceFile() const { return source_file_; }
int BoundaryDefinition::SourceLine() const { return source_line_; }
int BoundaryDefinition::Step() const { return step_; }
ConditionOperation BoundaryDefinition::Operation() const { return operation_; }

StatusCode BoundaryDefinition::SetOriginalTarget(const char* text,
                                                 Diagnostic* d) {
  if (!text)
    return model_storage::Invalid("BoundaryDefinition::SetOriginalTarget", d);
  std::unique_ptr<char[]> copy(model_storage::CopyText(text));
  if (!copy)
    return model_storage::AllocationFailure(
        "BoundaryDefinition::SetOriginalTarget", d);
  delete[] original_target_;
  original_target_ = copy.release();
  return StatusCode::kOk;
}
const char* BoundaryDefinition::OriginalTarget() const {
  return original_target_ ? original_target_ : "";
}

}  // namespace fesa
