#include "fesa/model/step_definition.h"

#include <cmath>
#include <memory>
#include <new>

#include "model/model_storage.h"

namespace fesa {

StepDefinition::StepDefinition(int id, char* name, double period,
    std::vector<const LoadDefinition*>* loads,
    std::vector<const BoundaryDefinition*>* boundaries)
    : id_(id), name_(name), period_(period), loads_(loads),
      boundaries_(boundaries) {}
StepDefinition::~StepDefinition() {
  delete heat_transfer_;
  delete thermal_conditions_;
  delete nonlinear_;
  delete loads_;
  delete boundaries_;
  delete[] name_;
}
StatusCode StepDefinition::Create(int id, const char* name, double period,
                                  StepDefinition** out, Diagnostic* d) {
  if (out == nullptr) {
    return model_storage::Invalid("StepDefinition::Create", d);
  }
  *out = nullptr;
  try {
    if (!model_storage::ValidId(id) || name == nullptr ||
        !std::isfinite(period) || period <= 0) {
      return model_storage::Invalid("StepDefinition::Create", d);
    }
    std::unique_ptr<char[]> copy(model_storage::CopyText(name));
    // Explicit size construction can propagate Debug iterator-proxy allocation
    // failure; MSVC's noexcept default vector constructor cannot.
    std::unique_ptr<std::vector<const LoadDefinition*>> loads(
        new (std::nothrow) std::vector<const LoadDefinition*>(0));
    std::unique_ptr<std::vector<const BoundaryDefinition*>> boundaries(
        new (std::nothrow) std::vector<const BoundaryDefinition*>(0));
    if (copy && loads && boundaries) {
      *out = new (std::nothrow) StepDefinition(id, copy.get(), period,
                                              loads.get(), boundaries.get());
    }
    if (*out == nullptr) {
      return model_storage::AllocationFailure("StepDefinition::Create", d);
    }
    copy.release();
    loads.release();
    boundaries.release();
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return model_storage::AllocationFailure("StepDefinition::Create", d);
  }
}
StatusCode StepDefinition::CreateDynamic(int id, const char* name,
                                         double period, StepDefinition** out,
                                         Diagnostic* d) {
  const auto status = Create(id, name, period, out, d);
  if (status == StatusCode::kOk) (*out)->dynamic_ = true;
  return status;
}
StatusCode StepDefinition::CreateHeatTransfer(
    int id, const char* name, const HeatTransferStepParameters* parameters,
    StepDefinition** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("Thermal step output", d);
  *out = nullptr;
  HeatTransferStepParameters* settings = nullptr;
  auto status = CloneHeatTransferStepParameters(parameters, &settings, d);
  std::unique_ptr<HeatTransferStepParameters> settings_owner(settings);
  if (status != StatusCode::kOk) return status;
  StepDefinition* raw = nullptr;
  status = Create(id, name, settings->period, &raw, d);
  std::unique_ptr<StepDefinition> step(raw);
  if (status != StatusCode::kOk) return status;
  try {
    step->thermal_conditions_ =
        new (std::nothrow) std::vector<const ThermalConditionDefinition*>(0);
  } catch (const std::bad_alloc&) {
    return model_storage::AllocationFailure("Thermal step condition list", d);
  }
  if (!step->thermal_conditions_)
    return model_storage::AllocationFailure("Thermal step condition list", d);
  step->heat_transfer_ = settings_owner.release();
  *out = step.release();
  return StatusCode::kOk;
}
bool StepDefinition::IsHeatTransfer() const {
  return heat_transfer_ != nullptr;
}
const HeatTransferStepParameters* StepDefinition::HeatTransfer() const {
  return heat_transfer_;
}
StatusCode StepDefinition::AddThermalCondition(
    const ThermalConditionDefinition* condition, Diagnostic* d) {
  if (!IsHeatTransfer() || !condition ||
      condition->specification()->step != id_)
    return model_storage::Invalid("Thermal step condition identity", d);
  try {
    thermal_conditions_->push_back(condition);
  } catch (const std::bad_alloc&) {
    return model_storage::AllocationFailure("Thermal step condition append", d);
  }
  return StatusCode::kOk;
}
const std::vector<const ThermalConditionDefinition*>*
StepDefinition::ThermalConditions() const {
  return thermal_conditions_;
}
StatusCode StepDefinition::CreateNonlinearDynamic(
    int id, const char* name, double period, StepDefinition** out, Diagnostic* d) {
  const auto status = CreateDynamic(id, name, period, out, d);
  if (status == StatusCode::kOk) (*out)->finite_dynamic_ = true;
  return status;
}
bool StepDefinition::HasFiniteGeometry() const {
  return finite_dynamic_ || nonlinear_ != nullptr;
}
bool StepDefinition::IsDynamic() const { return dynamic_; }
StatusCode StepDefinition::CreateNonlinearStatic(
    int id, const char* name, const NonlinearStaticParameters* parameters,
    StepDefinition** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("Nonlinear step output", d);
  *out = nullptr;
  auto status = ValidateNonlinearStaticParameters(parameters, d);
  if (status != StatusCode::kOk) return status;
  StepDefinition* raw = nullptr;
  status = Create(id, name, parameters->period, &raw, d);
  std::unique_ptr<StepDefinition> step(raw);
  if (status != StatusCode::kOk) return status;
  step->nonlinear_ = new (std::nothrow) NonlinearStaticParameters;
  if (!step->nonlinear_)
    return model_storage::AllocationFailure("Nonlinear step parameters", d);
  auto* copy = step->nonlinear_;
  copy->period = parameters->period;
  copy->initial_increment = parameters->initial_increment;
  copy->minimum_increment = parameters->minimum_increment;
  copy->has_maximum = parameters->has_maximum;
  copy->maximum_increment = parameters->maximum_increment;
  copy->increment_limit = parameters->increment_limit;
  copy->extrapolate = parameters->extrapolate;
  *out = step.release();
  return StatusCode::kOk;
}
const NonlinearStaticParameters* StepDefinition::NonlinearParameters() const {
  return nonlinear_;
}
StatusCode StepDefinition::AddLoad(const LoadDefinition* load, Diagnostic* d) {
  if (load == nullptr || IsHeatTransfer()) {
    return model_storage::Invalid("StepDefinition::AddLoad", d);
  }
  try {
    loads_->push_back(load);
  } catch (const std::bad_alloc&) {
    return model_storage::AllocationFailure("StepDefinition::AddLoad", d);
  }
  return StatusCode::kOk;
}
StatusCode StepDefinition::AddBoundary(const BoundaryDefinition* boundary,
                                      Diagnostic* d) {
  if (boundary == nullptr || IsHeatTransfer()) {
    return model_storage::Invalid("StepDefinition::AddBoundary", d);
  }
  try {
    boundaries_->push_back(boundary);
  } catch (const std::bad_alloc&) {
    return model_storage::AllocationFailure("StepDefinition::AddBoundary", d);
  }
  return StatusCode::kOk;
}
int StepDefinition::Id() const { return id_; }
const char* StepDefinition::Name() const { return name_; }
double StepDefinition::Period() const { return period_; }
const std::vector<const LoadDefinition*>* StepDefinition::Loads() const {
  return loads_;
}
const std::vector<const BoundaryDefinition*>*
StepDefinition::Boundaries() const {
  return boundaries_;
}

StatusCode StepDefinition::RecordConditionBlock(bool load,
                                                ConditionOperation operation,
                                                int line, Diagnostic* d) {
  int* stored_line = load ? &load_block_line_ : &boundary_block_line_;
  ConditionOperation* stored_op =
      load ? &load_operation_ : &boundary_operation_;
  if (line <= 0 || (*stored_line != 0 && *stored_op != operation))
    return model_storage::Invalid("StepDefinition mixed condition operations",
                                  d);
  if (*stored_line == 0) *stored_line = line;
  *stored_op = operation;
  return StatusCode::kOk;
}
bool StepDefinition::HasLoadBlock() const { return load_block_line_ != 0; }
bool StepDefinition::HasBoundaryBlock() const {
  return boundary_block_line_ != 0;
}
ConditionOperation StepDefinition::LoadOperation() const {
  return load_operation_;
}
ConditionOperation StepDefinition::BoundaryOperation() const {
  return boundary_operation_;
}
int StepDefinition::LoadBlockLine() const { return load_block_line_; }
int StepDefinition::BoundaryBlockLine() const { return boundary_block_line_; }
StatusCode StepDefinition::RecordThermalLoadBlock(ThermalLoadFamily family,
                                                  ConditionOperation operation,
                                                  int line, Diagnostic* d) {
  const int index = static_cast<int>(family);
  if (!IsHeatTransfer() || index < 0 || index >= 6 || line <= 0 ||
      (operation != ConditionOperation::kMod &&
       operation != ConditionOperation::kNew) ||
      (thermal_load_lines_[index] &&
       thermal_load_operations_[index] != operation))
    return model_storage::Invalid("Thermal keyword block operation", d);
  if (!thermal_load_lines_[index]) thermal_load_lines_[index] = line;
  thermal_load_operations_[index] = operation;
  return StatusCode::kOk;
}
ConditionOperation StepDefinition::ThermalLoadOperation(
    ThermalLoadFamily family) const {
  const int index = static_cast<int>(family);
  return index >= 0 && index < 6 ? thermal_load_operations_[index]
                                 : ConditionOperation::kMod;
}

}  // namespace fesa
