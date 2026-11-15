#include "analysis/nonlinear_condition_validator.h"

#include <cmath>
#include <limits>
#include <memory>
#include <new>

#include "fesa/elements/beam3d.h"
#include "fesa/math/vector.h"
#include "model/model_storage.h"

namespace fesa {
namespace {
using PhysicalStates = std::map<std::pair<int, int>, int>;
using PhysicalDefinitions =
    std::map<std::pair<int, int>, const BoundaryDefinition*>;

StatusCode Unsupported(const BoundaryDefinition* definition, int node, int step,
                       const char* operation, Diagnostic* d) {
  if (d) {
    d->source_file = definition->SourceFile();
    d->path_kind = DiagnosticPathKind::kLiteralFile;
    d->source_line = definition->SourceLine();
    d->first_source_line = 0;
    d->duplicate_boundary = false;
    d->step = step;
    d->node_id = node;
    d->element_id = 0;
    d->condition_dof = definition->Dof();
  }
  return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                       operation, d);
}

using RotationAxes = std::map<int, Vector*>;
void DeleteAxes(RotationAxes* axes) {
  if (!axes) return;
  for (auto& item : *axes) delete item.second;
  delete axes;
}
StatusCode CloneAxes(const RotationAxes* source, RotationAxes* target,
                     Diagnostic* d) {
  for (const auto& item : *source) {
    Vector* raw = nullptr;
    auto status = Vector::Create(3, &raw, d);
    std::unique_ptr<Vector> axis(raw);
    if (status != StatusCode::kOk) return status;
    status = axis->CopyFrom(item.second, d);
    if (status != StatusCode::kOk) return status;
    target->emplace(item.first, axis.get());
    axis.release();
  }
  return StatusCode::kOk;
}
StatusCode ValidateRotations(const PhysicalDefinitions* definitions, int step,
                             bool beam, RotationAxes* axes, Diagnostic* d) {
  std::map<int, int> counts;
  for (const auto& item : *definitions)
    if (item.first.second >= 4) ++counts[item.first.first];
  for (const auto& item : *definitions) {
    if (item.first.second < 4) continue;
    const int node = item.first.first;
    if (counts.at(node) != 3)
      return Unsupported(item.second, node, step,
                         "Nonlinear rotations require UR1 UR2 UR3 together", d);
    if (item.first.second != 4) continue;
    Vector* raw = nullptr;
    auto status = Vector::Create(3, &raw, d);
    std::unique_ptr<Vector> axis(raw);
    if (status != StatusCode::kOk) return status;
    int nonzero = 0;
    for (int component = 0; component < 3; ++component) {
      const double value = definitions->at({node, 4 + component})->Value();
      if (value != 0) ++nonzero;
      status = axis->Set(component, value, d);
      if (status != StatusCode::kOk) return status;
    }
    if (!beam && nonzero > 1)
      return Unsupported(item.second, node, step,
                         "Nonlinear prescribed rotation axis changed", d);
    double norm = 0;
    status = axis->Norm2(&norm, d);
    if (status != StatusCode::kOk) return status;
    if (norm == 0) continue;  // Zero does not establish or erase an axis.
    status = axis->Scale(1 / norm, d);
    if (status != StatusCode::kOk) return status;
    const auto prior = axes->find(node);
    if (prior == axes->end()) {
      axes->emplace(node, axis.get());
      axis.release();
      continue;
    }
    double projection = 0;
    status = axis->Dot(prior->second, &projection, d);
    if (status == StatusCode::kOk)
      status = axis->AddScaled(-projection, prior->second, d);
    if (status == StatusCode::kOk) status = axis->Norm2(&norm, d);
    if (status != StatusCode::kOk) return status;
    // Dimensionless unit-axis orthogonal residual: rounding tolerance only.
    if (norm > 64 * std::numeric_limits<double>::epsilon())
      return Unsupported(item.second, node, step,
                         "Nonlinear prescribed rotation axis changed", d);
  }
  return StatusCode::kOk;
}
}  // namespace

StatusCode NonlinearConditionValidator::Create(
    const Domain* domain, NonlinearConditionValidator** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("Condition validator output", d);
  *out = nullptr;
  if (!domain) return model_storage::Invalid("Condition validator Domain", d);
  try {
    std::unique_ptr<NonlinearConditionValidator> result(
        new NonlinearConditionValidator);
    result->states_ = new PhysicalStates;
    result->axes_ = new RotationAxes;
    result->beam_ = !domain->Elements()->empty();
    for (const auto& item : *domain->Elements())
      result->beam_ =
          result->beam_ && dynamic_cast<const Beam3D*>(item.second) != nullptr;
    for (const auto& item : *domain->Steps())
      result->beam_ = result->beam_ && !item.second->IsDynamic();
    result->domain_ = domain;
    *out = result.release();
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return model_storage::AllocationFailure("Condition validator storage", d);
  }
}
NonlinearConditionValidator::~NonlinearConditionValidator() {
  delete states_;
  DeleteAxes(axes_);
}

StatusCode NonlinearConditionValidator::ValidateStep(
    const EffectiveConditions* conditions, Diagnostic* d) {
  if (!conditions || conditions->SourceDomain() != domain_ ||
      conditions->StepNumber() != next_step_ ||
      !conditions->ReleasedBoundaries())
    return model_storage::Invalid("Condition validator snapshot/order", d);
  try {
    // Clone owned input history; numerical axis copies use the common API.
    PhysicalStates states(*states_);
    std::unique_ptr<RotationAxes, decltype(&DeleteAxes)> axes(new RotationAxes,
                                                              DeleteAxes);
    auto clone_status = CloneAxes(axes_, axes.get(), d);
    if (clone_status != StatusCode::kOk) return clone_status;
    PhysicalDefinitions current;
    for (const auto* row : *conditions->Boundaries())
      current.emplace(std::make_pair(row->node_id, row->definition->Dof()),
                      row->definition);
    if (next_step_ == 1) {
      PhysicalDefinitions initial;
      for (const auto* row : *conditions->Boundaries())
        if (row->physical_previous_definition)
          initial.emplace(std::make_pair(row->node_id, row->definition->Dof()),
                          row->physical_previous_definition);
      for (const auto* row : *conditions->ReleasedBoundaries())
        initial.emplace(std::make_pair(row->node_id, row->definition->Dof()),
                        row->definition);
      const auto status =
          ValidateRotations(&initial, next_step_, beam_, axes.get(), d);
      if (status != StatusCode::kOk) return status;
      for (const auto& item : initial) states[item.first] = 1;
    }
    for (const auto* row : *conditions->ReleasedBoundaries())
      states[std::make_pair(row->node_id, row->definition->Dof())] = 2;
    for (const auto& item : current) {
      const auto prior = states.find(item.first);
      if (prior != states.end() && prior->second == 2)
        return Unsupported(item.second, item.first.first, next_step_,
                           "Nonlinear released constraint reapplication", d);
      states[item.first] = 1;
    }
    const auto status =
        ValidateRotations(&current, next_step_, beam_, axes.get(), d);
    if (status != StatusCode::kOk) return status;
    states_->swap(states);
    axes_->swap(*axes.get());
    ++next_step_;
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return model_storage::AllocationFailure("Condition validator transition",
                                            d);
  }
}
}  // namespace fesa
