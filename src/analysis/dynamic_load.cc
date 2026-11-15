#include "analysis/dynamic_load.h"

#include <climits>
#include <cmath>
#include <memory>
#include <new>
#include <vector>

namespace fesa {
namespace {
void LoadContext(const EffectiveLoad* load, int step, Diagnostic* d) {
  if (!d) return;
  d->source_file = load->definition->SourceFile();
  d->source_line = load->definition->SourceLine();
  d->step = step;
  d->node_id = load->node_id;
  d->condition_dof = load->definition->Dof();
}
StatusCode EvaluateClLoads(bool derivative,
                           const EffectiveConditions* conditions,
                           const DofManager* dofs, double step_time,
                           double total_time, Vector** external,
                           Vector** origins, Diagnostic* d) {
  if (external) *external = nullptr;
  if (origins) *origins = nullptr;
  if (!conditions || !dofs || !external || !origins || external == origins ||
      !std::isfinite(step_time) || !std::isfinite(total_time) ||
      step_time < 0 || total_time < step_time)
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation, "Dynamic CLOAD inputs",
                         d);
  const auto* domain = conditions->SourceDomain();
  const StepDefinition* step = nullptr;
  auto status = domain->FindStepDefinition(conditions->StepNumber(), &step, d);
  if (status != StatusCode::kOk) return status;
  if (!step->IsDynamic() || domain->LineLoadProgram())
    return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                         "Dynamic CLOAD procedure/load support", d);
  if (step_time > step->Period() ||
      domain->Nodes()->size() != static_cast<size_t>(dofs->NodeCount()))
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation, "Dynamic CLOAD time/map",
                         d);
  for (const auto& node : *domain->Nodes()) {
    if (dofs->NodeDofCount(node.first) != 3 &&
        dofs->NodeDofCount(node.first) != 6)
      return ReportFailure(StatusCode::kUnsupported,
                           DiagnosticStage::kValidation,
                           "Dynamic CLOAD 3N/6N map", d);
  }
  const auto count = conditions->Loads()->size();
  if (count > INT_MAX)
    return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                         "Dynamic CLOAD origin count", d);
  try {
    Vector* raw = nullptr;
    status = Vector::Create(dofs->GlobalCount(), &raw, d);
    if (status != StatusCode::kOk) return status;
    std::unique_ptr<Vector> assembled(raw);
    std::unique_ptr<Vector> values;
    if (count > 0) {
      status = Vector::Create(static_cast<int>(count), &raw, d);
      if (status != StatusCode::kOk) return status;
      values.reset(raw);
    }
    std::vector<int> indices(count);
    int index = 0;
    for (const auto* load : *conditions->Loads()) {
      const auto* definition = load->definition;
      const int global = dofs->GlobalDof(load->node_id, definition->Dof());
      if (global < 0) {
        LoadContext(load, step->Id(), d);
        return ReportFailure(StatusCode::kInvalidArgument,
                             DiagnosticStage::kValidation,
                             "Dynamic CLOAD target mapping", d);
      }
      double factor = derivative ? 0 : 1;
      const auto* amplitude = definition->DynamicAmplitude();
      if (load->active && amplitude) {
        const StepDefinition* defining_step = nullptr;
        status =
            domain->FindStepDefinition(definition->Step(), &defining_step, d);
        if (status == StatusCode::kOk) {
          AmplitudeEvaluationContext context;
          context.definition_step = definition->Step();
          context.current_step = step->Id();
          context.definition_period = defining_step->Period();
          context.step_time = step_time;
          context.total_time = total_time;
          status = derivative ? amplitude->EvaluateRightDerivativeInStep(
                                    &context, &factor, d)
                              : amplitude->EvaluateInStep(&context, &factor, d);
        }
        if (status != StatusCode::kOk) {
          LoadContext(load, step->Id(), d);
          return status;
        }
      }
      // A single load's physical scalar coefficient, before vector assembly.
      const double value = load->active ? load->value * factor : 0;
      if (!std::isfinite(value)) {
        LoadContext(load, step->Id(), d);
        return ReportFailure(StatusCode::kLibraryFailure,
                             DiagnosticStage::kComputation,
                             derivative ? "Dynamic CLOAD rate overflow"
                                        : "Dynamic CLOAD value overflow",
                             d);
      }
      status = values->Set(index, value, d);
      if (status != StatusCode::kOk) {
        LoadContext(load, step->Id(), d);
        return status;
      }
      indices[index++] = global;
    }
    if (count > 0) {
      status = assembled->ScatterAdd(values.get(), indices.data(),
                                     static_cast<int>(count), d);
      if (status != StatusCode::kOk) return status;
    }
    *external = assembled.release();
    *origins = values.release();
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return ReportFailure(StatusCode::kAllocationFailure,
                         DiagnosticStage::kAllocation,
                         "Dynamic CLOAD allocation", d);
  }
}
}  // namespace
StatusCode EvaluateDynamicClLoads(const EffectiveConditions* conditions,
                                  const DofManager* dofs, double step_time,
                                  double total_time, Vector** external,
                                  Vector** origins, Diagnostic* d) {
  return EvaluateClLoads(false, conditions, dofs, step_time, total_time,
                         external, origins, d);
}
StatusCode EvaluateDynamicClLoadRates(const EffectiveConditions* conditions,
                                      const DofManager* dofs, double step_time,
                                      double total_time, Vector** external_rate,
                                      Vector** origin_rates, Diagnostic* d) {
  return EvaluateClLoads(true, conditions, dofs, step_time, total_time,
                         external_rate, origin_rates, d);
}
}  // namespace fesa
