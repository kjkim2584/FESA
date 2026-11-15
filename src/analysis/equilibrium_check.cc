#include "analysis/equilibrium_check.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <new>

#include "adapters/openmp/thread_runtime.h"
#include "analysis/field_decision.h"
#include "analysis/nonlinear_convergence.h"
#include "analysis/numerical_support.h"
#include "model/model_storage.h"
namespace fesa {
namespace {
using numerical_support::Buffers;
bool Nonnegative(double value) { return std::isfinite(value) && value >= 0; }
bool ValidControls(const FieldControls* controls) {
  return controls && Nonnegative(controls->rn) && Nonnegative(controls->cn) &&
         Nonnegative(controls->epsilon) && Nonnegative(controls->cepsilon) &&
         Nonnegative(controls->rl) && Nonnegative(controls->cf) &&
         Nonnegative(controls->epsilon_l) && Nonnegative(controls->epsilon_d) &&
         (controls->qu_explicit
              ? std::isfinite(controls->qu) && controls->qu > 0
              : std::isfinite(controls->q0) && controls->q0 > 0);
}
StatusCode NumericalFailure(Diagnostic* d) {
  return ReportFailure(StatusCode::kLibraryFailure,
                       DiagnosticStage::kComputation,
                       "Equilibrium threshold arithmetic", d);
}
void ResetField(FieldResidual* field) {
  field->dof_count = field->contribution_count = 0;
  field->maximum_node = field->maximum_dof = 0;
  field->maximum_element = 0;
  field->maximum_internal_end = -1;
  field->maximum = field->l2 = field->rms = 0;
  field->spatial_mean = field->raw_mean = field->used_mean = 0;
  field->maximum_flux = field->phi = 0;
  field->increment_max = field->correction_max = 0;
  field->residual_limit = field->correction_limit = 0;
  field->converted = false;
  field->passed = true;
  field->branch = ResidualBranch::kNoDofs;
}
StatusCode MakeVector(int size, std::unique_ptr<Vector>* out, Diagnostic* d) {
  Vector* raw = nullptr;
  const StatusCode status = Vector::Create(size, &raw, d);
  out->reset(raw);
  return status;
}
StatusCode FluxCount(const AnalysisModel* model,
                     const EffectiveConditions* conditions, int field, int* out,
                     Diagnostic* d) {
  const auto maximum =
      static_cast<std::size_t>(std::numeric_limits<int>::max());
  std::size_t count = 0;
  for (const auto& entry : *model->Elements()) {
    const Element* element = entry.second;
    if (!element || element->node_count() <= 0 ||
        element->node_count() > Buffers::kMaxElementDofs / 3 ||
        (element->dofs_per_node() != 3 && element->dofs_per_node() != 6) ||
        element->dof_count() !=
            element->node_count() * element->dofs_per_node() ||
        element->dof_count() > Buffers::kMaxElementDofs)
      return model_storage::Invalid("Equilibrium element layout", d);
    const int extra = field == 0 || element->dofs_per_node() == 6
                          ? 3 * element->node_count()
                          : 0;
    if (extra < 0 || static_cast<std::size_t>(extra) > maximum - count)
      return model_storage::Invalid("Equilibrium contribution count", d);
    count += extra;
  }
  for (const auto* load : *conditions->Loads())
    if (load->contributes_flux && (load->definition->Dof() - 1) / 3 == field) {
      if (count == maximum)
        return model_storage::Invalid("Equilibrium contribution count", d);
      ++count;
    }
  *out = static_cast<int>(count);
  return StatusCode::kOk;
}
struct ForceTaskContext {
  Buffers* buffers = nullptr;
  const Vector* candidate = nullptr;
  Vector* force_flux = nullptr;
  Vector* moment_flux = nullptr;
  Vector* lengths = nullptr;
  int first_element = 0;
  const int* force_offsets = nullptr;
  const int* moment_offsets = nullptr;
  const BeamLineLoadState* line_loads = nullptr;
};
StatusCode ComputeForce(int item, int worker, void* context, Diagnostic* d) {
  auto* task = static_cast<ForceTaskContext*>(context);
  Buffers* b = task->buffers;
  auto* local = b->workers[worker];
  const Element* element = b->elements[item];
  d->element_id = element->Id();
  const int size = element->dof_count();
  std::unique_ptr<numerical_support::WorkerBuffers> owned;
  if (local->u->Size() != size) {
    numerical_support::WorkerBuffers* raw = nullptr;
    const auto status =
        numerical_support::WorkerBuffers::Create(local->element, &raw, d, size);
    owned.reset(raw);
    if (status != StatusCode::kOk) return status;
    local = raw;
  }
  const int* global = b->global_indices + b->element_dofs * item;
  auto status = local->u->GatherFrom(task->candidate, global, size, 0, d);
  if (status != StatusCode::kOk) return status;
  if (task->line_loads) {
    status = numerical_support::EvaluateLineLoad(task->line_loads, element,
                                                 local, d);
    if (status != StatusCode::kOk) return status;
    status = element->ComputeLoadedInternalForce(
        local->u, local->line_global, local->element, local->force, d);
  } else {
    status = element->ComputeInternalForce(local->u, local->element,
                                           local->force, d);
  }
  if (status != StatusCode::kOk) return status;
  int indices[Buffers::kMaxElementDofs] = {};
  for (int field = 0; field < 2; ++field) {
    if (field == 1 && element->dofs_per_node() == 3) continue;
    int count = 0;
    for (int node = 0; node < element->node_count(); ++node)
      for (int component = 0; component < 3; ++component)
        indices[count++] =
            node * element->dofs_per_node() + 3 * field + component;
    Vector* flux = field == 0 ? task->force_flux : task->moment_flux;
    const int offset =
        field == 0 ? task->force_offsets[item] : task->moment_offsets[item];
    status = flux->GatherFrom(local->force, indices, count, offset, d);
    if (status != StatusCode::kOk) return status;
  }
  const int index = task->first_element + item;
  double length = 0;
  status = element->CharacteristicLength(&length, d);
  if (status != StatusCode::kOk) return status;
  status = task->lengths->Set(index, length, d);
  if (status != StatusCode::kOk) return status;
  for (int i = 0; i < size; ++i) indices[i] = i;
  return b->vector_values->GatherFrom(local->force, indices, size,
                                      b->vector_offsets[item], d);
}
StatusCode FieldNorms(const DofManager* dofs, const Vector* residual, int field,
                      FieldResidual* result, Diagnostic* d) {
  int count = 0;
  for (int row = 0; row < dofs->NodeCount(); ++row) {
    const int id = dofs->NodeId(row);
    for (int component = 3 * field + 1;
         component <= std::min(3 * field + 3, dofs->NodeDofCount(id));
         ++component)
      if (dofs->Equation(dofs->GlobalDof(id, component)) >= 0) ++count;
  }
  result->dof_count = count;
  if (count == 0) return StatusCode::kOk;
  std::unique_ptr<int[]> indices(new (std::nothrow) int[count]);
  if (!indices) return model_storage::AllocationFailure("Field indices", d);
  int j = 0;
  for (int row = 0; row < dofs->NodeCount(); ++row) {
    const int id = dofs->NodeId(row);
    for (int component = 3 * field + 1;
         component <= std::min(3 * field + 3, dofs->NodeDofCount(id));
         ++component) {
      const int global = dofs->GlobalDof(id, component);
      if (dofs->Equation(global) >= 0) indices[j++] = global;
    }
  }
  std::unique_ptr<Vector> values;
  StatusCode status = MakeVector(count, &values, d);
  if (status != StatusCode::kOk) return status;
  status = values->GatherFrom(residual, indices.get(), count, 0, d);
  if (status != StatusCode::kOk) return status;
  int maximum_index = 0;
  status = values->MaximumAbs(&result->maximum, &maximum_index, d);
  if (status != StatusCode::kOk) return status;
  for (int row = 0; row < dofs->NodeCount(); ++row) {
    const int id = dofs->NodeId(row);
    const int component = indices[maximum_index] - dofs->GlobalDof(id, 1) + 1;
    if (component > 0 && component <= dofs->NodeDofCount(id)) {
      result->maximum_node = id;
      result->maximum_dof = component;
      break;
    }
  }
  status = values->Norm2(&result->l2, d);
  if (status != StatusCode::kOk) return status;
  result->rms = result->l2 / std::sqrt(static_cast<double>(count));
  return StatusCode::kOk;
}
StatusCode IncrementMaximum(const DofManager* dofs, const Vector* increment,
                            int field, FieldResidual* result, Diagnostic* d) {
  int count = 0;
  for (int row = 0; row < dofs->NodeCount(); ++row)
    if (dofs->Participates(dofs->NodeId(row)) &&
        dofs->NodeDofCount(dofs->NodeId(row)) >= 3 * field + 3)
      count += 3;
  if (count == 0) return StatusCode::kOk;
  std::unique_ptr<int[]> indices(new (std::nothrow) int[count]);
  if (!indices) return model_storage::AllocationFailure("Increment indices", d);
  int j = 0;
  for (int row = 0; row < dofs->NodeCount(); ++row)
    if (dofs->Participates(dofs->NodeId(row)) &&
        dofs->NodeDofCount(dofs->NodeId(row)) >= 3 * field + 3)
      for (int i = 0; i < 3; ++i)
        indices[j++] = dofs->GlobalDof(dofs->NodeId(row), 3 * field + i + 1);
  std::unique_ptr<Vector> values;
  StatusCode status = MakeVector(count, &values, d);
  if (status != StatusCode::kOk) return status;
  status = values->GatherFrom(increment, indices.get(), count, 0, d);
  if (status != StatusCode::kOk) return status;
  int maximum_index = 0;
  status = values->MaximumAbs(&result->increment_max, &maximum_index, d);
  result->correction_max = result->increment_max;  // Actual first solve.
  return status;
}
StatusCode Average(const Vector* flux, const FieldControls* controls,
                   const AnalysisState* history, int field_number,
                   FieldResidual* field, Diagnostic* d) {
  const double previous_max =
      history ? history->CommittedMaximumFlux(field_number) : 0;
  const double initial = controls->q0_explicit || !history
                             ? controls->q0
                             : history->CommittedRawMean(field_number);
  if (field->contribution_count > 0) {
    int index = 0;
    StatusCode status = flux->MaximumAbs(&field->maximum_flux, &index, d);
    if (status != StatusCode::kOk) return status;
    const double threshold = controls->epsilon_l * previous_max;
    if (!std::isfinite(threshold)) return NumericalFailure(d);
    const bool full =
        controls->epsilon_l == 0 || field->maximum_flux < previous_max / 10;
    double sum = 0;
    if (full) {
      status = flux->AbsoluteSum(&sum, d);
      if (status != StatusCode::kOk) return status;
      field->spatial_mean = sum / field->contribution_count;
    } else {
      std::unique_ptr<int[]> indices(
          new (std::nothrow) int[field->contribution_count]);
      if (!indices)
        return model_storage::AllocationFailure("Active flux indices", d);
      int count = 0;
      for (int i = 0; i < field->contribution_count; ++i)
        if (std::abs(flux->Data()[i]) >= threshold) indices[count++] = i;
      if (count == 0) {
        const StatusCode failure = ReportFailure(
            StatusCode::kLibraryFailure, DiagnosticStage::kComputation,
            "Empty selected active flux mean", d);
        if (d) {
          d->active_flux_failure = true;
          d->flux_field = field_number;
          d->flux_epsilon_l = controls->epsilon_l;
          d->flux_history_max = previous_max;
          d->flux_threshold = threshold;
          d->flux_current_max = field->maximum_flux;
          d->flux_count = field->contribution_count;
          d->active_flux_count = count;
        }
        return failure;
      }
      std::unique_ptr<Vector> active;
      status = MakeVector(count, &active, d);
      if (status != StatusCode::kOk) return status;
      status = active->GatherFrom(flux, indices.get(), count, 0, d);
      if (status != StatusCode::kOk) return status;
      status = active->AbsoluteSum(&sum, d);
      if (status != StatusCode::kOk) return status;
      field->spatial_mean = sum / count;
    }
  }
  if (controls->qu_explicit) {
    field->raw_mean = controls->qu;
  } else if (field->contribution_count == 0) {
    field->raw_mean = initial;
  } else {
    const double update_threshold = controls->epsilon * initial;
    if (!std::isfinite(update_threshold)) return NumericalFailure(d);
    field->raw_mean =
        field->spatial_mean > update_threshold ? field->spatial_mean : initial;
  }
  field->used_mean = field->raw_mean;
  if (!std::isfinite(field->raw_mean) || field->raw_mean <= 0)
    return model_storage::Invalid("Equilibrium raw flux scale", d);
  return StatusCode::kOk;
}
StatusCode Convert(const FieldControls* controls, double factor,
                   const FieldResidual* other, FieldResidual* field,
                   Diagnostic* d) {
  if (controls->qu_explicit) return StatusCode::kOk;
  const double converted = factor * controls->cf * other->raw_mean;
  if (!std::isfinite(converted))
    return ReportFailure(StatusCode::kLibraryFailure,
                         DiagnosticStage::kComputation,
                         "Equilibrium field conversion", d);
  if (field->spatial_mean <= field->raw_mean && field->raw_mean < converted) {
    field->used_mean = converted;
    field->converted = true;
  }
  return StatusCode::kOk;
}
}  // namespace

StatusCode ValidateEquilibriumCounts(const AnalysisModel* model,
                                     const EffectiveConditions* conditions,
                                     Diagnostic* d) {
  if (!model || !conditions)
    return model_storage::Invalid("Equilibrium count input", d);
  int count = 0;
  auto status = FluxCount(model, conditions, 0, &count, d);
  if (status != StatusCode::kOk) return status;
  return FluxCount(model, conditions, 1, &count, d);
}
namespace {
StatusCode EvaluateFieldCriteriaImpl(const FieldControls* controls,
                                     bool displacement, double length,
                                     bool allow_linear,
                                     const FieldConvergenceOptions* options,
                                     FieldResidual* field, Diagnostic* d) {
  if (!field) return model_storage::Invalid("Field criteria", d);
  ScalarFieldDecisionInput input;
  input.dof_count = field->dof_count;
  input.maximum = field->maximum;
  input.spatial_mean = field->spatial_mean;
  input.used_mean = field->used_mean;
  input.increment_max = field->increment_max;
  input.correction_max = field->correction_max;
  input.allow_linear = allow_linear;
  ScalarFieldDecision decision;
  decision.phi = field->phi;
  decision.residual_limit = field->residual_limit;
  decision.correction_limit = field->correction_limit;
  // Both enums deliberately use the same four scalar branches; explicit
  // mapping keeps physical field identities out of the common decision seam.
  switch (field->branch) {
    case ResidualBranch::kNoDofs:
      decision.branch = ScalarFieldBranch::kNoDofs;
      break;
    case ResidualBranch::kZeroFlux:
      decision.branch = ScalarFieldBranch::kZeroFlux;
      break;
    case ResidualBranch::kLinear:
      decision.branch = ScalarFieldBranch::kLinear;
      break;
    case ResidualBranch::kGeneral:
      decision.branch = ScalarFieldBranch::kGeneral;
      break;
  }
  decision.passed = field->passed;
  const auto status = EvaluateScalarFieldDecision(
      controls, options, displacement, length, &input, &decision, d);
  field->phi = decision.phi;
  field->residual_limit = decision.residual_limit;
  field->correction_limit = decision.correction_limit;
  field->passed = decision.passed;
  switch (decision.branch) {
    case ScalarFieldBranch::kNoDofs:
      field->branch = ResidualBranch::kNoDofs;
      break;
    case ScalarFieldBranch::kZeroFlux:
      field->branch = ResidualBranch::kZeroFlux;
      break;
    case ScalarFieldBranch::kLinear:
      field->branch = ResidualBranch::kLinear;
      break;
    case ScalarFieldBranch::kGeneral:
      field->branch = ResidualBranch::kGeneral;
      break;
  }
  return status;
}

}  // namespace

StatusCode EvaluateFieldCriteria(const FieldControls* controls,
                                 bool displacement, double length,
                                 FieldResidual* field, Diagnostic* d) {
  return EvaluateFieldCriteriaImpl(controls, displacement, length, true,
                                   nullptr, field, d);
}

StatusCode EvaluateJointFieldCriteriaWithOptions(
    const ControlSettings* controls, double length,
    const FieldConvergenceOptions* force_options,
    const FieldConvergenceOptions* moment_options, FieldResidual* force,
    FieldResidual* moment, Diagnostic* d) {
  if (!controls || !force || !moment || force == moment)
    return model_storage::Invalid("Joint field criteria", d);
  StatusCode force_status = EvaluateFieldCriteriaImpl(
      controls->displacement, true, length, true, force_options, force, d);
  const bool force_valid = force_status == StatusCode::kOk ||
                           force_status == StatusCode::kEquilibriumFailure;
  // Complete the other field even after rejection, preserving the first error.
  StatusCode moment_status = EvaluateFieldCriteriaImpl(
      controls->rotation, false, length, true, moment_options, moment,
      force_valid ? d : nullptr);
  if (!force_valid) return force_status;
  if (moment_status != StatusCode::kOk &&
      moment_status != StatusCode::kEquilibriumFailure)
    return moment_status;
  const bool all_linear =
      (force->dof_count == 0 || force->phi <= controls->displacement->rl) &&
      (moment->dof_count == 0 || moment->phi <= controls->rotation->rl);
  // The public Linear Increments shortcut requires every active field.
  // Reclassify only provisional linear branches; zero-flux keeps priority.
  if (!all_linear) {
    if (force->branch == ResidualBranch::kLinear)
      force_status = EvaluateFieldCriteriaImpl(
          controls->displacement, true, length, false, force_options, force, d);
    if (moment->branch == ResidualBranch::kLinear)
      moment_status = EvaluateFieldCriteriaImpl(
          controls->rotation, false, length, false, moment_options, moment, d);
  }
  if (force_status != StatusCode::kOk &&
      force_status != StatusCode::kEquilibriumFailure)
    return force_status;
  if (moment_status != StatusCode::kOk &&
      moment_status != StatusCode::kEquilibriumFailure)
    return moment_status;
  if (force_status != StatusCode::kOk || moment_status != StatusCode::kOk) {
    const FieldResidual* failed = !force->passed ? force : moment;
    if (d) d->node_id = failed->maximum_node;
    return ReportFailure(StatusCode::kEquilibriumFailure,
                         DiagnosticStage::kComputation,
                         !force->passed ? "Equilibrium displacement field"
                                        : "Equilibrium rotation field",
                         d);
  }
  return StatusCode::kOk;
}

StatusCode EvaluateJointFieldCriteria(const ControlSettings* controls,
                                      double length, FieldResidual* force,
                                      FieldResidual* moment, Diagnostic* d) {
  return EvaluateJointFieldCriteriaWithOptions(controls, length, nullptr,
                                               nullptr, force, moment, d);
}

class EquilibriumEvaluator {
 public:
  static StatusCode Run(const AnalysisModel* model, const DofManager* dofs,
                        const EffectiveConditions* conditions,
                        const ControlSettings* controls,
                        const Vector* candidate, const AnalysisState* history,
                        StaticWorkspace* workspace, ResidualReport* report,
                        Diagnostic* d, ElementWorkspace** workers,
                        int worker_count) {
    StatusCode status =
        Contributions(model, dofs, conditions, controls, candidate, history,
                      workspace, report, d, workers, worker_count);
    if (status != StatusCode::kOk) return status;
    status = Residuals(dofs, candidate, workspace, report, d);
    return status == StatusCode::kOk ? Decide(controls, report, d) : status;
  }

 private:
  static StatusCode Contributions(
      const AnalysisModel* model, const DofManager* dofs,
      const EffectiveConditions* conditions, const ControlSettings* controls,
      const Vector* candidate, const AnalysisState* history,
      StaticWorkspace* workspace, ResidualReport* report, Diagnostic* d,
      ElementWorkspace** workers, int worker_count) {
    // Reports own summaries; reset scalar records without numerical data copy.
    ResetField(report->displacement_);
    ResetField(report->rotation_);
    report->length_ = 0;
    StatusCode status = workspace->Increment()->CopyFrom(candidate, d);
    if (status != StatusCode::kOk) return status;
    if (history) {
      status = workspace->Increment()->AddScaled(
          -1, history->CommittedDisplacement(), d);
      if (status != StatusCode::kOk) return status;
    }
    status = numerical_support::ExternalForce(dofs, conditions,
                                              workspace->ExternalForce(), d);
    if (status != StatusCode::kOk) return status;
    status = workspace->InternalForce()->Scale(0, d);
    if (status != StatusCode::kOk) return status;
    status = report->reaction_->Scale(0, d);
    if (status != StatusCode::kOk) return status;
    int counts[2] = {};
    for (int field = 0; field < 2; ++field) {
      status = FluxCount(model, conditions, field, &counts[field], d);
      if (status != StatusCode::kOk) return status;
    }
    std::unique_ptr<Vector> force_flux, moment_flux, lengths;
    status = MakeVector(std::max(1, counts[0]), &force_flux, d);
    if (status != StatusCode::kOk) return status;
    status = MakeVector(std::max(1, counts[1]), &moment_flux, d);
    if (status != StatusCode::kOk) return status;
    const int element_count = static_cast<int>(model->Elements()->size());
    status = MakeVector(std::max(1, element_count), &lengths, d);
    if (status != StatusCode::kOk) return status;
    Buffers* raw = nullptr;
    const int capacity =
        std::max(1, std::min(element_count, Buffers::kElements));
    int maximum_dofs = 12;
    bool beam_only = true;
    for (const auto& entry : *model->Elements()) {
      maximum_dofs = std::max(maximum_dofs, entry.second->dof_count());
      beam_only = beam_only && entry.second->family() == ElementFamily::kBeam;
    }
    status = Buffers::Create(false, capacity, &raw, d, workers, worker_count,
                             maximum_dofs, beam_only);
    std::unique_ptr<Buffers> buffers(raw);
    if (status != StatusCode::kOk) return status;
    status = numerical_support::PrepareLineLoads(conditions->LineLoadState(),
                                                 raw, d);
    if (status != StatusCode::kOk) return status;
    ForceTaskContext task;
    task.buffers = raw;
    task.candidate = candidate;
    task.force_flux = force_flux.get();
    task.moment_flux = moment_flux.get();
    task.lengths = lengths.get();
    task.line_loads = conditions->LineLoadState();
    std::unique_ptr<int[]> force_offsets(new (std::nothrow) int[capacity]);
    std::unique_ptr<int[]> moment_offsets(new (std::nothrow) int[capacity]);
    if (!force_offsets || !moment_offsets)
      return model_storage::AllocationFailure("Flux prefix indices", d);
    task.force_offsets = force_offsets.get();
    task.moment_offsets = moment_offsets.get();
    int force_offset = 0, moment_offset = 0;
    auto next = model->Elements()->begin();
    while (next != model->Elements()->end()) {
      int count = 0;
      raw->vector_offsets[0] = 0;
      while (count < capacity && next != model->Elements()->end()) {
        raw->elements[count] = next->second;
        const auto* element = next->second;
        int* global = raw->global_indices + raw->element_dofs * count;
        status =
            numerical_support::ElementIndices(next->second, dofs, global, d);
        if (status != StatusCode::kOk) {
          if (d) d->element_id = next->first;
          return status;
        }
        const int offset = raw->vector_offsets[count];
        for (int i = 0; i < element->dof_count(); ++i)
          raw->vector_indices[offset + i] = global[i];
        raw->vector_offsets[count + 1] = offset + element->dof_count();
        force_offsets[count] = force_offset;
        moment_offsets[count] = moment_offset;
        force_offset += 3 * element->node_count();
        if (element->dofs_per_node() == 6)
          moment_offset += 3 * element->node_count();
        ++count;
        ++next;
      }
      status = ExecuteElementTasks(count, worker_count, ComputeForce, &task,
                                   raw->diagnostics, nullptr, d);
      if (status != StatusCode::kOk) return status;
      status = workspace->InternalForce()->ScatterAdd(
          raw->vector_values, raw->vector_indices, raw->vector_offsets[count],
          d);
      if (status != StatusCode::kOk) return status;
      task.first_element += count;
    }
    if (element_count > 0) {
      status = lengths->AbsoluteSum(&report->length_, d);
      if (status != StatusCode::kOk) return status;
      report->length_ /= element_count;
    }
    for (const auto* load : *conditions->Loads()) {
      if (!load->contributes_flux) continue;
      status =
          load->definition->Dof() <= 3
              ? force_flux->Set(force_offset++, load->definition->Value(), d)
              : moment_flux->Set(moment_offset++, load->definition->Value(), d);
      if (status != StatusCode::kOk) return status;
    }
    FieldResidual* force = report->displacement_;
    FieldResidual* moment = report->rotation_;
    force->contribution_count = counts[0];
    moment->contribution_count = counts[1];
    status =
        Average(force_flux.get(), controls->displacement, history, 0, force, d);
    if (status != StatusCode::kOk) return status;
    bool has_rotation = false;
    for (int row = 0; row < dofs->NodeCount(); ++row)
      has_rotation = has_rotation || dofs->NodeDofCount(dofs->NodeId(row)) == 6;
    if (has_rotation) {
      status =
          Average(moment_flux.get(), controls->rotation, history, 1, moment, d);
      if (status != StatusCode::kOk) return status;
    }
    if (element_count > 0 && has_rotation) {
      status = Convert(controls->displacement, 1 / report->length_, moment,
                       force, d);
      if (status != StatusCode::kOk) return status;
      status = Convert(controls->rotation, report->length_, force, moment, d);
      if (status != StatusCode::kOk) return status;
    }
    return StatusCode::kOk;
  }
  static StatusCode Residuals(const DofManager* dofs, const Vector* candidate,
                              StaticWorkspace* workspace,
                              ResidualReport* report, Diagnostic* d) {
    FieldResidual* force = report->displacement_;
    FieldResidual* moment = report->rotation_;
    StatusCode status = StatusCode::kOk;
    std::unique_ptr<Vector> residual;
    status = MakeVector(candidate->Size(), &residual, d);
    if (status != StatusCode::kOk) return status;
    status = residual->CopyFrom(workspace->ExternalForce(), d);
    if (status != StatusCode::kOk) return status;
    status = residual->AddScaled(-1, workspace->InternalForce(), d);
    if (status != StatusCode::kOk) return status;
    status = FieldNorms(dofs, residual.get(), 0, force, d);
    if (status != StatusCode::kOk) return status;
    status = FieldNorms(dofs, residual.get(), 1, moment, d);
    if (status != StatusCode::kOk) return status;
    status = IncrementMaximum(dofs, workspace->Increment(), 0, force, d);
    if (status != StatusCode::kOk) return status;
    status = IncrementMaximum(dofs, workspace->Increment(), 1, moment, d);
    if (status != StatusCode::kOk) return status;
    // Gather constrained residuals, negate, and selector-scatter reactions.
    int boundary_count = 0;
    for (int i = 0; i < candidate->Size(); ++i)
      if (dofs->Equation(i) == -1) ++boundary_count;
    if (boundary_count > 0) {
      std::unique_ptr<int[]> indices(new (std::nothrow) int[boundary_count]);
      if (!indices)
        return model_storage::AllocationFailure("Reaction indices", d);
      int j = 0;
      for (int i = 0; i < candidate->Size(); ++i)
        if (dofs->Equation(i) == -1) indices[j++] = i;
      std::unique_ptr<Vector> constrained;
      status = MakeVector(boundary_count, &constrained, d);
      if (status != StatusCode::kOk) return status;
      status = constrained->GatherFrom(residual.get(), indices.get(),
                                       boundary_count, 0, d);
      if (status != StatusCode::kOk) return status;
      status = constrained->Scale(-1, d);
      if (status != StatusCode::kOk) return status;
      status = report->reaction_->ScatterAdd(constrained.get(), indices.get(),
                                             boundary_count, d);
      if (status != StatusCode::kOk) return status;
    }
    return StatusCode::kOk;
  }
  static StatusCode Decide(const ControlSettings* controls,
                           ResidualReport* report, Diagnostic* d) {
    return EvaluateJointFieldCriteria(
        controls, report->length_, report->displacement_, report->rotation_, d);
  }
};
namespace {
StatusCode CheckEquilibriumImpl(
    const AnalysisModel* model, const DofManager* dofs,
    const EffectiveConditions* conditions, const ControlSettings* controls,
    const Vector* candidate_u, const AnalysisState* history,
    StaticWorkspace* workspace, ResidualReport* report, Diagnostic* d,
    ElementWorkspace** workers, int worker_count) {
  StatusCode status =
      numerical_support::Validate(model, dofs, conditions, workspace, d);
  if (status != StatusCode::kOk) return status;
  if (!controls || !ValidControls(controls->displacement) ||
      !ValidControls(controls->rotation) || !candidate_u || !report ||
      candidate_u->Size() != dofs->GlobalCount() ||
      report->Reaction()->Size() != candidate_u->Size() ||
      candidate_u == workspace->ExternalForce() ||
      candidate_u == workspace->InternalForce() ||
      candidate_u == report->Reaction() ||
      (history &&
       history->CommittedDisplacement()->Size() != candidate_u->Size()))
    return model_storage::Invalid("CheckEquilibrium", d);
  if (!history && conditions->StepNumber() != 1)
    return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                         "CheckEquilibrium requires committed history", d);
  for (const auto* boundary : *conditions->Boundaries()) {
    const int index =
        dofs->GlobalDof(boundary->node_id, boundary->definition->Dof());
    if (index < 0 || dofs->Equation(index) != -1 ||
        candidate_u->Data()[index] != boundary->definition->Value())
      return model_storage::Invalid("Candidate prescribed displacement", d);
  }
  return EquilibriumEvaluator::Run(model, dofs, conditions, controls,
                                   candidate_u, history, workspace, report, d,
                                   workers, worker_count);
}
}  // namespace
StatusCode CheckEquilibrium(const AnalysisModel* model, const DofManager* dofs,
                            const EffectiveConditions* conditions,
                            const ControlSettings* controls,
                            const Vector* candidate_u,
                            StaticWorkspace* workspace, ResidualReport* report,
                            Diagnostic* d) {
  return CheckEquilibriumImpl(model, dofs, conditions, controls, candidate_u,
                              nullptr, workspace, report, d, nullptr, 1);
}
StatusCode CheckEquilibrium(const AnalysisModel* model, const DofManager* dofs,
                            const EffectiveConditions* conditions,
                            const ControlSettings* controls,
                            const Vector* candidate_u,
                            const AnalysisState* history,
                            StaticWorkspace* workspace, ResidualReport* report,
                            Diagnostic* d, ElementWorkspace** workers,
                            int worker_count) {
  if (!history) return model_storage::Invalid("CheckEquilibrium history", d);
  return CheckEquilibriumImpl(model, dofs, conditions, controls, candidate_u,
                              history, workspace, report, d, workers,
                              worker_count);
}
}  // namespace fesa
