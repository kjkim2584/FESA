#include "analysis/field_measurement.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <new>
#include <vector>

#include "analysis/linear_dynamic_element.h"
#include "fesa/analysis/analysis_model.h"
#include "fesa/analysis/analysis_state.h"
#include "fesa/analysis/beam_internal_state.h"
#include "fesa/analysis/rotation_state.h"
#include "fesa/math/rotation.h"
#include "model/model_storage.h"

namespace fesa {
namespace {
bool FiniteVector(const Vector* values) {
  if (!values) return false;
  for (int i = 0; i < values->Size(); ++i)
    if (!std::isfinite(values->Data()[i])) return false;
  return true;
}
bool BeamMapping(const DofManager* dofs) {
  if (!dofs || dofs->InternalDofCount() <= 0 ||
      dofs->InternalDofCount() % 2 != 0)
    return false;
  for (int row = 0; row < dofs->InternalDofCount() / 2; ++row)
    for (int end = 0; end < 2; ++end)
      if (dofs->InternalDof(dofs->InternalElementId(row), end) != 2 * row + end)
        return false;
  return true;
}
}  // namespace
StatusCode MeasureBeamResidualField(const Vector* physical_residual,
                                    const Vector* internal_residual,
                                    const DofManager* dofs,
                                    FieldResidual* output, Diagnostic* d) {
  if (!BeamMapping(dofs) || !output || !FiniteVector(physical_residual) ||
      !FiniteVector(internal_residual) ||
      physical_residual->Size() != dofs->GlobalCount() ||
      internal_residual->Size() != dofs->InternalDofCount())
    return model_storage::Invalid("Beam force field input", d);
  try {
    std::vector<int> indices, nodes, components;
    for (int row = 0; row < dofs->NodeCount(); ++row) {
      const int id = dofs->NodeId(row);
      for (int component = 1; component <= 3; ++component) {
        const int global = dofs->GlobalDof(id, component);
        if (dofs->Equation(global) < 0) continue;
        indices.push_back(global);
        nodes.push_back(id);
        components.push_back(component);
      }
    }
    const int physical = static_cast<int>(indices.size());
    const int count = physical + dofs->InternalDofCount();
    Vector* raw = nullptr;
    auto status = Vector::Create(count, &raw, d);
    std::unique_ptr<Vector> values(raw);
    if (status != StatusCode::kOk) return status;
    if (physical > 0)
      status =
          values->GatherFrom(physical_residual, indices.data(), physical, 0, d);
    if (status != StatusCode::kOk) return status;
    indices.resize(dofs->InternalDofCount());
    for (int i = 0; i < dofs->InternalDofCount(); ++i) indices[i] = i;
    status = values->GatherFrom(internal_residual, indices.data(),
                                dofs->InternalDofCount(), physical, d);
    if (status != StatusCode::kOk) return status;
    int maximum_index = 0;
    output->dof_count = count;
    status = values->MaximumAbs(&output->maximum, &maximum_index, d);
    if (status == StatusCode::kOk) status = values->Norm2(&output->l2, d);
    if (status != StatusCode::kOk) return status;
    output->rms = output->l2 / std::sqrt(static_cast<double>(count));
    output->maximum_node = output->maximum_dof = output->maximum_element = 0;
    output->maximum_internal_end = -1;
    if (maximum_index < physical) {
      output->maximum_node = nodes[maximum_index];
      output->maximum_dof = components[maximum_index];
    } else {
      const int internal = maximum_index - physical;
      output->maximum_element = dofs->InternalElementId(internal / 2);
      output->maximum_internal_end = internal % 2;
    }
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return model_storage::AllocationFailure("Beam force field mapping", d);
  }
}
StatusCode MeasureBeamInternalMotion(const DofManager* dofs,
                                     const AnalysisState* state,
                                     const Vector* actual_correction,
                                     FieldResidual* force, Diagnostic* d) {
  const auto* internal = state ? state->BeamInternals() : nullptr;
  if (!BeamMapping(dofs) || !state || !state->HasValidCandidate() || !force ||
      state->GlobalCount() != dofs->GlobalCount() || !internal ||
      internal->ElementCount() != dofs->InternalDofCount() / 2 ||
      !FiniteVector(internal->CandidatePosition()) ||
      !FiniteVector(internal->CommittedPosition()) ||
      internal->CandidatePosition()->Size() != dofs->InternalDofCount() ||
      (actual_correction &&
       (actual_correction->Size() != dofs->InternalDofCount() ||
        !FiniteVector(actual_correction))))
    return model_storage::Invalid("Beam motion field input", d);
  for (int row = 0; row < internal->ElementCount(); ++row)
    if (internal->ElementId(row) != dofs->InternalElementId(row) ||
        internal->Offset(internal->ElementId(row)) != 2 * row)
      return model_storage::Invalid("Beam motion field mapping", d);
  Vector* raw = nullptr;
  auto status = Vector::Create(dofs->InternalDofCount(), &raw, d);
  std::unique_ptr<Vector> increment(raw);
  if (status == StatusCode::kOk)
    status = increment->CopyFrom(internal->CandidatePosition(), d);
  if (status == StatusCode::kOk)
    status = increment->AddScaled(-1, internal->CommittedPosition(), d);
  double maximum = 0;
  int index = 0;
  if (status == StatusCode::kOk)
    status = increment->MaximumAbs(&maximum, &index, d);
  if (status != StatusCode::kOk) return status;
  force->increment_max = std::max(force->increment_max, maximum);
  if (!actual_correction) return StatusCode::kOk;
  status = actual_correction->MaximumAbs(&maximum, &index, d);
  if (status == StatusCode::kOk)
    force->correction_max = std::max(force->correction_max, maximum);
  return status;
}
StatusCode CreateElementFieldFlux(const AnalysisModel* model,
                                  const Vector* element_values, int block_count,
                                  int field, const Vector* external_origins,
                                  Vector** out, Diagnostic* d,
                                  bool material_damping, bool beam_internal) {
  if (out) *out = nullptr;
  if (!out || !model || block_count <= 0 || field < 0 || field > 1)
    return model_storage::Invalid("Element field flux inputs", d);
  const int limit = std::numeric_limits<int>::max();
  int required = 0;
  for (const auto& entry : *model->Elements()) {
    const auto* element = entry.second;
    const int width = element->dofs_per_node();
    const int count = element->dof_count() +
                      (beam_internal ? element->internal_dof_count() : 0);
    const bool extra = material_damping && HasMaterialRayleighDamping(element);
    if (extra && block_count == limit)
      return model_storage::Invalid("Element field flux block count", d);
    const int blocks = block_count + (extra ? 1 : 0);
    if ((beam_internal && (element->family() != ElementFamily::kBeam ||
                           element->dof_count() != 12 ||
                           element->internal_dof_count() != 2 || width != 6)) ||
        (width != 3 && width != 6) || count <= 0 ||
        (!beam_internal && count % width != 0) ||
        count > (limit - required) / blocks)
      return model_storage::Invalid("Element field flux layout/size", d);
    required += count * blocks;
  }
  if ((required > 0 && !element_values) ||
      (element_values && element_values->Size() != required))
    return model_storage::Invalid("Element field flux shape", d);
  try {
    std::vector<int> indices;
    int base = 0;
    for (const auto& entry : *model->Elements()) {
      const auto* element = entry.second;
      const int blocks = block_count +
          (material_damping && HasMaterialRayleighDamping(element) ? 1 : 0);
      const int physical = element->dof_count();
      const int count = physical +
          (beam_internal ? element->internal_dof_count() : 0);
      for (int block = 0; block < blocks; ++block) {
        for (int local = 0; local < count; ++local)
          if ((local < physical &&
               local % element->dofs_per_node() / 3 == field) ||
              (local >= physical && field == 0))
            indices.push_back(base + local);
        base += count;
      }
    }
    const int selected = static_cast<int>(indices.size());
    const int external = external_origins ? external_origins->Size() : 0;
    if (external > limit - selected)
      return model_storage::Invalid("Element field flux count", d);
    if (selected + external == 0) return StatusCode::kOk;
    Vector* raw = nullptr;
    auto status = Vector::Create(selected + external, &raw, d);
    std::unique_ptr<Vector> result(raw);
    if (status != StatusCode::kOk) return status;
    if (selected > 0) {
      status = raw->GatherFrom(element_values, indices.data(), selected, 0, d);
      if (status != StatusCode::kOk) return status;
    }
    if (external > 0) {
      indices.resize(external);
      for (int i = 0; i < external; ++i) indices[i] = i;
      status = raw->GatherFrom(external_origins, indices.data(), external,
                               selected, d);
      if (status != StatusCode::kOk) return status;
    }
    *out = result.release();
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return model_storage::AllocationFailure("Element field flux indices", d);
  }
}
StatusCode SplitLinearDynamicFluxOrigins(const AnalysisModel* model,
                                         const EffectiveConditions* conditions,
                                         const Vector* origins, Vector** force,
                                         Vector** moment, Diagnostic* d) {
  if (force) *force = nullptr;
  if (moment) *moment = nullptr;
  if (!force || !moment || force == moment || !model || !conditions ||
      !origins || model->Conditions() != conditions ||
      model->Elements()->empty())
    return model_storage::Invalid("Dynamic origin split inputs", d);
  const auto external = conditions->Loads()->size();
  const int limit = std::numeric_limits<int>::max();
  if (external > static_cast<size_t>(limit))
    return model_storage::Invalid("Dynamic origin split size", d);
  int element_entries = 0;
  int counts[2] = {};
  for (const auto& entry : *model->Elements()) {
    const auto* e = entry.second;
    const bool shell = e->family() == ElementFamily::kShell;
    const int size = e->dof_count(), components = e->dofs_per_node();
    if (size <= 0 || components != (shell ? 6 : 3) ||
        size / components != e->node_count() || size % components != 0 ||
        (!shell && e->family() != ElementFamily::kSolid))
      return model_storage::Invalid("Dynamic origin split element layout", d);
    const int blocks = HasMaterialRayleighDamping(e) ? 3 : 2;
    if (size > (limit - static_cast<int>(external) - element_entries) / blocks)
      return model_storage::Invalid("Dynamic origin split size", d);
    element_entries += size * blocks;
    counts[0] += (shell ? size / 2 : size) * blocks;
    counts[1] += shell ? size / 2 * blocks : 0;
  }
  if (origins->Size() != element_entries + static_cast<int>(external))
    return model_storage::Invalid("Dynamic origin split shape", d);
  for (const auto* load : *conditions->Loads()) {
    const int component = load->definition->Dof();
    if (component < 1 || component > 6)
      return model_storage::Invalid("Dynamic origin split load component", d);
    ++counts[component > 3 ? 1 : 0];
  }
  std::unique_ptr<int[]> indices[2];
  for (int field = 0; field < 2; ++field) {
    if (counts[field] == 0) continue;
    indices[field].reset(new (std::nothrow) int[counts[field]]);
    if (!indices[field])
      return model_storage::AllocationFailure("Dynamic origin split indices",
                                              d);
  }
  int offsets[2] = {}, base = 0;
  for (const auto& entry : *model->Elements()) {
    const int blocks = HasMaterialRayleighDamping(entry.second) ? 3 : 2;
    const int size = entry.second->dof_count();
    for (int kind = 0; kind < blocks; ++kind)
      for (int component = 0; component < size; ++component) {
        const int field =
            component % entry.second->dofs_per_node() >= 3 ? 1 : 0;
        indices[field][offsets[field]++] = base + size * kind + component;
      }
    base += size * blocks;
  }
  for (const auto* load : *conditions->Loads()) {
    const int field = load->definition->Dof() > 3 ? 1 : 0;
    indices[field][offsets[field]++] = base++;
  }
  std::unique_ptr<Vector> fields[2];
  for (int field = 0; field < 2; ++field) {
    if (counts[field] == 0) continue;
    Vector* raw = nullptr;
    auto status = Vector::Create(counts[field], &raw, d);
    fields[field].reset(raw);
    if (status == StatusCode::kOk)
      status =
          raw->GatherFrom(origins, indices[field].get(), counts[field], 0, d);
    if (status != StatusCode::kOk) return status;
  }
  *force = fields[0].release();
  *moment = fields[1].release();
  return StatusCode::kOk;
}
StatusCode ComputeSpatialFluxMean(const Vector* values,
                                  const FieldControls* controls,
                                  double previous, int field,
                                  FieldResidual* output, Diagnostic* d) {
  if (!output) return model_storage::Invalid("Spatial flux output", d);
  if (!values) return StatusCode::kOk;
  if (!controls || !std::isfinite(previous) || previous < 0 ||
      !std::isfinite(controls->epsilon_l) || controls->epsilon_l < 0 ||
      field < 0 || field > 1)
    return model_storage::Invalid("Spatial flux inputs", d);
  try {
    output->contribution_count = values->Size();
    int maximum_index = 0;
    auto status = values->MaximumAbs(&output->maximum_flux, &maximum_index, d);
    if (status != StatusCode::kOk) return status;
    const double threshold = controls->epsilon_l * previous;
    const bool all = output->maximum_flux < previous / 10 || threshold == 0;
    if (!all && !std::isfinite(threshold))
      return ReportFailure(StatusCode::kLibraryFailure,
                           DiagnosticStage::kComputation,
                           "Spatial active threshold overflow", d);
    std::vector<int> selected;
    for (int i = 0; i < values->Size(); ++i)
      if (all || std::abs(values->Data()[i]) >= threshold)
        selected.push_back(i);
    if (selected.empty()) {
      status = ReportFailure(StatusCode::kLibraryFailure,
                             DiagnosticStage::kComputation,
                             "Spatial active flux selection", d);
      if (d) {
        d->active_flux_failure = true;
        d->flux_field = field;
        d->flux_count = values->Size();
        d->active_flux_count = 0;
        d->flux_epsilon_l = controls->epsilon_l;
        d->flux_history_max = previous;
        d->flux_threshold = threshold;
        d->flux_current_max = output->maximum_flux;
      }
      return status;
    }
    Vector* raw = nullptr;
    status = Vector::Create(static_cast<int>(selected.size()), &raw, d);
    std::unique_ptr<Vector> active(raw);
    if (status != StatusCode::kOk) return status;
    status = active->GatherFrom(values, selected.data(), active->Size(), 0, d);
    if (status != StatusCode::kOk) return status;
    status = active->Scale(1.0 / active->Size(), d);
    if (status != StatusCode::kOk) return status;
    return active->AbsoluteSum(&output->spatial_mean, d);
  } catch (const std::bad_alloc&) {
    return model_storage::AllocationFailure("Spatial flux selection", d);
  }
}

StatusCode MeasureResidualField(const Vector* residual, const DofManager* dofs,
                                int field, FieldResidual* output,
                                Diagnostic* d) {
  if (!residual || !dofs || !output ||
      residual->Size() != dofs->GlobalCount() || field < 0 || field > 1)
    return model_storage::Invalid("Residual field input", d);
  output->maximum = output->l2 = output->rms = 0;
  output->maximum_node = output->maximum_dof = 0;
  output->maximum_element = 0;
  output->maximum_internal_end = -1;
  try {
    std::vector<int> indices;
    std::vector<int> node_ids, components;
    for (int row = 0; row < dofs->NodeCount(); ++row) {
      const int id = dofs->NodeId(row);
      for (int component = 1; component <= dofs->NodeDofCount(id);
           ++component) {
        const int global = dofs->GlobalDof(id, component);
        if ((component - 1) / 3 == field && dofs->Equation(global) >= 0) {
          indices.push_back(global);
          node_ids.push_back(id);
          components.push_back(component);
        }
      }
    }
    output->dof_count = static_cast<int>(indices.size());
    if (indices.empty()) return StatusCode::kOk;
    Vector* raw = nullptr;
    auto status = Vector::Create(output->dof_count, &raw, d);
    std::unique_ptr<Vector> values(raw);
    if (status != StatusCode::kOk) return status;
    status = values->GatherFrom(residual, indices.data(), values->Size(), 0, d);
    if (status != StatusCode::kOk) return status;
    int maximum_index = 0;
    status = values->MaximumAbs(&output->maximum, &maximum_index, d);
    if (status != StatusCode::kOk) return status;
    output->maximum_node = node_ids[maximum_index];
    output->maximum_dof = components[maximum_index];
    status = values->Norm2(&output->l2, d);
    if (status != StatusCode::kOk) return status;
    output->rms = output->l2 / std::sqrt(static_cast<double>(values->Size()));
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return model_storage::AllocationFailure("Residual field mapping", d);
  }
}

StatusCode MeasureLinearDynamicField(const DofManager* dofs,
                                     const Vector* residual,
                                     const Vector* committed,
                                     const Vector* candidate,
                                     const Vector* predictor,
                                     FieldResidual* output, Diagnostic* d) {
  if (dofs)
    for (int row = 0; row < dofs->NodeCount(); ++row)
      if (dofs->NodeDofCount(dofs->NodeId(row)) != 3)
        return model_storage::Invalid("Dynamic field requires 3N", d);
  return MeasureLinearDynamicField(dofs, 0, residual, committed, candidate,
                                   predictor, output, d);
}
StatusCode MeasureLinearDynamicField(const DofManager* dofs, int field,
                                     const Vector* residual,
                                     const Vector* committed,
                                     const Vector* candidate,
                                     const Vector* predictor,
                                     FieldResidual* output, Diagnostic* d) {
  if (!dofs || !output || field < 0 || field > 1)
    return model_storage::Invalid("Dynamic field context", d);
  for (const auto* input : {residual, committed, candidate, predictor}) {
    if (!input || input->Size() != dofs->GlobalCount())
      return model_storage::Invalid("Dynamic field vector size", d);
  }
  for (int row = 0; row < dofs->NodeCount(); ++row) {
    const int count = dofs->NodeDofCount(dofs->NodeId(row));
    if (count != 3 && count != 6)
      return model_storage::Invalid("Dynamic field node mapping", d);
  }
  auto status = MeasureResidualField(residual, dofs, field, output, d);
  if (status != StatusCode::kOk) return status;
  output->increment_max = output->correction_max = 0;
  try {
    std::vector<int> indices(0);
    for (int row = 0; row < dofs->NodeCount(); ++row) {
      const int id = dofs->NodeId(row);
      if (!dofs->Participates(id)) continue;
      for (int component = 1; component <= dofs->NodeDofCount(id); ++component)
        if ((component - 1) / 3 == field)
          indices.push_back(dofs->GlobalDof(id, component));
    }
    if (indices.empty()) return StatusCode::kOk;
    Vector *raw_current = nullptr, *raw_base = nullptr;
    status = Vector::Create(static_cast<int>(indices.size()), &raw_current, d);
    std::unique_ptr<Vector> current(raw_current);
    if (status != StatusCode::kOk) return status;
    status = Vector::Create(current->Size(), &raw_base, d);
    std::unique_ptr<Vector> base(raw_base);
    if (status != StatusCode::kOk) return status;
    for (bool correction : {false, true}) {
      status =
          current->GatherFrom(candidate, indices.data(), current->Size(), 0, d);
      if (status != StatusCode::kOk) return status;
      status = base->GatherFrom(correction ? predictor : committed,
                                indices.data(), base->Size(), 0, d);
      if (status != StatusCode::kOk) return status;
      status = current->AddScaled(-1, base.get(), d);
      if (status != StatusCode::kOk) return status;
      int index = 0;
      status = current->MaximumAbs(
          correction ? &output->correction_max : &output->increment_max, &index,
          d);
      if (status != StatusCode::kOk) return status;
    }
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return model_storage::AllocationFailure("Dynamic field mapping", d);
  }
}
StatusCode MeasureNonlinearMotionFields(const DofManager* dofs,
                                        const AnalysisState* state,
                                        const Vector* actual_correction,
                                        FieldResidual* force,
                                        FieldResidual* moment,
                                        bool* rotation_branch, Diagnostic* d) {
  if (rotation_branch) *rotation_branch = false;
  if (!dofs || !state || !force || !moment || force == moment ||
      !rotation_branch || !state->FieldHistory() ||
      !state->HasValidCandidate() || dofs->NodeCount() <= 0 ||
      dofs->GlobalCount() != state->GlobalCount() ||
      (actual_correction && actual_correction->Size() != state->GlobalCount()))
    return model_storage::Invalid("Nonlinear motion field inputs", d);
  if (actual_correction)
    for (int i = 0; i < actual_correction->Size(); ++i)
      if (!std::isfinite(actual_correction->Data()[i]))
        return model_storage::Invalid("Nonlinear motion correction values", d);
  const auto* rotations = state->Rotations();
  if (rotations && (rotations->NodeCount() != dofs->NodeCount() ||
                    !rotations->CandidateOrientations() ||
                    !rotations->CommittedOrientations()))
    return model_storage::Invalid("Nonlinear motion orientations", d);
  if (!rotations) {
    const auto status = state->ValidateMapping(dofs, d);
    if (status != StatusCode::kOk) return status;
  }
  for (int row = 0; row < dofs->NodeCount(); ++row) {
    const int components = dofs->NodeDofCount(dofs->NodeId(row));
    if ((components != 3 && components != 6) || (components == 6 && !rotations))
      return model_storage::Invalid("Nonlinear motion physical layout", d);
  }
  try {
    std::unique_ptr<Vector> scratch[6];
    const int sizes[] = {
        state->GlobalCount(), 3, 3, rotations ? 4 : 0, rotations ? 4 : 0,
        rotations ? 3 : 0};
    for (int i = 0; i < 6; ++i) {
      if (sizes[i] == 0) continue;
      Vector* raw = nullptr;
      const auto status = Vector::Create(sizes[i], &raw, d);
      scratch[i].reset(raw);
      if (status != StatusCode::kOk) return status;
    }
    auto* increments = scratch[0].get();
    auto* current_u = scratch[1].get();
    auto* committed_u = scratch[2].get();
    auto* current_q = scratch[3].get();
    auto* committed_q = scratch[4].get();
    auto* rotation = scratch[5].get();
    StatusCode status = StatusCode::kOk;
    std::vector<int> force_indices, moment_indices;
    for (int row = 0; row < dofs->NodeCount(); ++row) {
      if (!dofs->Participates(dofs->NodeId(row))) continue;
      const int id = dofs->NodeId(row);
      const int offset = dofs->GlobalDof(id, 1);
      const int u_indices[3] = {offset, offset + 1, offset + 2};
      const int r_indices[3] = {offset + 3, offset + 4, offset + 5};
      const int q_indices[4] = {4 * row, 4 * row + 1, 4 * row + 2, 4 * row + 3};
      status = state->GatherCandidateTranslation(row, current_u, d);
      if (status != StatusCode::kOk) return status;
      status = committed_u->GatherFrom(state->CommittedDisplacement(),
                                       u_indices, 3, 0, d);
      if (status != StatusCode::kOk) return status;
      status = current_u->AddScaled(-1, committed_u, d);
      if (status != StatusCode::kOk) return status;
      if (dofs->NodeDofCount(id) == 6) {
        status = state->Rotations()->CandidateOrientations()->GatherTo(
            q_indices, 4, current_q, 0, d);
        if (status != StatusCode::kOk) return status;
        status = state->Rotations()->CommittedOrientations()->GatherTo(
            q_indices, 4, committed_q, 0, d);
        if (status != StatusCode::kOk) return status;
        bool branch = false;
        status = SpatialRotationIncrement(committed_q, current_q, rotation,
                                          &branch, d);
        if (status != StatusCode::kOk) {
          if (status == StatusCode::kEquilibriumFailure && branch)
            *rotation_branch = true;
          return status;
        }
        status = increments->ScatterAdd(rotation, r_indices, 3, d);
        if (status != StatusCode::kOk) return status;
        moment_indices.insert(moment_indices.end(), r_indices, r_indices + 3);
      }
      status = increments->ScatterAdd(current_u, u_indices, 3, d);
      if (status != StatusCode::kOk) return status;
      force_indices.insert(force_indices.end(), u_indices, u_indices + 3);
    }
    const std::vector<int>* fields[2] = {&force_indices, &moment_indices};
    FieldResidual* outputs[2] = {force, moment};
    for (int field = 0; field < 2; ++field) {
      if (fields[field]->empty()) continue;
      Vector* raw_values = nullptr;
      status = Vector::Create(static_cast<int>(fields[field]->size()),
                              &raw_values, d);
      std::unique_ptr<Vector> values(raw_values);
      if (status != StatusCode::kOk) return status;
      status = values->GatherFrom(increments, fields[field]->data(),
                                  values->Size(), 0, d);
      if (status != StatusCode::kOk) return status;
      int maximum_index = 0;
      status =
          values->MaximumAbs(&outputs[field]->increment_max, &maximum_index, d);
      if (status != StatusCode::kOk) return status;
      if (!actual_correction) continue;
      status = values->GatherFrom(actual_correction, fields[field]->data(),
                                  values->Size(), 0, d);
      if (status != StatusCode::kOk) return status;
      status = values->MaximumAbs(&outputs[field]->correction_max,
                                  &maximum_index, d);
      if (status != StatusCode::kOk) return status;
    }
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return model_storage::AllocationFailure("Nonlinear motion fields", d);
  }
}
}  // namespace fesa
