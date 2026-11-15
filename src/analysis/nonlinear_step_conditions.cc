#include "analysis/nonlinear_step_conditions.h"

#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <new>
#include <vector>

#include "fesa/analysis/rotation_prescription.h"
#include "fesa/analysis/rotation_state.h"
#include "fesa/math/rotation.h"
#include "model/model_storage.h"

namespace fesa {
namespace {
using BoundaryRows = std::map<int, const EffectiveBoundary*>;
StatusCode Validate(const NonlinearStepConditionInput* in, Diagnostic* d) {
  if (!in || !in->conditions || !in->conditions->ReleasedBoundaries() ||
      !in->dofs || !in->start_displacement || !in->start_reaction ||
      in->dofs->NodeCount() <= 0 ||
      in->dofs->NodeCount() > std::numeric_limits<int>::max() / 6 ||
      in->start_displacement->Size() != in->dofs->GlobalCount() ||
      in->start_reaction->Size() != in->dofs->GlobalCount() ||
      !std::isfinite(in->lambda) ||
      !std::isfinite(in->committed_lambda) || in->committed_lambda < 0 ||
      in->lambda < in->committed_lambda || in->lambda > 1)
    return model_storage::Invalid("Nonlinear step condition snapshots", d);
  bool rotations = false;
  for (int row = 0; row < in->dofs->NodeCount(); ++row) {
    const int count = in->dofs->NodeDofCount(in->dofs->NodeId(row));
    if (count != 3 && count != 6)
      return model_storage::Invalid("Nonlinear condition physical layout", d);
    rotations = rotations || count == 6;
  }
  if ((rotations && (!in->start_orientations ||
                    in->start_orientations->Rows() != in->dofs->NodeCount() ||
                    in->start_orientations->Columns() != 4)) ||
      (!rotations && in->start_orientations))
    return model_storage::Invalid("Nonlinear condition orientation snapshot", d);
  return StatusCode::kOk;
}
void Context(const EffectiveBoundary* row, int step, Diagnostic* d) {
  if (!d) return;
  d->source_file = row->definition->SourceFile();
  d->path_kind = DiagnosticPathKind::kLiteralFile;
  d->source_line = row->definition->SourceLine();
  d->step = step;
  d->node_id = row->node_id;
  d->element_id = 0;
  d->condition_dof = row->definition->Dof();
  d->duplicate_boundary = false;
  d->first_source_line = 0;
}
StatusCode PackOrigins(const Vector* values, const std::vector<int>* components,
                       int field, Vector** out, Diagnostic* d) {
  std::vector<int> selected;
  for (int i = 0; i < values->Size(); ++i)
    if (((components->at(i) > 3) ? 1 : 0) == field && values->Data()[i] != 0)
      selected.push_back(i);
  if (selected.empty()) return StatusCode::kOk;
  Vector* raw = nullptr;
  auto status = Vector::Create(static_cast<int>(selected.size()), &raw, d);
  std::unique_ptr<Vector> result(raw);
  if (status != StatusCode::kOk) return status;
  status = result->GatherFrom(values, selected.data(),
                              static_cast<int>(selected.size()), 0, d);
  if (status == StatusCode::kOk) *out = result.release();
  return status;
}
StatusCode Interpolate(const Vector* start, const Vector* end, double lambda,
                       Vector* out, Diagnostic* d) {
  auto status = out->CopyFrom(start, d);
  if (status == StatusCode::kOk) status = out->Scale(1 - lambda, d);
  if (status == StatusCode::kOk) status = out->AddScaled(lambda, end, d);
  return status;
}
StatusCode TranslateNode(const NonlinearStepConditionInput* in, int node,
                         const BoundaryRows* rows, AnalysisState* state,
                         Diagnostic* d) {
  Vector *raw_start = nullptr, *raw_end = nullptr, *raw_value = nullptr,
         *raw_current = nullptr;
  auto status = Vector::Create(3, &raw_start, d);
  std::unique_ptr<Vector> start(raw_start);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(3, &raw_end, d);
  std::unique_ptr<Vector> end(raw_end);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(3, &raw_value, d);
  std::unique_ptr<Vector> value(raw_value);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(3, &raw_current, d);
  std::unique_ptr<Vector> current(raw_current);
  if (status != StatusCode::kOk) return status;
  const int offset = in->dofs->GlobalDof(in->dofs->NodeId(node), 1);
  const int indices[3] = {offset, offset + 1, offset + 2};
  status = start->GatherFrom(in->start_displacement, indices, 3, 0, d);
  if (status == StatusCode::kOk) status = end->CopyFrom(start.get(), d);
  if (status != StatusCode::kOk) return status;
  for (int component = 0; component < 3; ++component) {
    const auto found = rows->find(offset + component);
    if (found == rows->end()) continue;
    Context(found->second, in->conditions->StepNumber(), d);
    status = end->Set(component, found->second->definition->Value(), d);
    if (status != StatusCode::kOk) return status;
  }
  status = Interpolate(start.get(), end.get(), in->lambda, value.get(), d);
  if (status != StatusCode::kOk) return status;
  status = state->GatherCandidateTranslation(node, current.get(), d);
  if (status != StatusCode::kOk) return status;
  for (int component = 0; component < 3; ++component) {
    if (rows->find(offset + component) == rows->end()) continue;
    status = current->GatherFrom(value.get(), &component, 1, component, d);
    if (status != StatusCode::kOk) return status;
  }
  return state->SetCandidateTranslation(node, current.get(), d);
}
StatusCode RotateNode(const NonlinearStepConditionInput* in, int node,
                      const BoundaryRows* rows, AnalysisState* state,
                      NonlinearFailure* failure, Diagnostic* d) {
  Vector *raw_q = nullptr, *raw_start = nullptr, *raw_end = nullptr,
         *raw_target = nullptr;
  auto status = Vector::Create(4, &raw_q, d);
  std::unique_ptr<Vector> q(raw_q);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(3, &raw_start, d);
  std::unique_ptr<Vector> start(raw_start);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(3, &raw_end, d);
  std::unique_ptr<Vector> end(raw_end);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(4, &raw_target, d);
  std::unique_ptr<Vector> target(raw_target);
  if (status != StatusCode::kOk) return status;
  const int indices[4] = {4 * node, 4 * node + 1, 4 * node + 2, 4 * node + 3};
  status = in->start_orientations->GatherTo(indices, 4, q.get(), 0, d);
  if (status != StatusCode::kOk) return status;
  const int offset = in->dofs->GlobalDof(in->dofs->NodeId(node), 1);
  int prior_count = 0;
  for (int component = 0; component < 3; ++component) {
    const auto* row = rows->at(offset + 3 + component);
    Context(row, in->conditions->StepNumber(), d);
    status = end->Set(component, row->definition->Value(), d);
    if (status != StatusCode::kOk) return status;
    if (row->physical_previous_definition) {
      ++prior_count;
      status =
          start->Set(component, row->physical_previous_definition->Value(), d);
      if (status != StatusCode::kOk) return status;
    }
  }
  if (prior_count != 0 && prior_count != 3)
    return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                         "Partial prior nonlinear rotation prescription", d);
  if (prior_count == 0) {
    status = PrincipalRotationVector(q.get(), start.get(), d);
    if (status != StatusCode::kOk) return status;
  }
  bool branch = false;
  status = EvaluateRotationPrescription(q.get(), start.get(), end.get(),
                                        in->lambda, in->committed_lambda,
                                        target.get(), &branch, d);
  if (status == StatusCode::kOk)
    status = state->SetCandidateOrientation(node, target.get(), &branch, d);
  if (status == StatusCode::kEquilibriumFailure && branch)
    *failure = NonlinearFailure::kRotationBranch;
  return status;
}
StatusCode ApplyPrescriptions(const NonlinearStepConditionInput* in,
                              AnalysisState* state, NonlinearFailure* failure,
                              Diagnostic* d) {
  auto status = Validate(in, d);
  if (status != StatusCode::kOk) return status;
  if (!state || !failure || state->GlobalCount() != in->dofs->GlobalCount() ||
      !state->HasValidCandidate() ||
      (in->start_orientations && (!state->Rotations() ||
                                 !state->Rotations()->CandidateOrientations())))
    return model_storage::Invalid("Prepared nonlinear prescription state", d);
  if (!state->Rotations()) {
    status = state->ValidateMapping(in->dofs, d);
    if (status != StatusCode::kOk) return status;
  }
  BoundaryRows rows;
  for (const auto* row : *in->conditions->Boundaries()) {
    const int index = in->dofs->GlobalDof(row->node_id, row->definition->Dof());
    if (index < 0) return model_storage::Invalid("Prescription mapping", d);
    rows.emplace(index, row);
  }
  for (int node = 0; node < in->dofs->NodeCount(); ++node) {
    bool translation = false;
    int rotations = 0;
    const int id = in->dofs->NodeId(node);
    const int offset = in->dofs->GlobalDof(id, 1);
    for (int c = 0; c < in->dofs->NodeDofCount(id); ++c) {
      if (rows.find(offset + c) == rows.end()) continue;
      if (c < 3)
        translation = true;
      else
        ++rotations;
    }
    if (rotations != 0 && rotations != 3)
      return ReportFailure(StatusCode::kUnsupported,
                           DiagnosticStage::kValidation,
                           "Partial nonlinear rotation prescription", d);
    if (translation) {
      status = TranslateNode(in, node, &rows, state, d);
      if (status != StatusCode::kOk) return status;
    }
    if (rotations) {
      status = RotateNode(in, node, &rows, state, failure, d);
      if (status != StatusCode::kOk) return status;
    }
  }
  return StatusCode::kOk;
}
}  // namespace

StatusCode EvaluateNonlinearStepLoads(const NonlinearStepConditionInput* in,
                                      Vector** external, Vector** force_origins,
                                      Vector** moment_origins, Diagnostic* d) {
  if (external) *external = nullptr;
  if (force_origins) *force_origins = nullptr;
  if (moment_origins) *moment_origins = nullptr;
  if (!external || !force_origins || !moment_origins ||
      external == force_origins || external == moment_origins ||
      force_origins == moment_origins)
    return model_storage::Invalid("Nonlinear step load outputs", d);
  auto status = Validate(in, d);
  if (status != StatusCode::kOk) return status;
  try {
    const auto* loads = in->conditions->Loads();
    const auto* released = in->conditions->ReleasedBoundaries();
    const auto limit =
        static_cast<std::size_t>(std::numeric_limits<int>::max());
    if (loads->size() > limit || released->size() > limit - loads->size())
      return model_storage::Invalid("Nonlinear step load count overflow", d);
    const int count = static_cast<int>(loads->size() + released->size());
    Vector* raw = nullptr;
    status = Vector::Create(in->dofs->GlobalCount(), &raw, d);
    std::unique_ptr<Vector> global(raw);
    if (status != StatusCode::kOk) return status;
    if (count == 0) {
      *external = global.release();
      return StatusCode::kOk;
    }
    status = Vector::Create(count, &raw, d);
    std::unique_ptr<Vector> start(raw);
    if (status != StatusCode::kOk) return status;
    status = Vector::Create(count, &raw, d);
    std::unique_ptr<Vector> end(raw);
    if (status != StatusCode::kOk) return status;
    status = Vector::Create(count, &raw, d);
    std::unique_ptr<Vector> change(raw);
    if (status != StatusCode::kOk) return status;
    status = Vector::Create(count, &raw, d);
    std::unique_ptr<Vector> value(raw);
    if (status != StatusCode::kOk) return status;
    std::vector<int> globals, components;
    globals.reserve(count);
    int i = 0;
    for (const auto* load : *loads) {
      const int index =
          in->dofs->GlobalDof(load->node_id, load->definition->Dof());
      if (index < 0) return model_storage::Invalid("Step load mapping", d);
      globals.push_back(index);
      components.push_back(load->definition->Dof());
      status = end->Set(i, load->value, d);
      if (status == StatusCode::kOk) status = change->Set(i, load->change, d);
      if (status != StatusCode::kOk) return status;
      ++i;
    }
    status = start->CopyFrom(end.get(), d);
    if (status == StatusCode::kOk)
      status = start->AddScaled(-1, change.get(), d);
    if (status != StatusCode::kOk) return status;
    for (const auto* row : *released) {
      const int index =
          in->dofs->GlobalDof(row->node_id, row->definition->Dof());
      if (index < 0)
        return model_storage::Invalid("Released reaction mapping", d);
      globals.push_back(index);
      components.push_back(row->definition->Dof());
      status = start->GatherFrom(in->start_reaction, &index, 1, i++, d);
      if (status != StatusCode::kOk) return status;
    }
    status = Interpolate(start.get(), end.get(), in->lambda, value.get(), d);
    if (status != StatusCode::kOk) return status;
    status = global->ScatterAdd(value.get(), globals.data(), count, d);
    if (status != StatusCode::kOk) return status;
    Vector *raw_force = nullptr, *raw_moment = nullptr;
    status = PackOrigins(value.get(), &components, 0, &raw_force, d);
    std::unique_ptr<Vector> force(raw_force);
    if (status != StatusCode::kOk) return status;
    status = PackOrigins(value.get(), &components, 1, &raw_moment, d);
    std::unique_ptr<Vector> moment(raw_moment);
    if (status != StatusCode::kOk) return status;
    *external = global.release();
    *force_origins = force.release();
    *moment_origins = moment.release();
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return model_storage::AllocationFailure("Nonlinear step load storage", d);
  }
}

StatusCode ApplyNonlinearStepPrescriptions(
    const NonlinearStepConditionInput* in, AnalysisState* state,
    NonlinearFailure* failure, Diagnostic* d) {
  if (failure) *failure = NonlinearFailure::kNone;
  StatusCode status;
  try {
    status = ApplyPrescriptions(in, state, failure, d);
  } catch (const std::bad_alloc&) {
    status =
        model_storage::AllocationFailure("Nonlinear prescription storage", d);
  }
  if (status != StatusCode::kOk && state) state->DiscardCandidate();
  return status;
}
}  // namespace fesa
