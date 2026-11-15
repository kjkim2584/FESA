#include "analysis/nonlinear_assembler.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <limits>
#include <memory>
#include <new>
#include <vector>

#include "adapters/mkl/dense_ops.h"
#include "adapters/mkl/sparse_ops.h"
#include "adapters/openmp/thread_runtime.h"
#include "analysis/beam_dynamic_assembly.h"
#include "analysis/field_measurement.h"
#include "analysis/finite_hht_tangent.h"
#include "analysis/finite_shell_dynamic_forces.h"
#include "analysis/linear_dynamic_element.h"
#include "analysis/nonlinear_element_pose.h"
#include "analysis/nonlinear_task_failure.h"
#include "analysis/static_assembler.h"
#include "fesa/analysis/beam_internal_state.h"
#include "fesa/analysis/incremental_solid_history.h"
#include "fesa/analysis/rotation_state.h"
#include "fesa/elements/beam3d.h"
#include "fesa/elements/incremental_solid_element.h"
#include "fesa/elements/mitc4_shell.h"

namespace fesa {
namespace {
thread_local void (*beam_tangent_observer)(ElementEvaluationRequest,
                                           void*) = nullptr;
thread_local void* beam_tangent_observer_context = nullptr;
}  // namespace
void SetBeamTangentObserverForTest(void (*callback)(ElementEvaluationRequest,
                                                    void*),
                                   void* context) {
  beam_tangent_observer = callback;
  beam_tangent_observer_context = callback ? context : nullptr;
}
StatusCode ValidateNonlinearBeamAssemblyCounts(const AnalysisModel* model,
                                              const DofManager* dofs,
                                              Diagnostic* d) {
  const auto invalid = [d]() {
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation,
                         "Nonlinear Beam assembly count layout/range", d);
  };
  if (!model || !dofs || dofs->NodeCount() <= 0 ||
      dofs->NodeCount() > INT_MAX / 6 ||
      dofs->GlobalCount() != 6 * dofs->NodeCount() ||
      dofs->FreeCount() < 0 || dofs->InternalDofCount() <= 0 ||
      dofs->InternalDofCount() % 2 != 0 ||
      dofs->InternalDofCount() > INT_MAX - dofs->FreeCount() ||
      dofs->EquationCount() != dofs->FreeCount() + dofs->InternalDofCount() ||
      model->Nodes()->size() != static_cast<std::size_t>(dofs->NodeCount()) ||
      model->Elements()->size() !=
          static_cast<std::size_t>(dofs->InternalDofCount() / 2))
    return invalid();
  // Preserve the existing physical connectivity and LP64 size checks, then
  // validate the actual full-general retained-internal contribution count.
  const auto status = ValidateAssemblyCounts(model, dofs, d);
  if (status != StatusCode::kOk) return status;
  std::size_t count = 0;
  int element_row = 0;
  for (const auto& entry : *model->Elements()) {
    if (d) d->element_id = entry.first;
    const auto* beam = dynamic_cast<const Beam3D*>(entry.second);
    if (!beam || beam->node_count() != 2 || beam->dof_count() != 12 ||
        beam->internal_dof_count() != 2 ||
        dofs->InternalElementId(element_row) != entry.first)
      return invalid();
    int local_free = 2;
    for (int end = 0; end < 2; ++end) {
      if (dofs->InternalDof(entry.first, end) != 2 * element_row + end ||
          dofs->InternalEquation(entry.first, end) !=
              dofs->FreeCount() + 2 * element_row + end)
        return invalid();
      const auto* node = beam->GetNode(end);
      if (!node) return invalid();
      const auto found = model->Nodes()->find(node->Id());
      if (found == model->Nodes()->end() || found->second != node ||
          dofs->NodeDofCount(node->Id()) != 6)
        return invalid();
      for (int component = 1; component <= 6; ++component)
        if (dofs->Equation(dofs->GlobalDof(node->Id(), component)) >= 0)
          ++local_free;
    }
    const auto extra = static_cast<std::size_t>(local_free * local_free);
    if (extra > static_cast<std::size_t>(INT_MAX) - count ||
        count + extra > std::numeric_limits<std::size_t>::max() /
                            (sizeof(double) + 2 * sizeof(int)) ||
        !sparse_ops::ValidSizes(dofs->EquationCount(), count + extra))
      return invalid();
    count += extra;
    ++element_row;
  }
  return StatusCode::kOk;
}
namespace {

StatusCode Invalid(const char* operation, Diagnostic* d) {
  return ReportFailure(StatusCode::kInvalidArgument,
                       DiagnosticStage::kValidation, operation, d);
}
StatusCode Allocation(Diagnostic* d) {
  return ReportFailure(StatusCode::kAllocationFailure,
                       DiagnosticStage::kAllocation, "Nonlinear assembly", d);
}
StatusCode MakeVector(int size, std::unique_ptr<Vector>* out, Diagnostic* d) {
  Vector* raw = nullptr;
  const auto status = Vector::Create(size, &raw, d);
  out->reset(raw);
  return status;
}
StatusCode MakeMatrix(int rows, int columns, std::unique_ptr<Matrix>* out,
                      Diagnostic* d) {
  Matrix* raw = nullptr;
  const auto status = Matrix::Create(rows, columns, &raw, d);
  out->reset(raw);
  return status;
}
bool Finite(const Vector* vector) {
  return vector && dense_ops::AllFinite(vector->Size(), vector->Data());
}
bool ValidOrigins(const Vector* origins, bool allow_zero) {
  if (!origins) return true;
  if (!Finite(origins)) return false;
  for (int i = 0; i < origins->Size(); ++i)
    if (!allow_zero && origins->Data()[i] == 0) return false;
  return true;
}

int FiniteSize(const Element* element) {
  return element->dof_count() + element->internal_dof_count();
}
int LocalEquation(const Element* element, const DofManager* dofs, int local) {
  if (local >= element->dof_count())
    return dofs->InternalEquation(element->Id(), local - element->dof_count());
  return dofs->Equation(
      dofs->GlobalDof(element->GetNode(local / element->dofs_per_node())->Id(),
                      local % element->dofs_per_node() + 1));
}

// Only metadata and dimensions are examined here. No trial force or tangent is
// requested until every element's optional interface has passed this preflight.
StatusCode Validate(const AnalysisModel* model, const DofManager* dofs,
                    const AnalysisState* state,
                    const ElementWorkspace* const* workspaces,
                    int workspace_count, Diagnostic* d,
                    bool geometry_only = false) {
  if (!model || !dofs || !state || InElementExecution() ||
      dofs->NodeCount() <= 0 || dofs->NodeCount() > INT_MAX / 6 ||
      dofs->GlobalCount() != state->GlobalCount() ||
      !state->HasValidCandidate() || !state->CandidateReaction() ||
      (state->Rotations() &&
       (state->Rotations()->NodeCount() != dofs->NodeCount() ||
        !state->Rotations()->CandidateOrientations() ||
        !state->Rotations()->CandidateBodySpin())) ||
      workspace_count < 0 || (!workspaces && workspace_count != 0) ||
      (workspaces && workspace_count < ConfiguredThreadCount()) ||
      model->Nodes()->size() != static_cast<std::size_t>(dofs->NodeCount()) ||
      model->Elements()->size() > static_cast<std::size_t>(INT_MAX - 1))
    return Invalid("Nonlinear assembly layout", d);
  int offset = 0;
  for (int row = 0; row < dofs->NodeCount(); ++row) {
    const int id = dofs->NodeId(row);
    const auto found = model->Nodes()->find(id);
    if (found == model->Nodes()->end() || !found->second ||
        (found->second->dof_count() != 3 && found->second->dof_count() != 6) ||
        found->second->dof_count() != dofs->NodeDofCount(id) ||
        (found->second->dof_count() == 6 && !state->Rotations()) ||
        dofs->GlobalDof(id, 1) != offset)
      return Invalid("Nonlinear node layout", d);
    offset += dofs->NodeDofCount(id);
  }
  if (offset != state->GlobalCount()) return Invalid("Nonlinear offsets", d);
  const bool beam_role = dofs->InternalDofCount() > 0;
  if (beam_role) {
    const auto* internal = state->BeamInternals();
    if (!internal || dofs->GlobalCount() != 6 * dofs->NodeCount() ||
        model->Elements()->size() !=
            static_cast<std::size_t>(internal->ElementCount()) ||
        internal->ElementCount() > INT_MAX / 2 ||
        dofs->InternalDofCount() != 2 * internal->ElementCount() ||
        !internal->CandidatePosition() || !internal->CommittedPosition() ||
        internal->CandidatePosition()->Size() != dofs->InternalDofCount() ||
        internal->CommittedPosition()->Size() != dofs->InternalDofCount() ||
        dofs->EquationCount() != dofs->FreeCount() + dofs->InternalDofCount())
      return Invalid("Nonlinear Beam state layout", d);
    const auto status = state->ValidateMapping(dofs, d);
    if (status != StatusCode::kOk) return status;
    for (int row = 0; row < internal->ElementCount(); ++row) {
      const int id = internal->ElementId(row);
      if (dofs->InternalElementId(row) != id ||
          internal->Offset(id) != 2 * row ||
          dofs->InternalDof(id, 0) != 2 * row ||
          dofs->InternalDof(id, 1) != 2 * row + 1 ||
          dofs->InternalEquation(id, 0) != dofs->FreeCount() + 2 * row ||
          dofs->InternalEquation(id, 1) != dofs->FreeCount() + 2 * row + 1)
        return Invalid("Nonlinear Beam internal offsets", d);
    }
  }
  if (!state->Rotations()) {
    const auto status = state->ValidateMapping(dofs, d);
    if (status != StatusCode::kOk) return status;
  }
  std::size_t force_count = 0;
  for (const auto& entry : *model->Elements()) {
    const Element* element = entry.second;
    if (d) d->element_id = entry.first;
    if (!element || element->node_count() <= 0 ||
        element->node_count() > INT_MAX / 6 || (element->dofs_per_node() != 3 && element->dofs_per_node() != 6) ||
        element->dof_count() != element->dofs_per_node() * element->node_count() ||
        !dense_ops::ValidDimensions(element->dof_count(), element->dof_count()))
      return Invalid("Nonlinear element layout", d);
    if (beam_role ? !dynamic_cast<const Beam3D*>(element) ||
                        element->internal_dof_count() != 2
                  : element->internal_dof_count() != 0)
      return Invalid("Nonlinear internal element role", d);
    if (element->dofs_per_node() == 3) {
      if (!workspaces) return Invalid("Solid worker contexts required", d);
      for (int worker = 0; worker < ConfiguredThreadCount(); ++worker)
        if (!workspaces[worker]) return Invalid("Null Solid worker context", d);
    }
    if (!geometry_only && !element->SupportsTrialEvaluation())
      return ReportFailure(StatusCode::kUnsupported,
                           DiagnosticStage::kValidation,
                           "Element trial evaluation unavailable", d);
    if (static_cast<std::size_t>(FiniteSize(element)) > INT_MAX - force_count)
      return Invalid("Nonlinear contribution count", d);
    force_count += FiniteSize(element);
    for (int node = 0; node < element->node_count(); ++node) {
      const Node* connected = element->GetNode(node);
      if (!connected) return Invalid("Nonlinear connectivity", d);
      const auto found = model->Nodes()->find(connected->Id());
      if (found == model->Nodes()->end() || found->second != connected ||
          dofs->GlobalDof(connected->Id(), 1) < 0 ||
          dofs->NodeDofCount(connected->Id()) != element->dofs_per_node())
        return Invalid("Nonlinear connectivity mapping", d);
    }
  }
  return StatusCode::kOk;
}

bool OwnedPoseVector(const AnalysisState* state, const Vector* vector) {
  const auto* internal = state->BeamInternals();
  if (internal && vector &&
      (vector == internal->CandidatePosition() ||
       vector == internal->CommittedPosition() ||
       vector == internal->CandidateVelocity() ||
       vector == internal->CommittedVelocity() ||
       vector == internal->CandidateAcceleration() ||
       vector == internal->CommittedAcceleration()))
    return true;
  const auto* motion = state->Dynamics();
  if (motion && vector &&
      (vector == motion->CandidateVelocity() ||
       vector == motion->CommittedVelocity() ||
       vector == motion->CandidateAcceleration() ||
       vector == motion->CommittedAcceleration()))
    return true;
  return vector == state->CommittedDisplacement() ||
         vector == state->CommittedReaction() ||
         vector == state->CandidateDisplacement() ||
         (state->Rotations() &&
          (vector == state->Rotations()->CandidateBodySpin() ||
           vector == state->Rotations()->CommittedBodySpin()));
}
StatusCode ValidateResidual(const NonlinearResidualInput* input,
                            const NonlinearResidualOutput* output,
                            bool allow_zero_origins, Diagnostic* d) {
  if (!input || !output || !input->controls || !input->controls->displacement ||
      !input->controls->rotation || !output->residual || !output->reaction ||
      !output->force || !output->moment || output->force == output->moment ||
      output->residual == output->reaction || !Finite(input->external) ||
      !ValidOrigins(input->force_origins, allow_zero_origins) ||
      !ValidOrigins(input->moment_origins, allow_zero_origins) ||
      !std::isfinite(input->previous_force_max_mean) ||
      !std::isfinite(input->previous_moment_max_mean) ||
      input->previous_force_max_mean < 0 || input->previous_moment_max_mean < 0)
    return Invalid("Nonlinear residual input", d);
  auto status = Validate(input->model, input->dofs, input->state,
                         input->workspaces, input->workspace_count, d);
  if (status != StatusCode::kOk) return status;
  const int size = input->state->GlobalCount();
  const int internal_count = input->dofs->InternalDofCount();
  if (internal_count > 0
          ? !output->internal_residual ||
                output->internal_residual->Size() != internal_count
          : output->internal_residual || input->internal_actual_correction)
    return Invalid("Nonlinear internal residual role/shape", d);
  if (input->internal_actual_correction &&
      (!Finite(input->internal_actual_correction) ||
       input->internal_actual_correction->Size() != internal_count))
    return Invalid("Nonlinear internal correction shape", d);
  if (output->internal_residual &&
      (output->internal_residual == output->residual ||
       output->internal_residual == output->reaction ||
       output->internal_residual == input->state->CandidateReaction() ||
       OwnedPoseVector(input->state, output->internal_residual)))
    return Invalid("Nonlinear internal residual output alias", d);
  if (input->external->Size() != size || output->residual->Size() != size ||
      output->reaction->Size() != size ||
      (input->actual_correction && (input->actual_correction->Size() != size ||
                                    !Finite(input->actual_correction))))
    return Invalid("Nonlinear residual dimensions", d);
  for (const Vector* vector :
       {input->external, input->force_origins, input->moment_origins,
        input->actual_correction, input->internal_actual_correction})
    if (vector && (vector == output->residual || vector == output->reaction ||
                   vector == output->internal_residual))
      return Invalid("Nonlinear residual alias", d);
  if (OwnedPoseVector(input->state, output->residual) ||
      OwnedPoseVector(input->state, output->reaction) ||
      output->residual == input->state->CandidateReaction())
    return Invalid("Nonlinear state output alias", d);
  for (const FieldControls* controls :
       {input->controls->displacement, input->controls->rotation})
    if (!std::isfinite(controls->epsilon_l) || controls->epsilon_l < 0)
      return Invalid("Nonlinear active threshold", d);
  return StatusCode::kOk;
}

StatusCode ValidateFiniteDynamic(const AnalysisModel* model,
                                 const DofManager* dofs,
                                 const AnalysisState* state, MassMode mode,
                                 Diagnostic* d) {
  if (!state->Rotations() || !state->Dynamics() ||
      !state->Dynamics()->CandidateVelocity() ||
      !state->Dynamics()->CandidateAcceleration() ||
      model->Elements()->empty() ||
      (mode != MassMode::kConsistent && mode != MassMode::kRowSum))
    return Invalid("Finite dynamic context", d);
  auto status = state->ValidateMapping(dofs, d);
  if (status != StatusCode::kOk) return status;
  for (int row = 0; row < dofs->NodeCount(); ++row)
    if (dofs->NodeDofCount(dofs->NodeId(row)) != 6)
      return Invalid("Finite dynamic requires six physical DOFs", d);
  const bool beam_role = dofs->InternalDofCount() > 0;
  if (beam_role) {
    status = ValidateNonlinearBeamAssemblyCounts(model, dofs, d);
    if (status != StatusCode::kOk) return status;
    const auto* internal = state->BeamInternals();
    if (!internal || !internal->CandidateVelocity() ||
        !internal->CandidateAcceleration())
      return Invalid("Finite Beam joint motion", d);
    if (mode != MassMode::kConsistent)
      return ReportFailure(StatusCode::kUnsupported,
                           DiagnosticStage::kValidation,
                           "Finite Beam requires consistent mass", d);
  }
  for (const auto& entry : *model->Elements()) {
    if (d) d->element_id = entry.first;
    const auto* shell = dynamic_cast<const Mitc4Shell*>(entry.second);
    const auto* beam = dynamic_cast<const Beam3D*>(entry.second);
    if ((!beam_role && !shell) || (beam_role && !beam))
      return ReportFailure(StatusCode::kUnsupported,
                           DiagnosticStage::kValidation,
                           "Finite HHT element family", d);
    double density = 0;
    status = (beam_role ? beam->Section()->GetMaterial() : shell->material())
                 ->Density(&density, d);
    if (status != StatusCode::kOk) return status;
    if ((!beam_role && !shell->SupportsMass(mode)) ||
        !entry.second->SupportsDynamicEvaluation(
            DynamicKinematics::kFiniteRotation))
      return ReportFailure(StatusCode::kUnsupported,
                           DiagnosticStage::kValidation,
                           "Finite HHT dynamic capability", d);
  }
  return StatusCode::kOk;
}

// H6 composition shares the structural callback's copied physical pose. All
// numerical buffers live only for this item; no new sparse assembler or cache.
StatusCode ComposeLocalFiniteTangent(
    const Element* element, const AnalysisState* state,
    const NonlinearElementPose* pose,
    const FiniteHhtLinearization* linearization,
    const ElementWorkspace* workspace, Matrix* tangent,
    TrialEvaluationFailure* failure, Diagnostic* d) {
  const int physical_size = element->dof_count();
  const int size = FiniteSize(element);
  const int nodes = element->node_count();
  const bool beam = element->internal_dof_count() == 2;
  std::unique_ptr<Vector> vectors[6];
  for (int i = 0; i < (beam ? 6 : 4); ++i) {
    auto status = MakeVector(i < 2   ? physical_size
                             : i < 4 ? size
                                     : 2,
                             &vectors[i], d);
    if (status != StatusCode::kOk) return status;
  }
  auto status = vectors[0]->GatherFrom(state->Dynamics()->CandidateVelocity(),
                                       pose->GlobalDofs(), physical_size, 0, d);
  if (status == StatusCode::kOk)
    status = vectors[1]->GatherFrom(state->Dynamics()->CandidateAcceleration(),
                                    pose->GlobalDofs(), physical_size, 0, d);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<Matrix> matrices[8];
  for (int i = 0; i < 8; ++i) {
    status =
        MakeMatrix(i < 6 ? size : nodes, i < 6 ? size : 9, &matrices[i], d);
    if (status != StatusCode::kOk) return status;
  }
  for (int node = 0; node < nodes; ++node) {
    const int row = pose->GlobalDofs()[6 * node] / 6;
    status = matrices[6]->CopyBlockFrom(
        linearization->rotation->velocity_per_rotation, row, 0, 1, 9, node, 0,
        d);
    if (status == StatusCode::kOk)
      status = matrices[7]->CopyBlockFrom(
          linearization->rotation->acceleration_per_rotation, row, 0, 1, 9,
          node, 0, d);
    if (status != StatusCode::kOk) return status;
  }
  ElementMotionInput input;
  input.kinematics = DynamicKinematics::kFiniteRotation;
  input.finite_state = pose->Trial();
  input.velocity = vectors[0].get();
  input.acceleration = vectors[1].get();
  if (beam) {
    const auto* internal = state->BeamInternals();
    const int offset = internal->Offset(element->Id());
    const int indices[] = {offset, offset + 1};
    status =
        vectors[4]->GatherFrom(internal->CandidateVelocity(), indices, 2, 0, d);
    if (status == StatusCode::kOk)
      status = vectors[5]->GatherFrom(internal->CandidateAcceleration(),
                                      indices, 2, 0, d);
    if (status != StatusCode::kOk) return status;
    input.internal_velocity = vectors[4].get();
    input.internal_acceleration = vectors[5].get();
  }
  ElementDynamicResult physical;
  physical.inertia_force = vectors[2].get();
  physical.damping_force = vectors[3].get();
  physical.inertia_displacement = matrices[0].get();
  physical.inertia_velocity = matrices[1].get();
  physical.inertia_acceleration = matrices[2].get();
  physical.damping_displacement = matrices[3].get();
  physical.damping_velocity = matrices[4].get();
  status = element->EvaluateDynamic(
      &input, linearization->mode, workspace,
      DynamicEvaluationRequest::kForcesAndLinearization, &physical, failure, d);
  if (status != StatusCode::kOk) return status;
  RotationalLinearization rotation;
  rotation.velocity_per_rotation = matrices[6].get();
  rotation.acceleration_per_rotation = matrices[7].get();
  const auto compose =
      beam ? ComposeBeamFiniteHhtTangent : ComposeFiniteHhtTangent;
  status = compose(linearization->alpha, linearization->translation, &rotation,
                   tangent, &physical, matrices[5].get(), d);
  return status == StatusCode::kOk ? tangent->CopyFrom(matrices[5].get(), d)
                                   : status;
}

struct TaskContext {
  const DofManager* dofs = nullptr;
  const AnalysisState* state = nullptr;
  const std::vector<const Element*>* elements = nullptr;
  const std::vector<int>* offsets = nullptr;
  const std::vector<int>* matrix_offsets = nullptr;
  std::vector<TrialEvaluationFailure>* failures = nullptr;
  Vector* values = nullptr;  // Disjoint prefix ranges until the joined barrier.
  const FiniteHhtLinearization* dynamic = nullptr;  // Borrowed through join.
  bool tangent = false;
  bool geometry_only = false;
  bool committed_pose = false;  // Physical source only, not a state mutation.
  bool use_committed_geometric_stress = false;  // Read-only through join.
  ElementEvaluationRequest tangent_request = ElementEvaluationRequest::kTangent;
  MatrixType type = MatrixType::kGeneral;
  const ElementWorkspace* const* workspaces = nullptr;
  const AnalysisModel* recovery_model = nullptr;  // Borrowed mutable payloads.
};

StatusCode StageIncrementalStress(const Element* element,
    const ElementTrialState* trial, const AnalysisState* state, Diagnostic* d) {
  if (!trial->incremental_solid) return StatusCode::kOk;
  auto* history = state->SolidHistory();
  if (!history) return Invalid("Incremental assembly history", d);
  return history->StageStress(element->Id(),
                              trial->incremental_solid->candidate_stress, d);
}

StatusCode EvaluateElement(int item, int worker, void* context, Diagnostic* d) {
  auto* task = static_cast<TaskContext*>(context);
  const Element* element = (*task->elements)[item];
  d->element_id = element->Id();
  try {
    const int size = FiniteSize(element);
    const auto* workspace =
        task->workspaces ? task->workspaces[worker] : nullptr;
    IncrementalElasticWorkspace* incremental_workspace = nullptr;
    const auto* incremental =
        dynamic_cast<const IncrementalSolidElement*>(element);
    if (workspace && incremental) {
      incremental_workspace =
          workspace->incremental_solid_workspace(incremental->formulation());
    }
    NonlinearElementPose* raw_pose = nullptr;
    auto status = NonlinearElementPose::Create(
        element, task->dofs, task->state, task->committed_pose, &raw_pose, d,
        incremental_workspace,
        task->use_committed_geometric_stress);
    std::unique_ptr<NonlinearElementPose> pose(raw_pose);
    if (status != StatusCode::kOk) return status;
    std::unique_ptr<Vector> local_force;
    std::unique_ptr<Matrix> tangent;
    if (!task->recovery_model && !task->geometry_only) {
      if (task->tangent)
        status = MakeMatrix(size, size, &tangent, d);
      else
        status = MakeVector(size, &local_force, d);
    }
    if (status != StatusCode::kOk) return status;
    const auto* trial = pose->Trial();
    auto* typed = &(*task->failures)[item];
    if (task->geometry_only)
      return element->CheckTrialGeometry(trial, workspace, typed, d);
    if (task->recovery_model) {
      status = task->recovery_model->Elements()
          ->at(element->Id())
          ->UpdateTrialResults(trial, workspace, typed, d);
      if (status == StatusCode::kOk)
        status = StageIncrementalStress(element, trial, task->state, d);
      return status;
    }
    status = element->EvaluateTrial(
        trial, workspace,
        task->tangent ? task->tangent_request
                      : ElementEvaluationRequest::kInternalForce,
        local_force.get(), tangent.get(), typed, d);
    if (status != StatusCode::kOk) return status;
    if (!task->committed_pose) {
      status = StageIncrementalStress(element, trial, task->state, d);
      if (status != StatusCode::kOk) return status;
    }
    if (task->dynamic) {
      status = ComposeLocalFiniteTangent(element, task->state, pose.get(),
                                          task->dynamic, workspace,
                                          tangent.get(), typed, d);
      if (status != StatusCode::kOk) return status;
    }
    std::vector<int> indices;
    if (task->tangent) {
      for (int row = 0; row < size; ++row) {
        const int row_equation = LocalEquation(element, task->dofs, row);
        if (row_equation < 0) continue;
        for (int column = 0; column < size; ++column) {
          const int column_equation =
              LocalEquation(element, task->dofs, column);
          if (column_equation >= 0 && (task->type == MatrixType::kGeneral ||
                                       column_equation >= row_equation))
            indices.push_back(size * row + column);
        }
      }
      return tangent->GatherTo(indices.data(), static_cast<int>(indices.size()),
                               task->values, (*task->matrix_offsets)[item], d);
    }
    indices.resize(size);
    for (int i = 0; i < size; ++i) indices[i] = i;
    return task->values->GatherFrom(local_force.get(), indices.data(), size,
                                    (*task->offsets)[item], d);
  } catch (const std::bad_alloc&) {
    return Allocation(d);
  }
}

StatusCode Dispatch(TaskContext* task, NonlinearFailure* failure,
                    Diagnostic* d) {
  const int workers = ConfiguredThreadCount();
  std::unique_ptr<Diagnostic[]> diagnostics(new (std::nothrow)
                                                Diagnostic[workers]);
  if (!diagnostics) return Allocation(d);
  Diagnostic aggregate;
  Diagnostic* selected = d ? d : &aggregate;
  const auto status = ExecuteElementTasks(
      static_cast<int>(task->elements->size()), workers, EvaluateElement, task,
      diagnostics.get(), nullptr, selected);
  if (status != StatusCode::kOk && task->state->SolidHistory())
    task->state->SolidHistory()->DiscardCandidate();
  return SelectNonlinearTaskFailure(
      status, task->elements->data(), static_cast<int>(task->elements->size()),
      task->failures->data(), diagnostics.get(), workers, selected, failure);
}

void CollectElements(const AnalysisModel* model,
                     std::vector<const Element*>* elements,
                     std::vector<int>* offsets) {
  offsets->push_back(0);
  for (const auto& entry : *model->Elements()) {
    elements->push_back(entry.second);
    offsets->push_back(offsets->back() + FiniteSize(entry.second));
  }
}

StatusCode SetConstrainedReaction(const DofManager* dofs,
                                  const Vector* residual, Vector* reaction,
                                  Diagnostic* d) {
  StatusCode status = StatusCode::kOk;
  std::vector<int> constrained;
  for (int global = 0; global < dofs->GlobalCount(); ++global)
    if (dofs->Equation(global) == -1) constrained.push_back(global);
  if (!constrained.empty()) {
    std::unique_ptr<Vector> values;
    status = MakeVector(static_cast<int>(constrained.size()), &values, d);
    if (status != StatusCode::kOk) return status;
    status =
        values->GatherFrom(residual, constrained.data(), values->Size(), 0, d);
    if (status != StatusCode::kOk) return status;
    status = values->Scale(-1, d);
    if (status != StatusCode::kOk) return status;
    status = reaction->ScatterAdd(values.get(), constrained.data(),
                                  values->Size(), d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}
StatusCode MeasureNonlinearResidual(
    const NonlinearResidualInput* input, const Vector* residual,
    const Vector* packed, int block_count, FieldResidual* force_out,
    FieldResidual* moment_out, NonlinearFailure* failure, Diagnostic* d,
    bool material_damping = false, const Vector* internal_residual = nullptr) {
  FieldResidual force, moment;
  std::unique_ptr<Vector> force_flux, moment_flux;
  Vector* raw_flux = nullptr;
  StatusCode status = StatusCode::kOk;
  status = CreateElementFieldFlux(
      input->model, packed, block_count, 0, input->force_origins, &raw_flux, d,
      material_damping, input->dofs->InternalDofCount() > 0);
  force_flux.reset(raw_flux);
  if (status != StatusCode::kOk) return status;
  status = CreateElementFieldFlux(
      input->model, packed, block_count, 1, input->moment_origins, &raw_flux, d,
      material_damping, input->dofs->InternalDofCount() > 0);
  moment_flux.reset(raw_flux);
  if (status != StatusCode::kOk) return status;
  status =
      ComputeSpatialFluxMean(force_flux.get(), input->controls->displacement,
                             input->previous_force_max_mean, 0, &force, d);
  if (status != StatusCode::kOk) return status;
  status =
      ComputeSpatialFluxMean(moment_flux.get(), input->controls->rotation,
                             input->previous_moment_max_mean, 1, &moment, d);
  if (status != StatusCode::kOk) return status;
  status = internal_residual
               ? MeasureBeamResidualField(residual, internal_residual,
                                          input->dofs, &force, d)
               : MeasureResidualField(residual, input->dofs, 0, &force, d);
  if (status != StatusCode::kOk) return status;
  status = MeasureResidualField(residual, input->dofs, 1, &moment, d);
  if (status != StatusCode::kOk) return status;
  bool rotation_branch = false;
  status = MeasureNonlinearMotionFields(input->dofs, input->state,
                                        input->actual_correction, &force,
                                        &moment, &rotation_branch, d);
  if (status != StatusCode::kOk) {
    if (status == StatusCode::kEquilibriumFailure && rotation_branch)
      *failure = NonlinearFailure::kRotationBranch;
    return status;
  }
  if (internal_residual) {
    status =
        MeasureBeamInternalMotion(input->dofs, input->state,
                                  input->internal_actual_correction, &force, d);
    if (status != StatusCode::kOk) return status;
  }
  *force_out = force;
  *moment_out = moment;
  return StatusCode::kOk;
}
}  // namespace

StatusCode EvaluateNonlinearStructuralForce(
    const AnalysisModel* model, const DofManager* dofs,
    const AnalysisState* state, bool committed,
    const ElementWorkspace* const* workspaces, int workspace_count,
    Vector** global_force, Vector** element_force,
    double* characteristic_length, NonlinearFailure* failure, Diagnostic* d,
    Vector** global_internal_force) {
  if (global_internal_force) *global_internal_force = nullptr;
  if (global_force) *global_force = nullptr;
  if (element_force) *element_force = nullptr;
  if (failure) *failure = NonlinearFailure::kNone;
  if (!global_force || !element_force || global_force == element_force ||
      (global_internal_force && (global_internal_force == global_force ||
                                 global_internal_force == element_force)) ||
      !characteristic_length || !failure)
    return Invalid("Nonlinear structural force outputs", d);
  auto status = Validate(model, dofs, state, workspaces, workspace_count, d);
  if (status != StatusCode::kOk) return status;
  const bool beam_role = dofs->InternalDofCount() > 0;
  if (beam_role && !global_internal_force)
    return Invalid("Nonlinear Beam internal force output required", d);
  try {
    std::vector<const Element*> elements;
    std::vector<int> offsets, global_indices, physical_indices,
        internal_indices;
    CollectElements(model, &elements, &offsets);
    std::unique_ptr<Vector> packed, global, lengths, internal, physical;
    if (offsets.back() > 0) {
      status = MakeVector(offsets.back(), &packed, d);
      if (status != StatusCode::kOk) return status;
    }
    status = MakeVector(state->GlobalCount(), &global, d);
    if (status != StatusCode::kOk) return status;
    if (beam_role) {
      status = MakeVector(dofs->InternalDofCount(), &internal, d);
      if (status != StatusCode::kOk) return status;
    }
    std::vector<TrialEvaluationFailure> failures(elements.size(),
                                                 TrialEvaluationFailure::kNone);
    TaskContext task;
    task.dofs = dofs;
    task.state = state;
    task.workspaces = workspaces;
    task.elements = &elements;
    task.offsets = &offsets;
    task.failures = &failures;
    task.values = packed.get();
    task.committed_pose = committed;
    status = Dispatch(&task, failure, d);
    if (status != StatusCode::kOk) return status;
    double length_mean = 0;
    if (!elements.empty()) {
      status = MakeVector(static_cast<int>(elements.size()), &lengths, d);
      if (status != StatusCode::kOk) return status;
    }
    for (std::size_t item = 0; item < elements.size(); ++item) {
      const Element* element = elements[item];
      if (d) d->element_id = element->Id();
      double length = 0;
      status = element->CharacteristicLength(&length, d);
      if (status != StatusCode::kOk) return status;
      if (!std::isfinite(length) || length <= 0)
        return Invalid("Nonlinear characteristic length", d);
      status = lengths->Set(static_cast<int>(item), length, d);
      if (status != StatusCode::kOk) return status;
      for (int local = 0; local < element->dof_count(); ++local) {
        global_indices.push_back(dofs->GlobalDof(
            element->GetNode(local / element->dofs_per_node())->Id(),
            local % element->dofs_per_node() + 1));
        physical_indices.push_back(offsets[item] + local);
      }
      if (beam_role)
        for (int end = 0; end < 2; ++end)
          internal_indices.push_back(offsets[item] + element->dof_count() +
                                     end);
    }
    if (lengths) {
      status = lengths->Scale(1.0 / lengths->Size(), d);
      if (status == StatusCode::kOk)
        status = lengths->AbsoluteSum(&length_mean, d);
      if (status != StatusCode::kOk) return status;
    }
    if (packed) {
      const Vector* physical_source = packed.get();
      if (beam_role) {
        status =
            MakeVector(static_cast<int>(physical_indices.size()), &physical, d);
        if (status == StatusCode::kOk)
          status = physical->GatherFrom(packed.get(), physical_indices.data(),
                                        physical->Size(), 0, d);
        if (status == StatusCode::kOk)
          status = internal->GatherFrom(packed.get(), internal_indices.data(),
                                        internal->Size(), 0, d);
        if (status != StatusCode::kOk) return status;
        physical_source = physical.get();
      }
      status = global->ScatterAdd(physical_source, global_indices.data(),
                                  physical_source->Size(), d);
      if (status != StatusCode::kOk) return status;
    }
    *global_force = global.release();
    *element_force = packed.release();
    if (global_internal_force) *global_internal_force = internal.release();
    *characteristic_length = length_mean;
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return Allocation(d);
  }
}

StatusCode EvaluateNonlinearResidual(const NonlinearResidualInput* input,
                                     NonlinearResidualOutput* output,
                                     NonlinearFailure* failure, Diagnostic* d) {
  if (!failure) return Invalid("Nonlinear residual failure output", d);
  *failure = NonlinearFailure::kNone;
  auto status = ValidateResidual(input, output, false, d);
  if (status != StatusCode::kOk) return status;
  try {
    Vector *raw_packed = nullptr, *raw_residual = nullptr,
           *raw_internal = nullptr;
    double characteristic_length = 0;
    status = EvaluateNonlinearStructuralForce(
        input->model, input->dofs, input->state, false, input->workspaces,
        input->workspace_count, &raw_residual, &raw_packed,
        &characteristic_length, failure, d,
        input->dofs->InternalDofCount() > 0 ? &raw_internal : nullptr);
    std::unique_ptr<Vector> packed(raw_packed), residual(raw_residual),
        reaction, internal(raw_internal);
    if (status != StatusCode::kOk) return status;
    if (internal) {
      status = internal->Scale(-1, d);
      if (status != StatusCode::kOk) return status;
    }
    status = MakeVector(input->state->GlobalCount(), &reaction, d);
    if (status != StatusCode::kOk) return status;
    status = residual->Scale(-1, d);
    if (status != StatusCode::kOk) return status;
    status = residual->AddScaled(1, input->external, d);
    if (status != StatusCode::kOk) return status;
    status =
        SetConstrainedReaction(input->dofs, residual.get(), reaction.get(), d);
    if (status != StatusCode::kOk) return status;
    FieldResidual force, moment;
    status =
        MeasureNonlinearResidual(input, residual.get(), packed.get(), 1, &force,
                                 &moment, failure, d, false, internal.get());
    if (status != StatusCode::kOk) return status;
    status = output->residual->CopyFrom(residual.get(), d);
    if (status != StatusCode::kOk) return status;
    status = output->reaction->CopyFrom(reaction.get(), d);
    if (status != StatusCode::kOk) return status;
    if (internal) {
      status = output->internal_residual->CopyFrom(internal.get(), d);
      if (status != StatusCode::kOk) return status;
    }
    *output->force = force;
    *output->moment = moment;
    output->characteristic_length = characteristic_length;
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return Allocation(d);
  }
}

StatusCode EvaluateFiniteHhtResidual(const NonlinearResidualInput* input,
                                     const Vector* previous_external,
                                     double alpha, MassMode mode,
                                     NonlinearResidualOutput* output,
                                     NonlinearFailure* failure, Diagnostic* d) {
  if (!failure) return Invalid("Finite HHT failure output", d);
  *failure = NonlinearFailure::kNone;
  auto status = ValidateResidual(input, output, true, d);
  if (status != StatusCode::kOk) return status;
  const auto* motion = input->state->Dynamics();
  if (!std::isfinite(alpha) || alpha < -1.0 / 3 || alpha > 0 ||
      !Finite(previous_external) ||
      previous_external->Size() != input->state->GlobalCount() || !motion ||
      (alpha != 0 &&
       (!motion->CommittedVelocity() || !motion->CommittedAcceleration())) ||
      !motion->CandidateVelocity() || !motion->CandidateAcceleration() ||
      !input->state->Rotations() || input->model->Elements()->empty() ||
      (mode != MassMode::kConsistent && mode != MassMode::kRowSum))
    return Invalid("Finite HHT residual context", d);
  for (const Vector* vector :
       {previous_external, motion->CommittedVelocity(),
        motion->CommittedAcceleration(), motion->CandidateVelocity(),
        motion->CandidateAcceleration()})
    if (vector && (vector == output->residual || vector == output->reaction ||
                   vector == output->internal_residual))
      return Invalid("Finite HHT residual alias", d);
  status =
      ValidateFiniteDynamic(input->model, input->dofs, input->state, mode, d);
  if (status != StatusCode::kOk) return status;
  try {
    const bool beam = input->dofs->InternalDofCount() > 0;
    Vector *raw_current = nullptr, *raw_packed = nullptr,
           *raw_internal = nullptr;
    double length = 0;
    status = EvaluateNonlinearStructuralForce(
        input->model, input->dofs, input->state, false, input->workspaces,
        input->workspace_count, &raw_current, &raw_packed, &length, failure, d,
        beam ? &raw_internal : nullptr);
    std::unique_ptr<Vector> current(raw_current), packed(raw_packed),
        internal(raw_internal);
    if (status != StatusCode::kOk) return status;
    FiniteShellDynamicForces* raw_dynamic = nullptr;
    BeamDynamicForces* raw_beam_dynamic = nullptr;
    if (beam)
      status =
          BeamDynamicForces::Create(input->model, input->dofs, input->state,
                                    mode, false, &raw_beam_dynamic, failure, d,
                                    input->workspaces, input->workspace_count);
    else
      status = FiniteShellDynamicForces::Create(input->model, input->dofs,
                                                input->state, mode, false,
                                                &raw_dynamic, failure, d);
    std::unique_ptr<FiniteShellDynamicForces> dynamic(raw_dynamic);
    std::unique_ptr<BeamDynamicForces> beam_dynamic(raw_beam_dynamic);
    if (status != StatusCode::kOk) return status;
    const Vector* inertia = beam ? beam_dynamic->Inertia() : dynamic->Inertia();
    const Vector* damping = beam ? beam_dynamic->Damping() : dynamic->Damping();
    const Vector* inertia_origins =
        beam ? beam_dynamic->InertiaOrigins() : dynamic->InertiaOrigins();
    const Vector* damping_origins =
        beam ? beam_dynamic->DampingOrigins() : dynamic->DampingOrigins();
    bool damped = false;
    for (const auto& entry : *input->model->Elements())
      damped = damped || HasMaterialRayleighDamping(entry.second);
    std::unique_ptr<Vector> residual, reaction, physical_residual;
    status = MakeVector(current->Size(), &residual, d);
    if (status != StatusCode::kOk) return status;
    status = residual->CopyFrom(input->external, d);
    if (status == StatusCode::kOk)
      status = residual->AddScaled(-1, current.get(), d);
    if (status == StatusCode::kOk && damped)
      status = residual->AddScaled(-1, damping, d);
    if (status == StatusCode::kOk && beam) {
      status = internal->Scale(-1, d);
      if (status == StatusCode::kOk && damped)
        status = internal->AddScaled(-1, beam_dynamic->InternalDamping(), d);
    }
    if (status != StatusCode::kOk) return status;
    status = MakeVector(current->Size(), &physical_residual, d);
    if (status == StatusCode::kOk)
      status = physical_residual->CopyFrom(residual.get(), d);
    if (status == StatusCode::kOk)
      status = physical_residual->AddScaled(-1, inertia, d);
    if (status != StatusCode::kOk) return status;
    status = residual->Scale(1 + alpha, d);
    if (status == StatusCode::kOk) status = residual->AddScaled(-1, inertia, d);
    if (status == StatusCode::kOk && beam) {
      status = internal->Scale(1 + alpha, d);
      if (status == StatusCode::kOk)
        status = internal->AddScaled(-1, beam_dynamic->InternalInertia(), d);
    }
    if (status != StatusCode::kOk) return status;
    if (alpha != 0) {
      Vector *raw_old = nullptr, *raw_old_origins = nullptr,
             *raw_old_internal = nullptr;
      double old_length = 0;
      status = EvaluateNonlinearStructuralForce(
          input->model, input->dofs, input->state, true, input->workspaces,
          input->workspace_count, &raw_old, &raw_old_origins, &old_length,
          failure, d, beam ? &raw_old_internal : nullptr);
      std::unique_ptr<Vector> old(raw_old), old_origins(raw_old_origins),
          old_internal(raw_old_internal);
      if (status == StatusCode::kOk)
        status = residual->AddScaled(alpha, old.get(), d);
      if (status == StatusCode::kOk)
        status = residual->AddScaled(-alpha, previous_external, d);
      if (status == StatusCode::kOk && beam)
        status = internal->AddScaled(alpha, old_internal.get(), d);
      if (status == StatusCode::kOk && damped) {
        FiniteShellDynamicForces* raw_old_dynamic = nullptr;
        BeamDynamicForces* raw_old_beam = nullptr;
        if (beam)
          status =
              BeamDynamicForces::Create(input->model, input->dofs, input->state,
                                        mode, true, &raw_old_beam, failure, d,
                                        input->workspaces, input->workspace_count);
        else
          status = FiniteShellDynamicForces::Create(
              input->model, input->dofs, input->state, mode, true,
              &raw_old_dynamic, failure, d);
        std::unique_ptr<FiniteShellDynamicForces> old_dynamic(raw_old_dynamic);
        std::unique_ptr<BeamDynamicForces> old_beam(raw_old_beam);
        if (status == StatusCode::kOk)
          status = residual->AddScaled(
              alpha, beam ? old_beam->Damping() : old_dynamic->Damping(), d);
        if (status == StatusCode::kOk && beam)
          status = internal->AddScaled(alpha, old_beam->InternalDamping(), d);
      }
      if (status != StatusCode::kOk) return status;
    }
    status = MakeVector(current->Size(), &reaction, d);
    if (status == StatusCode::kOk)
      status = SetConstrainedReaction(input->dofs, physical_residual.get(),
                                      reaction.get(), d);
    if (status != StatusCode::kOk) return status;
    if (packed->Size() > INT_MAX / 2 ||
        packed->Size() != inertia_origins->Size() ||
        packed->Size() != damping_origins->Size())
      return Invalid("Finite HHT origin shape", d);
    int origin_size = 2 * packed->Size();
    for (const auto& entry : *input->model->Elements()) {
      if (!HasMaterialRayleighDamping(entry.second)) continue;
      if (origin_size > INT_MAX - FiniteSize(entry.second))
        return Invalid("Finite HHT damping origin size", d);
      origin_size += FiniteSize(entry.second);
    }
    std::unique_ptr<Vector> origins;
    status = MakeVector(origin_size, &origins, d);
    if (status != StatusCode::kOk) return status;
    std::vector<int> indices;
    int source = 0, target = 0;
    for (const auto& entry : *input->model->Elements()) {
      const int size = FiniteSize(entry.second);
      indices.resize(size);
      for (int i = 0; i < size; ++i) indices[i] = source + i;
      const Vector* blocks[] = {packed.get(), inertia_origins, damping_origins};
      for (const Vector* values : blocks) {
        if (values == damping_origins &&
            !HasMaterialRayleighDamping(entry.second))
          continue;
        status = origins->GatherFrom(values, indices.data(), size, target, d);
        if (status != StatusCode::kOk) return status;
        target += size;
      }
      source += size;
    }
    FieldResidual force, moment;
    status = MeasureNonlinearResidual(input, residual.get(), origins.get(), 2,
                                      &force, &moment, failure, d, true,
                                      internal.get());
    if (status == StatusCode::kOk)
      status = output->residual->CopyFrom(residual.get(), d);
    if (status == StatusCode::kOk)
      status = output->reaction->CopyFrom(reaction.get(), d);
    if (status == StatusCode::kOk && beam)
      status = output->internal_residual->CopyFrom(internal.get(), d);
    if (status != StatusCode::kOk) return status;
    *output->force = force;
    *output->moment = moment;
    output->characteristic_length = length;
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return Allocation(d);
  }
}

namespace {
StatusCode AssembleTangent(const AnalysisModel* model, const DofManager* dofs,
                           const AnalysisState* state,
                           const Vector* current_residual,
                           const ElementWorkspace* const* workspaces,
                           int workspace_count, MatrixType type,
                           const FiniteHhtLinearization* linearization,
                           EquationSystem** out, NonlinearFailure* failure,
                           Diagnostic* d, bool use_committed_geometric_stress,
                           const Vector* internal_residual,
                           bool beam_material_tangent) {
  if (out) *out = nullptr;
  if (!out || !failure) return Invalid("Nonlinear tangent output", d);
  *failure = NonlinearFailure::kNone;
  auto status = Validate(model, dofs, state, workspaces, workspace_count, d);
  if (status != StatusCode::kOk) return status;
  if (beam_material_tangent) {
    if (linearization || use_committed_geometric_stress ||
        model->Elements()->empty())
      return Invalid("B33 material tangent requires static B33 model", d);
    for (const auto& entry : *model->Elements())
      if (!dynamic_cast<const Beam3D*>(entry.second))
        return ReportFailure(StatusCode::kUnsupported,
                             DiagnosticStage::kValidation,
                             "B33 material tangent requires pure B33 model", d);
  }
  if (!Finite(current_residual) ||
      current_residual->Size() != state->GlobalCount())
    return Invalid("Nonlinear tangent residual", d);
  const int internal_count = dofs->InternalDofCount();
  const int equation_count = dofs->EquationCount();
  if (internal_count > 0 ? !Finite(internal_residual) ||
                               internal_residual->Size() != internal_count ||
                               internal_residual == current_residual ||
                               type != MatrixType::kGeneral
                         : internal_residual != nullptr)
    return Invalid("Nonlinear tangent internal residual/type", d);
  if (type != MatrixType::kGeneral && type != MatrixType::kSymmetricIndefinite)
    return Invalid("Nonlinear tangent matrix type", d);
  if (linearization) {
    const auto* translation = linearization->translation;
    const auto* rotation = linearization->rotation;
    if (!std::isfinite(linearization->alpha) ||
        linearization->alpha < -1.0 / 3 || linearization->alpha > 0 ||
        !translation || !rotation ||
        !std::isfinite(translation->velocity_per_displacement) ||
        !std::isfinite(translation->acceleration_per_displacement) ||
        translation->velocity_per_displacement <= 0 ||
        translation->acceleration_per_displacement <= 0 ||
        dofs->NodeCount() > INT_MAX / 9)
      return Invalid("Finite HHT tangent coefficients", d);
    for (const auto* matrix :
         {rotation->velocity_per_rotation, rotation->acceleration_per_rotation})
      if (!matrix || matrix->Rows() != dofs->NodeCount() ||
          matrix->Columns() != 9 ||
          !dense_ops::AllFinite(9 * dofs->NodeCount(), matrix->Data()))
        return Invalid("Finite HHT tangent angular blocks", d);
    status = ValidateFiniteDynamic(model, dofs, state, linearization->mode, d);
    if (status != StatusCode::kOk) return status;
  }
  if (equation_count == 0) return StatusCode::kOk;
  try {
    std::vector<const Element*> elements;
    std::vector<int> offsets, matrix_offsets(1, 0), rows, columns;
    CollectElements(model, &elements, &offsets);
    std::size_t count = 0;
    for (const Element* element : elements) {
      std::vector<int> free;
      for (int local = 0; local < FiniteSize(element); ++local) {
        const int equation = LocalEquation(element, dofs, local);
        if (equation >= 0) free.push_back(equation);
      }
      const std::size_t extra = type == MatrixType::kGeneral
          ? free.size() * free.size() : free.size() * (free.size() + 1) / 2;
      if (extra > static_cast<std::size_t>(INT_MAX) - count ||
          count + extra > std::numeric_limits<std::size_t>::max() /
                              (sizeof(double) + 2 * sizeof(int)) ||
          !sparse_ops::ValidSizes(equation_count, count + extra))
        return Invalid("Nonlinear general COO count", d);
      count += extra;
      matrix_offsets.push_back(static_cast<int>(count));
      for (int row : free)
        for (int column : free) {
          if (type != MatrixType::kGeneral && column < row) continue;
          rows.push_back(row);
          columns.push_back(column);
        }
    }
    std::unique_ptr<Vector> values, rhs;
    // A zero-entry matrix still provides a nonnull unused pack for tasks whose
    // free-free selection is empty; SparseMatrix receives the true count0.
    status = MakeVector(std::max(1, static_cast<int>(count)), &values, d);
    if (status != StatusCode::kOk) return status;
    std::vector<TrialEvaluationFailure> failures(elements.size(),
                                                 TrialEvaluationFailure::kNone);
    TaskContext task;
    task.dofs = dofs;
    task.state = state;
    task.workspaces = workspaces;
    task.elements = &elements;
    task.offsets = &offsets;
    task.matrix_offsets = &matrix_offsets;
    task.failures = &failures;
    task.values = values.get();
    task.dynamic = linearization;
    task.tangent = true;
    task.use_committed_geometric_stress = use_committed_geometric_stress;
    task.tangent_request = beam_material_tangent
                               ? ElementEvaluationRequest::kMaterialTangent
                               : ElementEvaluationRequest::kTangent;
    task.type = type;
    if (internal_count > 0 && beam_tangent_observer)
      beam_tangent_observer(task.tangent_request,
                            beam_tangent_observer_context);
    status = Dispatch(&task, failure, d);
    if (status != StatusCode::kOk) return status;
    SparseMatrix* raw_matrix = nullptr;
    status = type == MatrixType::kGeneral
                 ? SparseMatrix::CreateGeneral(
                       equation_count, static_cast<int>(count), rows.data(),
                       columns.data(), values.get(), &raw_matrix, d)
                 : SparseMatrix::Create(equation_count, static_cast<int>(count),
                                        rows.data(), columns.data(),
                                        values.get(), &raw_matrix, d);
    std::unique_ptr<SparseMatrix> matrix(raw_matrix);
    if (status != StatusCode::kOk) return status;
    status = MakeVector(equation_count, &rhs, d);
    if (status != StatusCode::kOk) return status;
    std::vector<int> indices(dofs->FreeCount());
    for (int i = 0; i < dofs->FreeCount(); ++i)
      indices[i] = dofs->GlobalDofForEquation(i);
    status = rhs->GatherFrom(current_residual, indices.data(),
                             dofs->FreeCount(), 0, d);
    if (status != StatusCode::kOk) return status;
    if (internal_count > 0) {
      indices.resize(internal_count);
      for (int i = 0; i < internal_count; ++i) indices[i] = i;
      status = rhs->GatherFrom(internal_residual, indices.data(),
                               internal_count, dofs->FreeCount(), d);
      if (status != StatusCode::kOk) return status;
    }
    status = EquationSystem::Create(matrix.get(), rhs.get(), type,
                                    out, d);
    if (status == StatusCode::kOk) {
      matrix.release();
      rhs.release();
    }
    return status;
  } catch (const std::bad_alloc&) {
    return Allocation(d);
  }
}

}  // namespace

StatusCode AssembleNonlinearTangent(
    const AnalysisModel* model, const DofManager* dofs,
    const AnalysisState* state, const Vector* current_residual,
    const ElementWorkspace* const* workspaces, int workspace_count,
    MatrixType type, EquationSystem** out, NonlinearFailure* failure,
    Diagnostic* d, bool use_committed_geometric_stress,
    const Vector* internal_residual, bool beam_material_tangent) {
  return AssembleTangent(model, dofs, state, current_residual, workspaces,
                         workspace_count, type, nullptr, out, failure, d,
                         use_committed_geometric_stress, internal_residual,
                         beam_material_tangent);
}
StatusCode AssembleFiniteHhtTangent(
    const AnalysisModel* model, const DofManager* dofs,
    const AnalysisState* state, const Vector* current_residual,
    const FiniteHhtLinearization* linearization, EquationSystem** out,
    NonlinearFailure* failure, Diagnostic* d, const Vector* internal_residual,
    const ElementWorkspace* const* workspaces, int workspace_count) {
  if (out) *out = nullptr;
  if (failure) *failure = NonlinearFailure::kNone;
  if (!linearization) return Invalid("Finite HHT tangent linearization", d);
  return AssembleTangent(model, dofs, state, current_residual, workspaces,
                         workspace_count, MatrixType::kGeneral, linearization,
                         out, failure, d, false, internal_residual, false);
}

StatusCode CheckNonlinearTrialGeometry(
    const AnalysisModel* model, const DofManager* dofs,
    const AnalysisState* state, const ElementWorkspace* const* workspaces,
    int workspace_count, NonlinearFailure* failure, Diagnostic* d) {
  if (!failure) return Invalid("Nonlinear geometry failure output", d);
  *failure = NonlinearFailure::kNone;
  const auto status =
      Validate(model, dofs, state, workspaces, workspace_count, d, true);
  if (status != StatusCode::kOk) return status;
  try {
    std::vector<const Element*> elements;
    for (const auto& entry : *model->Elements())
      if (entry.second->SupportsTrialGeometryCheck())
        elements.push_back(entry.second);
    if (elements.empty()) return StatusCode::kOk;
    std::vector<TrialEvaluationFailure> failures(elements.size(),
                                                 TrialEvaluationFailure::kNone);
    TaskContext task;
    task.dofs = dofs;
    task.state = state;
    task.workspaces = workspaces;
    task.elements = &elements;
    task.failures = &failures;
    task.geometry_only = true;
    return Dispatch(&task, failure, d);
  } catch (const std::bad_alloc&) {
    return Allocation(d);
  }
}

StatusCode RecoverNonlinearResults(const AnalysisModel* model,
                                   const DofManager* dofs,
                                   const AnalysisState* state,
                                   const ElementWorkspace* const* workspaces,
                                   int workspace_count,
                                   NonlinearFailure* failure, Diagnostic* d) {
  if (!failure) return Invalid("Nonlinear recovery failure output", d);
  *failure = NonlinearFailure::kNone;
  auto status = Validate(model, dofs, state, workspaces, workspace_count, d);
  if (status != StatusCode::kOk) return status;
  if (!state->CandidateDisplacement() ||
      state->Results() != ResultState::kUnconfirmed)
    return Invalid("Nonlinear recovery requires projected unconfirmed state",
                   d);
  try {
    std::unique_ptr<Vector> values;
    for (const auto& entry : *model->Nodes()) {
      const int components = dofs->NodeDofCount(entry.first);
      status = MakeVector(2 * components, &values, d);
      if (status != StatusCode::kOk) return status;
      if (d) {
        d->node_id = entry.first;
        d->element_id = 0;
      }
      int indices[6] = {};
      for (int c = 0; c < components; ++c)
        indices[c] = dofs->GlobalDof(entry.first, c + 1);
      status =
          values->GatherFrom(state->CandidateDisplacement(), indices, components, 0, d);
      if (status == StatusCode::kOk)
        status =
            values->GatherFrom(state->CandidateReaction(), indices, components, components, d);
      if (status != StatusCode::kOk) return status;
      status = entry.second->UpdateResults(values.get(), d);
      if (status != StatusCode::kOk) return status;
    }
    std::vector<const Element*> elements;
    std::vector<int> offsets;
    CollectElements(model, &elements, &offsets);
    std::vector<TrialEvaluationFailure> failures(elements.size(),
                                                 TrialEvaluationFailure::kNone);
    TaskContext task;
    task.dofs = dofs;
    task.state = state;
    task.workspaces = workspaces;
    task.elements = &elements;
    task.offsets = &offsets;
    task.failures = &failures;
    task.recovery_model = model;
    return Dispatch(&task, failure, d);
  } catch (const std::bad_alloc&) {
    return Allocation(d);
  }
}
}  // namespace fesa
