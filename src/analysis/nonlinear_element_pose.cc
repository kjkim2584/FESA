#include "analysis/nonlinear_element_pose.h"

#include <climits>
#include <memory>
#include <new>

#include "fesa/analysis/rotation_state.h"
#include "fesa/analysis/beam_internal_state.h"
#include "fesa/elements/beam3d.h"
#include "fesa/analysis/incremental_solid_history.h"
#include "fesa/elements/incremental_solid_element.h"
#include "model/model_storage.h"

namespace fesa {
NonlinearElementPose::~NonlinearElementPose() {
  delete trial_;
  delete[] global_dofs_;
  delete translations_;
  delete internal_position_;
  delete spin_;
  delete orientations_;
  delete committed_;
  delete committed_translations_;
  delete candidate_stress_;
  delete incremental_;
}
StatusCode NonlinearElementPose::Create(
    const Element* element, const DofManager* dofs, const AnalysisState* state,
    bool committed, NonlinearElementPose** out, Diagnostic* d,
    IncrementalElasticWorkspace* scratch,
    bool use_committed_geometric_stress) {
  if (!out) return model_storage::Invalid("Nonlinear element pose output", d);
  *out = nullptr;
  std::unique_ptr<NonlinearElementPose> pose(new (std::nothrow)
                                                 NonlinearElementPose);
  if (!pose)
    return model_storage::AllocationFailure("Nonlinear element pose", d);
  const auto status = pose->Initialize(element, dofs, state, committed, d, scratch,
                                       use_committed_geometric_stress);
  if (status != StatusCode::kOk) return status;
  *out = pose.release();
  return StatusCode::kOk;
}
StatusCode NonlinearElementPose::Initialize(const Element* element,
                                            const DofManager* dofs,
                                            const AnalysisState* state,
                                            bool committed, Diagnostic* d,
                                            IncrementalElasticWorkspace* scratch,
                                            bool use_committed_geometric_stress) {
  if (!element || !dofs || !state || !state->FieldHistory() ||
      element->node_count() <= 0 || element->node_count() > INT_MAX / 6 ||
      (element->dofs_per_node() != 3 && element->dofs_per_node() != 6) ||
      element->dof_count() !=
          element->node_count() * element->dofs_per_node() ||
      state->GlobalCount() != dofs->GlobalCount() || dofs->NodeCount() <= 0 ||
      (!committed && !state->HasValidCandidate()))
    return model_storage::Invalid("Nonlinear element pose inputs", d);
  const int count = element->node_count(),
            components = element->dofs_per_node();
  const auto* beam = dynamic_cast<const Beam3D*>(element);
  const auto* internal = state->BeamInternals();
  if (element->internal_dof_count() != 0 || dofs->InternalDofCount() != 0) {
    if (!beam || element->internal_dof_count() != 2 || !internal ||
        dofs->InternalDofCount() != 2 * internal->ElementCount() ||
        dofs->GlobalCount() != 6 * dofs->NodeCount() ||
        dofs->InternalDof(beam->Id(), 0) != internal->Offset(beam->Id()) ||
        internal->Offset(beam->Id()) < 0 ||
        dofs->InternalDof(beam->Id(), 1) != internal->Offset(beam->Id()) + 1)
      return model_storage::Invalid("Nonlinear Beam internal pose mapping", d);
    for (int row = 0; row < internal->ElementCount(); ++row)
      if (dofs->InternalElementId(row) != internal->ElementId(row))
        return model_storage::Invalid("Nonlinear Beam internal element IDs", d);
    const auto* position = committed ? internal->CommittedPosition()
                                     : internal->CandidatePosition();
    if (!position || position->Size() != dofs->InternalDofCount())
      return model_storage::Invalid("Nonlinear Beam internal pose source", d);
    auto status = state->ValidateMapping(dofs, d);
    if (status != StatusCode::kOk) return status;
    status = Vector::Create(2, &internal_position_, d);
    if (status != StatusCode::kOk) return status;
    const int indices[] = {internal->Offset(beam->Id()),
                           internal->Offset(beam->Id()) + 1};
    status = internal_position_->GatherFrom(position, indices, 2, 0, d);
    if (status != StatusCode::kOk) return status;
  }
  const bool rotations = components == 6;
  const auto* rotation = state->Rotations();
  if (rotations && (!rotation || rotation->NodeCount() != dofs->NodeCount() ||
                    (!committed && (!rotation->CandidateOrientations() ||
                                    !rotation->CandidateBodySpin()))))
    return model_storage::Invalid("Nonlinear element pose rotation state", d);
  if (!rotation) {
    const auto status = state->ValidateMapping(dofs, d);
    if (status != StatusCode::kOk) return status;
  }
  dof_count_ = element->dof_count();
  global_dofs_ = new (std::nothrow) int[dof_count_];
  trial_ = new (std::nothrow) ElementTrialState;
  if (!global_dofs_ || !trial_)
    return model_storage::AllocationFailure("Nonlinear element pose metadata",
                                            d);
  auto status = Vector::Create(3 * count, &translations_, d);
  if (status != StatusCode::kOk) return status;
  Vector* raw_row = nullptr;
  status = Vector::Create(3, &raw_row, d);
  std::unique_ptr<Vector> row_u(raw_row);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<int[]> spin_indices;
  if (rotations) {
    spin_indices.reset(new (std::nothrow) int[3 * count]);
    if (!spin_indices)
      return model_storage::AllocationFailure("Nonlinear element spin indices",
                                              d);
    status = Vector::Create(3 * count, &spin_, d);
    if (status == StatusCode::kOk)
      status = Matrix::Create(count, 4, &orientations_, d);
    if (status == StatusCode::kOk)
      status = Matrix::Create(count, 4, &committed_, d);
    if (status != StatusCode::kOk) return status;
  }
  const int axes[] = {0, 1, 2};
  for (int node = 0; node < count; ++node) {
    const auto* connected = element->GetNode(node);
    if (!connected || dofs->NodeDofCount(connected->Id()) != components)
      return model_storage::Invalid("Nonlinear element pose connectivity", d);
    int row = 0, end = dofs->NodeCount();
    while (row < end) {
      const int middle = row + (end - row) / 2;
      if (dofs->NodeId(middle) < connected->Id())
        row = middle + 1;
      else
        end = middle;
    }
    if (row == dofs->NodeCount() || dofs->NodeId(row) != connected->Id())
      return model_storage::Invalid("Nonlinear element pose node row", d);
    for (int c = 0; c < components; ++c) {
      const int index = dofs->GlobalDof(connected->Id(), c + 1);
      if (index < 0)
        return model_storage::Invalid("Nonlinear element pose DOF", d);
      global_dofs_[components * node + c] = index;
    }
    if (committed) {
      status = translations_->GatherFrom(state->CommittedDisplacement(),
                                         global_dofs_ + components * node, 3,
                                         3 * node, d);
    } else {
      status = state->GatherCandidateTranslation(row, row_u.get(), d);
      if (status == StatusCode::kOk)
        status = translations_->GatherFrom(row_u.get(), axes, 3, 3 * node, d);
    }
    if (status != StatusCode::kOk) return status;
    if (!rotations) continue;
    for (int c = 0; c < 3; ++c) spin_indices[3 * node + c] = 3 * row + c;
    status = orientations_->CopyBlockFrom(
        committed ? rotation->CommittedOrientations()
                  : rotation->CandidateOrientations(),
        row, 0, 1, 4, node, 0, d);
    if (status == StatusCode::kOk)
      status = committed_->CopyBlockFrom(rotation->CommittedOrientations(), row,
                                         0, 1, 4, node, 0, d);
    if (status != StatusCode::kOk) return status;
  }
  if (rotations) {
    status = spin_->GatherFrom(committed ? rotation->CommittedBodySpin()
                                         : rotation->CandidateBodySpin(),
                               spin_indices.get(), 3 * count, 0, d);
    if (status != StatusCode::kOk) return status;
  }
  trial_->translations = translations_;
  trial_->orientations = orientations_;
  trial_->committed_orientations = committed_;
  trial_->body_spin = spin_;
  trial_->beam_internal_position = internal_position_;
  if (dynamic_cast<const IncrementalSolidElement*>(element)) {
    const auto* history = state->SolidHistory();
    const auto* stress = history ? history->CommittedStress(element->Id())
                                 : nullptr;
    if (!stress)
      return model_storage::Invalid("Incremental pose history layout", d);
    incremental_ = new (std::nothrow) IncrementalSolidTrialContext;
    if (!incremental_)
      return model_storage::AllocationFailure("Incremental pose context", d);
    status = Vector::Create(3 * count, &committed_translations_, d);
    if (status == StatusCode::kOk)
      status = Matrix::Create(stress->Rows(), 6, &candidate_stress_, d);
    if (status == StatusCode::kOk)
      status = committed_translations_->GatherFrom(
          state->CommittedDisplacement(), global_dofs_, 3 * count, 0, d);
    if (status != StatusCode::kOk) return status;
    incremental_->workspace = scratch;
    incremental_->use_committed_geometric_stress = use_committed_geometric_stress;
    incremental_->committed_translations = committed_translations_;
    incremental_->committed_stress = stress;
    incremental_->candidate_stress = candidate_stress_;
    trial_->incremental_solid = incremental_;
  }
  return StatusCode::kOk;
}
const ElementTrialState* NonlinearElementPose::Trial() const { return trial_; }
const int* NonlinearElementPose::GlobalDofs() const { return global_dofs_; }
int NonlinearElementPose::DofCount() const { return dof_count_; }
}  // namespace fesa
