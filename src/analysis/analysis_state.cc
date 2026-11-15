#include "fesa/analysis/analysis_state.h"

#include <cassert>
#include <cmath>
#include <limits>
#include <memory>
#include <new>
#include <utility>

#include "fesa/analysis/beam_internal_state.h"
#include "fesa/analysis/dof_manager.h"
#include "fesa/analysis/incremental_solid_history.h"
#include "fesa/elements/beam3d.h"
#include "fesa/elements/incremental_solid_element.h"
#include "fesa/analysis/nonlinear_field_history.h"
#include "fesa/analysis/rotation_state.h"
#include "fesa/model/domain.h"
#include "model/model_storage.h"

namespace fesa {
const BeamInternalState* AnalysisState::BeamInternals() const {
  return beam_internals_;
}
IncrementalSolidHistory* AnalysisState::SolidHistory() const {
  return solid_history_;
}
namespace {
// Validation only; numerical copying/arithmetic remains in the math adapters.
bool FiniteStateInput(const Vector* input) {
  if (!input) return false;
  for (int i = 0; i < input->Size(); ++i)
    if (!std::isfinite(input->Data()[i])) return false;
  return true;
}
}  // namespace
AnalysisState::~AnalysisState() {
  delete beam_internals_;
  delete solid_history_;
  delete dynamics_;
  delete[] node_ids_;
  delete[] node_offsets_;
  delete field_history_;
  delete rotations_;
  delete committed_u_;
  delete candidate_u_;
  delete committed_reaction_;
  delete candidate_reaction_;
  delete zero_;
}
StatusCode AnalysisState::Create(const DofManager* dofs, AnalysisState** out,
                                 Diagnostic* d) {
  if (!dofs || dofs->GlobalCount() <= 0) {
    if (out) *out = nullptr;
    return model_storage::Invalid("AnalysisState mapping", d);
  }
  bool has_rotation = false;
  for (int row = 0; row < dofs->NodeCount(); ++row)
    if (dofs->NodeDofCount(dofs->NodeId(row)) == 6) has_rotation = true;
  return CreateLayout(dofs->GlobalCount(), has_rotation, out, d);
}
StatusCode AnalysisState::Create(int node_count, AnalysisState** out,
                                 Diagnostic* d) {
  if (node_count <= 0 || node_count > std::numeric_limits<int>::max() / 6) {
    if (out) *out = nullptr;
    return model_storage::Invalid("AnalysisState::Create", d);
  }
  return CreateLayout(6 * node_count, true, out, d);
}
StatusCode AnalysisState::CreateForDomain(const Domain* domain,
                                          AnalysisState** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("AnalysisState output", d);
  *out = nullptr;
  if (!domain || domain->Nodes()->empty())
    return model_storage::Invalid("AnalysisState domain layout", d);
  int global_count = 0;
  bool has_rotation = false;
  for (const auto& entry : *domain->Nodes()) {
    const int count = entry.second->dof_count();
    if (count == 0)
      return model_storage::Invalid("Thermal node in mechanical state", d);
    if ((count != 3 && count != 6) ||
        global_count > std::numeric_limits<int>::max() - count)
      return model_storage::Invalid("AnalysisState domain layout", d);
    global_count += count;
    has_rotation = has_rotation || count == 6;
  }
  AnalysisState* raw = nullptr;
  auto status = CreateLayout(global_count, has_rotation, &raw, d);
  std::unique_ptr<AnalysisState> state(raw);
  if (status != StatusCode::kOk) return status;
  if (domain->Nodes()->size() >= std::numeric_limits<int>::max())
    return model_storage::Invalid("AnalysisState node count", d);
  state->node_count_ = static_cast<int>(domain->Nodes()->size());
  state->node_ids_ = new (std::nothrow) int[state->node_count_];
  state->node_offsets_ = new (std::nothrow) int[state->node_count_ + 1];
  if (!state->node_ids_ || !state->node_offsets_)
    return model_storage::AllocationFailure("AnalysisState layout", d);
  int row = 0, offset = 0;
  for (const auto& entry : *domain->Nodes()) {
    state->node_ids_[row] = entry.first;
    state->node_offsets_[row++] = offset;
    offset += entry.second->dof_count();
  }
  state->node_offsets_[row] = offset;
  *out = state.release();
  return StatusCode::kOk;
}
StatusCode AnalysisState::CreateLayout(int global_count, bool has_rotation,
                                       AnalysisState** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("AnalysisState::Create", d);
  *out = nullptr;
  std::unique_ptr<AnalysisState> result(new (std::nothrow) AnalysisState);
  if (!result)
    return model_storage::AllocationFailure("AnalysisState::Create", d);
  Vector** stores[] = {&result->committed_u_, &result->candidate_u_,
                       &result->committed_reaction_,
                       &result->candidate_reaction_, &result->zero_};
  for (auto* store : stores) {
    const StatusCode status = Vector::Create(global_count, store, d);
    if (status != StatusCode::kOk) return status;
  }
  result->node_count_ = global_count / (has_rotation ? 6 : 3);
  result->has_rotation_ = has_rotation;
  if (!has_rotation) {
    result->committed_raw_[1] = 0;
    result->candidate_raw_[1] = 0;
  }
  *out = result.release();
  return StatusCode::kOk;
}
StatusCode AnalysisState::PrepareCandidate(Diagnostic* d) {
  DiscardCandidate();
  StatusCode status = candidate_u_->CopyFrom(committed_u_, d);
  if (status != StatusCode::kOk) return status;
  status = candidate_reaction_->CopyFrom(zero_, d);
  if (status != StatusCode::kOk) return status;
  for (int field = 0; field < 2; ++field) {
    candidate_raw_[field] = committed_raw_[field];
    candidate_max_[field] = committed_max_[field];
  }
  if (rotations_) {
    status = rotations_->PrepareCandidate(d);
    if (status != StatusCode::kOk) return status;
  }
  if (solid_history_) {
    status = solid_history_->PrepareCandidate(d);
    if (status != StatusCode::kOk) return status;
  }
  if (beam_internals_) {
    status = beam_internals_->PrepareCandidate(d);
    if (status != StatusCode::kOk) {
      DiscardCandidate();
      return status;
    }
  }
  candidate_ready_ = true;
  return StatusCode::kOk;
}
StatusCode AnalysisState::PrepareNonlinearCandidate(
    int step, const ControlSettings* controls, Diagnostic* d) {
  if (!field_history_ || dynamics_) {
    DiscardCandidate();
    return model_storage::Invalid(
        "Static nonlinear history requires static state", d);
  }
  if (field_history_->CommittedStepNumber() != StepNumber()) {
    DiscardCandidate();
    return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                         "Incomplete committed nonlinear field history", d);
  }
  auto status = PrepareCandidate(d);
  if (status == StatusCode::kOk)
    status = field_history_->PrepareCandidate(step, controls, d);
  if (status != StatusCode::kOk) {
    DiscardCandidate();
    return status;
  }
  nonlinear_history_step_ = step;
  nonlinear_history_required_ = true;
  return StatusCode::kOk;
}
const NonlinearFieldHistory* AnalysisState::FieldHistory() const {
  return field_history_;
}
StatusCode AnalysisState::PreviewNonlinearFluxHistory(
    double force_mean, double moment_mean, double force_max, double moment_max,
    double* force_raw, double* moment_raw, Diagnostic* d) {
  if (!candidate_ready_ || !nonlinear_history_required_ || !field_history_ ||
      !force_raw || !moment_raw || force_raw == moment_raw) {
    DiscardCandidate();
    return model_storage::Invalid("Nonlinear field preview state/outputs", d);
  }
  double raw_force = 0, raw_moment = 0;
  auto status = field_history_->Preview(force_mean, moment_mean, force_max,
                                        moment_max, &raw_force, &raw_moment, d);
  if (status == StatusCode::kOk)
    status = StageFluxHistory(raw_force, raw_moment, force_max, moment_max, d);
  if (status != StatusCode::kOk) {
    DiscardCandidate();
    return status;
  }
  *force_raw = raw_force;
  *moment_raw = raw_moment;
  return StatusCode::kOk;
}

StatusCode AnalysisState::EnableTranslationalDynamics(Diagnostic* d) {
  if (candidate_ready_)
    return model_storage::Invalid("Dynamic storage requires inactive candidate",
                                  d);
  if (has_rotation_ || field_history_)
    return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                         "Dynamic storage currently requires linear 3N state",
                         d);
  return EnableLinearDynamics(d);
}
StatusCode AnalysisState::EnableLinearDynamics(Diagnostic* d) {
  if (candidate_ready_)
    return model_storage::Invalid("Dynamic storage requires inactive candidate",
                                  d);
  if (field_history_)
    return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                         "Linear motion storage requires linear state", d);
  if (dynamics_) return StatusCode::kOk;
  return DynamicState::Create(GlobalCount(), &dynamics_, d);
}
StatusCode AnalysisState::EnableFiniteDynamics(Diagnostic* d) {
  if (candidate_ready_)
    return model_storage::Invalid("Dynamic storage requires inactive candidate",
                                  d);
  if (!field_history_ || !rotations_ || !has_rotation_ ||
      (GlobalCount() % 6 != 0 || GlobalCount() / 6 != rotations_->NodeCount()))
    return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                         "Finite dynamics requires nonlinear all-sixDOF state",
                         d);
  if (dynamics_) return StatusCode::kOk;
  if (beam_internals_) {
    const auto status = beam_internals_->EnableDynamics(d);
    if (status != StatusCode::kOk) return status;
  }
  return DynamicState::Create(GlobalCount(), &dynamics_, d);
}
StatusCode AnalysisState::EnableDynamicConvergenceHistory(Diagnostic* d) {
  if (!dynamics_ || candidate_ready_ ||
      (!field_history_ && step_number_ != 0 && increment_number_ != 1) ||
      dynamics_->CommittedVelocity())
    return model_storage::Invalid(
        "Dynamic history requires unused motion state", d);
  if (dynamics_->convergence_history_) return StatusCode::kOk;
  if (field_history_) {
    if (field_history_->CommittedStepNumber() != step_number_)
      return model_storage::Invalid("Finite dynamic history missing preload",
                                    d);
    dynamics_->convergence_history_ = field_history_;
    dynamics_->owns_history_ = false;
    return StatusCode::kOk;
  }
  return NonlinearFieldHistory::Create(has_rotation_,
                                       &dynamics_->convergence_history_, d);
}
StatusCode AnalysisState::PrepareDynamicConvergenceHistory(
    int step, const ControlSettings* controls, bool initial,
    double initial_maximum, Diagnostic* d) {
  if (has_rotation_) {
    DiscardCandidate();
    return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                         "Force-only dynamic history on rotational state", d);
  }
  return PrepareDynamicConvergenceHistory(step, controls, initial,
                                          initial_maximum, 0, d);
}
StatusCode AnalysisState::PrepareDynamicConvergenceHistory(
    int step, const ControlSettings* controls, bool initial,
    double initial_force_maximum, double initial_moment_maximum,
    Diagnostic* d) {
  if (!candidate_ready_ || !dynamics_ || !dynamics_->CandidateVelocity() ||
      !dynamics_->convergence_history_ ||
      (beam_internals_ && !beam_internals_->CanCommit(true))) {
    DiscardCandidate();
    return model_storage::Invalid("Dynamic history requires staged motion", d);
  }
  auto* history = dynamics_->convergence_history_;
  const bool finite_static_handoff =
      initial && field_history_ && step_number_ > 0 &&
      !dynamics_->CommittedVelocity() &&
      history->CommittedStepNumber() == step_number_;
  const bool static_handoff = initial && !field_history_ && step_number_ > 0 &&
                              increment_number_ == 1 &&
                              !dynamics_->CommittedVelocity() &&
                              history->CommittedStepNumber() == 0;
  const bool valid_initial =
      initial && step_number_ < std::numeric_limits<int>::max() &&
      step == step_number_ + 1 && dynamics_->CandidateStepTime() == 0 &&
      (static_handoff || finite_static_handoff ||
       dynamics_->CandidateTotalTime() == dynamics_->TotalTime());
  const bool valid_end =
      !initial && step > 0 && step == step_number_ &&
      dynamics_->CommittedVelocity() &&
      dynamics_->CandidateStepTime() > dynamics_->StepTime() &&
      dynamics_->CandidateTotalTime() > dynamics_->TotalTime();
  if ((!static_handoff && history->CommittedStepNumber() != step_number_) ||
      (!valid_initial && !valid_end)) {
    DiscardCandidate();
    return model_storage::Invalid("Dynamic history step/time mismatch", d);
  }
  const double maxima[2] = {initial_force_maximum, initial_moment_maximum};
  const auto status =
      static_handoff ? history->PrepareDynamicInitialAfterLinearStatic(
                           step_number_, controls, committed_raw_,
                           committed_max_, maxima, d)
      : initial      ? history->PrepareDynamicInitial(step, controls,
                                                      initial_force_maximum,
                                                      initial_moment_maximum, d)
                     : history->PrepareCandidate(step, controls, d);
  if (status != StatusCode::kOk) {
    DiscardCandidate();
    return status;
  }
  dynamics_->history_step_ = step;
  dynamics_->history_initial_ = initial;
  if (initial) {
    candidate_max_[0] = initial_force_maximum;
    candidate_max_[1] = has_rotation_ ? initial_moment_maximum : 0;
  }
  return StatusCode::kOk;
}
StatusCode AnalysisState::PreviewDynamicFluxHistory(double force_mean,
                                                    double force_max,
                                                    double* force_raw,
                                                    Diagnostic* d) {
  if (has_rotation_) {
    DiscardCandidate();
    return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                         "Force-only dynamic preview on rotational state", d);
  }
  double inactive = 0;
  return PreviewDynamicFluxHistory(force_mean, 0, force_max, 0, force_raw,
                                   &inactive, d);
}
StatusCode AnalysisState::PreviewDynamicFluxHistory(
    double force_mean, double moment_mean, double force_max, double moment_max,
    double* force_raw, double* moment_raw, Diagnostic* d) {
  if (!candidate_ready_ || !dynamics_ || !dynamics_->convergence_history_ ||
      dynamics_->history_step_ == 0 || dynamics_->history_initial_ ||
      !force_raw || !moment_raw || force_raw == moment_raw) {
    DiscardCandidate();
    return model_storage::Invalid("Dynamic field preview state/output", d);
  }
  double raw = 0, raw_moment = 0;
  auto status = dynamics_->convergence_history_->Preview(
      force_mean, moment_mean, force_max, moment_max, &raw, &raw_moment, d);
  if (status == StatusCode::kOk)
    status = StageFluxHistory(raw, raw_moment, force_max, moment_max, d);
  if (status != StatusCode::kOk) {
    DiscardCandidate();
    return status;
  }
  *force_raw = raw;
  *moment_raw = raw_moment;
  return StatusCode::kOk;
}
const DynamicState* AnalysisState::Dynamics() const { return dynamics_; }
StatusCode AnalysisState::StageDynamicKinematics(const Vector* velocity,
                                                 const Vector* acceleration,
                                                 double step_time,
                                                 double total_time,
                                                 Diagnostic* d) {
  if (beam_internals_) beam_internals_->DiscardMotion();
  if (!candidate_ready_ || !dynamics_) {
    DiscardCandidate();
    return model_storage::Invalid("Dynamic staging requires prepared state", d);
  }
  const auto status =
      dynamics_->Stage(velocity, acceleration, step_time, total_time, d);
  if (status != StatusCode::kOk) DiscardCandidate();
  return status;
}

StatusCode AnalysisState::StageBeamInternalKinematics(
    const Vector* velocity, const Vector* acceleration, Diagnostic* d) {
  if (!candidate_ready_ || !beam_internals_ || !dynamics_ ||
      !dynamics_->CandidateVelocity()) {
    DiscardCandidate();
    return model_storage::Invalid("Beam motion requires staged physical state",
                                  d);
  }
  // Internal restaging retains physical motion but its R11 history must be
  // prepared again for the new complete kinematic candidate.
  dynamics_->history_step_ = 0;
  dynamics_->history_initial_ = false;
  if (dynamics_->convergence_history_)
    dynamics_->convergence_history_->DiscardCandidate();
  const auto status = beam_internals_->StageMotion(velocity, acceleration, d);
  if (status != StatusCode::kOk) DiscardCandidate();
  return status;
}

void AnalysisState::MarkResultsUnconfirmed() {
  assert(step_number_ < std::numeric_limits<int>::max());
  MarkResultsUnconfirmed(step_number_ + 1, 1);
}
void AnalysisState::MarkResultsUnconfirmed(int step_number,
                                           int increment_number) {
  assert(!beam_internals_ || beam_internals_->CanCommit(dynamics_ != nullptr));
  assert(candidate_ready_ && step_number > 0 &&
         (increment_number > 0 || (increment_number == 0 && dynamics_ &&
                                   dynamics_->CandidateVelocity())));
  assert(!nonlinear_history_required_ ||
         step_number == nonlinear_history_step_);
  assert(!dynamics_ ||
         dynamics_->HistoryMatches(step_number, increment_number));
  result_step_number_ = step_number;
  result_increment_number_ = increment_number;
  result_state_ = ResultState::kUnconfirmed;
}
void AnalysisState::DiscardCandidate() noexcept {
  if (beam_internals_) beam_internals_->DiscardCandidate();
  if (solid_history_) solid_history_->DiscardCandidate();
  if (dynamics_) dynamics_->Discard();
  candidate_ready_ = false;
  nonlinear_output_ready_ = false;
  if (rotations_) rotations_->DiscardCandidate();
  if (field_history_) field_history_->DiscardCandidate();
  nonlinear_history_required_ = false;
  nonlinear_history_step_ = 0;
}
void AnalysisState::Commit(int step_number, int increment_number) noexcept {
  assert(!beam_internals_ || beam_internals_->CanCommit(dynamics_ != nullptr));
  assert(!solid_history_ || solid_history_->CanCommit());
  assert(candidate_ready_ && result_state_ == ResultState::kUnconfirmed &&
         step_number > 0 &&
         (increment_number > 0 || (increment_number == 0 && dynamics_)) &&
         result_step_number_ == step_number &&
         result_increment_number_ == increment_number);
  assert(!field_history_ || nonlinear_output_ready_);
  assert(!nonlinear_history_required_ ||
         (field_history_ && field_history_->CanCommit() &&
          step_number == nonlinear_history_step_));
  assert(!dynamics_ || dynamics_->CandidateVelocity());
  assert(!dynamics_ ||
         dynamics_->HistoryMatches(step_number, increment_number));
  if (dynamics_) {
    dynamics_->Commit();
    if (dynamics_->ConvergenceHistory()) {
      candidate_raw_[0] = dynamics_->ConvergenceHistory()->CommittedRawMean(0);
      candidate_raw_[1] = dynamics_->ConvergenceHistory()->CommittedRawMean(1);
    }
  }
  if (rotations_) rotations_->Commit();
  if (nonlinear_history_required_) field_history_->Commit();
  if (solid_history_) solid_history_->CommitCandidate();
  if (beam_internals_) beam_internals_->CommitCandidate(dynamics_ != nullptr);
  std::swap(candidate_u_, committed_u_);
  std::swap(candidate_reaction_, committed_reaction_);
  for (int field = 0; field < 2; ++field) {
    committed_raw_[field] = candidate_raw_[field];
    committed_max_[field] = candidate_max_[field];
  }
  result_step_number_ = step_number;
  result_increment_number_ = increment_number;
  step_number_ = step_number;
  increment_number_ = increment_number;
  result_state_ = ResultState::kConfirmed;
  candidate_ready_ = false;
  nonlinear_history_required_ = false;
  nonlinear_history_step_ = 0;
}
Vector* AnalysisState::CandidateDisplacement() {
  return candidate_ready_ && !field_history_ ? candidate_u_ : nullptr;
}
Vector* AnalysisState::CandidateReaction() {
  return candidate_ready_ ? candidate_reaction_ : nullptr;
}
const Vector* AnalysisState::CandidateDisplacement() const {
  return candidate_ready_ && (!field_history_ || nonlinear_output_ready_)
             ? candidate_u_
             : nullptr;
}
const Vector* AnalysisState::CandidateReaction() const {
  return candidate_ready_ ? candidate_reaction_ : nullptr;
}
const Vector* AnalysisState::CommittedDisplacement() const {
  return committed_u_;
}
const Vector* AnalysisState::CommittedReaction() const {
  return committed_reaction_;
}
StatusCode AnalysisState::StageFluxHistory(double force_raw, double moment_raw,
                                           double force_max, double moment_max,
                                           Diagnostic* d) {
  if (!candidate_ready_ || !std::isfinite(force_raw) || force_raw <= 0 ||
      !std::isfinite(force_max) || force_max < 0 ||
      (has_rotation_ && (!std::isfinite(moment_raw) || moment_raw <= 0 ||
                         !std::isfinite(moment_max) || moment_max < 0)))
    return model_storage::Invalid("AnalysisState::StageFluxHistory", d);
  candidate_raw_[0] = force_raw;
  candidate_raw_[1] = has_rotation_ ? moment_raw : 0;
  candidate_max_[0] = force_max;
  candidate_max_[1] = has_rotation_ ? moment_max : 0;
  return StatusCode::kOk;
}

StatusCode AnalysisState::CreateNonlinear(int node_count, AnalysisState** out,
                                          Diagnostic* d) {
  if (!out) return model_storage::Invalid("Nonlinear state output", d);
  *out = nullptr;
  AnalysisState* raw = nullptr;
  auto status = Create(node_count, &raw, d);
  std::unique_ptr<AnalysisState> state(raw);
  if (status != StatusCode::kOk) return status;
  status = RotationState::Create(node_count, &state->rotations_, d);
  if (status != StatusCode::kOk) return status;
  status = NonlinearFieldHistory::Create(&state->field_history_, d);
  if (status != StatusCode::kOk) return status;
  *out = state.release();
  return StatusCode::kOk;
}
StatusCode AnalysisState::CreateNonlinearForDomain(
    const Domain* domain, AnalysisState** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("Nonlinear state output", d);
  *out = nullptr;
  AnalysisState* raw = nullptr;
  auto status = CreateForDomain(domain, &raw, d);
  std::unique_ptr<AnalysisState> state(raw);
  if (status != StatusCode::kOk) return status;
  if (state->has_rotation_) {
    status = RotationState::Create(state->node_count_, &state->rotations_, d);
    if (status != StatusCode::kOk) return status;
  }
  status = NonlinearFieldHistory::Create(state->has_rotation_,
                                        &state->field_history_, d);
  if (status != StatusCode::kOk) return status;
  for (const auto& entry : *domain->Elements()) {
    if (!dynamic_cast<const IncrementalSolidElement*>(entry.second)) continue;
    status = IncrementalSolidHistory::Create(domain, &state->solid_history_, d);
    if (status != StatusCode::kOk) return status;
    break;
  }
  for (const auto& entry : *domain->Elements()) {
    if (!dynamic_cast<const Beam3D*>(entry.second)) continue;
    status = BeamInternalState::Create(domain, &state->beam_internals_, d);
    if (status != StatusCode::kOk) return status;
    break;
  }
  *out = state.release();
  return StatusCode::kOk;
}
bool AnalysisState::HasValidCandidate() const { return candidate_ready_; }
StatusCode AnalysisState::ApplyBeamInternalCorrection(const Vector* correction,
                                                      Diagnostic* d) {
  nonlinear_output_ready_ = false;
  if (dynamics_) dynamics_->Discard();
  if (!candidate_ready_ || !beam_internals_) {
    DiscardCandidate();
    return model_storage::Invalid("Beam correction requires prepared state", d);
  }
  const auto status = beam_internals_->ApplyCorrection(correction, d);
  if (status != StatusCode::kOk) DiscardCandidate();
  return status;
}
StatusCode AnalysisState::ValidateMapping(const DofManager* dofs,
                                          Diagnostic* d) const {
  if (!dofs || !node_ids_ || !node_offsets_ ||
      dofs->NodeCount() != node_count_ || dofs->GlobalCount() != GlobalCount())
    return model_storage::Invalid("AnalysisState physical mapping", d);
  for (int row = 0; row < node_count_; ++row) {
    const int id = node_ids_[row];
    if (dofs->NodeId(row) != id ||
        dofs->NodeDofCount(id) != node_offsets_[row + 1] - node_offsets_[row] ||
        dofs->GlobalDof(id, 1) != node_offsets_[row])
      return model_storage::Invalid("AnalysisState physical mapping", d);
  }
  return StatusCode::kOk;
}
int AnalysisState::PhysicalOffset(int node) const {
  return node_offsets_ ? node_offsets_[node] : 6 * node;
}
const RotationState* AnalysisState::Rotations() const { return rotations_; }
StatusCode AnalysisState::GatherCandidateTranslation(int node, Vector* out,
                                                     Diagnostic* d) const {
  if (!candidate_ready_ || !field_history_ || node < 0 ||
      node >= node_count_ || !out || out->Size() != 3)
    return model_storage::Invalid("Nonlinear translation gather", d);
  const int offset = PhysicalOffset(node);
  const int indices[3] = {offset, offset + 1, offset + 2};
  return out->GatherFrom(candidate_u_, indices, 3, 0, d);
}
StatusCode AnalysisState::SetCandidateTranslation(int node,
                                                  const Vector* translation,
                                                  Diagnostic* d) {
  nonlinear_output_ready_ = false;
  if (dynamics_) dynamics_->Discard();
  if (beam_internals_) beam_internals_->DiscardMotion();
  if (!candidate_ready_ || !field_history_ || node < 0 ||
      node >= node_count_ || !translation ||
      translation->Size() != 3 || !FiniteStateInput(translation)) {
    DiscardCandidate();
    return model_storage::Invalid("Nonlinear translation target", d);
  }
  const int components[3] = {0, 1, 2};
  auto status = solid_history_ ? solid_history_->PrepareCandidate(d)
                               : StatusCode::kOk;
  if (status == StatusCode::kOk)
    status = candidate_u_->GatherFrom(translation, components, 3,
                                       PhysicalOffset(node), d);
  if (status != StatusCode::kOk) DiscardCandidate();
  return status;
}
StatusCode AnalysisState::SetCandidateOrientation(int node, const Vector* q,
                                                  bool* branch_exceeded,
                                                  Diagnostic* d) {
  if (branch_exceeded) *branch_exceeded = false;
  nonlinear_output_ready_ = false;
  if (dynamics_) dynamics_->Discard();
  if (beam_internals_) beam_internals_->DiscardMotion();
  if (!candidate_ready_ || !rotations_ || node < 0 || node >= node_count_ ||
      PhysicalOffset(node + 1) - PhysicalOffset(node) != 6) {
    DiscardCandidate();
    return model_storage::Invalid("Nonlinear orientation target", d);
  }
  const auto status =
      rotations_->SetCandidateOrientation(node, q, branch_exceeded, d);
  if (status != StatusCode::kOk) DiscardCandidate();
  return status;
}
StatusCode AnalysisState::ApplyNonlinearCorrection(const Vector* correction,
                                                   bool* branch_exceeded,
                                                   Diagnostic* d) {
  if (branch_exceeded) *branch_exceeded = false;
  nonlinear_output_ready_ = false;
  if (dynamics_) dynamics_->Discard();
  if (beam_internals_) beam_internals_->DiscardMotion();
  if (!candidate_ready_ || !field_history_ || !branch_exceeded || !correction ||
      correction->Size() != GlobalCount() || correction == candidate_u_ ||
      correction == committed_u_ || correction == candidate_reaction_ ||
      correction == committed_reaction_ || correction == zero_ ||
      !FiniteStateInput(correction)) {
    DiscardCandidate();
    return model_storage::Invalid("Nonlinear correction inputs", d);
  }
  std::unique_ptr<Vector> values[3];
  for (auto& value : values) {
    Vector* raw = nullptr;
    const auto status = Vector::Create(3, &raw, d);
    value.reset(raw);
    if (status != StatusCode::kOk) {
      DiscardCandidate();
      return status;
    }
  }
  StatusCode status = solid_history_ ? solid_history_->PrepareCandidate(d)
                                     : StatusCode::kOk;
  for (int node = 0; node < node_count_ && status == StatusCode::kOk; ++node) {
    const int offset = PhysicalOffset(node);
    const bool rotates = PhysicalOffset(node + 1) - offset == 6;
    const int translations[3] = {offset, offset + 1, offset + 2};
    const int rotations[3] = {offset + 3, offset + 4, offset + 5};
    status = values[0]->GatherFrom(correction, translations, 3, 0, d);
    if (status != StatusCode::kOk) break;
    if (rotates)
      status = values[1]->GatherFrom(correction, rotations, 3, 0, d);
    if (status != StatusCode::kOk) break;
    status = values[2]->GatherFrom(candidate_u_, translations, 3, 0, d);
    if (status != StatusCode::kOk) break;
    status = values[2]->AddScaled(1, values[0].get(), d);
    if (status != StatusCode::kOk) break;
    if (rotates)
      status = rotations_->ApplyCorrection(node, values[1].get(),
                                          branch_exceeded, d);
    if (status != StatusCode::kOk) break;
    const int components[3] = {0, 1, 2};
    status = candidate_u_->GatherFrom(values[2].get(), components, 3,
                                       offset, d);
    if (status != StatusCode::kOk) break;
  }
  if (status != StatusCode::kOk) DiscardCandidate();
  return status;
}
StatusCode AnalysisState::PrepareNonlinearOutput(Diagnostic* d) {
  nonlinear_output_ready_ = false;
  if (!candidate_ready_ || !field_history_) {
    DiscardCandidate();
    return model_storage::Invalid("Nonlinear output preparation", d);
  }
  if (!rotations_) {
    nonlinear_output_ready_ = true;
    return StatusCode::kOk;
  }
  Vector* raw = nullptr;
  auto status = Vector::Create(3 * rotations_->NodeCount(), &raw, d);
  std::unique_ptr<Vector> principal(raw);
  if (status == StatusCode::kOk)
    status = rotations_->PrincipalOutput(principal.get(), d);
  if (status == StatusCode::kOk) {
    for (int node = 0; node < rotations_->NodeCount(); ++node) {
      if (PhysicalOffset(node + 1) - PhysicalOffset(node) != 6) continue;
      const int indices[3] = {3 * node, 3 * node + 1, 3 * node + 2};
      status = candidate_u_->GatherFrom(principal.get(), indices, 3,
                                        PhysicalOffset(node) + 3, d);
      if (status != StatusCode::kOk) break;
    }
  }
  if (status != StatusCode::kOk)
    DiscardCandidate();
  else
    nonlinear_output_ready_ = true;
  return status;
}

double AnalysisState::CommittedRawMean(int field) const {
  assert(field == 0 || field == 1);
  return committed_raw_[field];
}
double AnalysisState::CommittedMaximumFlux(int field) const {
  assert(field == 0 || field == 1);
  return committed_max_[field];
}
int AnalysisState::GlobalCount() const { return committed_u_->Size(); }
ResultState AnalysisState::Results() const { return result_state_; }
int AnalysisState::ResultStepNumber() const { return result_step_number_; }
int AnalysisState::ResultIncrementNumber() const {
  return result_increment_number_;
}
int AnalysisState::StepNumber() const { return step_number_; }
int AnalysisState::IncrementNumber() const { return increment_number_; }
}  // namespace fesa
