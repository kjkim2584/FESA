#include "analysis/finite_dynamic_kinematics.h"

#include <climits>
#include <memory>
#include <new>
#include <vector>

#include "fesa/analysis/beam_internal_state.h"
#include "fesa/analysis/rotation_state.h"
#include "model/model_storage.h"

namespace fesa {
namespace {
enum class Operation { kPredict, kUpdate, kLinearize };
/// Validates both immutable roles before any candidate mutation. All pointers
/// are borrowed; argument failures preserve state and linearization outputs.
StatusCode Validate(const DofManager* dofs, const TimeIntegrator* integrator,
                    Operation operation, const AnalysisState* state,
                    TranslationalLinearization* translation,
                    RotationalLinearization* rotation, bool* branch,
                    Diagnostic* d) {
  if (branch) *branch = false;
  if (!dofs || !integrator || !state || !branch ||
      !state->HasValidCandidate() || !state->Rotations() ||
      !state->Dynamics() || !state->Dynamics()->CommittedVelocity() ||
      !state->Dynamics()->CommittedAcceleration() || dofs->NodeCount() <= 0 ||
      dofs->NodeCount() > INT_MAX / 6 ||
      (operation == Operation::kLinearize &&
       (!translation || !rotation || !state->Dynamics()->CandidateVelocity() ||
        !state->Dynamics()->CandidateAcceleration())))
    return model_storage::Invalid("Finite kinematic state", d);
  auto status = state->ValidateMapping(dofs, d);
  if (status != StatusCode::kOk) return status;
  for (int row = 0; row < dofs->NodeCount(); ++row)
    if (dofs->NodeDofCount(dofs->NodeId(row)) != 6)
      return model_storage::Invalid("Finite kinematic six-DOF mapping", d);
  const auto* internal = state->BeamInternals();
  if ((internal != nullptr) != (dofs->InternalDofCount() > 0))
    return model_storage::Invalid("Finite kinematic Beam role", d);
  if (!internal) return StatusCode::kOk;
  if (internal->ElementCount() <= 0 || internal->ElementCount() > INT_MAX / 2 ||
      dofs->InternalDofCount() != 2 * internal->ElementCount() ||
      !internal->CandidatePosition() || !internal->CommittedVelocity() ||
      !internal->CommittedAcceleration() ||
      (operation == Operation::kLinearize &&
       (!internal->CandidateVelocity() || !internal->CandidateAcceleration())))
    return model_storage::Invalid("Finite kinematic Beam internal motion", d);
  for (int row = 0; row < internal->ElementCount(); ++row) {
    const int id = internal->ElementId(row);
    if (dofs->InternalElementId(row) != id || internal->Offset(id) != 2 * row ||
        dofs->InternalDof(id, 0) != 2 * row ||
        dofs->InternalDof(id, 1) != 2 * row + 1)
      return model_storage::Invalid("Finite kinematic Beam internal mapping",
                                    d);
  }
  return StatusCode::kOk;
}
/// Borrows validated joint state for this call. Predict changes only eta;
/// update stages eta rates after physical staging. Scratch never escapes.
StatusCode EvaluateBeamInternals(const TimeIntegrator* integrator, double dt,
                                 bool predict, AnalysisState* state,
                                 Diagnostic* d) {
  const auto* internal = state->BeamInternals();
  if (!internal) return StatusCode::kOk;
  std::unique_ptr<Vector> trial[3];
  for (auto& buffer : trial) {
    Vector* raw = nullptr;
    const auto status =
        Vector::Create(internal->CommittedPosition()->Size(), &raw, d);
    buffer.reset(raw);
    if (status != StatusCode::kOk) return status;
  }
  TranslationalStateView committed;
  committed.displacement = internal->CommittedPosition();
  committed.velocity = internal->CommittedVelocity();
  committed.acceleration = internal->CommittedAcceleration();
  TranslationalCandidate candidate;
  candidate.displacement = trial[0].get();
  candidate.velocity = trial[1].get();
  candidate.acceleration = trial[2].get();
  if (predict) {
    auto status = integrator->PrepareKinematics(&committed, dt, &candidate, d);
    // Reuse a distinct call-owned vector for predictor minus current eta.
    if (status == StatusCode::kOk)
      status = trial[1]->CopyFrom(candidate.displacement, d);
    if (status == StatusCode::kOk)
      status = trial[1]->AddScaled(-1, internal->CandidatePosition(), d);
    return status == StatusCode::kOk
               ? state->ApplyBeamInternalCorrection(trial[1].get(), d)
               : status;
  }
  auto status = trial[0]->CopyFrom(internal->CandidatePosition(), d);
  if (status == StatusCode::kOk)
    status = integrator->UpdateKinematics(&committed, dt, &candidate, d);
  return status == StatusCode::kOk ? state->StageBeamInternalKinematics(
                                         trial[1].get(), trial[2].get(), d)
                                   : status;
}
StatusCode Evaluate(const DofManager* dofs, const TimeIntegrator* integrator,
                    double dt, double time, double total, Operation operation,
                    const AnalysisState* state, AnalysisState* mutable_state,
                    TranslationalLinearization* translation,
                    RotationalLinearization* rotation, bool* branch,
                    Diagnostic* d) {
  auto status = Validate(dofs, integrator, operation, state, translation,
                         rotation, branch, d);
  if (status != StatusCode::kOk) return status;
  try {
    const int nodes = dofs->NodeCount(), size = 3 * nodes;
    std::vector<int> translations(size), rotations(size);
    for (int row = 0; row < nodes; ++row)
      for (int c = 0; c < 3; ++c) {
        translations[3 * row + c] = dofs->GlobalDof(dofs->NodeId(row), c + 1);
        rotations[3 * row + c] = dofs->GlobalDof(dofs->NodeId(row), c + 4);
      }
    // Five committed/trial roles: Cartesian u/v/a and spatial omega/eta.
    std::unique_ptr<Vector> old[5], trial[5];
    for (int i = 0; i < 5; ++i) {
      Vector* raw = nullptr;
      status = Vector::Create(size, &raw, d);
      old[i].reset(raw);
      if (status != StatusCode::kOk) return status;
      status = Vector::Create(size, &raw, d);
      trial[i].reset(raw);
      if (status != StatusCode::kOk) return status;
    }
    const auto* dynamics = state->Dynamics();
    const Vector* old_values[] = {
        state->CommittedDisplacement(), dynamics->CommittedVelocity(),
        dynamics->CommittedAcceleration(), dynamics->CommittedVelocity(),
        dynamics->CommittedAcceleration()};
    for (int i = 0; i < 5; ++i) {
      status = old[i]->GatherFrom(
          old_values[i], i < 3 ? translations.data() : rotations.data(), size,
          0, d);
      if (status != StatusCode::kOk) return status;
    }
    Matrix* raw_q = nullptr;
    status = Matrix::Create(nodes, 4, &raw_q, d);
    std::unique_ptr<Matrix> q(raw_q);
    if (status == StatusCode::kOk)
      status = q->CopyFrom(state->Rotations()->CandidateOrientations(), d);
    if (status != StatusCode::kOk) return status;
    Vector *raw_u = nullptr, *raw_node_q = nullptr;
    status = Vector::Create(3, &raw_u, d);
    std::unique_ptr<Vector> node_u(raw_u);
    if (status == StatusCode::kOk) status = Vector::Create(4, &raw_node_q, d);
    std::unique_ptr<Vector> node_q(raw_node_q);
    if (status != StatusCode::kOk) return status;
    if (operation != Operation::kPredict)
      for (int row = 0; row < nodes; ++row) {
        const int indices[] = {3 * row, 3 * row + 1, 3 * row + 2};
        status = state->GatherCandidateTranslation(row, node_u.get(), d);
        if (status == StatusCode::kOk)
          status = trial[0]->ScatterAdd(node_u.get(), indices, 3, d);
        if (status != StatusCode::kOk) return status;
      }
    TranslationalStateView committed;
    committed.displacement = old[0].get();
    committed.velocity = old[1].get();
    committed.acceleration = old[2].get();
    TranslationalCandidate candidate;
    candidate.displacement = trial[0].get();
    candidate.velocity = trial[1].get();
    candidate.acceleration = trial[2].get();
    RotationalStateView committed_rotation;
    committed_rotation.orientations =
        state->Rotations()->CommittedOrientations();
    committed_rotation.velocity = old[3].get();
    committed_rotation.acceleration = old[4].get();
    RotationalCandidate candidate_rotation;
    candidate_rotation.orientations = q.get();
    candidate_rotation.velocity = trial[3].get();
    candidate_rotation.acceleration = trial[4].get();
    if (operation == Operation::kLinearize) {
      status = trial[3]->GatherFrom(dynamics->CandidateVelocity(),
                                    rotations.data(), size, 0, d);
      if (status == StatusCode::kOk)
        status = trial[4]->GatherFrom(dynamics->CandidateAcceleration(),
                                      rotations.data(), size, 0, d);
      TranslationalLinearization local;
      if (status == StatusCode::kOk)
        status = integrator->EvaluateKinematicLinearization(dt, &local, d);
      if (status == StatusCode::kOk)
        status = integrator->EvaluateRotationalLinearization(
            &committed_rotation, &candidate_rotation, dt, rotation, branch, d);
      if (status == StatusCode::kOk) {
        translation->velocity_per_displacement =
            local.velocity_per_displacement;
        translation->acceleration_per_displacement =
            local.acceleration_per_displacement;
      }
      return status;
    }
    status = operation == Operation::kPredict
                 ? integrator->PrepareKinematics(&committed, dt, &candidate, d)
                 : integrator->UpdateKinematics(&committed, dt, &candidate, d);
    if (status != StatusCode::kOk) return status;
    status = operation == Operation::kPredict
                 ? integrator->PrepareRotationalKinematics(
                       &committed_rotation, dt, &candidate_rotation, branch, d)
                 : integrator->UpdateRotationalKinematics(
                       &committed_rotation, dt, &candidate_rotation, branch, d);
    if (status != StatusCode::kOk) return status;
    if (operation == Operation::kPredict)
      for (int row = 0; row < nodes; ++row) {
        const int u_indices[] = {3 * row, 3 * row + 1, 3 * row + 2};
        const int q_indices[] = {4 * row, 4 * row + 1, 4 * row + 2,
                                 4 * row + 3};
        status = node_u->GatherFrom(trial[0].get(), u_indices, 3, 0, d);
        if (status == StatusCode::kOk)
          status = q->GatherTo(q_indices, 4, node_q.get(), 0, d);
        if (status == StatusCode::kOk)
          status = mutable_state->SetCandidateTranslation(row, node_u.get(), d);
        if (status == StatusCode::kOk)
          status = mutable_state->SetCandidateOrientation(row, node_q.get(),
                                                          branch, d);
        if (status != StatusCode::kOk) return status;
      }
    std::unique_ptr<Vector> velocity, acceleration;
    Vector* raw = nullptr;
    status = Vector::Create(dofs->GlobalCount(), &raw, d);
    velocity.reset(raw);
    if (status != StatusCode::kOk) return status;
    status = Vector::Create(dofs->GlobalCount(), &raw, d);
    acceleration.reset(raw);
    if (status != StatusCode::kOk) return status;
    status = velocity->ScatterAdd(trial[1].get(), translations.data(), size, d);
    if (status == StatusCode::kOk)
      status = velocity->ScatterAdd(trial[3].get(), rotations.data(), size, d);
    if (status == StatusCode::kOk)
      status = acceleration->ScatterAdd(trial[2].get(), translations.data(),
                                        size, d);
    if (status == StatusCode::kOk)
      status =
          acceleration->ScatterAdd(trial[4].get(), rotations.data(), size, d);
    return status == StatusCode::kOk
               ? mutable_state->StageDynamicKinematics(
                     velocity.get(), acceleration.get(), time, total, d)
               : status;
  } catch (const std::bad_alloc&) {
    return model_storage::AllocationFailure("Finite kinematic scratch", d);
  }
}
}  // namespace
StatusCode StageFiniteHhtKinematics(const DofManager* dofs,
                                    const TimeIntegrator* integrator, double dt,
                                    double step_time, double total_time,
                                    bool predict, AnalysisState* state,
                                    bool* branch, Diagnostic* d) {
  const auto operation = predict ? Operation::kPredict : Operation::kUpdate;
  auto status =
      Validate(dofs, integrator, operation, state, nullptr, nullptr, branch, d);
  if (status == StatusCode::kOk && predict)
    status = EvaluateBeamInternals(integrator, dt, true, state, d);
  if (status == StatusCode::kOk)
    status = Evaluate(dofs, integrator, dt, step_time, total_time, operation,
                      state, state, nullptr, nullptr, branch, d);
  if (status == StatusCode::kOk)
    status = EvaluateBeamInternals(integrator, dt, false, state, d);
  if (status != StatusCode::kOk && state) state->DiscardCandidate();
  return status;
}
StatusCode EvaluateFiniteHhtKinematicLinearization(
    const DofManager* dofs, const TimeIntegrator* integrator, double dt,
    const AnalysisState* state, TranslationalLinearization* translation,
    RotationalLinearization* rotation, bool* branch, Diagnostic* d) {
  return Evaluate(dofs, integrator, dt, 0, 0, Operation::kLinearize, state,
                  nullptr, translation, rotation, branch, d);
}
}  // namespace fesa
