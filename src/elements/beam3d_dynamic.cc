#include "fesa/elements/beam3d.h"

#include <memory>

#include "elements/beam3d_dynamic_metrics.h"
#include "elements/beam3d_inertia.h"
#include "elements/beam3d_nonlinear_ops.h"
#include "fesa/elements/element_workspace.h"
#include "model/model_storage.h"

namespace fesa {
namespace {
StatusCode NewMatrix(int rows, int columns, std::unique_ptr<Matrix>* out,
                     Diagnostic* d) {
  Matrix* raw = nullptr;
  const auto status = Matrix::Create(rows, columns, &raw, d);
  out->reset(raw);
  return status;
}
StatusCode NewVector(int size, std::unique_ptr<Vector>* out, Diagnostic* d) {
  Vector* raw = nullptr;
  const auto status = Vector::Create(size, &raw, d);
  out->reset(raw);
  return status;
}
#define FESA_BEAM_DYNAMIC_TRY(expression)         \
  do {                                            \
    const auto status = (expression);             \
    if (status != StatusCode::kOk) return status; \
  } while (false)
}  // namespace

bool Beam3D::SupportsDynamicEvaluation(DynamicKinematics kinematics) const {
  return kinematics == DynamicKinematics::kFiniteRotation;
}

StatusCode Beam3D::EvaluateDynamic(const ElementMotionInput* input,
                                   MassMode mass_mode,
                                   const ElementWorkspace* workspace,
                                   DynamicEvaluationRequest request,
                                   ElementDynamicResult* out,
                                   TrialEvaluationFailure* failure,
                                   Diagnostic* d) const {
  if (failure) *failure = TrialEvaluationFailure::kNone;
  if (!failure) return model_storage::Invalid("Beam dynamic failure output", d);
  // Legacy linear Beam dynamics remain unsupported. Within the supported
  // finite role, common preflight is first and never reads previous outputs.
  if (input && input->kinematics == DynamicKinematics::kLinear)
    return Element::EvaluateDynamic(input, mass_mode, workspace, request, out,
                                    failure, d);
  FESA_BEAM_DYNAMIC_TRY(
      ValidateBeamElementDynamicEvaluation(input, request, out, d));
  if (mass_mode != MassMode::kConsistent)
    return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                         "Beam dynamic consistent mass required", d);
  const bool full =
      request == DynamicEvaluationRequest::kForcesAndLinearization;
  std::unique_ptr<Matrix> positions, relative_positions, tangents, normals;
  for (auto* owner : {&positions, &relative_positions, &tangents, &normals})
    FESA_BEAM_DYNAMIC_TRY(NewMatrix(2, 3, owner, d));
  std::unique_ptr<Vector> velocity, acceleration, direction, fi, fd;
  for (auto* owner : {&velocity, &acceleration, &direction, &fi, &fd})
    FESA_BEAM_DYNAMIC_TRY(NewVector(14, owner, d));
  std::unique_ptr<Matrix> gi, ci, mi, gd, cd;
  if (full)
    for (auto* owner : {&gi, &ci, &mi, &gd, &cd})
      FESA_BEAM_DYNAMIC_TRY(NewMatrix(14, 14, owner, d));
  BeamKinematicWorkspace* current =
      workspace ? workspace->beam_kinematics_ : nullptr;
  if (!current) {
    FESA_BEAM_DYNAMIC_TRY(BeamKinematicWorkspace::Create(&current, d));
    if (workspace) workspace->beam_kinematics_ = current;
  }
  std::unique_ptr<BeamKinematicWorkspace> current_owner(
      workspace ? nullptr : current);
  BeamTimeKinematicWorkspace* time =
      workspace ? workspace->beam_time_kinematics_ : nullptr;
  if (!time) {
    FESA_BEAM_DYNAMIC_TRY(BeamTimeKinematicWorkspace::Create(&time, d));
    if (workspace) workspace->beam_time_kinematics_ = time;
  }
  std::unique_ptr<BeamTimeKinematicWorkspace> time_owner(workspace ? nullptr : time);
  const int physical[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
  const int internal[] = {0, 1};
  FESA_BEAM_DYNAMIC_TRY(
      velocity->GatherFrom(input->velocity, physical, 12, 0, d));
  FESA_BEAM_DYNAMIC_TRY(
      velocity->GatherFrom(input->internal_velocity, internal, 2, 12, d));
  FESA_BEAM_DYNAMIC_TRY(
      acceleration->GatherFrom(input->acceleration, physical, 12, 0, d));
  FESA_BEAM_DYNAMIC_TRY(acceleration->GatherFrom(input->internal_acceleration,
                                                 internal, 2, 12, d));
  FESA_BEAM_DYNAMIC_TRY(
      PrepareBeamNonlinearPose(this, input->finite_state, positions.get(),
                               tangents.get(), normals.get(), d,
                               relative_positions.get()));
  BeamKinematicInput pose;
  pose.length = length_;
  pose.positions = positions.get();
  pose.relative_positions = relative_positions.get();
  pose.tangents = tangents.get();
  pose.normals = normals.get();
  pose.internal_position = input->finite_state->beam_internal_position;
  pose.direction = direction.get();
  FESA_BEAM_DYNAMIC_TRY(EvaluateBeamInertia(
      &pose, section_, velocity.get(), acceleration.get(), current,
      time, request, fi.get(), gi.get(), ci.get(), mi.get(), failure, d));
  const auto* material = section_->GetMaterial();
  if (material->rayleigh_alpha() != 0 || material->rayleigh_beta() != 0) {
    BeamRayleighWorkspace* rayleigh =
        workspace ? workspace->beam_rayleigh_ : nullptr;
    if (!rayleigh) {
      FESA_BEAM_DYNAMIC_TRY(BeamRayleighWorkspace::Create(&rayleigh, d));
      if (workspace) workspace->beam_rayleigh_ = rayleigh;
    }
    std::unique_ptr<BeamRayleighWorkspace> rayleigh_owner(workspace ? nullptr
                                                                    : rayleigh);
    std::unique_ptr<Matrix> damping, dc, row, basis;
    FESA_BEAM_DYNAMIC_TRY(NewMatrix(14, 14, &damping, d));
    FESA_BEAM_DYNAMIC_TRY(rayleigh->EvaluateBase(
        &pose, section_, full ? mi.get() : nullptr, damping.get(), failure, d));
    FESA_BEAM_DYNAMIC_TRY(Matrix::MultiplyVector(damping.get(), false,
                                                 velocity.get(), fd.get(), d));
    if (full) {
      FESA_BEAM_DYNAMIC_TRY(cd->CopyFrom(damping.get(), d));
      FESA_BEAM_DYNAMIC_TRY(NewMatrix(14, 14, &dc, d));
      FESA_BEAM_DYNAMIC_TRY(NewMatrix(1, 14, &row, d));
      FESA_BEAM_DYNAMIC_TRY(NewMatrix(1, 14, &basis, d));
      std::unique_ptr<Vector> column;
      FESA_BEAM_DYNAMIC_TRY(NewVector(14, &column, d));
      for (int j = 0; j < 14; ++j) {
        if (j > 0) FESA_BEAM_DYNAMIC_TRY(direction->Set(j - 1, 0, d));
        FESA_BEAM_DYNAMIC_TRY(direction->Set(j, 1, d));
        FESA_BEAM_DYNAMIC_TRY(
            rayleigh->EvaluateDirection(&pose, dc.get(), failure, d));
        FESA_BEAM_DYNAMIC_TRY(Matrix::MultiplyVector(
            dc.get(), false, velocity.get(), column.get(), d));
        FESA_BEAM_DYNAMIC_TRY(row->SetRow(0, column.get(), d));
        FESA_BEAM_DYNAMIC_TRY(basis->SetRow(0, direction.get(), d));
        // Library outer product inserts one full spatial derivative column.
        FESA_BEAM_DYNAMIC_TRY(Matrix::MultiplyAdd(row.get(), true, basis.get(),
                                                  false, 1, 1, gd.get(), d));
      }
    }
  }
  // Fresh zero fd/gd/cd also overwrite nonfinite values from an old attempt.
  // No caller-visible output is touched before all requested calculations pass.
  FESA_BEAM_DYNAMIC_TRY(out->inertia_force->CopyFrom(fi.get(), d));
  FESA_BEAM_DYNAMIC_TRY(out->damping_force->CopyFrom(fd.get(), d));
  if (full) {
    FESA_BEAM_DYNAMIC_TRY(out->inertia_displacement->CopyFrom(gi.get(), d));
    FESA_BEAM_DYNAMIC_TRY(out->inertia_velocity->CopyFrom(ci.get(), d));
    FESA_BEAM_DYNAMIC_TRY(out->inertia_acceleration->CopyFrom(mi.get(), d));
    FESA_BEAM_DYNAMIC_TRY(out->damping_displacement->CopyFrom(gd.get(), d));
    FESA_BEAM_DYNAMIC_TRY(out->damping_velocity->CopyFrom(cd.get(), d));
  }
  return StatusCode::kOk;
}
#undef FESA_BEAM_DYNAMIC_TRY
}  // namespace fesa
