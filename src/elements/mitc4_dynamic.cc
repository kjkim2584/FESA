#include <cmath>
#include <memory>

#include "elements/mitc4_geometry.h"
#include "elements/mitc4_initial.h"
#include "elements/mitc4_rayleigh.h"
#include "fesa/elements/mitc4_shell.h"
#include "fesa/math/rotation.h"
#include "model/model_storage.h"

namespace fesa {

namespace {
StatusCode NewMatrix(int rows, int columns, std::unique_ptr<Matrix>* owner,
                     Diagnostic* d) {
  Matrix* raw = nullptr;
  const auto status = Matrix::Create(rows, columns, &raw, d);
  owner->reset(raw);
  return status;
}
StatusCode NewVector(int count, std::unique_ptr<Vector>* owner, Diagnostic* d) {
  Vector* raw = nullptr;
  const auto status = Vector::Create(count, &raw, d);
  owner->reset(raw);
  return status;
}
StatusCode Skew(const Vector* v, Matrix* out, Diagnostic* d) {
  const double* x = v->Data();
  // Coefficient placement for S(x)y=x cross y; products use common BLAS.
  const double values[] = {0, -x[2], x[1], x[2], 0, -x[0], -x[1], x[0], 0};
  for (int i = 0; i < 9; ++i) {
    const auto status = out->Set(i / 3, i % 3, values[i], d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}
StatusCode Row(const Matrix* matrix, int row, Vector* out, Diagnostic* d) {
  int indices[4] = {};
  for (int c = 0; c < matrix->Columns(); ++c)
    indices[c] = row * matrix->Columns() + c;
  return matrix->GatherTo(indices, matrix->Columns(), out, 0, d);
}
StatusCode BuildFiniteMaps(const ElementMotionInput* input, const Matrix* frame,
                           Matrix* transform, Matrix* position_map,
                           Matrix* velocity_map, Vector* centrifugal,
                           Matrix* directors, Diagnostic* d) {
  std::unique_ptr<Matrix> rotation, sd, sw, sa, square, product, skew_rate;
  for (auto* owner :
       {&rotation, &sd, &sw, &sa, &square, &product, &skew_rate}) {
    const auto status = NewMatrix(3, 3, owner, d);
    if (status != StatusCode::kOk) return status;
  }
  std::unique_ptr<Vector> q, initial, director, omega, eta, rate, cross;
  auto status = NewVector(4, &q, d);
  if (status != StatusCode::kOk) return status;
  for (auto* owner : {&initial, &director, &omega, &eta, &rate, &cross}) {
    status = NewVector(3, owner, d);
    if (status != StatusCode::kOk) return status;
  }
  status = Row(frame, 2, initial.get(), d);
  if (status != StatusCode::kOk) return status;
  const int identity[] = {0, 1, 2};
  for (int node = 0; node < 4; ++node) {
    status = Row(input->finite_state->orientations, node, q.get(), d);
    if (status != StatusCode::kOk) return status;
    status = RotationMatrix(q.get(), rotation.get(), d);
    if (status != StatusCode::kOk) return status;
    status = Matrix::MultiplyVector(rotation.get(), false, initial.get(),
                                    director.get(), d);
    if (status != StatusCode::kOk) return status;
    status = directors->SetRow(node, director.get(), d);
    if (status != StatusCode::kOk) return status;
    status = Skew(director.get(), sd.get(), d);
    if (status != StatusCode::kOk) return status;
    status = product->CopyFrom(sd.get(), d);
    if (status != StatusCode::kOk) return status;
    status = product->Scale(-1, d);
    if (status != StatusCode::kOk) return status;
    const int offset = 6 * node + 3;
    status =
        transform->CopyBlockFrom(product.get(), 0, 0, 3, 3, offset, offset, d);
    if (status != StatusCode::kOk) return status;
    for (int c = 0; c < 3; ++c) {
      status = transform->Set(6 * node + c, 6 * node + c, 1, d);
      if (status != StatusCode::kOk) return status;
    }
    const int indices[] = {offset, offset + 1, offset + 2};
    status = omega->GatherFrom(input->velocity, indices, 3, 0, d);
    if (status != StatusCode::kOk) return status;
    status = Skew(omega.get(), sw.get(), d);
    if (status != StatusCode::kOk) return status;
    status =
        Matrix::Multiply(sw.get(), false, sw.get(), false, square.get(), d);
    if (status != StatusCode::kOk) return status;
    status = Matrix::MultiplyVector(square.get(), false, director.get(),
                                    cross.get(), d);
    if (status != StatusCode::kOk) return status;
    status = centrifugal->GatherFrom(cross.get(), identity, 3, offset, d);
    if (status != StatusCode::kOk) return status;
    if (!position_map) continue;
    status = eta->GatherFrom(input->acceleration, indices, 3, 0, d);
    if (status != StatusCode::kOk) return status;
    status = Skew(eta.get(), sa.get(), d);
    if (status != StatusCode::kOk) return status;
    status = square->AddScaled(1, sa.get(), d);
    if (status != StatusCode::kOk) return status;
    status = Matrix::Multiply(square.get(), false, sd.get(), false,
                              product.get(), d);
    if (status != StatusCode::kOk) return status;
    status = product->Scale(-1, d);
    if (status != StatusCode::kOk) return status;
    status = position_map->CopyBlockFrom(product.get(), 0, 0, 3, 3, offset,
                                         offset, d);
    if (status != StatusCode::kOk) return status;
    status =
        Matrix::MultiplyVector(sw.get(), false, director.get(), rate.get(), d);
    if (status != StatusCode::kOk) return status;
    status = Skew(rate.get(), skew_rate.get(), d);
    if (status != StatusCode::kOk) return status;
    status =
        Matrix::Multiply(sw.get(), false, sd.get(), false, product.get(), d);
    if (status != StatusCode::kOk) return status;
    status = product->AddScaled(1, skew_rate.get(), d);
    if (status != StatusCode::kOk) return status;
    status = product->Scale(-1, d);
    if (status != StatusCode::kOk) return status;
    status = velocity_map->CopyBlockFrom(product.get(), 0, 0, 3, 3, offset,
                                         offset, d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}
// The derivative of T^T times a director-coordinate resultant, shared by
// inertia and RD10 mass drag; no force or matrix survives this call.
StatusCode AddDirectorGeometric(const Matrix* directors, const Vector* weighted,
                                 Matrix* jacobian, Diagnostic* d) {
  auto status = StatusCode::kOk;
  std::unique_ptr<Matrix> sd, ss, product, diagonal;
  for (auto* owner : {&sd, &ss, &product, &diagonal}) {
    status = NewMatrix(3, 3, owner, d);
    if (status != StatusCode::kOk) return status;
  }
  std::unique_ptr<Vector> director, resultant;
  status = NewVector(3, &director, d);
  if (status != StatusCode::kOk) return status;
  status = NewVector(3, &resultant, d);
  if (status != StatusCode::kOk) return status;
  for (int node = 0; node < 4; ++node) {
    const int offset = 6 * node + 3;
    const int indices[] = {offset, offset + 1, offset + 2};
    status = Row(directors, node, director.get(), d);
    if (status != StatusCode::kOk) return status;
    status = resultant->GatherFrom(weighted, indices, 3, 0, d);
    if (status != StatusCode::kOk) return status;
    status = Skew(director.get(), sd.get(), d);
    if (status != StatusCode::kOk) return status;
    status = Skew(resultant.get(), ss.get(), d);
    if (status != StatusCode::kOk) return status;
    status =
        Matrix::Multiply(ss.get(), false, sd.get(), false, product.get(), d);
    if (status != StatusCode::kOk) return status;
    status = diagonal->CopyBlockFrom(jacobian, offset, offset, 3, 3, 0, 0, d);
    if (status != StatusCode::kOk) return status;
    status = diagonal->AddScaled(1, product.get(), d);
    if (status != StatusCode::kOk) return status;
    status = jacobian->CopyBlockFrom(diagonal.get(), 0, 0, 3, 3, offset,
                                      offset, d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}
StatusCode FiniteMassDamping(double alpha, const ElementMotionInput* input,
                             const Matrix* base, const Matrix* transform,
                             const Matrix* directors, Vector* force,
                             Matrix* gd, Matrix* cd, Diagnostic* d) {
  std::unique_ptr<Vector> rate, weighted;
  auto status = NewVector(24, &rate, d);
  if (status == StatusCode::kOk) status = NewVector(24, &weighted, d);
  if (status == StatusCode::kOk)
    status = Matrix::MultiplyVector(transform, false, input->velocity,
                                    rate.get(), d);
  if (status == StatusCode::kOk)
    status = Matrix::MultiplyVector(base, false, rate.get(), weighted.get(), d);
  if (status == StatusCode::kOk)
    status = Matrix::MultiplyVector(transform, true, weighted.get(), force, d);
  if (status == StatusCode::kOk) status = force->Scale(alpha, d);
  if (status != StatusCode::kOk || !gd) return status;
  std::unique_ptr<Matrix> map, work, sd, sw, product;
  for (auto* owner : {&map, &work}) {
    status = NewMatrix(24, 24, owner, d);
    if (status != StatusCode::kOk) return status;
  }
  for (auto* owner : {&sd, &sw, &product}) {
    status = NewMatrix(3, 3, owner, d);
    if (status != StatusCode::kOk) return status;
  }
  std::unique_ptr<Vector> director, omega;
  status = NewVector(3, &director, d);
  if (status == StatusCode::kOk) status = NewVector(3, &omega, d);
  if (status != StatusCode::kOk) return status;
  for (int node = 0; node < 4; ++node) {
    const int offset = 6 * node + 3;
    const int indices[] = {offset, offset + 1, offset + 2};
    status = Row(directors, node, director.get(), d);
    if (status == StatusCode::kOk)
      status = omega->GatherFrom(input->velocity, indices, 3, 0, d);
    if (status == StatusCode::kOk) status = Skew(director.get(), sd.get(), d);
    if (status == StatusCode::kOk) status = Skew(omega.get(), sw.get(), d);
    if (status == StatusCode::kOk)
      status = Matrix::Multiply(sw.get(), false, sd.get(), false,
                                product.get(), d);
    if (status == StatusCode::kOk) status = product->Scale(-1, d);
    if (status == StatusCode::kOk)
      status = map->CopyBlockFrom(product.get(), 0, 0, 3, 3, offset, offset, d);
    if (status != StatusCode::kOk) return status;
  }
  status = Matrix::Multiply(base, false, map.get(), false, work.get(), d);
  if (status == StatusCode::kOk)
    status = Matrix::Multiply(transform, true, work.get(), false, gd, d);
  if (status == StatusCode::kOk)
    status = AddDirectorGeometric(directors, weighted.get(), gd, d);
  if (status == StatusCode::kOk) status = gd->Scale(alpha, d);
  if (status == StatusCode::kOk)
    status = Matrix::Multiply(base, false, transform, false, work.get(), d);
  if (status == StatusCode::kOk)
    status = Matrix::Multiply(transform, true, work.get(), false, cd, d);
  return status == StatusCode::kOk ? cd->Scale(alpha, d) : status;
}
StatusCode FiniteInertia(const Mitc4Shell* element, double thickness,
                         const ElementMotionInput* input, MassMode mode,
                         DynamicEvaluationRequest request,
                         ElementDynamicResult* out, Diagnostic* d) {
  const bool derivatives =
      request == DynamicEvaluationRequest::kForcesAndLinearization;
  std::unique_ptr<Matrix> initial_mass, base, transform, frame, directors,
      position_map, velocity_map, work, gi, ci, mi, zero_matrix, gd, cd;
  for (auto* owner : {&initial_mass, &base, &transform}) {
    const auto status = NewMatrix(24, 24, owner, d);
    if (status != StatusCode::kOk) return status;
  }
  auto status = element->ComputeMass(mode, nullptr, initial_mass.get(), d);
  if (status != StatusCode::kOk) return status;
  status = NewMatrix(3, 3, &frame, d);
  if (status != StatusCode::kOk) return status;
  status = element->ComputeInitialFrame(frame.get(), d);
  if (status != StatusCode::kOk) return status;
  status = NewMatrix(4, 3, &directors, d);
  if (status != StatusCode::kOk) return status;
  if (derivatives) {
    for (auto* owner :
         {&position_map, &velocity_map, &work, &gi, &ci, &mi, &zero_matrix}) {
      status = NewMatrix(24, 24, owner, d);
      if (status != StatusCode::kOk) return status;
    }
  }
  // The reference scalar mass is isotropic in global translational components.
  // B's rotational coordinates are director accelerations, not angular DOFs.
  const double radius2 = thickness * thickness / 12;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      const double m = initial_mass->Data()[(6 * i) * 24 + 6 * j];
      const double rotary = radius2 * m;
      if (!std::isfinite(rotary) || (m > 0 && !(rotary > 0)))
        return ReportFailure(StatusCode::kLibraryFailure,
                             DiagnosticStage::kComputation,
                             "MITC4 director inertia coefficient range", d);
      for (int c = 0; c < 6; ++c) {
        status = base->Set(6 * i + c, 6 * j + c, c < 3 ? m : rotary, d);
        if (status != StatusCode::kOk) return status;
      }
    }
  }
  std::unique_ptr<Vector> centrifugal, drive, weighted, force, zero;
  for (auto* owner : {&centrifugal, &drive, &weighted, &force, &zero}) {
    status = NewVector(24, owner, d);
    if (status != StatusCode::kOk) return status;
  }
  status = BuildFiniteMaps(input, frame.get(), transform.get(),
                           position_map.get(), velocity_map.get(),
                           centrifugal.get(), directors.get(), d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::MultiplyVector(transform.get(), false, input->acceleration,
                                  drive.get(), d);
  if (status != StatusCode::kOk) return status;
  status = drive->AddScaled(1, centrifugal.get(), d);
  if (status != StatusCode::kOk) return status;
  status =
      Matrix::MultiplyVector(base.get(), false, drive.get(), weighted.get(), d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::MultiplyVector(transform.get(), true, weighted.get(),
                                  force.get(), d);
  if (status != StatusCode::kOk) return status;
  if (derivatives) {
    const Matrix* maps[] = {position_map.get(), velocity_map.get(),
                            transform.get()};
    Matrix* results[] = {gi.get(), ci.get(), mi.get()};
    for (int i = 0; i < 3; ++i) {
      status =
          Matrix::Multiply(base.get(), false, maps[i], false, work.get(), d);
      if (status != StatusCode::kOk) return status;
      status = Matrix::Multiply(transform.get(), true, work.get(), false,
                                results[i], d);
      if (status != StatusCode::kOk) return status;
    }
    status = AddDirectorGeometric(directors.get(), weighted.get(), gi.get(), d);
    if (status != StatusCode::kOk) return status;
  }
  const double alpha = element->material()->rayleigh_alpha();
  if (alpha != 0) {
    if (derivatives) {
      status = NewMatrix(24, 24, &gd, d);
      if (status == StatusCode::kOk) status = NewMatrix(24, 24, &cd, d);
      if (status != StatusCode::kOk) return status;
    }
    status = FiniteMassDamping(alpha, input, base.get(), transform.get(),
                               directors.get(), zero.get(), gd.get(), cd.get(),
                               d);
    if (status != StatusCode::kOk) return status;
  }
  status = out->inertia_force->CopyFrom(force.get(), d);
  if (status != StatusCode::kOk) return status;
  status = out->damping_force->CopyFrom(zero.get(), d);
  if (status != StatusCode::kOk || !derivatives) return status;
  const Matrix* values[] = {gi.get(), ci.get(), mi.get(),
                            gd ? gd.get() : zero_matrix.get(),
                            cd ? cd.get() : zero_matrix.get()};
  Matrix* outputs[] = {out->inertia_displacement, out->inertia_velocity,
                       out->inertia_acceleration, out->damping_displacement,
                       out->damping_velocity};
  for (int i = 0; i < 5; ++i) {
    status = outputs[i]->CopyFrom(values[i], d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}
}  // namespace

StatusCode Mitc4Shell::ComputeInitialFrame(Matrix* out, Diagnostic* d) const {
  if (!out || out->Rows() != 3 || out->Columns() != 3)
    return model_storage::Invalid("MITC4 initial frame output", d);
  Matrix* raw = nullptr;
  auto status = Matrix::Create(4, 3, &raw, d);
  std::unique_ptr<Matrix> coordinates(raw);
  if (status != StatusCode::kOk) return status;
  double diameter = 0, length = 0;
  return ValidateMitc4Geometry(nodes_, 4, coordinates.get(), out, &diameter,
                               &length, d);
}
bool Mitc4Shell::SupportsDynamicEvaluation(DynamicKinematics kinematics) const {
  return kinematics == DynamicKinematics::kLinear ||
         kinematics == DynamicKinematics::kFiniteRotation;
}

StatusCode Mitc4Shell::EvaluateDynamic(const ElementMotionInput* input,
                                       MassMode mode,
                                       const ElementWorkspace* workspace,
                                       DynamicEvaluationRequest request,
                                       ElementDynamicResult* out,
                                       TrialEvaluationFailure* failure,
                                       Diagnostic* d) const {
  if (failure) *failure = TrialEvaluationFailure::kNone;
  if (!SupportsDynamicEvaluation(input ? input->kinematics
                                       : DynamicKinematics::kLinear))
    return Element::EvaluateDynamic(input, mode, workspace, request, out,
                                    failure, d);
  if (!failure)
    return model_storage::Invalid("MITC4 dynamic failure output", d);
  auto status = ValidateElementDynamicEvaluation(4, 6, input, request, out, d);
  if (status != StatusCode::kOk) return status;
  if (!SupportsMass(mode))
    return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                         "MITC4 dynamic mass mode", d);
  if (input->kinematics == DynamicKinematics::kFiniteRotation) {
    status = FiniteInertia(this, section_->thickness(), input, mode, request,
                           out, d);
    if (status != StatusCode::kOk || material_->rayleigh_beta() == 0)
      return status;
    const bool derivatives =
        request == DynamicEvaluationRequest::kForcesAndLinearization;
    std::unique_ptr<Vector> force;
    std::unique_ptr<Matrix> gd, cd;
    status = NewVector(24, &force, d);
    if (derivatives && status == StatusCode::kOk) {
      status = NewMatrix(24, 24, &gd, d);
      if (status == StatusCode::kOk) status = NewMatrix(24, 24, &cd, d);
    }
    if (status == StatusCode::kOk)
      status = EvaluateMitc4RayleighBeta(nodes_, section_, input, request,
          force.get(), gd.get(), cd.get(), failure, d);
    if (status == StatusCode::kOk)
      status = out->damping_force->AddScaled(1, force.get(), d);
    if (derivatives && status == StatusCode::kOk)
      status = out->damping_displacement->AddScaled(1, gd.get(), d);
    if (derivatives && status == StatusCode::kOk)
      status = out->damping_velocity->AddScaled(1, cd.get(), d);
    return status;
  }
  Matrix* raw_mass = nullptr;
  status = Matrix::Create(24, 24, &raw_mass, d);
  std::unique_ptr<Matrix> mass(raw_mass);
  if (status != StatusCode::kOk) return status;
  // The initial shell mass ignores workspace and retains no numerical cache.
  status = ComputeMass(mode, nullptr, mass.get(), d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::MultiplyVector(mass.get(), false, input->acceleration,
                                  out->inertia_force, d);
  if (status != StatusCode::kOk) return status;
  Matrix* raw_damping = nullptr;
  status = Matrix::Create(24, 24, &raw_damping, d);
  std::unique_ptr<Matrix> damping(raw_damping);
  if (status != StatusCode::kOk) return status;
  if (material_->rayleigh_beta() != 0) {
    status = ComputeMitc4InitialElasticStiffness(nodes_, 4, section_,
                                                 damping.get(), d);
    if (status == StatusCode::kOk)
      status = damping->Scale(material_->rayleigh_beta(), d);
  }
  if (status == StatusCode::kOk && material_->rayleigh_alpha() != 0)
    status = damping->AddScaled(material_->rayleigh_alpha(), mass.get(), d);
  if (status == StatusCode::kOk)
    status = Matrix::MultiplyVector(damping.get(), false, input->velocity,
                                    out->damping_force, d);
  if (status != StatusCode::kOk || request == DynamicEvaluationRequest::kForces)
    return status;
  status = out->inertia_acceleration->CopyFrom(mass.get(), d);
  if (status != StatusCode::kOk) return status;
  Matrix* raw_zero = nullptr;
  status = Matrix::Create(24, 24, &raw_zero, d);
  std::unique_ptr<Matrix> zero(raw_zero);
  if (status != StatusCode::kOk) return status;
  for (auto* derivative : {out->inertia_displacement, out->inertia_velocity,
                           out->damping_displacement}) {
    status = derivative->CopyFrom(zero.get(), d);
    if (status != StatusCode::kOk) return status;
  }
  return out->damping_velocity->CopyFrom(damping.get(), d);
}
}  // namespace fesa
