#include "fesa/elements/beam3d.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <new>

#include "fesa/elements/element_workspace.h"
#include "fesa/math/matrix.h"
#include "fesa/math/rotation.h"
#include "elements/beam3d_nonlinear_ops.h"
#include "model/model_storage.h"

namespace fesa {
namespace {
StatusCode Invalid(Diagnostic* d) {
  return ReportFailure(StatusCode::kInputError, DiagnosticStage::kValidation,
                       "Beam3D invalid length/orientation/coefficients", d);
}
StatusCode InvalidCall(const char* operation, Diagnostic* d) {
  return ReportFailure(StatusCode::kInvalidArgument,
                       DiagnosticStage::kValidation, operation, d);
}
// Evaluate one B2/B3 scalar coefficient without overflow in intermediate
// products whose final result is representable. This is not vector arithmetic.
double Coefficient(double factor, double a, double b, double length,
                   int power) {
  int ea = 0, eb = 0, el = 0;
  const double ma = std::frexp(a, &ea);
  const double mb = std::frexp(b, &eb);
  const double ml = std::frexp(length, &el);
  double mantissa = factor * ma * mb;
  for (int i = 0; i < power; ++i) mantissa /= ml;
  return std::ldexp(mantissa, ea + eb - power * el);
}
// B5 scalar reciprocal of a product without intermediate overflow/underflow.
double ReciprocalProduct(double a, double b) {
  int ea = 0, eb = 0;
  const double ma = std::frexp(a, &ea);
  const double mb = std::frexp(b, &eb);
  return std::ldexp(1. / (ma * mb), -ea - eb);
}
// Finite scan selects a scale, while all vector arithmetic uses MKL wrappers.
StatusCode Normalize(Vector* v, bool require_length, double* length,
                     Diagnostic* d) {
  double maximum = 0;
  for (int i = 0; i < 3; ++i)
    maximum = std::max(maximum, std::abs(v->Data()[i]));
  if (maximum == 0 || !std::isfinite(maximum)) return Invalid(d);
  int exponent = 0;
  std::frexp(maximum, &exponent);
  StatusCode status = StatusCode::kOk;
  if (-exponent > 1023) {
    status = v->Scale(std::ldexp(1., 1023), d);
    if (status != StatusCode::kOk) return status;
    status = v->Scale(std::ldexp(1., -exponent - 1023), d);
  } else
    status = v->Scale(std::ldexp(1., -exponent), d);
  if (status != StatusCode::kOk) return status;
  double norm = 0;
  status = v->Norm2(&norm, d);
  if (status != StatusCode::kOk) return status;
  if (require_length) {
    *length = std::ldexp(norm, exponent);
    if (!std::isfinite(*length) || *length <= 0) return Invalid(d);
  }
  return v->Scale(1. / norm, d);
}
}  // namespace
Beam3D::Beam3D(int id, const Node* first, const Node* second,
               const BeamSection* section)
    : Element(id), first_(first), second_(second), section_(section) {}
Beam3D::~Beam3D() {
  delete local_axes_;
  delete internal_force_;
  delete end_results_;
}
ElementFamily Beam3D::family() const { return ElementFamily::kBeam; }
const char* Beam3D::type_name() const { return "B33"; }
int Beam3D::node_count() const { return 2; }
int Beam3D::dofs_per_node() const { return 6; }
int Beam3D::dof_count() const { return 12; }
int Beam3D::internal_dof_count() const { return 2; }
const Node* Beam3D::GetNode(int end) const {
  return end == 0 ? first_ : (end == 1 ? second_ : nullptr);
}
StatusCode Beam3D::CharacteristicLength(double* out, Diagnostic* d) const {
  if (!out) return InvalidCall("Beam3D characteristic length output", d);
  *out = length_;
  return StatusCode::kOk;
}
const BeamSection* Beam3D::Section() const { return section_; }
const Vector* Beam3D::InternalForce() const { return internal_force_; }
const Vector* Beam3D::EndResults() const { return end_results_; }
Vector* Beam3D::MutableInternalForce() { return internal_force_; }
Vector* Beam3D::MutableEndResults() { return end_results_; }
StatusCode Beam3D::CreateResults(Diagnostic* d) {
  auto status = Vector::Create(12, &internal_force_, d);
  if (status != StatusCode::kOk) return status;
  return Vector::Create(16, &end_results_, d);
}
const Matrix* Beam3D::LocalAxes() const { return local_axes_; }
bool Beam3D::SupportsTrialEvaluation() const { return true; }
bool Beam3D::SupportsTrialGeometryCheck() const { return true; }

StatusCode Beam3D::EvaluateTrial(const ElementTrialState* state,
                                 const ElementWorkspace* context,
                                 ElementEvaluationRequest request,
                                 Vector* force, Matrix* tangent,
                                 TrialEvaluationFailure* failure,
                                 Diagnostic* d) const {
  if (failure) *failure = TrialEvaluationFailure::kNone;
  if (!state || !failure) return InvalidCall("Beam3D trial input", d);
  BeamKinematicWorkspace* raw = context ? context->beam_kinematics_ : nullptr;
  auto status = raw ? StatusCode::kOk : BeamKinematicWorkspace::Create(&raw, d);
  std::unique_ptr<BeamKinematicWorkspace> owner(context ? nullptr : raw);
  if (context) context->beam_kinematics_ = raw;
  if (status != StatusCode::kOk) return status;
  return EvaluateBeamNonlinearResponse(this, state,
                                       state->beam_internal_position, raw,
                                       request, force, tangent, failure, d);
}

StatusCode Beam3D::CheckTrialGeometry(const ElementTrialState* state,
                                      const ElementWorkspace* context,
                                      TrialEvaluationFailure* failure,
                                      Diagnostic* d) const {
  if (failure) *failure = TrialEvaluationFailure::kNone;
  if (!state || !failure || !state->beam_internal_position ||
      state->beam_internal_position->Size() != 2 || state->incremental_solid ||
      (state->body_spin && state->body_spin->Size() != 6) ||
      (state->committed_orientations &&
       (state->committed_orientations->Rows() != 2 ||
        state->committed_orientations->Columns() != 4)))
    return InvalidCall("Beam3D geometry trial input", d);
  for (const auto* values : {state->body_spin, state->beam_internal_position})
    if (values && !std::all_of(values->Data(), values->Data() + values->Size(),
                               [](double x) { return std::isfinite(x); }))
      return InvalidCall("Beam3D finite geometry input", d);
  if (state->committed_orientations &&
      !std::all_of(state->committed_orientations->Data(),
                   state->committed_orientations->Data() + 8,
                   [](double x) { return std::isfinite(x); }))
    return InvalidCall("Beam3D finite geometry history", d);
  Matrix *positions = nullptr, *tangents = nullptr, *normals = nullptr;
  auto status = Matrix::Create(2, 3, &positions, d);
  std::unique_ptr<Matrix> x(positions);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Create(2, 3, &tangents, d);
  std::unique_ptr<Matrix> t(tangents);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Create(2, 3, &normals, d);
  std::unique_ptr<Matrix> n(normals);
  if (status != StatusCode::kOk) return status;
  Matrix* raw_relative = nullptr;
  status = Matrix::Create(2, 3, &raw_relative, d);
  std::unique_ptr<Matrix> relative(raw_relative);
  if (status != StatusCode::kOk) return status;
  status =
      PrepareBeamNonlinearPose(this, state, positions, tangents, normals, d,
                               relative.get());
  if (status != StatusCode::kOk) return status;
  Vector* raw_direction = nullptr;
  status = Vector::Create(14, &raw_direction, d);
  std::unique_ptr<Vector> direction(raw_direction);
  if (status != StatusCode::kOk) return status;
  BeamKinematicWorkspace* raw = context ? context->beam_kinematics_ : nullptr;
  status = raw ? StatusCode::kOk : BeamKinematicWorkspace::Create(&raw, d);
  std::unique_ptr<BeamKinematicWorkspace> owner(context ? nullptr : raw);
  if (context) context->beam_kinematics_ = raw;
  if (status != StatusCode::kOk) return status;
  BeamKinematicInput input;
  input.length = length_;
  input.positions = positions;
  input.relative_positions = relative.get();
  input.tangents = tangents;
  input.normals = normals;
  input.direction = direction.get();
  input.internal_position = state->beam_internal_position;
  const double inner = std::sqrt((3 - 2 * std::sqrt(6. / 5)) / 7);
  const double outer = std::sqrt((3 + 2 * std::sqrt(6. / 5)) / 7);
  for (double point : {-outer, -inner, inner, outer}) {
    status = raw->Evaluate(&input, (point + 1) / 2, failure, d,
                           BeamKinematicDerivative::kValue);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}

StatusCode Beam3D::UpdateTrialResults(const ElementTrialState* state,
                                      const ElementWorkspace* workspace,
                                      TrialEvaluationFailure* failure,
                                      Diagnostic* d) {
  if (failure) *failure = TrialEvaluationFailure::kNone;
  if (!state || !failure) return InvalidCall("Beam3D recovery trial input", d);
  Vector* raw_force = nullptr;
  auto status = Vector::Create(14, &raw_force, d);
  std::unique_ptr<Vector> force(raw_force);
  if (status != StatusCode::kOk) return status;
  status =
      EvaluateTrial(state, workspace, ElementEvaluationRequest::kInternalForce,
                    force.get(), nullptr, failure, d);
  if (status != StatusCode::kOk) return status;
  Vector *raw_local = nullptr, *raw_end = nullptr, *raw_q = nullptr;
  Vector *raw_global = nullptr, *raw_projected = nullptr;
  status = Vector::Create(12, &raw_local, d);
  std::unique_ptr<Vector> local(raw_local);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(16, &raw_end, d);
  std::unique_ptr<Vector> ends(raw_end);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(4, &raw_q, d);
  std::unique_ptr<Vector> q(raw_q);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(6, &raw_global, d);
  std::unique_ptr<Vector> global(raw_global);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(6, &raw_projected, d);
  std::unique_ptr<Vector> projected(raw_projected);
  if (status != StatusCode::kOk) return status;
  Matrix *raw_rotation = nullptr, *raw_axes = nullptr, *raw_transform = nullptr;
  status = Matrix::Create(3, 3, &raw_rotation, d);
  std::unique_ptr<Matrix> rotation(raw_rotation);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Create(3, 3, &raw_axes, d);
  std::unique_ptr<Matrix> axes(raw_axes);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Create(6, 6, &raw_transform, d);
  std::unique_ptr<Matrix> transform(raw_transform);
  if (status != StatusCode::kOk) return status;
  const int six[] = {0, 1, 2, 3, 4, 5};
  for (int end = 0; end < 2; ++end) {
    const int qi[] = {4 * end, 4 * end + 1, 4 * end + 2, 4 * end + 3};
    const int fi[] = {6 * end,     6 * end + 1, 6 * end + 2,
                      6 * end + 3, 6 * end + 4, 6 * end + 5};
    status = state->orientations->GatherTo(qi, 4, q.get(), 0, d);
    if (status == StatusCode::kOk)
      status = RotationMatrix(q.get(), rotation.get(), d);
    if (status == StatusCode::kOk)
      status = Matrix::Multiply(local_axes_, false, rotation.get(), true,
                                axes.get(), d);
    if (status == StatusCode::kOk)
      status = transform->CopyBlockFrom(axes.get(), 0, 0, 3, 3, 0, 0, d);
    if (status == StatusCode::kOk)
      status = transform->CopyBlockFrom(axes.get(), 0, 0, 3, 3, 3, 3, d);
    if (status == StatusCode::kOk)
      status = global->GatherFrom(force.get(), fi, 6, 0, d);
    if (status == StatusCode::kOk)
      status = Matrix::MultiplyVector(transform.get(), false, global.get(),
                                      projected.get(), d);
    if (status == StatusCode::kOk)
      status = local->GatherFrom(projected.get(), six, 6, 6 * end, d);
    if (status == StatusCode::kOk && end == 0) status = projected->Scale(-1, d);
    if (status == StatusCode::kOk)
      status = ends->GatherFrom(projected.get(), six, 6, 6 * end, d);
    if (status != StatusCode::kOk) return status;
    const double stretch =
        1 + state->beam_internal_position->Data()[end] / length_;
    const double green = (stretch * stretch - 1) / 2;
    status = ends->Set(14 + end, green, d);
    if (status == StatusCode::kOk)
      status = ends->Set(12 + end, section_->GetMaterial()->young() * green, d);
    if (status != StatusCode::kOk) return status;
  }
  status = internal_force_->CopyFrom(local.get(), d);
  return status == StatusCode::kOk ? end_results_->CopyFrom(ends.get(), d)
                                   : status;
}
StatusCode Beam3D::Create(int id, const Node* first, const Node* second,
                          const ElasticMaterial* material,
                          const BeamSection* section, Beam3D** out,
                          Diagnostic* d) {
  if (!out) return model_storage::Invalid("Beam3D output", d);
  *out = nullptr;
  if (!model_storage::ValidId(id) || !first || !second || !section ||
      first->Id() == second->Id() || material != section->GetMaterial())
    return Invalid(d);
  Vector *a = nullptr, *b = nullptr, *direction = nullptr, *cross = nullptr;
  StatusCode status = Vector::Create(3, &a, d);
  std::unique_ptr<Vector> a_owner(a);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(3, &b, d);
  std::unique_ptr<Vector> b_owner(b);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(3, &direction, d);
  std::unique_ptr<Vector> direction_owner(direction);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(3, &cross, d);
  std::unique_ptr<Vector> cross_owner(cross);
  if (status != StatusCode::kOk) return status;
  for (int i = 0; i < 3; ++i) {
    double value = 0;
    status = first->Coordinate(i, &value, d);
    if (status != StatusCode::kOk) return status;
    status = a->Set(i, value, d);
    if (status != StatusCode::kOk) return status;
    status = second->Coordinate(i, &value, d);
    if (status != StatusCode::kOk) return status;
    status = b->Set(i, value, d);
    if (status != StatusCode::kOk) return status;
  }
  status = b->AddScaled(-1, a, d);
  if (status != StatusCode::kOk) return status;
  double length = 0;
  status = Normalize(b, true, &length, d);
  if (status != StatusCode::kOk) return status;
  status = direction->CopyFrom(section->ApproximateN1(), d);
  if (status != StatusCode::kOk) return status;
  status = Normalize(direction, false, nullptr, d);
  if (status != StatusCode::kOk) return status;
  Matrix* skew = nullptr;
  status = Matrix::Create(3, 3, &skew, d);
  std::unique_ptr<Matrix> skew_owner(skew);
  if (status != StatusCode::kOk) return status;
  // B1 cross product expressed as the skew matrix times the input direction.
  const double* t = b->Data();
  const int row[6] = {0, 0, 1, 1, 2, 2};
  const int col[6] = {1, 2, 0, 2, 0, 1};
  const double coefficient[6] = {-t[2], t[1], t[2], -t[0], -t[1], t[0]};
  for (int i = 0; i < 6; ++i) {
    status = skew->Set(row[i], col[i], coefficient[i], d);
    if (status != StatusCode::kOk) return status;
  }
  status = Matrix::MultiplyVector(skew, false, direction, cross, d);
  if (status != StatusCode::kOk) return status;
  double sine = 0;
  status = cross->Norm2(&sine, d);
  if (status != StatusCode::kOk) return status;
  if (sine <= std::sqrt(std::numeric_limits<double>::epsilon()))
    return Invalid(d);
  const double young = section->GetMaterial()->young();
  const double shear =
      Coefficient(1, young, 1, 2 * (1 + section->GetMaterial()->poisson()), 1);
  // B2/B3 scalar definitions, checked before any solver/HDF5 construction.
  const double coefficients[] = {
      shear,
      Coefficient(1, young, section->Area(), length, 1),
      Coefficient(1, shear, section->TorsionJ(), length, 1),
      Coefficient(12, young, section->I11(), length, 3),
      Coefficient(12, young, section->I22(), length, 3),
      Coefficient(6, young, section->I11(), length, 2),
      Coefficient(6, young, section->I22(), length, 2),
      Coefficient(4, young, section->I11(), length, 1),
      Coefficient(4, young, section->I22(), length, 1),
      Coefficient(2, young, section->I11(), length, 1),
      Coefficient(2, young, section->I22(), length, 1)};
  for (double coefficient_value : coefficients)
    if (!std::isfinite(coefficient_value) || coefficient_value <= 0)
      return Invalid(d);
  status = cross->Scale(1. / sine, d);  // B1 n2.
  if (status != StatusCode::kOk) return status;
  // skew(t)^T n2 = n2 cross t = n1, using the common MKL operation.
  status = Matrix::MultiplyVector(skew, true, cross, a, d);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<Beam3D> candidate(
      new (std::nothrow) Beam3D(id, first, second, section));
  if (!candidate) return model_storage::AllocationFailure("Beam3D::Create", d);
  candidate->length_ = length;
  status = Matrix::Create(3, 3, &candidate->local_axes_, d);
  if (status != StatusCode::kOk) return status;
  // B1 coefficient placement: rows are t, n1, n2, not a matrix copy.
  for (int column = 0; column < 3; ++column) {
    status = candidate->local_axes_->Set(0, column, b->Data()[column], d);
    if (status != StatusCode::kOk) return status;
    status = candidate->local_axes_->Set(1, column, a->Data()[column], d);
    if (status != StatusCode::kOk) return status;
    status = candidate->local_axes_->Set(2, column, cross->Data()[column], d);
    if (status != StatusCode::kOk) return status;
  }
  status = candidate->CreateResults(d);
  if (status != StatusCode::kOk) return status;
  *out = candidate.release();
  return StatusCode::kOk;
}

StatusCode Beam3D::Prepare(ElementWorkspace* workspace, Diagnostic* d) const {
  if (!workspace) return InvalidCall("Beam3D workspace", d);
  const BeamSection* section = Section();
  const double young = section->GetMaterial()->young();
  const double shear =
      Coefficient(1, young, 1, 2 * (1 + section->GetMaterial()->poisson()), 1);
  const int axial_indices[] = {0, 6};
  const int torsion_indices[] = {3, 9};
  const double axial = Coefficient(1, young, section->Area(), length_, 1);
  const double torsion = Coefficient(1, shear, section->TorsionJ(), length_, 1);
  for (int row = 0; row < 2; ++row) {
    for (int column = 0; column < 2; ++column) {
      const double sign = row == column ? 1. : -1.;
      StatusCode status = workspace->local_k_->Set(
          axial_indices[row], axial_indices[column], sign * axial, d);
      if (status != StatusCode::kOk) return status;
      status = workspace->local_k_->Set(
          torsion_indices[row], torsion_indices[column], sign * torsion, d);
      if (status != StatusCode::kOk) return status;
    }
  }
  // B3 scalar coefficient tables, including theta_y=-w' signs. Only index
  // placement occurs in these loops; no matrix product is expanded by hand.
  const int bending_indices[2][4] = {{1, 5, 7, 11}, {2, 4, 8, 10}};
  for (int plane = 0; plane < 2; ++plane) {
    const double inertia = plane == 0 ? section->I22() : section->I11();
    const double a = Coefficient(12, young, inertia, length_, 3);
    const double b =
        Coefficient(plane == 0 ? 6 : -6, young, inertia, length_, 2);
    const double c = Coefficient(4, young, inertia, length_, 1);
    const double e = Coefficient(2, young, inertia, length_, 1);
    const double coefficients[16] = {a,  b,  -a, b,  b, c, -b, e,
                                     -a, -b, a,  -b, b, e, -b, c};
    for (int row = 0; row < 4; ++row) {
      for (int column = 0; column < 4; ++column) {
        const StatusCode status = workspace->local_k_->Set(
            bending_indices[plane][row], bending_indices[plane][column],
            coefficients[row * 4 + column], d);
        if (status != StatusCode::kOk) return status;
      }
    }
  }
  // B4 T=diag(R,R,R,R). All other entries remain factory-initialized zeros.
  for (int block = 0; block < 4; ++block) {
    for (int row = 0; row < 3; ++row) {
      for (int column = 0; column < 3; ++column) {
        const StatusCode status = workspace->transform_->Set(
            block * 3 + row, block * 3 + column,
            local_axes_->Data()[row * 3 + column], d);
        if (status != StatusCode::kOk) return status;
      }
    }
  }
  return StatusCode::kOk;
}

StatusCode Beam3D::ComputeStiffness(ElementWorkspace* workspace,
                                    Matrix* global_k, Diagnostic* d) const {
  if (!workspace || !global_k || global_k->Rows() != 12 ||
      global_k->Columns() != 12)
    return InvalidCall("Beam3D::ComputeStiffness", d);
  StatusCode status = Prepare(workspace, d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Multiply(workspace->local_k_, false, workspace->transform_,
                            false, workspace->product_, d);
  if (status != StatusCode::kOk) return status;
  return Matrix::Multiply(workspace->transform_, true, workspace->product_,
                          false, global_k, d);
}

StatusCode Beam3D::EvaluateLocal(const Vector* global_u,
                                 ElementWorkspace* workspace,
                                 Diagnostic* d) const {
  if (!workspace || !global_u || global_u->Size() != 12)
    return InvalidCall("Beam3D candidate requires 12 global DOFs", d);
  StatusCode status = Prepare(workspace, d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::MultiplyVector(workspace->transform_, false, global_u,
                                  workspace->local_u_, d);
  if (status != StatusCode::kOk) return status;
  return Matrix::MultiplyVector(workspace->local_k_, false, workspace->local_u_,
                                workspace->local_force_, d);
}

StatusCode Beam3D::ComputeInternalForce(const Vector* global_u,
                                        ElementWorkspace* workspace,
                                        Vector* global_force,
                                        Diagnostic* d) const {
  if (!global_force || global_force == global_u || global_force->Size() != 12)
    return InvalidCall("Beam3D::ComputeInternalForce", d);
  StatusCode status = EvaluateLocal(global_u, workspace, d);
  if (status != StatusCode::kOk) return status;
  return Matrix::MultiplyVector(workspace->transform_, true,
                                workspace->local_force_, global_force, d);
}

StatusCode Beam3D::UpdateResults(const Vector* global_u,
                                 ElementWorkspace* workspace, Diagnostic* d) {
  const auto status = StoreLocalForce(global_u, workspace, d);
  return status == StatusCode::kOk ? RecoverEndForces(workspace, d) : status;
}

StatusCode Beam3D::StoreLocalForce(const Vector* global_u,
                                   ElementWorkspace* workspace, Diagnostic* d) {
  StatusCode status = EvaluateLocal(global_u, workspace, d);
  if (status != StatusCode::kOk) return status;
  return MutableInternalForce()->CopyFrom(workspace->local_force_, d);
}

StatusCode Beam3D::RecoverEndForces(ElementWorkspace* workspace,
                                    Diagnostic* d) {
  if (!workspace) return InvalidCall("Beam3D::RecoverEndForces", d);
  // B5 linear map: signed end force plus axial-only stress and strain.
  // Avoid intermediate overflow/underflow when 1/(E*A) is representable.
  const double inverse_area = ReciprocalProduct(1, Section()->Area());
  const double inverse_ea =
      ReciprocalProduct(Section()->GetMaterial()->young(), Section()->Area());
  if (!std::isfinite(inverse_area) || !std::isfinite(inverse_ea) ||
      inverse_area <= 0 || inverse_ea <= 0)
    return ReportFailure(StatusCode::kLibraryFailure,
                         DiagnosticStage::kComputation,
                         "Beam3D B5 unrepresentable recovery coefficients", d);
  for (int index = 0; index < 12; ++index) {
    const StatusCode status =
        workspace->recovery_->Set(index, index, index < 6 ? -1. : 1., d);
    if (status != StatusCode::kOk) return status;
  }
  for (int end = 0; end < 2; ++end) {
    const double sign = end == 0 ? -1. : 1.;
    StatusCode status =
        workspace->recovery_->Set(12 + end, 6 * end, sign * inverse_area, d);
    if (status != StatusCode::kOk) return status;
    status = workspace->recovery_->Set(14 + end, 6 * end, sign * inverse_ea, d);
    if (status != StatusCode::kOk) return status;
  }
  StatusCode status = Matrix::MultiplyVector(
      workspace->recovery_, false, InternalForce(), workspace->end_results_, d);
  if (status != StatusCode::kOk) return status;
  return MutableEndResults()->CopyFrom(workspace->end_results_, d);
}

StatusCode Beam3D::ComputeLoadedInternalForce(
    const Vector* global_u, const Vector* element_load,
    ElementWorkspace* workspace, Vector* global_force, Diagnostic* d) const {
  if (!element_load)
    return ComputeInternalForce(global_u, workspace, global_force, d);
  if (element_load->Size() != 12 || element_load == global_force ||
      !std::all_of(element_load->Data(), element_load->Data() + 12,
                   [](double value) { return std::isfinite(value); }))
    return InvalidCall("Beam3D loaded global force input/alias", d);
  const StatusCode status =
      ComputeInternalForce(global_u, workspace, global_force, d);
  return status == StatusCode::kOk
             ? global_force->AddScaled(-1, element_load, d) : status;
}

StatusCode Beam3D::UpdateLoadedResults(const Vector* global_u,
                                      const Vector* local_element_load,
                                      ElementWorkspace* workspace,
                                      Diagnostic* d) {
  if (!local_element_load) return UpdateResults(global_u, workspace, d);
  if (local_element_load->Size() != 12 ||
      local_element_load == InternalForce() ||
      !std::all_of(local_element_load->Data(), local_element_load->Data() + 12,
                   [](double value) { return std::isfinite(value); }))
    return InvalidCall("Beam3D loaded local results input/alias", d);
  auto status = StoreLocalForce(global_u, workspace, d);
  if (status != StatusCode::kOk) return status;
  status = MutableInternalForce()->AddScaled(-1, local_element_load, d);
  return status == StatusCode::kOk ? RecoverEndForces(workspace, d) : status;
}
}  // namespace fesa
