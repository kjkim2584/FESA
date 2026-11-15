#include "elements/beam3d_nonlinear_ops.h"

#include <cmath>
#include <memory>

#include "fesa/math/rotation.h"

namespace fesa {
namespace {

#define FESA_BEAM_TRY(expression)                 \
  do {                                            \
    const StatusCode status = (expression);       \
    if (status != StatusCode::kOk) return status; \
  } while (false)

bool FiniteVector(const Vector* value, int size) {
  if (!value || value->Size() != size) return false;
  for (int i = 0; i < size; ++i)
    if (!std::isfinite(value->Data()[i])) return false;
  return true;
}

bool FiniteMatrix(const Matrix* value, int rows, int columns) {
  if (!value || value->Rows() != rows || value->Columns() != columns)
    return false;
  for (int i = 0; i < rows * columns; ++i)
    if (!std::isfinite(value->Data()[i])) return false;
  return true;
}

StatusCode Invalid(const char* operation, Diagnostic* d) {
  return ReportFailure(StatusCode::kInvalidArgument,
                       DiagnosticStage::kValidation, operation, d);
}

StatusCode CreateMatrix(int rows, int columns, std::unique_ptr<Matrix>* owner,
                        Diagnostic* d) {
  Matrix* raw = nullptr;
  const auto status = Matrix::Create(rows, columns, &raw, d);
  owner->reset(raw);
  return status;
}

StatusCode CreateVector(int size, std::unique_ptr<Vector>* owner,
                        Diagnostic* d) {
  Vector* raw = nullptr;
  const auto status = Vector::Create(size, &raw, d);
  owner->reset(raw);
  return status;
}

StatusCode Validate(const Beam3D* beam, const ElementTrialState* state,
                    const Vector* eta, const BeamKinematicWorkspace* workspace,
                    ElementEvaluationRequest request, const Vector* force,
                    const Matrix* tangent, TrialEvaluationFailure* failure,
                    Diagnostic* d) {
  if (!beam || !state || !workspace || !failure)
    return Invalid("BeamNonlinear required inputs", d);
  if (request != ElementEvaluationRequest::kInternalForce &&
      request != ElementEvaluationRequest::kTangent &&
      request != ElementEvaluationRequest::kBoth &&
      request != ElementEvaluationRequest::kMaterialTangent)
    return Invalid("BeamNonlinear request", d);
  if (((request != ElementEvaluationRequest::kTangent &&
        request != ElementEvaluationRequest::kMaterialTangent)
           ? !force || force->Size() != 14
           : force != nullptr) ||
      (request != ElementEvaluationRequest::kInternalForce
           ? !tangent || tangent->Rows() != 14 || tangent->Columns() != 14
           : tangent != nullptr))
    return Invalid("BeamNonlinear output shape/request", d);
  if (!FiniteVector(state->translations, 6) ||
      !FiniteMatrix(state->orientations, 2, 4) || !FiniteVector(eta, 2) ||
      (state->committed_orientations &&
       !FiniteMatrix(state->committed_orientations, 2, 4)) ||
      (state->body_spin && !FiniteVector(state->body_spin, 6)))
    return Invalid("BeamNonlinear finite pose/history shape", d);
  if (force && (force == eta || force == state->translations ||
                force == state->body_spin || force == beam->InternalForce() ||
                force == beam->EndResults()))
    return Invalid("BeamNonlinear force alias", d);
  if (tangent && (tangent == state->orientations ||
                  tangent == state->committed_orientations ||
                  tangent == beam->LocalAxes()))
    return Invalid("BeamNonlinear tangent alias", d);
  // A workspace getter is borrowed input storage and expires on Evaluate.
  constexpr BeamKinematicField kFields[] = {
      BeamKinematicField::kPosition, BeamKinematicField::kTangent,
      BeamKinematicField::kNormal1,  BeamKinematicField::kNormal2,
      BeamKinematicField::kStrain,   BeamKinematicField::kStretch};
  for (const auto field : kFields) {
    if (tangent && (tangent == workspace->Value(field) ||
                    tangent == workspace->Variation(field) ||
                    tangent == workspace->Directional(field) ||
                    tangent == workspace->Mixed(field)))
      return Invalid("BeamNonlinear workspace output alias", d);
  }
  if (state->incremental_solid)
    return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                         "BeamNonlinear incremental solid context", d);
  return StatusCode::kOk;
}

}  // namespace

StatusCode PrepareBeamNonlinearPose(const Beam3D* beam,
                                    const ElementTrialState* state,
                                    Matrix* positions, Matrix* tangents,
                                    Matrix* normals, Diagnostic* d,
                                    Matrix* relative_positions) {
  if (!beam || !state || !FiniteVector(state->translations, 6) ||
      !FiniteMatrix(state->orientations, 2, 4) || !positions || !tangents ||
      !normals || positions == tangents || positions == normals ||
      tangents == normals || positions->Rows() != 2 ||
      positions->Columns() != 3 || tangents->Rows() != 2 ||
      tangents->Columns() != 3 || normals->Rows() != 2 ||
      normals->Columns() != 3 || positions == state->orientations ||
      tangents == state->orientations || normals == state->orientations ||
      positions == state->committed_orientations ||
      tangents == state->committed_orientations ||
      normals == state->committed_orientations ||
      positions == beam->LocalAxes() || tangents == beam->LocalAxes() ||
      normals == beam->LocalAxes())
    return Invalid("BeamNonlinear pose preparation inputs", d);
  if (relative_positions &&
      (relative_positions->Rows() != 2 || relative_positions->Columns() != 3 ||
       relative_positions == positions || relative_positions == tangents ||
       relative_positions == normals ||
       relative_positions == state->orientations ||
       relative_positions == state->committed_orientations ||
       relative_positions == beam->LocalAxes()))
    return Invalid("BeamNonlinear relative pose output", d);
  std::unique_ptr<Vector> x, row, q;
  std::unique_ptr<Matrix> rotation, axes;
  FESA_BEAM_TRY(CreateVector(6, &x, d));
  FESA_BEAM_TRY(CreateVector(3, &row, d));
  FESA_BEAM_TRY(CreateVector(4, &q, d));
  FESA_BEAM_TRY(CreateMatrix(3, 3, &rotation, d));
  FESA_BEAM_TRY(CreateMatrix(3, 3, &axes, d));
  for (int end = 0; end < 2; ++end) {
    for (int axis = 0; axis < 3; ++axis) {
      double coordinate = 0;
      FESA_BEAM_TRY(beam->GetNode(end)->Coordinate(axis, &coordinate, d));
      FESA_BEAM_TRY(x->Set(3 * end + axis, coordinate, d));
    }
  }
  std::unique_ptr<Vector> relative, origin;
  if (relative_positions) {
    FESA_BEAM_TRY(CreateVector(6, &relative, d));
    FESA_BEAM_TRY(CreateVector(6, &origin, d));
    for (int axis = 0; axis < 3; ++axis) {
      FESA_BEAM_TRY(origin->Set(axis, x->Data()[axis], d));
      FESA_BEAM_TRY(origin->Set(3 + axis, x->Data()[axis], d));
    }
    FESA_BEAM_TRY(relative->CopyFrom(x.get(), d));
    FESA_BEAM_TRY(relative->AddScaled(-1, origin.get(), d));
    FESA_BEAM_TRY(relative->AddScaled(1, state->translations, d));
  }
  FESA_BEAM_TRY(x->AddScaled(1, state->translations, d));
  for (int end = 0; end < 2; ++end) {
    const int xyz[] = {3 * end, 3 * end + 1, 3 * end + 2};
    const int wxyz[] = {4 * end, 4 * end + 1, 4 * end + 2, 4 * end + 3};
    FESA_BEAM_TRY(row->GatherFrom(x.get(), xyz, 3, 0, d));
    FESA_BEAM_TRY(positions->SetRow(end, row.get(), d));
    if (relative_positions) {
      FESA_BEAM_TRY(row->GatherFrom(relative.get(), xyz, 3, 0, d));
      FESA_BEAM_TRY(relative_positions->SetRow(end, row.get(), d));
    }
    FESA_BEAM_TRY(state->orientations->GatherTo(wxyz, 4, q.get(), 0, d));
    FESA_BEAM_TRY(RotationMatrix(q.get(), rotation.get(), d));
    FESA_BEAM_TRY(Matrix::Multiply(beam->LocalAxes(), false, rotation.get(),
                                   true, axes.get(), d));
    FESA_BEAM_TRY(tangents->CopyBlockFrom(axes.get(), 0, 0, 1, 3, end, 0, d));
    FESA_BEAM_TRY(normals->CopyBlockFrom(axes.get(), 1, 0, 1, 3, end, 0, d));
  }
  return StatusCode::kOk;
}

StatusCode EvaluateBeamNonlinearResponse(
    const Beam3D* beam, const ElementTrialState* state,
    const Vector* internal_position, BeamKinematicWorkspace* workspace,
    ElementEvaluationRequest request, Vector* force, Matrix* tangent,
    TrialEvaluationFailure* failure, Diagnostic* d) {
  if (failure) *failure = TrialEvaluationFailure::kNone;
  FESA_BEAM_TRY(Validate(beam, state, internal_position, workspace, request,
                         force, tangent, failure, d));
  const BeamSection* section = beam->Section();
  const ElasticMaterial* material = section->GetMaterial();
  const double young = material->young();
  const double shear = young / (2 * (1 + material->poisson()));
  const double coefficients[] = {young * section->Area(),
                                 young * section->I11(), young * section->I22(),
                                 shear * section->TorsionJ()};
  for (double coefficient : coefficients) {
    if (!std::isfinite(coefficient) || coefficient <= 0)
      return ReportFailure(StatusCode::kInputError,
                           DiagnosticStage::kComputation,
                           "BeamNonlinear section coefficient range", d);
  }
  std::unique_ptr<Matrix> positions, relative_positions, tangents, normals;
  std::unique_ptr<Matrix> c, strain, b, stress;
  std::unique_ptr<Matrix> cb, kg, column, integrated_force, integrated_tangent;
  std::unique_ptr<Vector> direction;
  FESA_BEAM_TRY(CreateMatrix(2, 3, &positions, d));
  FESA_BEAM_TRY(CreateMatrix(2, 3, &relative_positions, d));
  FESA_BEAM_TRY(CreateMatrix(2, 3, &tangents, d));
  FESA_BEAM_TRY(CreateMatrix(2, 3, &normals, d));
  FESA_BEAM_TRY(PrepareBeamNonlinearPose(
      beam, state, positions.get(), tangents.get(), normals.get(), d,
      relative_positions.get()));
  FESA_BEAM_TRY(CreateMatrix(4, 4, &c, d));
  FESA_BEAM_TRY(CreateMatrix(4, 1, &strain, d));
  FESA_BEAM_TRY(CreateMatrix(4, 14, &b, d));
  FESA_BEAM_TRY(CreateMatrix(4, 1, &stress, d));
  FESA_BEAM_TRY(CreateVector(14, &direction, d));
  if (force) FESA_BEAM_TRY(CreateMatrix(14, 1, &integrated_force, d));
  if (tangent) {
    FESA_BEAM_TRY(CreateMatrix(4, 14, &cb, d));
    if (request != ElementEvaluationRequest::kMaterialTangent) {
      FESA_BEAM_TRY(CreateMatrix(14, 14, &kg, d));
      FESA_BEAM_TRY(CreateMatrix(14, 1, &column, d));
    }
    FESA_BEAM_TRY(CreateMatrix(14, 14, &integrated_tangent, d));
  }
  BeamKinematicInput input;
  FESA_BEAM_TRY(beam->CharacteristicLength(&input.length, d));
  input.positions = positions.get();
  input.relative_positions = relative_positions.get();
  input.tangents = tangents.get();
  input.normals = normals.get();
  input.internal_position = internal_position;
  input.direction = direction.get();
  // The same selected energy defines both residual and tangent: axial Gauss3,
  // bending/twist Gauss4. Map each rule once by dS=L/2 dxi; see P3 contract.
  const double inner = std::sqrt((3 - 2 * std::sqrt(6. / 5)) / 7);
  const double outer = std::sqrt((3 + 2 * std::sqrt(6. / 5)) / 7);
  const double nodes[] = {-std::sqrt(3. / 5), 0, std::sqrt(3. / 5),
                          -outer, -inner, inner, outer};
  const double weights[] = {
      5. / 9, 8. / 9, 5. / 9,
      (18 - std::sqrt(30.)) / 36, (18 + std::sqrt(30.)) / 36,
      (18 + std::sqrt(30.)) / 36, (18 - std::sqrt(30.)) / 36};
  for (int point = 0; point < 7; ++point) {
    const double xi = (nodes[point] + 1) / 2;
    const double measure = input.length * weights[point] / 2;
    const bool axial = point < 3;
    for (int i = 0; i < 4; ++i) {
      const bool active = (i == 0) == axial;
      FESA_BEAM_TRY(c->Set(i, i, active ? coefficients[i] : 0, d));
    }
    FESA_BEAM_TRY(direction->Scale(0, d));
    FESA_BEAM_TRY(workspace->Evaluate(&input, xi, failure, d,
                                      BeamKinematicDerivative::kFirst));
    // Getter storage is overwritten by every directional evaluation.
    FESA_BEAM_TRY(
        strain->CopyFrom(workspace->Value(BeamKinematicField::kStrain), d));
    FESA_BEAM_TRY(
        b->CopyFrom(workspace->Variation(BeamKinematicField::kStrain), d));
    FESA_BEAM_TRY(
        Matrix::Multiply(c.get(), false, strain.get(), false, stress.get(), d));
    if (force)
      FESA_BEAM_TRY(Matrix::MultiplyAdd(b.get(), true, stress.get(), false,
                                        measure, 1, integrated_force.get(), d));
    if (!tangent) continue;
    FESA_BEAM_TRY(
        Matrix::Multiply(c.get(), false, b.get(), false, cb.get(), d));
    FESA_BEAM_TRY(Matrix::MultiplyAdd(b.get(), true, cb.get(), false, measure,
                                      1, integrated_tangent.get(), d));
    // Approximate first correction only; the exact residual and default
    // consistent tangent retain the selected energy and all quadrature points.
    if (request == ElementEvaluationRequest::kMaterialTangent) continue;
    for (int j = 0; j < 14; ++j) {
      if (j > 0) FESA_BEAM_TRY(direction->Set(j - 1, 0, d));
      FESA_BEAM_TRY(direction->Set(j, 1, d));
      FESA_BEAM_TRY(workspace->Evaluate(&input, xi, failure, d));
      FESA_BEAM_TRY(
          Matrix::Multiply(workspace->Mixed(BeamKinematicField::kStrain), true,
                           stress.get(), false, column.get(), d));
      FESA_BEAM_TRY(kg->CopyBlockFrom(column.get(), 0, 0, 14, 1, 0, j, d));
    }
    FESA_BEAM_TRY(integrated_tangent->AddScaled(measure, kg.get(), d));
  }
  if (force) {
    const int indices[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
    FESA_BEAM_TRY(integrated_force->GatherTo(indices, 14, force, 0, d));
  }
  if (tangent) FESA_BEAM_TRY(tangent->CopyFrom(integrated_tangent.get(), d));
  return StatusCode::kOk;
}

#undef FESA_BEAM_TRY
}  // namespace fesa
