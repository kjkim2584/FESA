#include "elements/mitc4_trial.h"

#include <cmath>
#include <memory>

#include "elements/mitc4_geometry.h"
#include "elements/mitc4_linear_ops.h"
#include "elements/mitc4_nonlinear_ops.h"
#include "fesa/math/rotation.h"
#include "model/model_storage.h"

namespace fesa {
namespace {

constexpr int kRawPoint[7] = {0, 0, 0, 1, 2, 3, 4};
constexpr int kRawA[7] = {0, 1, 0, 0, 0, 1, 1};
constexpr int kRawB[7] = {0, 1, 1, 2, 2, 2, 2};

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

StatusCode NumericalFailure(const char* operation, Diagnostic* d) {
  return ReportFailure(StatusCode::kLibraryFailure,
                       DiagnosticStage::kComputation, operation, d);
}

StatusCode Skew(const Vector* vector, Matrix* out, Diagnostic* d) {
  const double* v = vector->Data();
  const double coefficients[9] = {0.0,   -v[2], v[1], v[2], 0.0,
                                  -v[0], -v[1], v[0], 0.0};
  for (int i = 0; i < 9; ++i) {
    const auto status = out->Set(i / 3, i % 3, coefficients[i], d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}

StatusCode MatrixRow(const Matrix* source, int row, Vector* out,
                     Diagnostic* d) {
  int indices[4] = {};
  for (int i = 0; i < source->Columns(); ++i)
    indices[i] = source->Columns() * row + i;
  return source->GatherTo(indices, source->Columns(), out, 0, d);
}

StatusCode ValidateInput(const Node* const* nodes, int count,
                         const ShellSection* section,
                         const ElementTrialState* state,
                         ElementEvaluationRequest request, double drilling,
                         Vector* force, Matrix* tangent, Diagnostic* d) {
  const bool want_force = request == ElementEvaluationRequest::kInternalForce ||
                          request == ElementEvaluationRequest::kBoth;
  const bool want_tangent = request == ElementEvaluationRequest::kTangent ||
                            request == ElementEvaluationRequest::kBoth;
  if ((!want_force && !want_tangent) || !nodes || count != 4 || !section ||
      !state || !std::isfinite(drilling) || drilling < 0.0 ||
      (want_force ? !force || force->Size() != 24 : force != nullptr) ||
      (want_tangent
           ? !tangent || tangent->Rows() != 24 || tangent->Columns() != 24
           : tangent != nullptr)) {
    return model_storage::Invalid("MITC4 trial request/shape", d);
  }
  if (!FiniteVector(state->translations, 12) ||
      !FiniteVector(state->body_spin, 12) ||
      !FiniteMatrix(state->orientations, 4, 4) ||
      !FiniteMatrix(state->committed_orientations, 4, 4) ||
      (force && (force == state->translations || force == state->body_spin)) ||
      (tangent && (tangent == state->orientations ||
                   tangent == state->committed_orientations))) {
    return model_storage::Invalid("MITC4 trial state/alias/nonfinite", d);
  }
  // A positive input thickness still needs a representable half-thickness.
  // This is an invalid initial definition, never a trial cutback condition.
  if (!(section->thickness() / 2.0 > 0.0))
    return ReportFailure(StatusCode::kInputError, DiagnosticStage::kValidation,
                         "MITC4 initial half-thickness range", d);
  for (int i = 0; i < 4; ++i) {
    if (!nodes[i] || nodes[i]->dof_count() != 6)
      return model_storage::Invalid("MITC4 trial node components", d);
    for (int j = 0; j < i; ++j)
      if (nodes[i] == nodes[j] || nodes[i]->Id() == nodes[j]->Id())
        return model_storage::Invalid("MITC4 trial duplicate node", d);
  }
  return StatusCode::kOk;
}

StatusCode RelativeCoordinates(Matrix* coordinates, Diagnostic* d) {
  std::unique_ptr<Matrix> anchor;
  std::unique_ptr<Matrix> ones;
  std::unique_ptr<Matrix> translated;
  auto status = CreateMatrix(1, 3, &anchor, d);
  if (status != StatusCode::kOk) return status;
  status = CreateMatrix(4, 1, &ones, d);
  if (status != StatusCode::kOk) return status;
  status = CreateMatrix(4, 3, &translated, d);
  if (status != StatusCode::kOk) return status;
  for (int i = 0; i < 4; ++i) {
    status = ones->Set(i, 0, 1.0, d);
    if (status != StatusCode::kOk) return status;
  }
  status = anchor->CopyBlockFrom(coordinates, 0, 0, 1, 3, 0, 0, d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Multiply(ones.get(), false, anchor.get(), false,
                            translated.get(), d);
  if (status != StatusCode::kOk) return status;
  return coordinates->AddScaled(-1.0, translated.get(), d);
}

StatusCode BuildPositions(const Node* const* nodes, const Vector* u, Vector* y0,
                          Vector* y, Diagnostic* d) {
  std::unique_ptr<Matrix> initial;
  std::unique_ptr<Matrix> current;
  std::unique_ptr<Matrix> displacement;
  std::unique_ptr<Vector> row;
  auto status = CreateMatrix(4, 3, &initial, d);
  if (status != StatusCode::kOk) return status;
  status = CreateMatrix(4, 3, &current, d);
  if (status != StatusCode::kOk) return status;
  status = CreateMatrix(4, 3, &displacement, d);
  if (status != StatusCode::kOk) return status;
  status = CreateVector(3, &row, d);
  if (status != StatusCode::kOk) return status;
  for (int node = 0; node < 4; ++node) {
    for (int axis = 0; axis < 3; ++axis) {
      double coordinate = 0.0;
      status = nodes[node]->Coordinate(axis, &coordinate, d);
      if (status != StatusCode::kOk) return status;
      status = initial->Set(node, axis, coordinate, d);
      if (status != StatusCode::kOk) return status;
    }
    const int indices[3] = {3 * node, 3 * node + 1, 3 * node + 2};
    status = row->GatherFrom(u, indices, 3, 0, d);
    if (status != StatusCode::kOk) return status;
    status = displacement->SetRow(node, row.get(), d);
    if (status != StatusCode::kOk) return status;
  }
  status = RelativeCoordinates(initial.get(), d);
  if (status != StatusCode::kOk) return status;
  status = RelativeCoordinates(displacement.get(), d);
  if (status != StatusCode::kOk) return status;
  status = current->CopyFrom(initial.get(), d);
  if (status != StatusCode::kOk) return status;
  status = current->AddScaled(1.0, displacement.get(), d);
  if (status != StatusCode::kOk) return status;
  // y positions are (Xi-X1)+(ui-u1). Partition of unity removes the
  // common translation without projecting any initial out-of-plane roundoff.
  for (int node = 0; node < 4; ++node) {
    const int indices[3] = {6 * node, 6 * node + 1, 6 * node + 2};
    status = MatrixRow(initial.get(), node, row.get(), d);
    if (status != StatusCode::kOk) return status;
    status = y0->ScatterAdd(row.get(), indices, 3, d);
    if (status != StatusCode::kOk) return status;
    status = MatrixRow(current.get(), node, row.get(), d);
    if (status != StatusCode::kOk) return status;
    status = y->ScatterAdd(row.get(), indices, 3, d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}

StatusCode BuildDirectors(const ElementTrialState* state, const Vector* normal,
                          Vector* y0, Vector* y, Matrix* directors,
                          Vector* angles, Matrix* angle_gradients,
                          TrialEvaluationFailure* failure, Diagnostic* d) {
  std::unique_ptr<Vector> vectors[7];
  const int sizes[7] = {4, 4, 3, 3, 3, 3, 3};
  for (int i = 0; i < 7; ++i) {
    if (!angle_gradients && i >= 5) continue;
    const auto status = CreateVector(sizes[i], &vectors[i], d);
    if (status != StatusCode::kOk) return status;
  }
  std::unique_ptr<Matrix> rotation;
  std::unique_ptr<Matrix> committed_rotation;
  std::unique_ptr<Matrix> inverse_jacobian;
  auto status = CreateMatrix(3, 3, &rotation, d);
  if (status != StatusCode::kOk) return status;
  if (angle_gradients) {
    status = CreateMatrix(3, 3, &committed_rotation, d);
    if (status != StatusCode::kOk) return status;
    status = CreateMatrix(3, 3, &inverse_jacobian, d);
    if (status != StatusCode::kOk) return status;
  }
  Vector* q = vectors[0].get();
  Vector* qc = vectors[1].get();
  Vector* body_increment = vectors[2].get();
  Vector* director = vectors[3].get();
  Vector* spin = vectors[4].get();
  Vector* temporary = vectors[5].get();
  Vector* gradient = vectors[6].get();
  for (int node = 0; node < 4; ++node) {
    status = MatrixRow(state->orientations, node, q, d);
    if (status != StatusCode::kOk) return status;
    status = MatrixRow(state->committed_orientations, node, qc, d);
    if (status != StatusCode::kOk) return status;
    bool branch = false;
    status = BodyRotationIncrement(qc, q, body_increment, &branch, d);
    if (status != StatusCode::kOk) {
      if (branch) *failure = TrialEvaluationFailure::kRotationBranch;
      return status;
    }
    status = RotationMatrix(q, rotation.get(), d);
    if (status != StatusCode::kOk) return status;
    status = Matrix::MultiplyVector(rotation.get(), false, normal, director, d);
    if (status != StatusCode::kOk) return status;
    status = directors->SetRow(node, director, d);
    if (status != StatusCode::kOk) return status;
    const int director_indices[3] = {6 * node + 3, 6 * node + 4, 6 * node + 5};
    status = y0->ScatterAdd(normal, director_indices, 3, d);
    if (status != StatusCode::kOk) return status;
    status = y->ScatterAdd(director, director_indices, 3, d);
    if (status != StatusCode::kOk) return status;
    const int spin_indices[3] = {3 * node, 3 * node + 1, 3 * node + 2};
    status = spin->GatherFrom(state->body_spin, spin_indices, 3, 0, d);
    if (status != StatusCode::kOk) return status;
    double angle = 0.0;
    status = normal->Dot(spin, &angle, d);
    if (status != StatusCode::kOk) return status;
    status = angles->Set(node, angle, d);
    if (status != StatusCode::kOk) return status;
    if (!angle_gradients) continue;
    // gpsi=Rc*JleftInverse(b)^T*n0: derivative of n0 dot (ac+b)
    // under the same spatial Exp(dw) correction used by the physical force.
    status =
        InverseLeftRotationJacobian(body_increment, inverse_jacobian.get(), d);
    if (status != StatusCode::kOk) return status;
    status = RotationMatrix(qc, committed_rotation.get(), d);
    if (status != StatusCode::kOk) return status;
    status = Matrix::MultiplyVector(inverse_jacobian.get(), true, normal,
                                    temporary, d);
    if (status != StatusCode::kOk) return status;
    status = Matrix::MultiplyVector(committed_rotation.get(), false, temporary,
                                    gradient, d);
    if (status != StatusCode::kOk) return status;
    status = angle_gradients->SetRow(node, gradient, d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}

StatusCode FillDerivatives(double r, double s, double zeta, double thickness,
                           Matrix* const* derivatives, Diagnostic* d) {
  constexpr double kR[4] = {-1.0, 1.0, 1.0, -1.0};
  constexpr double kS[4] = {-1.0, -1.0, 1.0, 1.0};
  const double z = thickness * zeta / 2.0;
  for (int node = 0; node < 4; ++node) {
    const double n = (1.0 + kR[node] * r) * (1.0 + kS[node] * s) / 4.0;
    const double nr = kR[node] * (1.0 + kS[node] * s) / 4.0;
    const double ns = kS[node] * (1.0 + kR[node] * r) / 4.0;
    for (int axis = 0; axis < 3; ++axis) {
      auto status = derivatives[0]->Set(axis, 6 * node + axis, nr, d);
      if (status != StatusCode::kOk) return status;
      status = derivatives[1]->Set(axis, 6 * node + axis, ns, d);
      if (status != StatusCode::kOk) return status;
      status = derivatives[0]->Set(axis, 6 * node + 3 + axis, z * nr, d);
      if (status != StatusCode::kOk) return status;
      status = derivatives[1]->Set(axis, 6 * node + 3 + axis, z * ns, d);
      if (status != StatusCode::kOk) return status;
      status = derivatives[2]->Set(axis, 6 * node + 3 + axis,
                                   (thickness / 2.0) * n, d);
      if (status != StatusCode::kOk) return status;
    }
  }
  return StatusCode::kOk;
}

StatusCode CheckCurrentGeometry(Vector* const* initial, Vector* const* current,
                                double* initial_determinant,
                                Mitc4VolumeWorkspace* workspace,
                                TrialEvaluationFailure* failure,
                                Diagnostic* d) {
  double jacobian = 0;
  return workspace->Validate(initial, current, initial_determinant,
                             &jacobian, failure, d);
}

StatusCode CovariantRow(const Matrix* da, const Matrix* db, const Vector* ga,
                        const Vector* gb, const Vector* initial_a,
                        const Vector* initial_b, Vector* by, double* strain,
                        Diagnostic* d) {
  std::unique_ptr<Vector> second;
  auto status = CreateVector(24, &second, d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::MultiplyVector(da, true, gb, by, d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::MultiplyVector(db, true, ga, second.get(), d);
  if (status != StatusCode::kOk) return status;
  status = by->AddScaled(1.0, second.get(), d);
  if (status != StatusCode::kOk) return status;
  status = by->Scale(0.5, d);
  if (status != StatusCode::kOk) return status;
  double metric = 0.0;
  double initial_metric = 0.0;
  status = ga->Dot(gb, &metric, d);
  if (status != StatusCode::kOk) return status;
  status = initial_a->Dot(initial_b, &initial_metric, d);
  if (status != StatusCode::kOk) return status;
  const double value = 0.5 * (metric - initial_metric);
  if (!std::isfinite(value))
    return NumericalFailure("MITC4 trial covariant strain range", d);
  *strain = value;
  return StatusCode::kOk;
}

StatusCode EvaluateCovariantPoints(
    double r, double s, double zeta, double thickness, const Vector* y0,
    const Vector* y, const Matrix* frame_transpose, Matrix* const* derivatives,
    Matrix* p, Vector* raw_strain, Matrix* raw_by, double* weight,
    TrialEvaluationFailure* failure, Diagnostic* d,
    Mitc4VolumeWorkspace* workspace, Matrix* cached_by, Vector* cached_strain,
    bool* cached_ready) {
  std::unique_ptr<Vector> initial_owners[3];
  std::unique_ptr<Vector> current_owners[3];
  Vector* initial[3] = {};
  Vector* current[3] = {};
  for (int axis = 0; axis < 3; ++axis) {
    auto status = CreateVector(3, &initial_owners[axis], d);
    if (status != StatusCode::kOk) return status;
    status = CreateVector(3, &current_owners[axis], d);
    if (status != StatusCode::kOk) return status;
    initial[axis] = initial_owners[axis].get();
    current[axis] = current_owners[axis].get();
  }
  std::unique_ptr<Matrix> basis;
  std::unique_ptr<Vector> by;
  auto status = CreateMatrix(3, 3, &basis, d);
  if (status != StatusCode::kOk) return status;
  status = CreateVector(24, &by, d);
  if (status != StatusCode::kOk) return status;
  const double point_r[5] = {r, 0.0, 0.0, -1.0, 1.0};
  const double point_s[5] = {s, -1.0, 1.0, 0.0, 0.0};
  for (int point = 0; point < 5; ++point) {
    Matrix* const* operators = derivatives + 3 * point;
    status = FillDerivatives(point_r[point], point_s[point], zeta, thickness,
                             operators, d);
    if (status != StatusCode::kOk) return status;
    // Keep FillDerivatives above: the later geometric Hessian still reads it.
    const int cache_slot = (zeta < 0 ? 0 : 4) + point - 1;
    if (point > 0 && cached_ready[cache_slot]) {
      int indices[24] = {};
      for (int i = 0; i < 24; ++i) indices[i] = 24 * cache_slot + i;
      status = cached_by->GatherTo(indices, 24, by.get(), 0, d);
      if (status == StatusCode::kOk)
        status =
            raw_strain->Set(point + 2, cached_strain->Data()[cache_slot], d);
      if (status == StatusCode::kOk)
        status = raw_by->SetRow(point + 2, by.get(), d);
      if (status != StatusCode::kOk) return status;
      continue;
    }
    for (int axis = 0; axis < 3; ++axis) {
      status =
          Matrix::MultiplyVector(operators[axis], false, y0, initial[axis], d);
      if (status != StatusCode::kOk) return status;
      status =
          Matrix::MultiplyVector(operators[axis], false, y, current[axis], d);
      if (status != StatusCode::kOk) return status;
      status = basis->SetRow(axis, initial[axis], d);
      if (status != StatusCode::kOk) return status;
    }
    double determinant = 0.0;
    status = CheckCurrentGeometry(initial, current, &determinant, workspace,
                                   failure, d);
    if (status != StatusCode::kOk) {
      // Diagnostic::Record clears point metadata: retain the actual basis
      // location, including tying points distinct from the parent Gauss point.
      if (d) {
        d->geometry_point_available = true;
        d->natural_xi = point_r[point];
        d->natural_eta = point_s[point];
        d->natural_zeta = zeta;
      }
      return status;
    }
    if (point == 0) {
      *weight = determinant;
      status = BuildMitc4NonlinearCartesianTransform(
          basis.get(), frame_transpose, p, d);
      if (status != StatusCode::kOk) return status;
    }
    for (int row = 0; row < 7; ++row) {
      if (kRawPoint[row] != point) continue;
      const int a = kRawA[row];
      const int b = kRawB[row];
      double strain = 0.0;
      status = CovariantRow(operators[a], operators[b], current[a], current[b],
                            initial[a], initial[b], by.get(), &strain, d);
      if (status != StatusCode::kOk) return status;
      status = raw_strain->Set(row, strain, d);
      if (status != StatusCode::kOk) return status;
      status = raw_by->SetRow(row, by.get(), d);
      if (status != StatusCode::kOk) return status;
      if (point > 0) {
        status = cached_strain->Set(cache_slot, strain, d);
        if (status == StatusCode::kOk)
          status = cached_by->SetRow(cache_slot, by.get(), d);
        if (status != StatusCode::kOk) return status;
        cached_ready[cache_slot] = true;
      }
    }
  }
  return StatusCode::kOk;
}

StatusCode Integrate(const Vector* y0, const Vector* y, const Matrix* frame,
                     const ShellSection* section, Vector* py, Matrix* hy,
                     TrialEvaluationFailure* failure, Diagnostic* d) {
  std::unique_ptr<Matrix> derivative_owners[15];
  Matrix* derivatives[15] = {};
  for (int i = 0; i < 15; ++i) {
    const auto status = CreateMatrix(3, 24, &derivative_owners[i], d);
    if (status != StatusCode::kOk) return status;
    derivatives[i] = derivative_owners[i].get();
  }
  std::unique_ptr<Matrix> matrices[9];
  const int rows[9] = {3, 3, 5, 5, 7, 5, 5, 5, 5};
  const int columns[9] = {3, 3, 5, 7, 24, 24, 24, 5, 24};
  // No matrix in the force-only integration path has 24x24 entries.
  for (int i = 0; i < 9; ++i) {
    if (!hy && i >= 8) continue;
    const auto status = CreateMatrix(rows[i], columns[i], &matrices[i], d);
    if (status != StatusCode::kOk) return status;
  }
  Matrix* identity = matrices[0].get();
  Matrix* frame_transpose = matrices[1].get();
  Matrix* p = matrices[2].get();
  Matrix* selector = matrices[3].get();
  Matrix* raw_by = matrices[4].get();
  Matrix* covariant_by = matrices[5].get();
  Matrix* by = matrices[6].get();
  Matrix* constitutive = matrices[7].get();
  Matrix* d_by = matrices[8].get();
  std::unique_ptr<Vector> vectors[7];
  const int sizes[7] = {7, 5, 5, 5, 24, 5, 7};
  for (int i = 0; i < 7; ++i) {
    if (!hy && i >= 5) continue;
    const auto status = CreateVector(sizes[i], &vectors[i], d);
    if (status != StatusCode::kOk) return status;
  }
  Vector* raw_strain = vectors[0].get();
  Vector* covariant_strain = vectors[1].get();
  Vector* strain = vectors[2].get();
  Vector* stress = vectors[3].get();
  Vector* point_force = vectors[4].get();
  Vector* covariant_stress = vectors[5].get();
  Vector* raw_stress = vectors[6].get();
  StatusCode status = StatusCode::kOk;
  for (int i = 0; i < 3; ++i) {
    status = identity->Set(i, i, 1.0, d);
    if (status != StatusCode::kOk) return status;
  }
  status = Matrix::Multiply(frame, true, identity, false, frame_transpose, d);
  if (status != StatusCode::kOk) return status;
  status = BuildMitc4LinearConstitutive(section->material(), constitutive, d);
  if (status != StatusCode::kOk) return status;
  // Fixed four tying points at each of two zeta layers. y0/y/frame/thickness
  // are invariant during this call; destroy cache on success or failure.
  // Store raw strain and its gradient with respect to 24 position/director
  // components, before selector/P. Never reuse across states or workers.
  std::unique_ptr<Matrix> cached_by;
  std::unique_ptr<Vector> cached_strain;
  bool cached_ready[8] = {};
  status = CreateMatrix(8, 24, &cached_by, d);
  if (status == StatusCode::kOk)
    status = CreateVector(8, &cached_strain, d);
  if (status != StatusCode::kOk) return status;
  Mitc4VolumeWorkspace volume_workspace;
  status = volume_workspace.Initialize(d);
  if (status != StatusCode::kOk) return status;
  const double gauss = 1.0 / std::sqrt(3.0);
  for (double r : {-gauss, gauss}) {
    for (double s : {-gauss, gauss}) {
      status = BuildMitc4NonlinearTyingSelector(r, s, selector, d);
      if (status != StatusCode::kOk) return status;
      for (double zeta : {-gauss, gauss}) {
        double weight = 0.0;
        status = EvaluateCovariantPoints(
            r, s, zeta, section->thickness(), y0, y, frame_transpose,
            derivatives, p, raw_strain, raw_by, &weight, failure, d,
            &volume_workspace, cached_by.get(), cached_strain.get(),
            cached_ready);
        if (status != StatusCode::kOk) return status;
        status = Matrix::MultiplyVector(selector, false, raw_strain,
                                        covariant_strain, d);
        if (status != StatusCode::kOk) return status;
        status = Matrix::MultiplyVector(p, false, covariant_strain, strain, d);
        if (status != StatusCode::kOk) return status;
        status =
            Matrix::Multiply(selector, false, raw_by, false, covariant_by, d);
        if (status != StatusCode::kOk) return status;
        status = Matrix::Multiply(p, false, covariant_by, false, by, d);
        if (status != StatusCode::kOk) return status;
        status = Matrix::MultiplyVector(constitutive, false, strain, stress, d);
        if (status != StatusCode::kOk) return status;
        status = Matrix::MultiplyVector(by, true, stress, point_force, d);
        if (status != StatusCode::kOk) return status;
        status = py->AddScaled(weight, point_force, d);
        if (status != StatusCode::kOk) return status;
        if (!hy) continue;
        status = Matrix::Multiply(constitutive, false, by, false, d_by, d);
        if (status != StatusCode::kOk) return status;
        status = Matrix::MultiplyAdd(by, true, d_by, false, weight, 1.0, hy, d);
        if (status != StatusCode::kOk) return status;
        // tau_raw=Selector^T P^T s. This contracts the same tying/P maps
        // into the covariant Hessians without storing five A_alpha arrays.
        status = Matrix::MultiplyVector(p, true, stress, covariant_stress, d);
        if (status != StatusCode::kOk) return status;
        status = Matrix::MultiplyVector(selector, true, covariant_stress,
                                        raw_stress, d);
        if (status != StatusCode::kOk) return status;
        for (int row = 0; row < 7; ++row) {
          Matrix* da = derivatives[3 * kRawPoint[row] + kRawA[row]];
          Matrix* db = derivatives[3 * kRawPoint[row] + kRawB[row]];
          const double coefficient = weight * raw_stress->Data()[row];
          if (!std::isfinite(coefficient))
            return NumericalFailure("MITC4 trial geometric weight range", d);
          if (da == db) {
            status = Matrix::MultiplyAdd(da, true, db, false, coefficient, 1.0,
                                         hy, d);
          } else {
            status = Matrix::MultiplyAdd(da, true, db, false, 0.5 * coefficient,
                                         1.0, hy, d);
            if (status == StatusCode::kOk)
              status = Matrix::MultiplyAdd(db, true, da, false,
                                           0.5 * coefficient, 1.0, hy, d);
          }
          if (status != StatusCode::kOk) return status;
        }
      }
    }
  }
  return StatusCode::kOk;
}

StatusCode AssembleForce(const Vector* py, const Matrix* directors,
                         const Vector* angles, double drilling, Vector* out,
                         Diagnostic* d) {
  std::unique_ptr<Vector> part;
  std::unique_ptr<Vector> director;
  std::unique_ptr<Vector> moment;
  auto status = CreateVector(3, &part, d);
  if (status != StatusCode::kOk) return status;
  status = CreateVector(3, &director, d);
  if (status != StatusCode::kOk) return status;
  status = CreateVector(3, &moment, d);
  if (status != StatusCode::kOk) return status;
  for (int node = 0; node < 4; ++node) {
    const int position_indices[3] = {6 * node, 6 * node + 1, 6 * node + 2};
    const int rotation_indices[3] = {6 * node + 3, 6 * node + 4, 6 * node + 5};
    status = part->GatherFrom(py, position_indices, 3, 0, d);
    if (status != StatusCode::kOk) return status;
    status = out->ScatterAdd(part.get(), position_indices, 3, d);
    if (status != StatusCode::kOk) return status;
    status = part->GatherFrom(py, rotation_indices, 3, 0, d);
    if (status != StatusCode::kOk) return status;
    status = MatrixRow(directors, node, director.get(), d);
    if (status != StatusCode::kOk) return status;
    status = CrossMitc4Vectors(director.get(), part.get(), moment.get(), d);
    if (status != StatusCode::kOk) return status;
    status =
        moment->AddScaled(drilling * angles->Data()[node], director.get(), d);
    if (status != StatusCode::kOk) return status;
    status = out->ScatterAdd(moment.get(), rotation_indices, 3, d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}

StatusCode AssembleTangent(const Vector* py, const Matrix* hy,
                           const Matrix* directors, const Vector* angles,
                           const Matrix* angle_gradients, double drilling,
                           Matrix* out, Diagnostic* d) {
  std::unique_ptr<Matrix> matrices[9];
  const int rows[9] = {24, 24, 24, 3, 3, 3, 3, 1, 1};
  const int columns[9] = {24, 24, 24, 3, 3, 3, 3, 3, 3};
  for (int i = 0; i < 9; ++i) {
    const auto status = CreateMatrix(rows[i], columns[i], &matrices[i], d);
    if (status != StatusCode::kOk) return status;
  }
  Matrix* l = matrices[0].get();
  Matrix* product = matrices[1].get();
  Matrix* connection = matrices[2].get();
  Matrix* skew_d = matrices[3].get();
  Matrix* skew_pd = matrices[4].get();
  Matrix* block = matrices[5].get();
  Matrix* drilling_block = matrices[6].get();
  Matrix* row_d = matrices[7].get();
  Matrix* row_g = matrices[8].get();
  std::unique_ptr<Vector> vectors[3];
  for (int i = 0; i < 3; ++i) {
    const auto status = CreateVector(3, &vectors[i], d);
    if (status != StatusCode::kOk) return status;
  }
  Vector* director = vectors[0].get();
  Vector* pd = vectors[1].get();
  Vector* gradient = vectors[2].get();
  for (int node = 0; node < 4; ++node) {
    auto status = MatrixRow(directors, node, director, d);
    if (status != StatusCode::kOk) return status;
    status = Skew(director, skew_d, d);
    if (status != StatusCode::kOk) return status;
    for (int axis = 0; axis < 3; ++axis) {
      status = l->Set(6 * node + axis, 6 * node + axis, 1.0, d);
      if (status != StatusCode::kOk) return status;
    }
    status = skew_d->Scale(-1.0, d);
    if (status != StatusCode::kOk) return status;
    status =
        l->CopyBlockFrom(skew_d, 0, 0, 3, 3, 6 * node + 3, 6 * node + 3, d);
    if (status != StatusCode::kOk) return status;
    status = skew_d->Scale(-1.0, d);
    if (status != StatusCode::kOk) return status;
    const int indices[3] = {6 * node + 3, 6 * node + 4, 6 * node + 5};
    status = pd->GatherFrom(py, indices, 3, 0, d);
    if (status != StatusCode::kOk) return status;
    status = Skew(pd, skew_pd, d);
    if (status != StatusCode::kOk) return status;
    // d x pd differentiated under spatial rotation contributes
    // [pd]x[d]x in addition to the mapped y-Hessian.
    status = Matrix::Multiply(skew_pd, false, skew_d, false, block, d);
    if (status != StatusCode::kOk) return status;
    status = MatrixRow(angle_gradients, node, gradient, d);
    if (status != StatusCode::kOk) return status;
    status = row_d->SetRow(0, director, d);
    if (status != StatusCode::kOk) return status;
    status = row_g->SetRow(0, gradient, d);
    if (status != StatusCode::kOk) return status;
    status = Matrix::Multiply(row_d, true, row_g, false, drilling_block, d);
    if (status != StatusCode::kOk) return status;
    status = drilling_block->AddScaled(-angles->Data()[node], skew_d, d);
    if (status != StatusCode::kOk) return status;
    status = block->AddScaled(drilling, drilling_block, d);
    if (status != StatusCode::kOk) return status;
    status = connection->CopyBlockFrom(block, 0, 0, 3, 3, 6 * node + 3,
                                       6 * node + 3, d);
    if (status != StatusCode::kOk) return status;
  }
  auto status = Matrix::Multiply(hy, false, l, false, product, d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Multiply(l, true, product, false, out, d);
  if (status != StatusCode::kOk) return status;
  return out->AddScaled(1.0, connection, d);
}

}  // namespace

StatusCode EvaluateMitc4Trial(const Node* const* nodes, int count,
                              const ShellSection* section,
                              const ElementTrialState* state,
                              ElementEvaluationRequest request,
                              double fixed_drilling, Vector* force,
                              Matrix* tangent, TrialEvaluationFailure* failure,
                              Diagnostic* d) {
  if (!failure) return model_storage::Invalid("MITC4 trial failure flag", d);
  *failure = TrialEvaluationFailure::kNone;
  auto status = ValidateInput(nodes, count, section, state, request,
                              fixed_drilling, force, tangent, d);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<Matrix> local;
  std::unique_ptr<Matrix> frame;
  std::unique_ptr<Matrix> directors;
  std::unique_ptr<Matrix> angle_gradients;
  std::unique_ptr<Matrix> hy;
  status = CreateMatrix(4, 3, &local, d);
  if (status != StatusCode::kOk) return status;
  status = CreateMatrix(3, 3, &frame, d);
  if (status != StatusCode::kOk) return status;
  double diameter = 0.0;
  double length = 0.0;
  status = ValidateMitc4Geometry(nodes, count, local.get(), frame.get(),
                                 &diameter, &length, d);
  if (status != StatusCode::kOk) return status;
  status = CreateMatrix(4, 3, &directors, d);
  if (status != StatusCode::kOk) return status;
  if (tangent) {
    status = CreateMatrix(4, 3, &angle_gradients, d);
    if (status != StatusCode::kOk) return status;
    status = CreateMatrix(24, 24, &hy, d);
    if (status != StatusCode::kOk) return status;
  }
  std::unique_ptr<Vector> normal;
  std::unique_ptr<Vector> y0;
  std::unique_ptr<Vector> y;
  std::unique_ptr<Vector> py;
  std::unique_ptr<Vector> angles;
  status = CreateVector(3, &normal, d);
  if (status != StatusCode::kOk) return status;
  status = CreateVector(24, &y0, d);
  if (status != StatusCode::kOk) return status;
  status = CreateVector(24, &y, d);
  if (status != StatusCode::kOk) return status;
  status = CreateVector(24, &py, d);
  if (status != StatusCode::kOk) return status;
  status = CreateVector(4, &angles, d);
  if (status != StatusCode::kOk) return status;
  status = MatrixRow(frame.get(), 2, normal.get(), d);
  if (status != StatusCode::kOk) return status;
  status = BuildPositions(nodes, state->translations, y0.get(), y.get(), d);
  if (status != StatusCode::kOk) return status;
  status =
      BuildDirectors(state, normal.get(), y0.get(), y.get(), directors.get(),
                     angles.get(), angle_gradients.get(), failure, d);
  if (status != StatusCode::kOk) return status;
  status = Integrate(y0.get(), y.get(), frame.get(), section, py.get(),
                     hy.get(), failure, d);
  if (status != StatusCode::kOk) return status;
  if (tangent) {
    status = AssembleTangent(py.get(), hy.get(), directors.get(), angles.get(),
                             angle_gradients.get(), fixed_drilling, tangent, d);
    if (status != StatusCode::kOk) return status;
  }
  if (!force) return StatusCode::kOk;
  std::unique_ptr<Vector> global_force;
  status = CreateVector(24, &global_force, d);
  if (status != StatusCode::kOk) return status;
  status = AssembleForce(py.get(), directors.get(), angles.get(),
                         fixed_drilling, global_force.get(), d);
  if (status != StatusCode::kOk) return status;
  return force->CopyFrom(global_force.get(), d);
}

}  // namespace fesa
