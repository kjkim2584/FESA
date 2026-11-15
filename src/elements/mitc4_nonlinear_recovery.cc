#include "elements/mitc4_nonlinear_recovery.h"

#include <cmath>
#include <memory>

#include "elements/mitc4_geometry.h"
#include "elements/mitc4_linear_ops.h"
#include "elements/mitc4_nonlinear_ops.h"
#include "fesa/math/dense_factorization.h"
#include "fesa/math/rotation.h"
#include "model/model_storage.h"
namespace fesa {
namespace {
StatusCode Make(int rows, int columns, std::unique_ptr<Matrix>* owner,
                Diagnostic* d) {
  Matrix* raw = nullptr;
  const auto status = Matrix::Create(rows, columns, &raw, d);
  owner->reset(raw);
  return status;
}
StatusCode MakeVector(int size, std::unique_ptr<Vector>* owner, Diagnostic* d) {
  Vector* raw = nullptr;
  const auto status = Vector::Create(size, &raw, d);
  owner->reset(raw);
  return status;
}
StatusCode Geometry(TrialEvaluationFailure* failure, Diagnostic* d) {
  *failure = TrialEvaluationFailure::kGeometry;
  return ReportFailure(StatusCode::kEquilibriumFailure,
                       DiagnosticStage::kComputation,
                       "MITC4 recovery current geometry/stretch", d);
}
StatusCode Located(StatusCode status, double r, double s, double zeta,
                   Diagnostic* d) {
  if (status != StatusCode::kOk && d && !d->geometry_point_available) {
    d->geometry_point_available = true;
    d->natural_xi = r;
    d->natural_eta = s;
    d->natural_zeta = zeta;
  }
  return status;
}
bool Shape(const Matrix* m, int rows, int columns) {
  return m && m->Rows() == rows && m->Columns() == columns;
}
StatusCode Preflight(const Node* const* nodes, int count,
                     const ShellSection* section, const Vector* u,
                     const Matrix* q, Matrix* stress, Matrix* strain,
                     Matrix* axes, Diagnostic* d) {
  if (!nodes || count != 4 || !section || !u || u->Size() != 12 ||
      !Shape(q, 4, 4) || !Shape(stress, 16, 6) || !Shape(strain, 16, 6) ||
      !Shape(axes, 3, 3) || stress == strain || stress == axes ||
      strain == axes || q == stress || q == strain || q == axes)
    return model_storage::Invalid("MITC4 finite recovery shape/alias", d);
  for (int i = 0; i < 12; ++i)
    if (!std::isfinite(u->Data()[i]))
      return model_storage::Invalid("MITC4 finite recovery translations", d);
  for (int i = 0; i < 16; ++i)
    if (!std::isfinite(q->Data()[i]))
      return model_storage::Invalid("MITC4 finite recovery orientations", d);
  if (!(section->thickness() / 2 > 0))
    return ReportFailure(StatusCode::kInputError, DiagnosticStage::kValidation,
                         "MITC4 initial half-thickness range", d);
  for (int i = 0; i < 4; ++i) {
    if (!nodes[i] || nodes[i]->dof_count() != 6)
      return model_storage::Invalid("MITC4 finite recovery node", d);
    for (int j = 0; j < i; ++j)
      if (nodes[i] == nodes[j] || nodes[i]->Id() == nodes[j]->Id())
        return model_storage::Invalid("MITC4 finite recovery duplicate node",
                                      d);
  }
  return StatusCode::kOk;
}
StatusCode RelativeRows(Matrix* values, Diagnostic* d) {
  std::unique_ptr<Matrix> anchor, ones, translated;
  auto status = Make(1, 3, &anchor, d);
  if (status == StatusCode::kOk) status = Make(4, 1, &ones, d);
  if (status == StatusCode::kOk) status = Make(4, 3, &translated, d);
  if (status != StatusCode::kOk) return status;
  for (int i = 0; i < 4; ++i) {
    status = ones->Set(i, 0, 1, d);
    if (status != StatusCode::kOk) return status;
  }
  status = anchor->CopyBlockFrom(values, 0, 0, 1, 3, 0, 0, d);
  if (status == StatusCode::kOk)
    status = Matrix::Multiply(ones.get(), false, anchor.get(), false,
                              translated.get(), d);
  return status == StatusCode::kOk ? values->AddScaled(-1, translated.get(), d)
                                   : status;
}
StatusCode BuildCoordinates(const Node* const* nodes, const Vector* u,
                            const Matrix* q, const Matrix* frame, Matrix* x0,
                            Matrix* x, Matrix* d0, Matrix* directors,
                            Diagnostic* d) {
  std::unique_ptr<Matrix> du, rotation;
  std::unique_ptr<Vector> normal, translation, quaternion, director;
  auto status = Make(4, 3, &du, d);
  if (status == StatusCode::kOk) status = Make(3, 3, &rotation, d);
  if (status == StatusCode::kOk) status = MakeVector(3, &normal, d);
  if (status == StatusCode::kOk) status = MakeVector(3, &translation, d);
  if (status == StatusCode::kOk) status = MakeVector(4, &quaternion, d);
  if (status == StatusCode::kOk) status = MakeVector(3, &director, d);
  if (status != StatusCode::kOk) return status;
  const int normal_indices[3] = {6, 7, 8};
  status = frame->GatherTo(normal_indices, 3, normal.get(), 0, d);
  if (status != StatusCode::kOk) return status;
  for (int node = 0; node < 4; ++node) {
    for (int axis = 0; axis < 3; ++axis) {
      double coordinate = 0;
      status = nodes[node]->Coordinate(axis, &coordinate, d);
      if (status == StatusCode::kOk)
        status = x0->Set(node, axis, coordinate, d);
      if (status != StatusCode::kOk) return status;
    }
    const int ui[3] = {3 * node, 3 * node + 1, 3 * node + 2};
    status = translation->GatherFrom(u, ui, 3, 0, d);
    if (status == StatusCode::kOk)
      status = du->SetRow(node, translation.get(), d);
    const int qi[4] = {4 * node, 4 * node + 1, 4 * node + 2, 4 * node + 3};
    if (status == StatusCode::kOk)
      status = q->GatherTo(qi, 4, quaternion.get(), 0, d);
    if (status == StatusCode::kOk)
      status = RotationMatrix(quaternion.get(), rotation.get(), d);
    if (status == StatusCode::kOk)
      status = Matrix::MultiplyVector(rotation.get(), false, normal.get(),
                                      director.get(), d);
    if (status == StatusCode::kOk) status = d0->SetRow(node, normal.get(), d);
    if (status == StatusCode::kOk)
      status = directors->SetRow(node, director.get(), d);
    if (status != StatusCode::kOk) return status;
  }
  status = RelativeRows(x0, d);
  if (status == StatusCode::kOk) status = RelativeRows(du.get(), d);
  if (status == StatusCode::kOk) status = x->CopyFrom(x0, d);
  return status == StatusCode::kOk ? x->AddScaled(1, du.get(), d) : status;
}
StatusCode ApplyBasis(const Matrix* wp, const Matrix* wd, const Matrix* x,
                      const Matrix* directors, Matrix* tmp, Matrix* basis,
                      Diagnostic* d) {
  auto status = Matrix::Multiply(wp, false, x, false, basis, d);
  if (status == StatusCode::kOk)
    status = Matrix::Multiply(wd, false, directors, false, tmp, d);
  return status == StatusCode::kOk ? basis->AddScaled(1, tmp, d) : status;
}
StatusCode Bases(const Matrix* x0, const Matrix* x, const Matrix* d0,
                 const Matrix* directors, double r, double s, double zeta,
                 double thickness, Matrix* initial, Matrix* current,
                 Diagnostic* d) {
  std::unique_ptr<Matrix> wp, wd, tmp;
  auto status = Make(3, 4, &wp, d);
  if (status == StatusCode::kOk) status = Make(3, 4, &wd, d);
  if (status == StatusCode::kOk) status = Make(3, 3, &tmp, d);
  if (status != StatusCode::kOk) return status;
  constexpr double kr[4] = {-1, 1, 1, -1}, ks[4] = {-1, -1, 1, 1};
  for (int node = 0; node < 4; ++node) {
    const double n = (1 + kr[node] * r) * (1 + ks[node] * s) / 4;
    const double nr = kr[node] * (1 + ks[node] * s) / 4,
                 ns = ks[node] * (1 + kr[node] * r) / 4;
    const double position[3] = {nr, ns, 0};
    const double director[3] = {thickness * zeta / 2 * nr,
                                thickness * zeta / 2 * ns, thickness / 2 * n};
    for (int axis = 0; axis < 3; ++axis) {
      status = wp->Set(axis, node, position[axis], d);
      if (status == StatusCode::kOk)
        status = wd->Set(axis, node, director[axis], d);
      if (status != StatusCode::kOk) return status;
    }
  }
  status = ApplyBasis(wp.get(), wd.get(), x0, d0, tmp.get(), initial, d);
  return status == StatusCode::kOk
             ? ApplyBasis(wp.get(), wd.get(), x, directors, tmp.get(), current,
                          d)
             : status;
}
StatusCode ReadBasis(const Matrix* basis, Vector* const* rows, Diagnostic* d) {
  for (int i = 0; i < 3; ++i) {
    const int indices[3] = {3 * i, 3 * i + 1, 3 * i + 2};
    const auto status = basis->GatherTo(indices, 3, rows[i], 0, d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}
StatusCode OutputFrame(const Matrix* x0, const Matrix* x, const Matrix* d0,
                       const Matrix* dirs, double thickness, Matrix* axes,
                       TrialEvaluationFailure* failure, Diagnostic* d) {
  std::unique_ptr<Matrix> initial, current;
  auto status = Make(3, 3, &initial, d);
  if (status == StatusCode::kOk) status = Make(3, 3, &current, d);
  if (status != StatusCode::kOk) return status;
  status = Bases(x0, x, d0, dirs, 0, 0, 0, thickness, initial.get(),
                 current.get(), d);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<Vector> values[9];
  Vector* g0[3] = {};
  Vector* g[3] = {};
  for (int i = 0; i < 9; ++i) {
    status = MakeVector(3, &values[i], d);
    if (status != StatusCode::kOk) return status;
  }
  for (int i = 0; i < 3; ++i) {
    g0[i] = values[i].get();
    g[i] = values[i + 3].get();
  }
  status = ReadBasis(initial.get(), g0, d);
  if (status == StatusCode::kOk) status = ReadBasis(current.get(), g, d);
  double determinant = 0, jacobian = 0;
  if (status == StatusCode::kOk)
    status =
        ValidateMitc4VolumeBasis(g0, g, &determinant, &jacobian, failure, d);
  if (status != StatusCode::kOk) return status;
  Vector* a = values[6].get();
  Vector* normal = values[7].get();
  Vector* b = values[8].get();
  status = a->CopyFrom(g[0], d);
  double norm = 0;
  if (status == StatusCode::kOk) status = a->Norm2(&norm, d);
  if (status == StatusCode::kOk) status = a->Scale(1 / norm, d);
  if (status == StatusCode::kOk)
    status = CrossMitc4Vectors(g[0], g[1], normal, d);
  if (status == StatusCode::kOk) status = normal->Norm2(&norm, d);
  if (status == StatusCode::kOk && (!(norm > 0) || !std::isfinite(1 / norm)))
    return Geometry(failure, d);
  if (status == StatusCode::kOk) status = normal->Scale(1 / norm, d);
  if (status == StatusCode::kOk) status = CrossMitc4Vectors(normal, a, b, d);
  double direction = 0;
  if (status == StatusCode::kOk) status = normal->Dot(g[2], &direction, d);
  if (status == StatusCode::kOk && !(direction > 0))
    return Geometry(failure, d);
  if (status == StatusCode::kOk) status = axes->SetRow(0, a, d);
  if (status == StatusCode::kOk) status = axes->SetRow(1, b, d);
  if (status == StatusCode::kOk) status = axes->SetRow(2, normal, d);
  return status;
}
StatusCode RotateTensor(const Matrix* axes, const Matrix* tensor,
                        Matrix* temporary, Matrix* out, Diagnostic* d) {
  auto status = Matrix::Multiply(axes, false, tensor, false, temporary, d);
  return status == StatusCode::kOk
             ? Matrix::Multiply(temporary, false, axes, true, out, d)
             : status;
}
StatusCode Point(const Matrix* x0, const Matrix* x, const Matrix* d0,
                 const Matrix* dirs, const Matrix* frame_transpose,
                 const Matrix* axes, const Matrix* constitutive,
                 const Matrix* engineering, const Matrix* pk2_selector,
                 double r, double s, double zeta, double thickness,
                 Vector* stress, Vector* strain,
                 TrialEvaluationFailure* failure, Diagnostic* d) {
  std::unique_ptr<Matrix> m[13];
  for (int i = 0; i < 13; ++i) {
    const int rows = i == 3 || i == 4 ? 5 : 3, columns = i == 3   ? 5
                                                         : i == 4 ? 7
                                                                  : 3;
    const auto status = Make(rows, columns, &m[i], d);
    if (status != StatusCode::kOk) return status;
  }
  Matrix *g0 = m[0].get(), *g = m[1].get(), *h = m[2].get(), *p = m[3].get(),
         *selector = m[4].get();
  Matrix *ft = m[5].get(), *pk = m[6].get(), *f_local = m[7].get(),
         *tmp = m[8].get();
  Matrix *sigma = m[9].get(), *left_metric = m[10].get(),
         *log_metric = m[11].get(), *rotated = m[12].get();
  std::unique_ptr<Vector> vectors[11];
  const int sizes[11] = {7, 5, 5, 5, 9, 3, 3, 3, 3, 3, 3};
  for (int i = 0; i < 11; ++i) {
    const auto status = MakeVector(sizes[i], &vectors[i], d);
    if (status != StatusCode::kOk) return status;
  }
  Vector *raw = vectors[0].get(), *covariant = vectors[1].get(),
         *e = vectors[2].get(), *pk5 = vectors[3].get(),
         *flat = vectors[4].get();
  Vector* initial[3] = {vectors[5].get(), vectors[6].get(), vectors[7].get()};
  Vector* current[3] = {vectors[8].get(), vectors[9].get(), vectors[10].get()};
  const double rr[5] = {r, 0, 0, -1, 1}, ss[5] = {s, -1, 1, 0, 0};
  constexpr int point_index[7] = {0, 0, 0, 1, 2, 3, 4};
  constexpr int aa[7] = {0, 1, 0, 0, 0, 1, 1}, bb[7] = {0, 1, 1, 2, 2, 2, 2};
  double jacobian = 0;
  for (int point = 0; point < 5; ++point) {
    auto status =
        Bases(x0, x, d0, dirs, rr[point], ss[point], zeta, thickness, g0, g, d);
    if (status == StatusCode::kOk) status = ReadBasis(g0, initial, d);
    if (status == StatusCode::kOk) status = ReadBasis(g, current, d);
    double determinant = 0, j = 0;
    if (status == StatusCode::kOk)
      status = ValidateMitc4VolumeBasis(initial, current, &determinant, &j,
                                        failure, d);
    if (status != StatusCode::kOk)
      return Located(status, rr[point], ss[point], zeta, d);
    if (point == 0) {
      // One initial G_rows factor serves H and physical F^T; the assumed
      // MITC strain is never used to reconstruct the deformation gradient.
      jacobian = j;
      DenseLuFactorization* raw_factor = nullptr;
      status = DenseLuFactorization::Create(g0, &raw_factor, d);
      std::unique_ptr<DenseLuFactorization> factor(raw_factor);
      if (status == StatusCode::kOk)
        status = factor->Solve(frame_transpose, true, h, d);
      if (status == StatusCode::kOk) status = factor->Solve(g, false, ft, d);
      if (status == StatusCode::kOk)
        status = BuildMitc4NonlinearTransformFromH(h, p, d);
      if (status != StatusCode::kOk) return status;
    }
    for (int row = 0; row < 7; ++row)
      if (point_index[row] == point) {
        double actual = 0, reference = 0;
        status = current[aa[row]]->Dot(current[bb[row]], &actual, d);
        if (status == StatusCode::kOk)
          status = initial[aa[row]]->Dot(initial[bb[row]], &reference, d);
        if (status != StatusCode::kOk) return status;
        const double value = (actual - reference) / 2;
        if (!std::isfinite(value))
          return ReportFailure(StatusCode::kLibraryFailure,
                               DiagnosticStage::kComputation,
                               "MITC4 recovery covariant strain range", d);
        status = raw->Set(row, value, d);
        if (status != StatusCode::kOk) return status;
      }
  }
  auto status = BuildMitc4NonlinearTyingSelector(r, s, selector, d);
  if (status == StatusCode::kOk)
    status = Matrix::MultiplyVector(selector, false, raw, covariant, d);
  if (status == StatusCode::kOk)
    status = Matrix::MultiplyVector(p, false, covariant, e, d);
  if (status == StatusCode::kOk)
    status = Matrix::MultiplyVector(constitutive, false, e, pk5, d);
  if (status != StatusCode::kOk) return status;
  // Expand constitutive5 into a symmetric tensor using the constant linear
  // selector, then explicit library gather/copies (no numerical copy loop).
  status = Matrix::MultiplyVector(pk2_selector, false, pk5, flat, d);
  if (status != StatusCode::kOk) return status;
  Vector* tensor_row = current[0];  // Volume rows are no longer needed.
  for (int row = 0; row < 3; ++row) {
    const int indices[3] = {3 * row, 3 * row + 1, 3 * row + 2};
    status = tensor_row->GatherFrom(flat, indices, 3, 0, d);
    if (status == StatusCode::kOk) status = pk->SetRow(row, tensor_row, d);
    if (status != StatusCode::kOk) return status;
  }
  // M=F Q0^T; sigma=M S_local M^T/J, then rotate all components to Q_E.
  status = Matrix::Multiply(ft, true, frame_transpose, false, f_local, d);
  if (status == StatusCode::kOk)
    status = Matrix::Multiply(f_local, false, pk, false, tmp, d);
  if (status == StatusCode::kOk)
    status = Matrix::Multiply(tmp, false, f_local, true, sigma, d);
  if (status == StatusCode::kOk && !std::isfinite(1 / jacobian))
    return Geometry(failure, d);
  if (status == StatusCode::kOk) status = sigma->Scale(1 / jacobian, d);
  if (status == StatusCode::kOk)
    status = RotateTensor(axes, sigma, tmp, rotated, d);
  constexpr int six[6] = {0, 4, 8, 1, 2, 5};
  if (status == StatusCode::kOk)
    status = rotated->GatherTo(six, 6, stress, 0, d);
  // LE=.5 log(F F^T), including physical thickness/shear components.
  if (status == StatusCode::kOk)
    status = Matrix::Multiply(ft, true, ft, false, left_metric, d);
  bool not_positive = false;
  if (status == StatusCode::kOk)
    status = Matrix::SymmetricPositiveDefiniteLog(left_metric, log_metric,
                                                  &not_positive, d);
  if (not_positive) return Geometry(failure, d);
  if (status == StatusCode::kOk) status = log_metric->Scale(.5, d);
  if (status == StatusCode::kOk)
    status = RotateTensor(axes, log_metric, tmp, rotated, d);
  constexpr int all[9] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  if (status == StatusCode::kOk) status = rotated->GatherTo(all, 9, flat, 0, d);
  return status == StatusCode::kOk
             ? Matrix::MultiplyVector(engineering, false, flat, strain, d)
             : status;
}
}  // namespace
StatusCode RecoverMitc4FiniteSamples(
    const Node* const* nodes, int count, const ShellSection* section,
    const Vector* translations, const Matrix* orientations,
    Matrix* stress_samples, Matrix* strain_samples, Matrix* output_axes,
    TrialEvaluationFailure* failure, Diagnostic* d) {
  if (!failure) return model_storage::Invalid("MITC4 recovery failure flag", d);
  *failure = TrialEvaluationFailure::kNone;
  auto status = Preflight(nodes, count, section, translations, orientations,
                          stress_samples, strain_samples, output_axes, d);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<Matrix> local, frame, frame_transpose, identity, x0, x, d0,
      dirs, constitutive, engineering, pk2_selector;
  status = Make(4, 3, &local, d);
  if (status == StatusCode::kOk) status = Make(3, 3, &frame, d);
  if (status == StatusCode::kOk) status = Make(3, 3, &frame_transpose, d);
  if (status == StatusCode::kOk) status = Make(3, 3, &identity, d);
  if (status == StatusCode::kOk) status = Make(4, 3, &x0, d);
  if (status == StatusCode::kOk) status = Make(4, 3, &x, d);
  if (status == StatusCode::kOk) status = Make(4, 3, &d0, d);
  if (status == StatusCode::kOk) status = Make(4, 3, &dirs, d);
  if (status == StatusCode::kOk) status = Make(5, 5, &constitutive, d);
  if (status == StatusCode::kOk) status = Make(6, 9, &engineering, d);
  if (status == StatusCode::kOk) status = Make(9, 5, &pk2_selector, d);
  if (status != StatusCode::kOk) return status;
  double diameter = 0, length = 0;
  status = ValidateMitc4Geometry(nodes, count, local.get(), frame.get(),
                                 &diameter, &length, d);
  if (status == StatusCode::kOk)
    status = BuildCoordinates(nodes, translations, orientations, frame.get(),
                              x0.get(), x.get(), d0.get(), dirs.get(), d);
  for (int i = 0; status == StatusCode::kOk && i < 3; ++i)
    status = identity->Set(i, i, 1, d);
  if (status == StatusCode::kOk)
    status = Matrix::Multiply(frame.get(), true, identity.get(), false,
                              frame_transpose.get(), d);
  if (status == StatusCode::kOk)
    status = BuildMitc4LinearConstitutive(section->material(),
                                          constitutive.get(), d);
  constexpr int indices[6] = {0, 4, 8, 1, 2, 5};
  for (int i = 0; status == StatusCode::kOk && i < 6; ++i)
    status = engineering->Set(i, indices[i], i < 3 ? 1 : 2, d);
  constexpr int pk2_components[8] = {0, 2, 4, 2, 1, 3, 4, 3};
  for (int row = 0; status == StatusCode::kOk && row < 8; ++row)
    status = pk2_selector->Set(row, pk2_components[row], 1, d);
  if (status != StatusCode::kOk) return status;
  status = OutputFrame(x0.get(), x.get(), d0.get(), dirs.get(),
                       section->thickness(), output_axes, failure, d);
  if (status != StatusCode::kOk) return Located(status, 0, 0, 0, d);
  std::unique_ptr<Vector> stress, strain;
  status = MakeVector(6, &stress, d);
  if (status == StatusCode::kOk) status = MakeVector(6, &strain, d);
  if (status != StatusCode::kOk) return status;
  const double a = 1 / std::sqrt(3.), rr[4] = {-a, a, a, -a},
               ss[4] = {-a, -a, a, a}, zeta[4] = {-a, a, -1, 1};
  for (int layer = 0; layer < 4; ++layer)
    for (int plane = 0; plane < 4; ++plane) {
      status =
          Point(x0.get(), x.get(), d0.get(), dirs.get(), frame_transpose.get(),
                output_axes, constitutive.get(), engineering.get(),
                pk2_selector.get(), rr[plane], ss[plane], zeta[layer],
                section->thickness(), stress.get(), strain.get(), failure, d);
      if (status != StatusCode::kOk)
        return Located(status, rr[plane], ss[plane], zeta[layer], d);
      status = stress_samples->SetRow(4 * layer + plane, stress.get(), d);
      if (status == StatusCode::kOk)
        status = strain_samples->SetRow(4 * layer + plane, strain.get(), d);
      if (status != StatusCode::kOk)
        return Located(status, rr[plane], ss[plane], zeta[layer], d);
    }
  return StatusCode::kOk;
}
}  // namespace fesa
