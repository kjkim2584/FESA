#include "elements/mitc4_rayleigh.h"

#include <cmath>
#include <memory>

#include "elements/mitc4_geometry.h"
#include "elements/mitc4_nonlinear_ops.h"
#include "fesa/math/dense_factorization.h"
#include "fesa/math/rotation.h"
#include "model/model_storage.h"

namespace fesa {
StatusCode CrossMitc4RayleighVectors(Matrix* scratch, const Vector* a,
                                    const Vector* b, Vector* out,
                                    Diagnostic* d) {
  if (!a || !b || !out || !scratch || a->Size() != 3 || b->Size() != 3 ||
      out->Size() != 3 || out == a || out == b || scratch->Rows() != 3 ||
      scratch->Columns() != 3)
    return model_storage::Invalid("MITC4 borrowed cross shape/alias", d);
  for (int i = 0; i < 3; ++i)
    if (!std::isfinite(a->Data()[i]) || !std::isfinite(b->Data()[i]))
      return model_storage::Invalid("MITC4 borrowed cross finite input", d);
  const double* v = a->Data();
  const double coefficients[9] = {0.0, -v[2], v[1], v[2], 0.0,
                                  -v[0], -v[1], v[0], 0.0};
  for (int i = 0; i < 9; ++i) {
    const auto status = scratch->Set(i / 3, i % 3, coefficients[i], d);
    if (status != StatusCode::kOk) return status;
  }
  return Matrix::MultiplyVector(scratch, false, b, out, d);
}

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
StatusCode Row(const Matrix* source, int row, Vector* out, Diagnostic* d) {
  int indices[4] = {};
  for (int i = 0; i < source->Columns(); ++i)
    indices[i] = row * source->Columns() + i;
  return source->GatherTo(indices, source->Columns(), out, 0, d);
}
StatusCode SymmetricProduct(const Matrix* a, const Matrix* b, Matrix* out,
                            Diagnostic* d) {
  auto status = Matrix::MultiplyAdd(a, true, b, false, .5, 0, out, d);
  return status == StatusCode::kOk
             ? Matrix::MultiplyAdd(b, true, a, false, .5, 1, out, d)
             : status;
}

/// Owns only call-local 3x3 tensors and scratch, reused across this call's
/// quadrature points/directions. RD4/RD7/RD8 are independent of HHT.
class RayleighPoint {
 public:
  RayleighPoint() = default;
  ~RayleighPoint() {
    for (auto* value : {f_, rate_, a_, e_, t_, s_, p_, inverse_, identity_})
      delete value;
    for (auto* value : work_) delete value;
    delete diagonal_;
    delete ones_;
  }
  RayleighPoint(const RayleighPoint&) = delete;
  RayleighPoint& operator=(const RayleighPoint&) = delete;
  RayleighPoint(RayleighPoint&&) = delete;
  RayleighPoint& operator=(RayleighPoint&&) = delete;
  StatusCode Initialize(Diagnostic* d) {
    for (auto** value : {&f_, &rate_, &a_, &e_, &t_, &s_, &p_, &inverse_,
                         &identity_}) {
      const auto status = Matrix::Create(3, 3, value, d);
      if (status != StatusCode::kOk) return status;
    }
    for (auto& value : work_) {
      const auto status = Matrix::Create(3, 3, &value, d);
      if (status != StatusCode::kOk) return status;
    }
    auto status = Vector::Create(3, &diagonal_, d);
    if (status == StatusCode::kOk) status = Vector::Create(3, &ones_, d);
    for (int i = 0; status == StatusCode::kOk && i < 3; ++i)
      status = ones_->Set(i, 1, d);
    return status == StatusCode::kOk
               ? Matrix::Diagonal(ones_, identity_, d) : status;
  }
  StatusCode SetPoint(const Matrix* f, const Matrix* rate, const Matrix* h,
                       double jacobian, double young, double nu, double beta,
                       Diagnostic* d) {
    j_ = jacobian;
    beta_ = beta;
    mu_ = young / (2 * (1 + nu));
    lambda_ = young * nu / (1 - nu * nu);
    auto status = f_->CopyFrom(f, d);
    if (status == StatusCode::kOk) status = rate_->CopyFrom(rate, d);
    DenseLuFactorization* raw = nullptr;
    if (status == StatusCode::kOk)
      status = DenseLuFactorization::Create(f_, &raw, d);
    std::unique_ptr<DenseLuFactorization> factor(raw);
    if (status == StatusCode::kOk)
      status = factor->Solve(identity_, false, inverse_, d);
    std::unique_ptr<Matrix> fh, metric, id2, metric_inverse, product;
    if (status == StatusCode::kOk) status = Make(3, 2, &fh, d);
    if (status == StatusCode::kOk) status = Make(2, 2, &metric, d);
    if (status == StatusCode::kOk) status = Make(2, 2, &id2, d);
    if (status == StatusCode::kOk) status = Make(2, 2, &metric_inverse, d);
    if (status == StatusCode::kOk) status = Make(3, 2, &product, d);
    for (int i = 0; status == StatusCode::kOk && i < 2; ++i)
      status = id2->Set(i, i, 1, d);
    if (status == StatusCode::kOk)
      status = Matrix::Multiply(f_, false, h, false, fh.get(), d);
    if (status == StatusCode::kOk)
      status = Matrix::Multiply(fh.get(), true, fh.get(), false,
                                  metric.get(), d);
    raw = nullptr;
    if (status == StatusCode::kOk)
      status = DenseLuFactorization::Create(metric.get(), &raw, d);
    std::unique_ptr<DenseLuFactorization> metric_factor(raw);
    if (status == StatusCode::kOk)
      status = metric_factor->Solve(id2.get(), false, metric_inverse.get(), d);
    if (status == StatusCode::kOk)
      status = Matrix::Multiply(h, false, metric_inverse.get(), false,
                                product.get(), d);
    if (status == StatusCode::kOk)
      status = Matrix::Multiply(product.get(), false, h, true, a_, d);
    if (status == StatusCode::kOk)
      status = SymmetricProduct(f_, rate_, e_, d);
    if (status == StatusCode::kOk) status = ElasticRate(e_, t_, &q_, d);
    if (status == StatusCode::kOk) status = s_->CopyFrom(t_, d);
    if (status == StatusCode::kOk) status = s_->Scale(beta_ * j_, d);
    return status == StatusCode::kOk
               ? Matrix::Multiply(f_, false, s_, false, p_, d) : status;
  }
  const Matrix* piola() const { return p_; }
  StatusCode Direction(const Matrix* df, const Matrix* drate, bool fixed_pose,
                        Matrix* out, Diagnostic* d) {
    Matrix *da = work_[0], *de = work_[1], *dt = work_[2], *tmp = work_[3],
           *product = work_[4], *ds = work_[5], *metric = work_[6];
    auto status = SymmetricProduct(f_, drate, de, d);
    if (!fixed_pose && status == StatusCode::kOk) {
      status = SymmetricProduct(df, rate_, tmp, d);
      if (status == StatusCode::kOk) status = de->AddScaled(1, tmp, d);
    }
    double dq = 0;
    if (status == StatusCode::kOk) status = ElasticRate(de, dt, &dq, d);
    double dj = 0;
    if (!fixed_pose && status == StatusCode::kOk) {
      // RD7: dA=-A*(dF^T F+F^T dF)*A; dJ=J tr(F^-1 dF).
      status = SymmetricProduct(df, f_, metric, d);
      if (status == StatusCode::kOk) status = metric->Scale(2, d);
      if (status == StatusCode::kOk)
        status = Matrix::Multiply(a_, false, metric, false, tmp, d);
      if (status == StatusCode::kOk)
        status = Matrix::Multiply(tmp, false, a_, false, da, d);
      if (status == StatusCode::kOk) status = da->Scale(-1, d);
      if (status == StatusCode::kOk)
        status = Matrix::Multiply(inverse_, false, df, false, tmp, d);
      if (status == StatusCode::kOk) status = Trace(tmp, &dj, d);
      dj *= j_;
      if (status == StatusCode::kOk)
        status = Matrix::Multiply(da, false, e_, false, tmp, d);
      if (status == StatusCode::kOk) status = Trace(tmp, &dq, d);
      if (status == StatusCode::kOk)
        status = Matrix::Multiply(tmp, false, a_, false, product, d);
      if (status == StatusCode::kOk)
        status = dt->AddScaled(2 * mu_, product, d);
      if (status == StatusCode::kOk)
        status = Matrix::Multiply(a_, false, e_, false, tmp, d);
      if (status == StatusCode::kOk)
        status = Matrix::Multiply(tmp, false, da, false, product, d);
      if (status == StatusCode::kOk)
        status = dt->AddScaled(2 * mu_, product, d);
      if (status == StatusCode::kOk)
        status = dt->AddScaled(lambda_ * dq, a_, d);
      if (status == StatusCode::kOk)
        status = dt->AddScaled(lambda_ * q_, da, d);
    }
    if (status == StatusCode::kOk) status = ds->CopyFrom(dt, d);
    if (status == StatusCode::kOk) status = ds->Scale(beta_ * j_, d);
    if (status == StatusCode::kOk) status = ds->AddScaled(beta_ * dj, t_, d);
    if (status == StatusCode::kOk)
      status = Matrix::Multiply(f_, false, ds, false, out, d);
    return !fixed_pose && status == StatusCode::kOk
               ? Matrix::MultiplyAdd(df, false, s_, false, 1, 1, out, d)
               : status;
  }

 private:
  StatusCode Trace(const Matrix* value, double* out, Diagnostic* d) {
    const int indices[] = {0, 4, 8};
    const auto status = value->GatherTo(indices, 3, diagonal_, 0, d);
    return status == StatusCode::kOk ? diagonal_->Dot(ones_, out, d) : status;
  }
  StatusCode ElasticRate(const Matrix* rate, Matrix* out, double* trace,
                          Diagnostic* d) {
    auto status = Matrix::Multiply(a_, false, rate, false, work_[3], d);
    if (status == StatusCode::kOk) status = Trace(work_[3], trace, d);
    if (status == StatusCode::kOk)
      status = Matrix::Multiply(work_[3], false, a_, false, out, d);
    if (status == StatusCode::kOk) status = out->Scale(2 * mu_, d);
    return status == StatusCode::kOk
               ? out->AddScaled(lambda_ * *trace, a_, d) : status;
  }
  Matrix *f_ = nullptr, *rate_ = nullptr, *a_ = nullptr, *e_ = nullptr,
         *t_ = nullptr, *s_ = nullptr, *p_ = nullptr, *inverse_ = nullptr,
         *identity_ = nullptr;  ///< Owned tensors, one element-call lifetime.
  Matrix* work_[7] = {};        ///< Owned disjoint directional scratch.
  Vector *diagonal_ = nullptr, *ones_ = nullptr;  ///< Owned trace scratch.
  double j_ = 0, beta_ = 0, mu_ = 0, lambda_ = 0, q_ = 0;
};

StatusCode PointGradients(const Matrix* x0, const Matrix* x,
                           const Matrix* initial_directors,
                           const Matrix* directors, double r, double s,
                           double zeta, double thickness, Matrix* a, Matrix* b,
                           double* weight, double* jacobian,
                           Mitc4VolumeWorkspace* workspace,
                           TrialEvaluationFailure* failure, Diagnostic* d) {
  std::unique_ptr<Matrix> wp, wd, initial, current, solved, identity;
  auto status = Make(3, 4, &wp, d);
  if (status == StatusCode::kOk) status = Make(3, 4, &wd, d);
  for (auto* owner : {&initial, &current, &identity}) {
    if (status == StatusCode::kOk) status = Make(3, 3, owner, d);
  }
  if (status == StatusCode::kOk) status = Make(3, 4, &solved, d);
  if (status != StatusCode::kOk) return status;
  constexpr double kr[] = {-1, 1, 1, -1}, ks[] = {-1, -1, 1, 1};
  for (int node = 0; node < 4; ++node) {
    const double n = (1 + kr[node] * r) * (1 + ks[node] * s) / 4;
    const double nr = kr[node] * (1 + ks[node] * s) / 4;
    const double ns = ks[node] * (1 + kr[node] * r) / 4;
    const double translation[] = {nr, ns, 0};
    const double rotation[] = {thickness * zeta / 2 * nr,
                               thickness * zeta / 2 * ns, thickness / 2 * n};
    for (int c = 0; c < 3; ++c) {
      status = wp->Set(c, node, translation[c], d);
      if (status == StatusCode::kOk) status = wd->Set(c, node, rotation[c], d);
      if (status != StatusCode::kOk) return status;
    }
  }
  status = Matrix::Multiply(wp.get(), false, x0, false, initial.get(), d);
  if (status == StatusCode::kOk)
    status = Matrix::MultiplyAdd(wd.get(), false, initial_directors, false,
                                 1, 1, initial.get(), d);
  if (status == StatusCode::kOk)
    status = Matrix::Multiply(wp.get(), false, x, false, current.get(), d);
  if (status == StatusCode::kOk)
    status = Matrix::MultiplyAdd(wd.get(), false, directors, false,
                                 1, 1, current.get(), d);
  std::unique_ptr<Vector> rows[6];
  const Vector* g0[3] = {};
  const Vector* g[3] = {};
  for (int i = 0; status == StatusCode::kOk && i < 6; ++i) {
    status = MakeVector(3, &rows[i], d);
    if (status == StatusCode::kOk)
      status = Row(i < 3 ? initial.get() : current.get(), i % 3,
                    rows[i].get(), d);
    if (i < 3) g0[i] = rows[i].get();
    else g[i - 3] = rows[i].get();
  }
  if (status == StatusCode::kOk)
    status = workspace->Validate(g0, g, weight, jacobian, failure, d);
  DenseLuFactorization* raw = nullptr;
  if (status == StatusCode::kOk)
    status = DenseLuFactorization::Create(initial.get(), &raw, d);
  std::unique_ptr<DenseLuFactorization> factor(raw);
  for (int i = 0; status == StatusCode::kOk && i < 3; ++i)
    status = identity->Set(i, i, 1, d);
  for (int i = 0; status == StatusCode::kOk && i < 2; ++i) {
    status = factor->Solve(i == 0 ? wp.get() : wd.get(), false,
                            solved.get(), d);
    if (status == StatusCode::kOk)
      status = Matrix::Multiply(solved.get(), true, identity.get(), false,
                                i == 0 ? a : b, d);
  }
  return status;
}

StatusCode Coordinates(const Node* const* nodes, const ElementMotionInput* in,
                         Matrix* frame, Matrix* x0, Matrix* x, Matrix* d0,
                         Matrix* directors, Matrix* velocity, Matrix* omega,
                         Matrix* director_rate, Matrix* cross_scratch,
                         Diagnostic* d) {
  std::unique_ptr<Matrix> local, du, rotation, anchor, ones, shift;
  auto status = Make(4, 3, &local, d);
  if (status == StatusCode::kOk) status = Make(4, 3, &du, d);
  if (status == StatusCode::kOk) status = Make(3, 3, &rotation, d);
  if (status == StatusCode::kOk) status = Make(1, 3, &anchor, d);
  if (status == StatusCode::kOk) status = Make(4, 1, &ones, d);
  if (status == StatusCode::kOk) status = Make(4, 3, &shift, d);
  double diameter = 0, length = 0;
  if (status == StatusCode::kOk)
    status = ValidateMitc4Geometry(nodes, 4, local.get(), frame,
                                    &diameter, &length, d);
  if (status == StatusCode::kOk)
    status = Matrix::Multiply(local.get(), false, frame, false, x0, d);
  std::unique_ptr<Vector> normal, q, row, director, w, rate;
  if (status == StatusCode::kOk) status = MakeVector(3, &normal, d);
  if (status == StatusCode::kOk) status = MakeVector(4, &q, d);
  for (auto* owner : {&row, &director, &w, &rate})
    if (status == StatusCode::kOk) status = MakeVector(3, owner, d);
  if (status == StatusCode::kOk) status = Row(frame, 2, normal.get(), d);
  for (int node = 0; status == StatusCode::kOk && node < 4; ++node) {
    int ti[] = {3 * node, 3 * node + 1, 3 * node + 2};
    int vi[] = {6 * node, 6 * node + 1, 6 * node + 2};
    int wi[] = {6 * node + 3, 6 * node + 4, 6 * node + 5};
    status = row->GatherFrom(in->finite_state->translations, ti, 3, 0, d);
    if (status == StatusCode::kOk) status = du->SetRow(node, row.get(), d);
    if (status == StatusCode::kOk)
      status = Row(in->finite_state->orientations, node, q.get(), d);
    if (status == StatusCode::kOk)
      status = RotationMatrix(q.get(), rotation.get(), d);
    if (status == StatusCode::kOk)
      status = Matrix::MultiplyVector(rotation.get(), false, normal.get(),
                                      director.get(), d);
    if (status == StatusCode::kOk)
      status = directors->SetRow(node, director.get(), d);
    if (status == StatusCode::kOk) status = d0->SetRow(node, normal.get(), d);
    if (status == StatusCode::kOk)
      status = row->GatherFrom(in->velocity, vi, 3, 0, d);
    if (status == StatusCode::kOk)
      status = velocity->SetRow(node, row.get(), d);
    if (status == StatusCode::kOk)
      status = w->GatherFrom(in->velocity, wi, 3, 0, d);
    if (status == StatusCode::kOk) status = omega->SetRow(node, w.get(), d);
    if (status == StatusCode::kOk)
      status = CrossMitc4RayleighVectors(cross_scratch, w.get(),
                                        director.get(), rate.get(), d);
    if (status == StatusCode::kOk)
      status = director_rate->SetRow(node, rate.get(), d);
    if (status == StatusCode::kOk) status = ones->Set(node, 0, 1, d);
  }
  if (status == StatusCode::kOk) status = x->CopyFrom(x0, d);
  if (status == StatusCode::kOk) status = x->AddScaled(1, du.get(), d);
  // Remove the common current translation before differentiating shape data.
  if (status == StatusCode::kOk)
    status = anchor->CopyBlockFrom(x, 0, 0, 1, 3, 0, 0, d);
  if (status == StatusCode::kOk)
    status = Matrix::Multiply(ones.get(), false, anchor.get(), false,
                              shift.get(), d);
  return status == StatusCode::kOk ? x->AddScaled(-1, shift.get(), d) : status;
}

StatusCode Gradient(const Matrix* translation, const Matrix* rotation,
                      const Matrix* a, const Matrix* b, Matrix* out,
                      Diagnostic* d) {
  const auto status = Matrix::Multiply(translation, true, a, false, out, d);
  return status == StatusCode::kOk
             ? Matrix::MultiplyAdd(rotation, true, b, false, 1, 1, out, d)
             : status;
}
StatusCode Outer(const Vector* left, const Matrix* right, int row, Matrix* out,
                   Matrix* left_row, Matrix* right_row, Vector* scratch,
                   Diagnostic* d) {
  auto status = left_row->SetRow(0, left, d);
  if (status == StatusCode::kOk) status = Row(right, row, scratch, d);
  if (status == StatusCode::kOk) status = right_row->SetRow(0, scratch, d);
  return status == StatusCode::kOk
             ? Matrix::Multiply(left_row, true, right_row, false, out, d)
             : status;
}
StatusCode NodalForce(const Matrix* piola, const Matrix* directors,
                       const Matrix* a, const Matrix* b, Vector* out,
                       Matrix* const* matrices, Vector* const* vectors,
                       int changed_node, const Vector* dd,
                       const Matrix* base_pb, Matrix* cross_scratch,
                       Diagnostic* d) {
  // Reuses exclusively borrowed call-local scratch; matrices[1] holds P*b_i.
  auto status = Matrix::Multiply(a, false, piola, true, matrices[0], d);
  if (status == StatusCode::kOk)
    status = Matrix::Multiply(b, false, piola, true, matrices[1], d);
  const int identity[] = {0, 1, 2};
  for (int node = 0; status == StatusCode::kOk && node < 4; ++node) {
    status = Row(directors, node, vectors[0], d);
    if (status == StatusCode::kOk)
      status = Row(matrices[1], node, vectors[1], d);
    if (status == StatusCode::kOk)
      status = CrossMitc4RayleighVectors(cross_scratch, vectors[0],
                                        vectors[1], vectors[2], d);
    if (status == StatusCode::kOk && node == changed_node) {
      status = Row(base_pb, node, vectors[1], d);
      if (status == StatusCode::kOk)
        status = CrossMitc4RayleighVectors(cross_scratch, dd, vectors[1],
                                          vectors[4], d);
      if (status == StatusCode::kOk)
        status = vectors[2]->AddScaled(1, vectors[4], d);
    }
    if (status == StatusCode::kOk)
      status = Row(matrices[0], node, vectors[3], d);
    if (status == StatusCode::kOk)
      status = out->GatherFrom(vectors[3], identity, 3, 6 * node, d);
    if (status == StatusCode::kOk)
      status = out->GatherFrom(vectors[2], identity, 3, 6 * node + 3, d);
  }
  return status;
}
}  // namespace

StatusCode EvaluateMitc4RayleighBeta(const Node* const* nodes,
                                     const ShellSection* section,
                                     const ElementMotionInput* input,
                                     DynamicEvaluationRequest request,
                                     Vector* force, Matrix* gd, Matrix* cd,
                                     TrialEvaluationFailure* failure,
                                     Diagnostic* d) {
  const bool derivatives =
      request == DynamicEvaluationRequest::kForcesAndLinearization;
  if (!nodes || !section || !input || !input->finite_state || !force ||
      !failure || force->Size() != 24 ||
      (derivatives ? !gd || !cd : gd || cd))
    return model_storage::Invalid("MITC4 beta damping context", d);
  *failure = TrialEvaluationFailure::kNone;
  // All matrices are exclusive per call. Row-oriented derivative storage
  // permits SetRow's explicit external copy, followed by one final transpose.
  std::unique_ptr<Matrix> frame, x0, x, d0, dirs, v, omega, rate, a, b, h,
      selector, f, fdot, df, drate, dp, left_row, right_row, base_pb, id24,
      gdt, cdt, nodal[2];
  auto status = Make(3, 3, &frame, d);
  for (auto* owner : {&x0, &x, &d0, &dirs, &v, &omega, &rate, &a, &b,
                      &base_pb, &nodal[0], &nodal[1]})
    if (status == StatusCode::kOk) status = Make(4, 3, owner, d);
  for (auto* owner : {&f, &fdot, &df, &drate, &dp})
    if (status == StatusCode::kOk) status = Make(3, 3, owner, d);
  if (status == StatusCode::kOk) status = Make(3, 2, &h, d);
  if (status == StatusCode::kOk) status = Make(3, 2, &selector, d);
  if (status == StatusCode::kOk) status = Make(1, 3, &left_row, d);
  if (status == StatusCode::kOk) status = Make(1, 3, &right_row, d);
  std::unique_ptr<Vector> total, point_force, direction, unit, director, dd,
      w, ddrate, scratch, mapping[5];
  for (auto* owner : {&total, &point_force, &direction})
    if (status == StatusCode::kOk) status = MakeVector(24, owner, d);
  for (auto* owner : {&unit, &director, &dd, &w, &ddrate, &scratch,
                      &mapping[0], &mapping[1], &mapping[2], &mapping[3],
                      &mapping[4]})
    if (status == StatusCode::kOk) status = MakeVector(3, owner, d);
  if (derivatives) {
    for (auto* owner : {&gdt, &cdt, &id24})
      if (status == StatusCode::kOk) status = Make(24, 24, owner, d);
    for (int i = 0; status == StatusCode::kOk && i < 24; ++i)
      status = id24->Set(i, i, 1, d);
  }
  // dp has no live value here; Direction overwrites it after each cross.
  Matrix* cross_scratch = dp.get();
  if (status == StatusCode::kOk)
    status = Coordinates(nodes, input, frame.get(), x0.get(), x.get(), d0.get(),
                           dirs.get(), v.get(), omega.get(), rate.get(),
                           df.get(), d);
  for (int i = 0; status == StatusCode::kOk && i < 2; ++i)
    status = selector->Set(i, i, 1, d);
  if (status == StatusCode::kOk)
    status = Matrix::Multiply(frame.get(), true, selector.get(), false,
                                h.get(), d);
  RayleighPoint point;
  if (status == StatusCode::kOk) status = point.Initialize(d);
  if (status != StatusCode::kOk) return status;
  Matrix* nodal_work[] = {nodal[0].get(), nodal[1].get()};
  Vector* mapping_work[5] = {};
  for (int i = 0; i < 5; ++i) mapping_work[i] = mapping[i].get();
  const auto* material = section->material();
  Mitc4VolumeWorkspace volume_workspace;
  status = volume_workspace.Initialize(d);
  if (status != StatusCode::kOk) return status;
  const double gauss = 1 / std::sqrt(3.);
  for (double r : {-gauss, gauss})
    for (double s : {-gauss, gauss})
      for (double zeta : {-gauss, gauss}) {
        double weight = 0, jacobian = 0;
        status = PointGradients(x0.get(), x.get(), d0.get(), dirs.get(),
            r, s, zeta, section->thickness(), a.get(), b.get(), &weight,
            &jacobian, &volume_workspace, failure, d);
        if (status != StatusCode::kOk) {
          if (d) {
            d->geometry_point_available = true;
            d->natural_xi = r;
            d->natural_eta = s;
            d->natural_zeta = zeta;
          }
          return status;
        }
        status = Gradient(x.get(), dirs.get(), a.get(), b.get(), f.get(), d);
        if (status == StatusCode::kOk)
          status = Gradient(v.get(), rate.get(), a.get(), b.get(),
                               fdot.get(), d);
        if (status == StatusCode::kOk)
          status = point.SetPoint(f.get(), fdot.get(), h.get(), jacobian,
              material->young(), material->poisson(),
              material->rayleigh_beta(), d);
        // df is not live until the following direction loop initializes it.
        if (status == StatusCode::kOk)
          status = NodalForce(point.piola(), dirs.get(), a.get(), b.get(),
              point_force.get(), nodal_work, mapping_work, -1, nullptr,
              nullptr, df.get(), d);
        if (status == StatusCode::kOk)
          status = base_pb->CopyFrom(nodal[1].get(), d);
        if (status == StatusCode::kOk)
          status = total->AddScaled(weight, point_force.get(), d);
        if (status != StatusCode::kOk) return status;
        if (!derivatives) continue;
        for (int variable = 0; variable < 2; ++variable)
          for (int column = 0; column < 24; ++column) {
            const int node = column / 6, c = column % 6;
            status = unit->Scale(0, d);
            if (status == StatusCode::kOk) status = unit->Set(c % 3, 1, d);
            if (status == StatusCode::kOk) status = df->Scale(0, d);
            if (status == StatusCode::kOk) status = drate->Scale(0, d);
            if (status == StatusCode::kOk) status = dd->Scale(0, d);
            if (c < 3 && status == StatusCode::kOk) {
              status = Outer(unit.get(), a.get(), node,
                  variable == 0 ? df.get() : drate.get(), left_row.get(),
                  right_row.get(), scratch.get(), d);
            } else if (status == StatusCode::kOk) {
              status = Row(dirs.get(), node, director.get(), d);
              if (status == StatusCode::kOk)
                status = CrossMitc4RayleighVectors(
                    cross_scratch, unit.get(), director.get(), dd.get(), d);
              if (status == StatusCode::kOk)
                status = Outer(dd.get(), b.get(), node,
                    variable == 0 ? df.get() : drate.get(), left_row.get(),
                    right_row.get(), scratch.get(), d);
              if (variable == 0 && status == StatusCode::kOk) {
                status = Row(omega.get(), node, w.get(), d);
                if (status == StatusCode::kOk)
                  status = CrossMitc4RayleighVectors(
                      cross_scratch, w.get(), dd.get(), ddrate.get(), d);
                if (status == StatusCode::kOk)
                  status = Outer(ddrate.get(), b.get(), node, drate.get(),
                      left_row.get(), right_row.get(), scratch.get(), d);
              }
            }
            if (status == StatusCode::kOk)
              status = point.Direction(df.get(), drate.get(), variable == 1,
                                         dp.get(), d);
            // Direction consumed df; NodalForce may now overwrite it.
            if (status == StatusCode::kOk)
              status = NodalForce(dp.get(), dirs.get(), a.get(), b.get(),
                  direction.get(), nodal_work, mapping_work,
                  variable == 0 && c >= 3 ? node : -1, dd.get(),
                  base_pb.get(), df.get(), d);
            if (status == StatusCode::kOk) status = direction->Scale(weight, d);
            // Add this point's column through an explicit row copy/BLAS sum.
            if (status != StatusCode::kOk) return status;
            Matrix* target = variable == 0 ? gdt.get() : cdt.get();
            int indices[24] = {};
            for (int row = 0; row < 24; ++row) indices[row] = 24 * column + row;
            status = target->GatherTo(indices, 24, point_force.get(), 0, d);
            if (status == StatusCode::kOk)
              status = point_force->AddScaled(1, direction.get(), d);
            if (status == StatusCode::kOk)
              status = target->SetRow(column, point_force.get(), d);
            if (status != StatusCode::kOk) return status;
          }
      }
  status = force->CopyFrom(total.get(), d);
  if (derivatives && status == StatusCode::kOk)
    status = Matrix::Multiply(gdt.get(), true, id24.get(), false, gd, d);
  if (derivatives && status == StatusCode::kOk)
    status = Matrix::Multiply(cdt.get(), true, id24.get(), false, cd, d);
  return status;
}
}  // namespace fesa
