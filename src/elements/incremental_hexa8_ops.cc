#include "elements/incremental_hexa8_ops.h"

#include <cmath>
#include <memory>

#include "elements/incremental_hexa8_workspace.h"
#include "fesa/math/dense_factorization.h"
#include "model/model_storage.h"

namespace fesa {
namespace {
using H = IncrementalHexa8Workspace;
struct Means {
  double volume = 0, initial_volume = 0, midpoint_volume = 0;
  double theta = 0, pressure = 0, committed_pressure = 0;
  double jbar = 0, mu = 0, bulk = 0;
};
bool Shape(const Matrix* x, int r, int c) {
  return x && x->Rows() == r && x->Columns() == c;
}
bool Finite(const Vector* x) {
  if (!x) return false;
  for (int i = 0; i < x->Size(); ++i)
    if (!std::isfinite(x->Data()[i])) return false;
  return true;
}
StatusCode Trace(const Matrix* x, IncrementalElasticWorkspace* w, double* out,
                 Diagnostic* d) {
  const int indices[] = {0, 4, 8};
  auto s = x->GatherTo(indices, 3, w->diagonal(), 0, d);
  if (s == StatusCode::kOk) s = w->diagonal()->Dot(w->ones(), out, d);
  return s;
}
StatusCode TensorRows(IncrementalElasticWorkspace* w, Matrix* out,
                      Diagnostic* d) {
  for (int row = 0; row < 3; ++row) {
    const int indices[] = {3 * row, 3 * row + 1, 3 * row + 2};
    auto s = w->diagonal()->GatherFrom(w->flat(), indices, 3, 0, d);
    if (s == StatusCode::kOk) s = out->SetRow(row, w->diagonal(), d);
    if (s != StatusCode::kOk) return s;
  }
  return StatusCode::kOk;
}
StatusCode Unpack(const Matrix* sn, int g, IncrementalElasticWorkspace* w,
                  Matrix* out, Diagnostic* d) {
  const int packed[] = {0, 3, 4, 3, 1, 5, 4, 5, 2};
  int indices[9] = {};
  for (int i = 0; i < 9; ++i) indices[i] = 6 * g + packed[i];
  auto s = sn->GatherTo(indices, 9, w->flat(), 0, d);
  if (s == StatusCode::kOk) s = TensorRows(w, out, d);
  return s;
}
StatusCode Map(const Matrix* gradient, Matrix* out, Diagnostic* d) {
  // dF_(iA)/du_(ak)=delta_ik H_(Aa), a scalar map coefficient definition.
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      for (int a = 0; a < 8; ++a)
        for (int k = 0; k < 3; ++k) {
          auto s = out->Set(3 * i + j, 3 * a + k,
                            i == k ? gradient->Data()[8 * j + a] : 0, d);
          if (s != StatusCode::kOk) return s;
        }
  return StatusCode::kOk;
}
StatusCode Deformation(const Vector* u, IncrementalElasticWorkspace* w,
                       Matrix* out, Diagnostic* d) {
  auto s = Matrix::MultiplyVector(w->hexa8_workspace()->buffer(H::kMap), false,
                                  u, w->flat(), d);
  if (s == StatusCode::kOk) s = TensorRows(w, out, d);
  if (s == StatusCode::kOk) s = out->AddScaled(1, w->matrix(0), d);
  return s;
}
StatusCode Inverse(const Matrix* f, const Matrix* identity, Matrix* inverse,
                   double* determinant, Diagnostic* d) {
  DenseLuFactorization* raw = nullptr;
  auto s = DenseLuFactorization::Create(f, &raw, d);
  std::unique_ptr<DenseLuFactorization> factor(raw);
  if (s == StatusCode::kOk) s = factor->Determinant(determinant, d);
  if (s == StatusCode::kOk) s = factor->Solve(identity, false, inverse, d);
  return s;
}
StatusCode Point(const Matrix* gradient, double q, const Vector* un,
                 const Vector* u, const Matrix* sn,
                 const ElasticMaterial* material,
                 IncrementalElasticWorkspace* w, int g, bool exact,
                 TrialEvaluationFailure* failure, Diagnostic* d) {
  auto* h = w->hexa8_workspace();
  auto* fn = h->point(g, H::kFn);
  auto* f = w->candidate_f(g);
  auto s = Map(gradient, h->buffer(H::kMap), d);
  if (s == StatusCode::kOk) s = Deformation(un, w, fn, d);
  if (s == StatusCode::kOk) s = Deformation(u, w, f, d);
  if (s == StatusCode::kOk) s = Unpack(sn, g, w, h->buffer(H::kStressTemp), d);
  double pressure = 0;
  if (s == StatusCode::kOk)
    s = Trace(h->buffer(H::kStressTemp), w, &pressure, d);
  if (s == StatusCode::kOk) s = h->coefficient(H::kPn)->Set(g, pressure / 3, d);
  if (s == StatusCode::kOk)
    s = EvaluateIncrementalElasticPoint(
        fn, f, h->buffer(H::kStressTemp), material, w,
        h->point(g, H::kOrdinaryStress), h->point(g, H::kOrdinaryPk1),
        exact ? h->point(g, H::kOrdinaryJacobian) : nullptr, failure, d);
  // Preserve the point API's committed/candidate/backend distinction before
  // any additional midpoint computation can label an error as geometry.
  if (s != StatusCode::kOk) return s;
  double j = 0;
  s = Inverse(f, w->matrix(0), h->point(g, H::kInverse), &j, d);
  if (s == StatusCode::kOk) s = h->coefficient(H::kJ)->Set(g, j, d);
  if (s == StatusCode::kOk) s = h->coefficient(H::kW)->Set(g, q * j, d);
  if (s == StatusCode::kOk) s = h->buffer(H::kMid)->CopyFrom(fn, d);
  if (s == StatusCode::kOk) s = h->buffer(H::kMid)->AddScaled(1, f, d);
  if (s == StatusCode::kOk) s = h->buffer(H::kMid)->Scale(.5, d);
  if (s != StatusCode::kOk) return s;
  double jm = 0;
  s = Inverse(h->buffer(H::kMid), w->matrix(0), h->point(g, H::kMidInverse),
              &jm, d);
  if (s == StatusCode::kSingular || (s == StatusCode::kOk && !(jm > 0))) {
    *failure = TrialEvaluationFailure::kGeometry;
    return ReportFailure(StatusCode::kSingular, DiagnosticStage::kComputation,
                         "Incremental HEX8 nonpositive midpoint mapping", d);
  }
  if (s == StatusCode::kOk) s = h->coefficient(H::kM)->Set(g, q * jm, d);
  if (s == StatusCode::kOk) s = h->buffer(H::kDelta)->CopyFrom(f, d);
  if (s == StatusCode::kOk) s = h->buffer(H::kDelta)->AddScaled(-1, fn, d);
  if (s == StatusCode::kOk)
    s = Matrix::Multiply(h->buffer(H::kDelta), false,
                         h->point(g, H::kMidInverse), false, h->point(g, H::kG),
                         d);
  double theta = 0;
  if (s == StatusCode::kOk) s = Trace(h->point(g, H::kG), w, &theta, d);
  if (s == StatusCode::kOk) s = h->coefficient(H::kTheta)->Set(g, theta, d);
  if (s == StatusCode::kOk)
    s = Matrix::Multiply(h->point(g, H::kInverse), true, gradient, false,
                         h->point(g, H::kSpatialGradient), d);
  return s;
}
StatusCode Average(const Vector* q, const ElasticMaterial* material, H* h,
                   Means* mean, Diagnostic* d) {
  auto s = q->Dot(h->coefficient(H::kOnes), &mean->initial_volume, d);
  if (s == StatusCode::kOk)
    s = h->coefficient(H::kM)->Dot(h->coefficient(H::kOnes),
                                   &mean->midpoint_volume, d);
  if (s == StatusCode::kOk)
    s = h->coefficient(H::kW)->Dot(h->coefficient(H::kOnes), &mean->volume, d);
  if (s == StatusCode::kOk)
    s = h->coefficient(H::kM)->Dot(h->coefficient(H::kTheta), &mean->theta, d);
  if (s == StatusCode::kOk)
    s = q->Dot(h->coefficient(H::kPn), &mean->committed_pressure, d);
  if (s != StatusCode::kOk) return s;
  mean->theta /= mean->midpoint_volume;
  mean->committed_pressure /= mean->initial_volume;
  mean->jbar = mean->volume / mean->initial_volume;
  mean->mu = material->young() / (2 * (1 + material->poisson()));
  mean->bulk = material->young() / (3 * (1 - 2 * material->poisson()));
  mean->pressure = mean->committed_pressure + mean->bulk * mean->theta;
  if (!std::isfinite(mean->pressure) || !std::isfinite(mean->jbar))
    return ReportFailure(StatusCode::kLibraryFailure,
                         DiagnosticStage::kComputation,
                         "Incremental HEX8 average range", d);
  return StatusCode::kOk;
}
StatusCode DerivativeRows(IncrementalElasticWorkspace* w, int g,
                          Diagnostic* d) {
  auto* h = w->hexa8_workspace();
  auto s = h->buffer(H::kTemp)->CopyFrom(w->matrix(0), d);
  if (s == StatusCode::kOk) s = h->buffer(H::kTemp)->Scale(2, d);
  if (s == StatusCode::kOk)
    s = h->buffer(H::kTemp)->AddScaled(-1, h->point(g, H::kG), d);
  if (s == StatusCode::kOk)
    s = Matrix::Multiply(h->point(g, H::kMidInverse), false,
                         h->buffer(H::kTemp), false, h->buffer(H::kTemp2), d);
  if (s != StatusCode::kOk) return s;
  // Scalar derivatives of trace G, det Fm and log J (transpose indices).
  for (int i = 0; i < 3; ++i)
    for (int a = 0; a < 3; ++a) {
      s = h->buffer(H::kThetaRow)
              ->Set(0, 3 * i + a, .5 * h->buffer(H::kTemp2)->Data()[3 * a + i],
                    d);
      if (s == StatusCode::kOk)
        s = h->buffer(H::kMidRow)
                ->Set(0, 3 * i + a,
                      .5 * h->coefficient(H::kM)->Data()[g] *
                          h->point(g, H::kMidInverse)->Data()[3 * a + i],
                      d);
      if (s == StatusCode::kOk)
        s = h->buffer(H::kDivRow)
                ->Set(0, 3 * i + a, h->point(g, H::kInverse)->Data()[3 * a + i],
                      d);
      if (s != StatusCode::kOk) return s;
    }
  return StatusCode::kOk;
}
StatusCode GlobalRows(const Matrix* const* gradients,
                      IncrementalElasticWorkspace* w, const Means* mean,
                      Diagnostic* d) {
  auto* h = w->hexa8_workspace();
  for (int g = 0; g < 8; ++g) {
    auto s = Map(gradients[g], h->buffer(H::kMap), d);
    if (s == StatusCode::kOk) s = DerivativeRows(w, g, d);
    if (s == StatusCode::kOk)
      s = Matrix::MultiplyAdd(h->buffer(H::kDivRow), false, h->buffer(H::kMap),
                              false,
                              h->coefficient(H::kW)->Data()[g] / mean->volume,
                              g ? 1 : 0, h->buffer(H::kGlobalDiv), d);
    if (s == StatusCode::kOk)
      s = Matrix::MultiplyAdd(
          h->buffer(H::kThetaRow), false, h->buffer(H::kMap), false,
          h->coefficient(H::kM)->Data()[g] / mean->midpoint_volume, g ? 1 : 0,
          h->buffer(H::kGlobalTheta), d);
    if (s == StatusCode::kOk)
      s = Matrix::MultiplyAdd(
          h->buffer(H::kMidRow), false, h->buffer(H::kMap), false,
          (h->coefficient(H::kTheta)->Data()[g] - mean->theta) /
              mean->midpoint_volume,
          1, h->buffer(H::kGlobalTheta), d);
    if (s != StatusCode::kOk) return s;
  }
  return StatusCode::kOk;
}
StatusCode FlattenRow(const Matrix* x, IncrementalElasticWorkspace* w,
                      Matrix* row, Diagnostic* d) {
  const int indices[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  auto s = x->GatherTo(indices, 9, w->flat(), 0, d);
  if (s == StatusCode::kOk) s = row->SetRow(0, w->flat(), d);
  return s;
}
StatusCode ExactPoint(const Matrix* gradient, double q, int g,
                      const Means* mean, IncrementalElasticWorkspace* w,
                      Matrix* exact, Diagnostic* d) {
  auto* h = w->hexa8_workspace();
  auto s = Map(gradient, h->buffer(H::kMap), d);
  if (s == StatusCode::kOk) s = DerivativeRows(w, g, d);
  if (s != StatusCode::kOk) return s;
  const double j = h->coefficient(H::kJ)->Data()[g], t = mean->jbar / j;
  const double po = h->coefficient(H::kPn)->Data()[g] +
                    mean->bulk * h->coefficient(H::kTheta)->Data()[g];
  const double eta = j * mean->pressure - mean->jbar * po;
  const double* inverse = h->point(g, H::kInverse)->Data();
  const double* p = h->point(g, H::kOrdinaryPk1)->Data();
  const double* aa = h->point(g, H::kOrdinaryJacobian)->Data();
  // Local dP/dF plus two global rank-one terms, E6 exact derivation.
  for (int i = 0; i < 3; ++i)
    for (int a = 0; a < 3; ++a)
      for (int k = 0; k < 3; ++k)
        for (int b = 0; b < 3; ++b) {
          const int row = 3 * i + a, col = 3 * k + b;
          const double qi = inverse[3 * a + i], qk = inverse[3 * b + k];
          const double value = t * aa[9 * row + col] +
                               (-t * p[row] + j * mean->pressure * qi) * qk -
                               mean->jbar * mean->bulk * qi *
                                   h->buffer(H::kThetaRow)->Data()[col] -
                               eta * inverse[3 * b + i] * inverse[3 * a + k];
          s = h->buffer(H::kLocal)->Set(row, col, value, d);
          if (s != StatusCode::kOk) return s;
        }
  s = Matrix::Multiply(h->buffer(H::kLocal), false, h->buffer(H::kMap), false,
                       h->buffer(H::kProduct), d);
  if (s == StatusCode::kOk)
    s = h->buffer(H::kTemp)->CopyFrom(h->point(g, H::kOrdinaryStress), d);
  if (s == StatusCode::kOk)
    s = h->buffer(H::kTemp)->AddScaled(-po, w->matrix(0), d);
  if (s == StatusCode::kOk) s = h->buffer(H::kTemp)->Scale(mean->jbar, d);
  if (s == StatusCode::kOk)
    s = Matrix::Multiply(h->buffer(H::kTemp), false, h->point(g, H::kInverse),
                         true, h->buffer(H::kTemp2), d);
  if (s == StatusCode::kOk)
    s = FlattenRow(h->buffer(H::kTemp2), w, h->buffer(H::kRankRow), d);
  if (s == StatusCode::kOk)
    s = Matrix::MultiplyAdd(h->buffer(H::kRankRow), true,
                            h->buffer(H::kGlobalDiv), false, 1, 1,
                            h->buffer(H::kProduct), d);
  if (s == StatusCode::kOk)
    s = Matrix::Multiply(w->matrix(0), false, h->point(g, H::kInverse), true,
                         h->buffer(H::kTemp2), d);
  if (s == StatusCode::kOk)
    s = FlattenRow(h->buffer(H::kTemp2), w, h->buffer(H::kRankRow), d);
  if (s == StatusCode::kOk)
    s = Matrix::MultiplyAdd(h->buffer(H::kRankRow), true,
                            h->buffer(H::kGlobalTheta), false, mean->bulk * j,
                            1, h->buffer(H::kProduct), d);
  if (s == StatusCode::kOk)
    s = Matrix::MultiplyAdd(h->buffer(H::kMap), true, h->buffer(H::kProduct),
                            false, q, g ? 1 : 0, exact, d);
  return s;
}
StatusCode RatePoint(int g, double q, const Matrix* sn, bool lagged,
                     const Means* mean, IncrementalElasticWorkspace* w,
                     Matrix* rate, Diagnostic* d) {
  auto* h = w->hexa8_workspace();
  const Matrix* source = h->point(g, H::kOrdinaryStress);
  auto s = StatusCode::kOk;
  if (lagged) {
    s = Unpack(sn, g, w, h->buffer(H::kStressTemp), d);
    source = h->buffer(H::kStressTemp);
  }
  double pressure = 0;
  if (s == StatusCode::kOk) s = Trace(source, w, &pressure, d);
  if (s != StatusCode::kOk) return s;
  pressure /= 3;
  const double avg_p = lagged ? mean->committed_pressure : mean->pressure;
  const double alpha = q * mean->jbar,
               weight = h->coefficient(H::kW)->Data()[g];
  // Scalar fourth-order coefficients: deviatoric work uses modified alpha;
  // pressure variation uses physical w. No hand matrix/vector operations.
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 3; ++k)
        for (int l = 0; l < 3; ++l) {
          const double sil = source->Data()[3 * i + l] - pressure * (i == l);
          const double sik = source->Data()[3 * i + k] - pressure * (i == k);
          const double sjl = source->Data()[3 * j + l] - pressure * (j == l);
          const double sjk = source->Data()[3 * j + k] - pressure * (j == k);
          const double value =
              alpha * ((-2 * mean->mu / 3) * (i == j) * (k == l) +
                       mean->mu * ((i == k) * (j == l) + (i == l) * (j == k)) +
                       (i == k) * sjl -
                       .5 * ((j == k) * sil + (j == l) * sik + (i == k) * sjl +
                             (i == l) * sjk)) +
              avg_p * weight * ((i == j) * (k == l) - (i == l) * (j == k));
          s = h->buffer(H::kLocal)->Set(3 * i + j, 3 * k + l, value, d);
          if (s != StatusCode::kOk) return s;
        }
  s = Map(h->point(g, H::kSpatialGradient), h->buffer(H::kMap), d);
  if (s == StatusCode::kOk)
    s = Matrix::Multiply(h->buffer(H::kLocal), false, h->buffer(H::kMap), false,
                         h->buffer(H::kProduct), d);
  if (s == StatusCode::kOk)
    s = Matrix::MultiplyAdd(h->buffer(H::kMap), true, h->buffer(H::kProduct),
                            false, 1, g ? 1 : 0, rate, d);
  return s;
}
}  // namespace
StatusCode EvaluateIncrementalHexa8Points(
    const Matrix* const* gradients, const Vector* q, const Vector* un,
    const Vector* u, const Matrix* sn, const ElasticMaterial* material,
    IncrementalElasticWorkspace* w, ElementEvaluationRequest request,
    Vector* force, Matrix* exact, Matrix* rate, Matrix* stress,
    TrialEvaluationFailure* failure, Diagnostic* d, bool lagged,
    int* failed_point) {
  if (failure) *failure = TrialEvaluationFailure::kNone;
  if (failed_point) *failed_point = -1;
  const bool want_force = request == ElementEvaluationRequest::kBoth ||
                          request == ElementEvaluationRequest::kInternalForce;
  const bool want_rate = request == ElementEvaluationRequest::kBoth ||
                         request == ElementEvaluationRequest::kTangent;
  if (!failure || !w || w->node_count() != 8 || w->point_count() != 8 ||
      !w->hexa8_workspace() || !gradients || !q || q->Size() != 8 || !un ||
      !u || un->Size() != 24 || u->Size() != 24 || !Finite(un) || !Finite(u) ||
      !material || !Shape(sn, 8, 6) || !Shape(stress, 8, 6) || sn == stress ||
      (!want_force && !want_rate) ||
      (want_force ? (!force || force->Size() != 24) : force != nullptr) ||
      (want_rate ? !Shape(rate, 24, 24) : rate != nullptr) ||
      (exact && (!want_rate || !Shape(exact, 24, 24) || exact == rate)) ||
      (force && (force == un || force == u || force == w->point_force())))
    return model_storage::Invalid("Incremental coupled HEX8 input", d);
  auto* h = w->hexa8_workspace();
  for (int i = 0; i < IncrementalElasticWorkspace::kElementMatrixCount; ++i)
    if ((exact && exact == w->element_matrix(i)) ||
        (rate && rate == w->element_matrix(i)))
      return model_storage::Invalid("Incremental HEX8 output alias", d);
  for (int i = 0; i < H::kCoefficientCount; ++i)
    if (q == h->coefficient(static_cast<H::Coefficient>(i)))
      return model_storage::Invalid("Incremental HEX8 weight scratch alias", d);
  for (int g = 0; g < 8; ++g) {
    if (!Shape(gradients[g], 3, 8) || !std::isfinite(q->Data()[g]) ||
        !(q->Data()[g] > 0)) {
      if (failed_point) *failed_point = g;
      return model_storage::Invalid("Incremental HEX8 initial H/q", d);
    }
    for (int p = 0; p < 8; ++p)
      if (gradients[g] == h->point(p, H::kSpatialGradient))
        return model_storage::Invalid("Incremental HEX8 gradient scratch alias",
                                      d);
    for (int i = 0; i < 24; ++i)
      if (!std::isfinite(gradients[g]->Data()[i])) {
        if (failed_point) *failed_point = g;
        return model_storage::Invalid("Incremental HEX8 initial gradient", d);
      }
    for (int i = 0; i < 6; ++i)
      if (!std::isfinite(sn->Data()[6 * g + i])) {
        if (failed_point) *failed_point = g;
        return model_storage::Invalid("Incremental HEX8 committed stress", d);
      }
  }
  for (int g = 0; g < 8; ++g) {
    const auto s = Point(gradients[g], q->Data()[g], un, u, sn, material, w, g,
                         exact != nullptr, failure, d);
    if (s != StatusCode::kOk) {
      if (failed_point) *failed_point = g;
      return s;
    }
  }
  Means mean;
  auto s = Average(q, material, h, &mean, d);
  if (s == StatusCode::kOk && want_rate) s = GlobalRows(gradients, w, &mean, d);
  if (s != StatusCode::kOk) return s;
  for (int g = 0; g < 8; ++g) {
    const double theta = h->coefficient(H::kTheta)->Data()[g];
    const double j = h->coefficient(H::kJ)->Data()[g];
    const double po = h->coefficient(H::kPn)->Data()[g] + mean.bulk * theta;
    s = h->buffer(H::kStressTemp)->CopyFrom(h->point(g, H::kOrdinaryStress), d);
    if (s == StatusCode::kOk)
      s = h->buffer(H::kStressTemp)
              ->AddScaled(mean.bulk * (mean.theta - theta), w->matrix(0), d);
    const int components[] = {0, 4, 8, 1, 2, 5};
    if (s == StatusCode::kOk)
      s = h->buffer(H::kStressTemp)->GatherTo(components, 6, w->tensor(), 0, d);
    if (s == StatusCode::kOk) s = stress->SetRow(g, w->tensor(), d);
    if (s == StatusCode::kOk)
      s = h->point(g, H::kModifiedF)->CopyFrom(w->candidate_f(g), d);
    if (s == StatusCode::kOk)
      s = h->point(g, H::kModifiedF)->Scale(std::cbrt(mean.jbar / j), d);
    if (s == StatusCode::kOk && want_force) {
      s = h->buffer(H::kWeightedStress)
              ->CopyFrom(h->point(g, H::kOrdinaryStress), d);
      if (s == StatusCode::kOk)
        s = h->buffer(H::kWeightedStress)->Scale(mean.jbar, d);
      if (s == StatusCode::kOk)
        s = h->buffer(H::kWeightedStress)
                ->AddScaled(j * mean.pressure - mean.jbar * po, w->matrix(0),
                            d);
      if (s == StatusCode::kOk)
        s = Matrix::Multiply(h->buffer(H::kWeightedStress), false,
                             h->point(g, H::kInverse), true,
                             h->buffer(H::kEffectivePk1), d);
      if (s == StatusCode::kOk) s = Map(gradients[g], h->buffer(H::kMap), d);
      const int indices[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
      if (s == StatusCode::kOk)
        s = h->buffer(H::kEffectivePk1)->GatherTo(indices, 9, w->flat(), 0, d);
      if (s == StatusCode::kOk)
        s = Matrix::MultiplyVector(h->buffer(H::kMap), true, w->flat(),
                                   w->point_force(), d);
      if (s == StatusCode::kOk) s = w->point_force()->Scale(q->Data()[g], d);
      if (s == StatusCode::kOk)
        s = g ? force->AddScaled(1, w->point_force(), d)
              : force->CopyFrom(w->point_force(), d);
    }
    if (s == StatusCode::kOk && exact)
      s = ExactPoint(gradients[g], q->Data()[g], g, &mean, w, exact, d);
    if (s == StatusCode::kOk && want_rate)
      s = RatePoint(g, q->Data()[g], sn, lagged, &mean, w, rate, d);
    if (s != StatusCode::kOk) {
      if (failed_point) *failed_point = g;
      return s;
    }
  }
  if (want_rate) {
    const double p = lagged ? mean.committed_pressure : mean.pressure;
    s = Matrix::MultiplyAdd(h->buffer(H::kGlobalDiv), true,
                            h->buffer(H::kGlobalDiv), false,
                            mean.volume * (mean.bulk - p), 1, rate, d);
  }
  return s;
}
}  // namespace fesa
