#include "elements/solid_nonlinear_ops.h"

#include <cmath>
#include <memory>

#include "elements/solid_access.h"
#include "elements/solid_stvk.h"
#include "elements/solid_nonlinear_recovery.h"
#include "elements/solid_recovery.h"
#include "fesa/math/dense_factorization.h"

namespace fesa {
namespace {
StatusCode NewMatrix(int rows, int columns, std::unique_ptr<Matrix>* owner,
                     Diagnostic* d) {
  Matrix* raw = nullptr;
  const auto status = Matrix::Create(rows, columns, &raw, d);
  owner->reset(raw);
  return status;
}
StatusCode NewVector(int size, std::unique_ptr<Vector>* owner, Diagnostic* d) {
  Vector* raw = nullptr;
  const auto status = Vector::Create(size, &raw, d);
  owner->reset(raw);
  return status;
}

// B_F maps node-major translations into row-major vec9(F-I), M1/M4.
StatusCode Deformation(const Matrix* gradient, const Vector* u, Matrix* b,
                       Matrix* f, Diagnostic* d) {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int a = 0; a < gradient->Columns(); ++a) {
        const auto status = b->Set(3 * i + j, 3 * a + i,
            gradient->Data()[gradient->Columns() * j + a], d);
        if (status != StatusCode::kOk) return status;
      }
    }
  }
  std::unique_ptr<Vector> values, row;
  auto status = NewVector(9, &values, d);
  if (status != StatusCode::kOk) return status;
  status = NewVector(3, &row, d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::MultiplyVector(b, false, u, values.get(), d);
  if (status != StatusCode::kOk) return status;
  for (int i = 0; i < 3; ++i) {
    const int indices[] = {3 * i, 3 * i + 1, 3 * i + 2};
    status = row->GatherFrom(values.get(), indices, 3, 0, d);
    if (status != StatusCode::kOk) return status;
    status = f->SetRow(i, row.get(), d);
    if (status != StatusCode::kOk) return status;
  }
  std::unique_ptr<Matrix> identity;
  status = NewMatrix(3, 3, &identity, d);
  if (status != StatusCode::kOk) return status;
  for (int i = 0; i < 3; ++i) {
    status = identity->Set(i, i, 1, d);
    if (status != StatusCode::kOk) return status;
  }
  return f->AddScaled(1, identity.get(), d);
}

// Call-local owner; no trial state or constitutive cache escapes evaluation.
struct FbarPoint {
  FbarPoint() = default;
  FbarPoint(const FbarPoint&) = delete;
  FbarPoint& operator=(const FbarPoint&) = delete;
  FbarPoint(FbarPoint&&) = delete;
  FbarPoint& operator=(FbarPoint&&) = delete;
  ~FbarPoint() {
    delete f;
    delete b;
    delete inverse_transpose;
    delete modified;
    delete p;
    delete a;
    delete ell;
    delete t;
  }
  Matrix* f = nullptr;
  Matrix* b = nullptr;
  Matrix* inverse_transpose = nullptr;
  Matrix* modified = nullptr;
  Matrix* p = nullptr;
  Matrix* a = nullptr;
  Matrix* ell = nullptr;
  Matrix* t = nullptr;
  double determinant = 0;
  double scale = 0;
  double contraction = 0;
};

StatusCode FbarNumericalFailure(Diagnostic* d) {
  return ReportFailure(StatusCode::kLibraryFailure,
                       DiagnosticStage::kComputation, "F-bar nonfinite", d);
}

StatusCode FlattenRow(const Matrix* source, Vector* entries, Matrix* row,
                       Diagnostic* d) {
  const int indices[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  auto status = source->GatherTo(indices, 9, entries, 0, d);
  if (status != StatusCode::kOk) return status;
  return row->SetRow(0, entries, d);
}

StatusCode PrepareFbarPoint(const Matrix* g, const Vector* u,
                            FbarPoint* point, Diagnostic* d) {
  auto status = Matrix::Create(3, 3, &point->f, d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Create(9, 24, &point->b, d);
  if (status != StatusCode::kOk) return status;
  status = Deformation(g, u, point->b, point->f, d);
  if (status != StatusCode::kOk) return status;
  DenseLuFactorization* raw = nullptr;
  status = DenseLuFactorization::Create(point->f, &raw, d);
  std::unique_ptr<DenseLuFactorization> factor(raw);
  if (status == StatusCode::kSingular) {
    status = ReportFailure(StatusCode::kEquilibriumFailure,
                         DiagnosticStage::kComputation, "F-bar singular F", d);
    if (d) {
      d->geometry_quantity = GeometryQuantity::kDeformationJacobian;
      d->geometry_limit = 0;
    }
    return status;
  }
  if (status != StatusCode::kOk) return status;
  status = factor->Determinant(&point->determinant, d);
  if (status != StatusCode::kOk) return status;
  if (!(point->determinant > 0)) {
    status = ReportFailure(StatusCode::kEquilibriumFailure,
                         DiagnosticStage::kComputation,
                         "F-bar nonpositive J", d);
    if (d) {
      d->geometry_quantity = GeometryQuantity::kDeformationJacobian;
      d->geometry_value_available = true;
      d->geometry_value = point->determinant;
      d->geometry_limit = 0;
    }
    return status;
  }
  std::unique_ptr<Matrix> identity;
  status = NewMatrix(3, 3, &identity, d);
  if (status != StatusCode::kOk) return status;
  for (int i = 0; i < 3; ++i) {
    status = identity->Set(i, i, 1, d);
    if (status != StatusCode::kOk) return status;
  }
  status = Matrix::Create(3, 3, &point->inverse_transpose, d);
  if (status != StatusCode::kOk) return status;
  return factor->Solve(identity.get(), true, point->inverse_transpose, d);
}

// M6/M8: H_p=e_i G_a^T implies ell_p=D_ia and t_pq=D_ja D_ib,
// where D=F^-T G. The outer products/reductions use common math operations.
StatusCode PrepareFbarDerivatives(FbarPoint* point, const Matrix* g,
                                  Diagnostic* d) {
  std::unique_ptr<Matrix> spatial;
  auto status = NewMatrix(3, 8, &spatial, d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Multiply(point->inverse_transpose, false, g, false,
                             spatial.get(), d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Create(1, 24, &point->ell, d);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<Vector> entries;
  status = NewVector(24, &entries, d);
  if (status != StatusCode::kOk) return status;
  int indices[24];
  for (int p = 0; p < 24; ++p) indices[p] = (p % 3) * 8 + p / 3;
  status = spatial->GatherTo(indices, 24, entries.get(), 0, d);
  if (status != StatusCode::kOk) return status;
  status = point->ell->SetRow(0, entries.get(), d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Create(24, 24, &point->t, d);
  if (status != StatusCode::kOk) return status;
  for (int p = 0; p < 24; ++p) {
    for (int q = 0; q < 24; ++q) {
      const double coefficient = spatial->Data()[(q % 3) * 8 + p / 3] *
                                  spatial->Data()[(p % 3) * 8 + q / 3];
      if (!std::isfinite(coefficient)) return FbarNumericalFailure(d);
      status = point->t->Set(p, q, coefficient, d);
      if (status != StatusCode::kOk) return status;
    }
  }
  return StatusCode::kOk;
}

StatusCode AccumulateFbarTangent(const FbarPoint* point, const Matrix* mean,
    const Matrix* centered_mean, double weight, Matrix* tangent,
    Diagnostic* d) {
  std::unique_ptr<Matrix> b, modified_b, f_row, p_row, h, second, ab;
  auto status = NewMatrix(1, 24, &b, d);
  if (status != StatusCode::kOk) return status;
  status = b->CopyFrom(mean, d);
  if (status != StatusCode::kOk) return status;
  status = b->AddScaled(-1, point->ell, d);
  if (status != StatusCode::kOk) return status;
  status = b->Scale(1. / 3, d);
  if (status != StatusCode::kOk) return status;
  status = NewMatrix(9, 24, &modified_b, d);
  if (status != StatusCode::kOk) return status;
  status = modified_b->CopyFrom(point->b, d);
  if (status != StatusCode::kOk) return status;
  status = NewMatrix(1, 9, &f_row, d);
  if (status != StatusCode::kOk) return status;
  status = NewMatrix(1, 9, &p_row, d);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<Vector> entries;
  status = NewVector(9, &entries, d);
  if (status != StatusCode::kOk) return status;
  status = FlattenRow(point->f, entries.get(), f_row.get(), d);
  if (status != StatusCode::kOk) return status;
  status = FlattenRow(point->p, entries.get(), p_row.get(), d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::MultiplyAdd(f_row.get(), true, b.get(), false, 1, 1,
                                modified_b.get(), d);
  if (status != StatusCode::kOk) return status;
  status = modified_b->Scale(point->scale, d);
  if (status != StatusCode::kOk) return status;
  status = NewMatrix(9, 24, &ab, d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Multiply(point->a, false, modified_b.get(), false,
                             ab.get(), d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::MultiplyAdd(modified_b.get(), true, ab.get(), false,
                                weight, 1, tangent, d);
  if (status != StatusCode::kOk) return status;
  // Full M8 second derivative; no symmetry averaging or finite differences.
  status = NewMatrix(24, 24, &second, d);
  if (status != StatusCode::kOk) return status;
  status = second->CopyFrom(centered_mean, d);
  if (status != StatusCode::kOk) return status;
  status = second->AddScaled(1, point->t, d);
  if (status != StatusCode::kOk) return status;
  status = second->Scale(1. / 3, d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::MultiplyAdd(b.get(), true, b.get(), false, 1, 1,
                                second.get(), d);
  if (status != StatusCode::kOk) return status;
  status = second->Scale(point->contraction, d);
  if (status != StatusCode::kOk) return status;
  status = NewMatrix(1, 24, &h, d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Multiply(p_row.get(), false, point->b, false, h.get(), d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::MultiplyAdd(h.get(), true, b.get(), false, point->scale, 1,
                                second.get(), d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::MultiplyAdd(b.get(), true, h.get(), false, point->scale, 1,
                                second.get(), d);
  if (status != StatusCode::kOk) return status;
  return tangent->AddScaled(weight, second.get(), d);
}

// Called only after the public internal boundary has checked all shapes,
// request/alias contracts and all positive initial weights.
StatusCode EvaluateFbar(const Matrix* const* gradients, const Vector* weights,
    const Vector* u, const ElasticMaterial* material, double* energy,
    Vector* force, Matrix* tangent, Matrix* const* material_gradients,
    TrialEvaluationFailure* failure, Diagnostic* d) {
  FbarPoint points[8];
  std::unique_ptr<Vector> determinants, ones, densities, contractions;
  for (auto* owner : {&determinants, &ones, &densities}) {
    const auto status = NewVector(8, owner, d);
    if (status != StatusCode::kOk) return status;
  }
  if (force) {
    const auto status = NewVector(8, &contractions, d);
    if (status != StatusCode::kOk) return status;
  }
  // Every actual J must pass before any average is formed (M5).
  for (int g = 0; g < 8; ++g) {
    auto status = PrepareFbarPoint(gradients[g], u, &points[g], d);
    if (status == StatusCode::kEquilibriumFailure) {
      *failure = TrialEvaluationFailure::kGeometry;
      if (d) {
        const double location = 1 / std::sqrt(3.);
        d->geometry_point_available = true;
        d->natural_xi = (g & 1) ? location : -location;
        d->natural_eta = (g & 2) ? location : -location;
        d->natural_zeta = (g & 4) ? location : -location;
      }
    }
    if (status != StatusCode::kOk) return status;
    status = determinants->Set(g, points[g].determinant, d);
    if (status != StatusCode::kOk) return status;
    status = ones->Set(g, 1, d);
    if (status != StatusCode::kOk) return status;
  }
  double initial_volume = 0, current_volume = 0;
  auto status = weights->Dot(ones.get(), &initial_volume, d);
  if (status != StatusCode::kOk) return status;
  status = weights->Dot(determinants.get(), &current_volume, d);
  if (status != StatusCode::kOk) return status;
  if (!(initial_volume > 0) || !(current_volume > 0))
    return FbarNumericalFailure(d);
  const double average = current_volume / initial_volume;
  if (!(average > 0) || !std::isfinite(average)) return FbarNumericalFailure(d);
  std::unique_ptr<Vector> p_entries, f_entries;
  status = NewVector(9, &p_entries, d);
  if (status != StatusCode::kOk) return status;
  status = NewVector(9, &f_entries, d);
  if (status != StatusCode::kOk) return status;
  const int indices[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  for (int g = 0; g < 8; ++g) {
    auto* point = &points[g];
    point->scale = std::cbrt(average / point->determinant);
    if (!(point->scale > 0) || !std::isfinite(point->scale))
      return FbarNumericalFailure(d);
    status = Matrix::Create(3, 3, &point->modified, d);
    if (status != StatusCode::kOk) return status;
    status = point->modified->CopyFrom(point->f, d);
    if (status != StatusCode::kOk) return status;
    status = point->modified->Scale(point->scale, d);
    if (status != StatusCode::kOk) return status;
    status = Matrix::Create(3, 3, &point->p, d);
    if (status != StatusCode::kOk) return status;
    if (tangent) {
      status = Matrix::Create(9, 9, &point->a, d);
      if (status != StatusCode::kOk) return status;
    }
    double density = 0;
    status = EvaluateStvkPoint(point->modified, material,
        energy ? &density : nullptr, nullptr, point->p, point->a, d);
    if (status == StatusCode::kEquilibriumFailure) {
      // All actual determinants already passed. Losing that property while
      // rescaling is arithmetic failure, not retryable physical geometry.
      return ReportFailure(StatusCode::kLibraryFailure,
                           DiagnosticStage::kComputation,
                           "F-bar modified determinant loss", d);
    }
    if (status != StatusCode::kOk) return status;
    status = densities->Set(g, density, d);
    if (status != StatusCode::kOk) return status;
    status = point->p->GatherTo(indices, 9, p_entries.get(), 0, d);
    if (status != StatusCode::kOk) return status;
    status = point->modified->GatherTo(indices, 9, f_entries.get(), 0, d);
    if (status != StatusCode::kOk) return status;
    status = p_entries->Dot(f_entries.get(), &point->contraction, d);
    if (status != StatusCode::kOk) return status;
    if (force) {
      status = contractions->Set(g, point->contraction, d);
      if (status != StatusCode::kOk) return status;
    }
  }
  if (force) {
    ObserveSolidEvaluation(SolidEvaluationEvent::kForce);
    double total_contraction = 0;
    status = weights->Dot(contractions.get(), &total_contraction, d);
    if (status != StatusCode::kOk) return status;
    std::unique_ptr<Vector> sum, contribution;
    status = NewVector(24, &sum, d);
    if (status != StatusCode::kOk) return status;
    status = NewVector(24, &contribution, d);
    if (status != StatusCode::kOk) return status;
    std::unique_ptr<Matrix> effective;
    status = NewMatrix(3, 3, &effective, d);
    if (status != StatusCode::kOk) return status;
    for (int g = 0; g < 8; ++g) {
      const auto* point = &points[g];
      const double r =
          (point->determinant * total_contraction / current_volume -
           point->contraction) / 3;
      if (!std::isfinite(r)) return FbarNumericalFailure(d);
      status = effective->CopyFrom(point->p, d);
      if (status != StatusCode::kOk) return status;
      status = effective->Scale(point->scale, d);
      if (status != StatusCode::kOk) return status;
      status = effective->AddScaled(r, point->inverse_transpose, d);
      if (status != StatusCode::kOk) return status;
      status = effective->GatherTo(indices, 9, p_entries.get(), 0, d);
      if (status != StatusCode::kOk) return status;
      status = Matrix::MultiplyVector(point->b, true, p_entries.get(),
                                       contribution.get(), d);
      if (status != StatusCode::kOk) return status;
      status = sum->AddScaled(weights->Data()[g], contribution.get(), d);
      if (status != StatusCode::kOk) return status;
    }
    status = force->CopyFrom(sum.get(), d);
    if (status != StatusCode::kOk) return status;
  }
  if (tangent) {
    ObserveSolidEvaluation(SolidEvaluationEvent::kStiffness);
    std::unique_ptr<Matrix> mean, centered_mean, sum;
    status = NewMatrix(1, 24, &mean, d);
    if (status != StatusCode::kOk) return status;
    status = NewMatrix(24, 24, &centered_mean, d);
    if (status != StatusCode::kOk) return status;
    status = NewMatrix(24, 24, &sum, d);
    if (status != StatusCode::kOk) return status;
    for (int g = 0; g < 8; ++g) {
      auto* point = &points[g];
      status = PrepareFbarDerivatives(point, gradients[g], d);
      if (status != StatusCode::kOk) return status;
      const double fraction = weights->Data()[g] * point->determinant /
                              current_volume;
      if (!std::isfinite(fraction)) return FbarNumericalFailure(d);
      status = mean->AddScaled(fraction, point->ell, d);
      if (status != StatusCode::kOk) return status;
      status = Matrix::MultiplyAdd(point->ell, true, point->ell, false,
                                    fraction, 1, centered_mean.get(), d);
      if (status != StatusCode::kOk) return status;
      status = centered_mean->AddScaled(-fraction, point->t, d);
      if (status != StatusCode::kOk) return status;
    }
    status = Matrix::MultiplyAdd(mean.get(), true, mean.get(), false, -1, 1,
                                  centered_mean.get(), d);
    if (status != StatusCode::kOk) return status;
    for (int g = 0; g < 8; ++g) {
      status = AccumulateFbarTangent(&points[g], mean.get(),
          centered_mean.get(), weights->Data()[g], sum.get(), d);
      if (status != StatusCode::kOk) return status;
    }
    status = tangent->CopyFrom(sum.get(), d);
    if (status != StatusCode::kOk) return status;
  }
  if (material_gradients) {
    for (int g = 0; g < 8; ++g) {
      status = material_gradients[g]->CopyFrom(points[g].modified, d);
      if (status != StatusCode::kOk) return status;
    }
  }
  if (energy) return weights->Dot(densities.get(), energy, d);
  return StatusCode::kOk;
}

StatusCode EvaluatePoints(SolidFormulation formulation,
    const Matrix* const* gradients, const Vector* weights,
    const Vector* translations, const ElasticMaterial* material,
    ElementEvaluationRequest request, double* energy, Vector* force,
    Matrix* tangent, Matrix* const* material_gradients,
    TrialEvaluationFailure* failure, Diagnostic* d) {
  if (failure) *failure = TrialEvaluationFailure::kNone;
  const bool fbar = formulation == SolidFormulation::kHexa8Full;
  const int nodes = fbar ? 8 : formulation == SolidFormulation::kTetra10Full
      ? 10 : formulation == SolidFormulation::kWedge6 ? 6
      : formulation == SolidFormulation::kPyramid5Full ? 5 : 4;
  const int count = fbar ? 8 : nodes == 10 ? 4 : nodes == 6 ? 2
      : nodes == 5 ? 5 : 1;
  const int dofs = 3 * nodes;
  if (!fbar && formulation != SolidFormulation::kTetra4Full &&
      formulation != SolidFormulation::kTetra10Full &&
      formulation != SolidFormulation::kWedge6 &&
      formulation != SolidFormulation::kPyramid5Full) {
    return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                         "nonlinear fixed-gradient formulation", d);
  }
  const bool need_force = request == ElementEvaluationRequest::kInternalForce ||
                          request == ElementEvaluationRequest::kBoth;
  const bool need_tangent = request == ElementEvaluationRequest::kTangent ||
                            request == ElementEvaluationRequest::kBoth;
  bool valid = failure && gradients && gradients[0] && weights &&
      translations && material && (need_force || need_tangent) &&
      ((need_force && force && force->Size() == dofs) ||
       (!need_force && !force)) &&
      ((need_tangent && tangent && tangent->Rows() == dofs &&
        tangent->Columns() == dofs) || (!need_tangent && !tangent));
  if (valid) {
    valid = weights->Size() == count && translations->Size() == dofs &&
            force != translations && force != weights;
    for (int g = 0; valid && g < count; ++g) {
      valid = gradients[g] && gradients[g]->Rows() == 3 &&
              gradients[g]->Columns() == nodes && tangent != gradients[g] &&
              std::isfinite(weights->Data()[g]) && weights->Data()[g] > 0;
      for (int i = 0; valid && i < dofs; ++i) {
        valid = std::isfinite(gradients[g]->Data()[i]) &&
                std::isfinite(translations->Data()[i]);
      }
    }
  }
  if (!valid) {
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation,
                         "EvaluateSolidNonlinearPoints", d);
  }
  if (fbar) {
    return EvaluateFbar(gradients, weights, translations, material, energy,
                         force, tangent, material_gradients, failure, d);
  }
  std::unique_ptr<Matrix> b, f, p, a, ab;
  std::unique_ptr<Vector> entries, point_force, densities;
  auto status = NewMatrix(9, dofs, &b, d);
  if (status != StatusCode::kOk) return status;
  status = NewMatrix(3, 3, &f, d);
  if (status != StatusCode::kOk) return status;
  if (need_force) {
    status = NewMatrix(3, 3, &p, d);
    if (status != StatusCode::kOk) return status;
    status = NewVector(9, &entries, d);
    if (status != StatusCode::kOk) return status;
    status = NewVector(dofs, &point_force, d);
    if (status != StatusCode::kOk) return status;
    ObserveSolidEvaluation(SolidEvaluationEvent::kForce);
  }
  if (need_tangent) {
    status = NewMatrix(9, 9, &a, d);
    if (status != StatusCode::kOk) return status;
    status = NewMatrix(9, dofs, &ab, d);
    if (status != StatusCode::kOk) return status;
    ObserveSolidEvaluation(SolidEvaluationEvent::kStiffness);
  }
  if (energy) {
    status = NewVector(count, &densities, d);
    if (status != StatusCode::kOk) return status;
  }
  // M4: integrate each actual-F StVK response with its initial measure.
  for (int g = 0; g < count; ++g) {
    status = Deformation(gradients[g], translations, b.get(), f.get(), d);
    if (status != StatusCode::kOk) return status;
    double density = 0;
    status = EvaluateStvkPoint(f.get(), material, energy ? &density : nullptr,
                               nullptr, p.get(), a.get(), d);
    if (status == StatusCode::kEquilibriumFailure) {
      *failure = TrialEvaluationFailure::kGeometry;
      if (d) {
        d->geometry_point_available = true;
        d->natural_xi = d->natural_eta = d->natural_zeta = .25;
        // Same ordered TET4/WEDGE2/PYR5_FELIPPA points as quadrature.md.
        if (nodes == 10) {
          const double high = (5 + 3 * std::sqrt(5.)) / 20;
          const double low = (5 - std::sqrt(5.)) / 20;
          d->natural_xi = g == 2 ? high : low;
          d->natural_eta = g == 1 ? high : low;
          d->natural_zeta = g == 0 ? high : low;
        } else if (nodes == 6) {
          d->natural_xi = d->natural_eta = 1. / 3;
          d->natural_zeta = (g == 0 ? -1 : 1) / std::sqrt(3.);
        } else if (nodes == 5) {
          const int signs[4][2] = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}};
          d->natural_xi = g == 4 ? 0 : signs[g][0] * std::sqrt(32. / 135);
          d->natural_eta = g == 4 ? 0 : signs[g][1] * std::sqrt(32. / 135);
          d->natural_zeta = g == 4 ? 7. / 10 : 1. / 6;
        }
      }
    }
    if (status != StatusCode::kOk) return status;
    const double weight = weights->Data()[g];
    if (need_force) {
      const int indices[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
      status = p->GatherTo(indices, 9, entries.get(), 0, d);
      if (status != StatusCode::kOk) return status;
      status = Matrix::MultiplyVector(b.get(), true, entries.get(),
                                       point_force.get(), d);
      if (status != StatusCode::kOk) return status;
      if (g == 0) {
        status = force->CopyFrom(point_force.get(), d);
        if (status == StatusCode::kOk) status = force->Scale(weight, d);
      } else {
        status = force->AddScaled(weight, point_force.get(), d);
      }
      if (status != StatusCode::kOk) return status;
    }
    if (need_tangent) {
      status = Matrix::Multiply(a.get(), false, b.get(), false, ab.get(), d);
      if (status != StatusCode::kOk) return status;
      status = Matrix::MultiplyAdd(b.get(), true, ab.get(), false, weight,
                                   g == 0 ? 0 : 1, tangent, d);
      if (status != StatusCode::kOk) return status;
    }
    if (material_gradients) {
      status = material_gradients[g]->CopyFrom(f.get(), d);
      if (status != StatusCode::kOk) return status;
    }
    if (energy) {
      status = densities->Set(g, density, d);
      if (status != StatusCode::kOk) return status;
    }
  }
  if (energy) return weights->Dot(densities.get(), energy, d);
  return StatusCode::kOk;
}
}  // namespace

StatusCode EvaluateSolidNonlinearPoints(SolidFormulation formulation,
    const Matrix* const* gradients, const Vector* weights,
    const Vector* translations, const ElasticMaterial* material,
    ElementEvaluationRequest request, double* energy, Vector* force,
    Matrix* tangent, TrialEvaluationFailure* failure, Diagnostic* d) {
  return EvaluatePoints(formulation, gradients, weights, translations, material,
      request, energy, force, tangent, nullptr, failure, d);
}

StatusCode RecoverSolidNonlinearPoints(const SolidReferenceData* data,
    const Matrix* const* gradients, const Vector* weights,
    const Vector* translations, const ElasticMaterial* material, Vector* force,
    Matrix* ip_stress, Matrix* ip_strain, Matrix* nodal_stress,
    Matrix* nodal_strain, TrialEvaluationFailure* failure, Diagnostic* d) {
  if (failure) *failure = TrialEvaluationFailure::kNone;
  auto status = ValidateSolidTensorRecovery(data, ip_stress, ip_strain,
                                              nodal_stress, nodal_strain, d);
  if (status != StatusCode::kOk) return status;
  if (data->formulation() != SolidFormulation::kTetra4Full &&
      data->formulation() != SolidFormulation::kHexa8Full &&
      data->formulation() != SolidFormulation::kTetra10Full &&
      data->formulation() != SolidFormulation::kWedge6 &&
      data->formulation() != SolidFormulation::kPyramid5Full)
    return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                         "Solid nonlinear recovery formulation", d);
  std::unique_ptr<Matrix> owners[8];
  Matrix* outputs[8] = {};
  const Matrix* inputs[8] = {};
  for (int g = 0; g < data->point_count(); ++g) {
    status = NewMatrix(3, 3, &owners[g], d);
    if (status != StatusCode::kOk) return status;
    outputs[g] = owners[g].get();
    inputs[g] = outputs[g];
  }
  status = EvaluatePoints(data->formulation(), gradients, weights, translations,
      material, ElementEvaluationRequest::kInternalForce, nullptr, force,
      nullptr, outputs, failure, d);
  if (status != StatusCode::kOk) return status;
  return RecoverSolidPointOutputs(data, inputs, material, ip_stress, ip_strain,
                                  nodal_stress, nodal_strain, failure, d);
}
}  // namespace fesa
