#include "elements/incremental_solid_ops.h"

#include <cmath>
#include <memory>

#include "fesa/math/dense_factorization.h"
#include "model/model_storage.h"

namespace fesa {
namespace {
enum ElementBuffer {
  kMap, kPointJacobian, kProduct, kExact, kSpatialGradient, kCommittedF,
  kCandidateF, kCommittedStress, kStress, kPk1
};
static_assert(kPk1 + 1 == IncrementalElasticWorkspace::kElementMatrixCount);
bool Shape(const Matrix* matrix, int rows, int columns) {
  return matrix && matrix->Rows() == rows && matrix->Columns() == columns;
}
bool Finite(const Vector* vector) {
  if (!vector) return false;
  for (int i = 0; i < vector->Size(); ++i)
    if (!std::isfinite(vector->Data()[i])) return false;
  return true;
}
StatusCode TensorRows(IncrementalElasticWorkspace* workspace, Matrix* out,
                       Diagnostic* d) {
  for (int row = 0; row < 3; ++row) {
    const int indices[] = {3 * row, 3 * row + 1, 3 * row + 2};
    auto status = workspace->diagonal()->GatherFrom(
        workspace->flat(), indices, 3, 0, d);
    if (status == StatusCode::kOk)
      status = out->SetRow(row, workspace->diagonal(), d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}
StatusCode Deformation(const Vector* u, IncrementalElasticWorkspace* w,
                        Matrix* out, Diagnostic* d) {
  auto status = Matrix::MultiplyVector(w->element_matrix(kMap), false, u,
                                       w->flat(), d);
  if (status == StatusCode::kOk) status = TensorRows(w, out, d);
  if (status == StatusCode::kOk) status = out->AddScaled(1, w->matrix(0), d);
  return status;
}
// IE5-R: current-volume rate correction, independent of the exact Jacobian.
// B maps global nodal increments to vec(grad_x du), row3*i+j.
// Scalar entries below define A, not a handwritten matrix multiplication.
StatusCode RateCorrection(const Matrix* initial_gradient, double q,
                          const ElasticMaterial* material,
                          IncrementalElasticWorkspace* w, Matrix* out,
                          Diagnostic* d, bool use_committed_geometric_stress) {
  Matrix* gradient = w->element_matrix(kSpatialGradient);
  DenseLuFactorization* raw_factor = nullptr;
  auto status = DenseLuFactorization::Create(w->element_matrix(kCandidateF),
                                             &raw_factor, d);
  std::unique_ptr<DenseLuFactorization> factor(raw_factor);
  double determinant = 0;
  if (status == StatusCode::kOk) status = factor->Determinant(&determinant, d);
  if (status == StatusCode::kOk)
    status = factor->Solve(initial_gradient, true, gradient, d);
  if (status != StatusCode::kOk) return status;
  const double mu = material->young() / (2 * (1 + material->poisson()));
  const double lambda =
      material->young() * material->poisson() /
      ((1 + material->poisson()) * (1 - 2 * material->poisson()));
  const double* sigma =
      w->element_matrix(use_committed_geometric_stress ? kCommittedStress
                                                       : kStress)
          ->Data();
  auto* coefficients = w->element_matrix(kPointJacobian);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        for (int l = 0; l < 3; ++l) {
          const double value =
              lambda * (i == j) * (k == l) +
              mu * ((i == k) * (j == l) + (i == l) * (j == k)) +
              (i == k) * sigma[3 * j + l] -
              .5 * ((j == k) * sigma[3 * i + l] + (j == l) * sigma[3 * i + k] +
                    (i == k) * sigma[3 * j + l] + (i == l) * sigma[3 * j + k]);
          status = coefficients->Set(3 * i + j, 3 * k + l, value, d);
          if (status != StatusCode::kOk) return status;
        }
      }
    }
  }
  auto* map = w->element_matrix(kMap);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int a = 0; a < initial_gradient->Columns(); ++a) {
        for (int k = 0; k < 3; ++k) {
          status = map->Set(3 * i + j, 3 * a + k,
                            i == k ? gradient->Data()[initial_gradient->Columns() * j + a] : 0, d);
          if (status != StatusCode::kOk) return status;
        }
      }
    }
  }
  status = Matrix::Multiply(coefficients, false, map, false,
                            w->element_matrix(kProduct), d);
  if (status == StatusCode::kOk)
    status =
        Matrix::Multiply(map, true, w->element_matrix(kProduct), false, out, d);
  if (status == StatusCode::kOk) status = out->Scale(q * determinant, d);
  return status;
}

}  // namespace
StatusCode EvaluateIncrementalSolidPoints(const Matrix* const* gradients,
    const Vector* weights, const Vector* un, const Vector* u, const Matrix* sn,
    const ElasticMaterial* material, IncrementalElasticWorkspace* w,
    ElementEvaluationRequest request, Vector* force, Matrix* exact,
    Matrix* symmetric, Matrix* stress, TrialEvaluationFailure* failure,
    Diagnostic* d, bool use_committed_geometric_stress, int* failed_point) {
  if (failure) *failure = TrialEvaluationFailure::kNone;
  if (failed_point) *failed_point = -1;
  const bool want_force = request == ElementEvaluationRequest::kInternalForce ||
                          request == ElementEvaluationRequest::kBoth;
  const bool want_tangent = request == ElementEvaluationRequest::kTangent ||
                            request == ElementEvaluationRequest::kBoth;
  const int n = u ? u->Size() / 3 : 0;
  const int nd = 3 * n, ng = weights ? weights->Size() : 0;
  if (w && w->node_count() == 8)
    return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                         "Incremental Hexa8 requires coupled operator", d);
  if (!failure || !w || !material || (!want_force && !want_tangent) ||
      !gradients || !weights || ng != w->point_count() ||
      n != w->node_count() || !un || un->Size() != nd || !u ||
      u->Size() != nd || !Finite(un) || !Finite(u) ||
      !Shape(sn, ng, 6) || !Shape(stress, ng, 6) || sn == stress ||
      (want_force ? (!force || force->Size() != nd) : force != nullptr) ||
      (want_tangent ? !Shape(symmetric, nd, nd) : symmetric != nullptr) ||
      (exact && (!want_tangent || !Shape(exact, nd, nd))) ||
      (exact && exact == symmetric) || force == un || force == u ||
      (force && force == w->point_force()))
    return model_storage::Invalid("Incremental ordinary operator input", d);
  auto invalid_point = [failed_point, d](int point, const char* operation) {
    if (failed_point) *failed_point = point;
    return model_storage::Invalid(operation, d);
  };
  for (int g = 0; g < ng; ++g) {
    if (!Shape(gradients[g], 3, n) ||
        !std::isfinite(weights->Data()[g]) || !(weights->Data()[g] > 0))
      return invalid_point(g, "Incremental initial point shape/weight");
    for (int i = 0; i < 3 * n; ++i)
      if (!std::isfinite(gradients[g]->Data()[i]))
        return invalid_point(g, "Incremental initial gradient finite");
    for (int i = 0; i < 6; ++i)
      if (!std::isfinite(sn->Data()[6 * g + i]))
        return invalid_point(g, "Incremental committed stress finite");
  }
  for (int i = 0; i < IncrementalElasticWorkspace::kElementMatrixCount; ++i)
    if ((exact && exact == w->element_matrix(i)) ||
        (symmetric && symmetric == w->element_matrix(i)))
      return model_storage::Invalid("Incremental tangent scratch alias", d);
  auto fail = [failed_point](StatusCode status, int point) {
    if (failed_point) *failed_point = point;
    return status;
  };
  for (int g = 0; g < ng; ++g) {
    const double q = weights->Data()[g];
    auto* b = w->element_matrix(kMap);
    // IE5 L_(iA,ak)=delta_ik H_aA,g. Rebuild from initial gradients at every
    // point: RateCorrection overwrites this map with current gradients.
    for (int i = 0; i < 3; ++i)
      for (int a = 0; a < n; ++a)
        for (int axis = 0; axis < 3; ++axis)
          for (int k = 0; k < 3; ++k) {
            const auto status = b->Set(3 * i + axis, 3 * a + k,
                i == k ? gradients[g]->Data()[n * axis + a] : 0, d);
            if (status != StatusCode::kOk) return fail(status, g);
          }
    auto status = Deformation(un, w, w->element_matrix(kCommittedF), d);
    if (status == StatusCode::kOk)
      status = Deformation(u, w, w->element_matrix(kCandidateF), d);
    const int packed[] = {0, 3, 4, 3, 1, 5, 4, 5, 2};
    int unpack[9] = {};
    for (int i = 0; i < 9; ++i) unpack[i] = 6 * g + packed[i];
    if (status == StatusCode::kOk)
      status = sn->GatherTo(unpack, 9, w->flat(), 0, d);
    if (status == StatusCode::kOk)
      status = TensorRows(w, w->element_matrix(kCommittedStress), d);
    if (status != StatusCode::kOk) return fail(status, g);
    status = EvaluateIncrementalElasticPoint(w->element_matrix(kCommittedF),
        w->element_matrix(kCandidateF), w->element_matrix(kCommittedStress),
        material, w, w->element_matrix(kStress), w->element_matrix(kPk1),
        want_tangent ? w->element_matrix(kPointJacobian) : nullptr, failure, d);
    if (status == StatusCode::kOk)
      status = w->candidate_f(g)->CopyFrom(w->element_matrix(kCandidateF), d);
    if (status != StatusCode::kOk) return fail(status, g);
    if (want_force) {
      const int indices[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
      status = w->element_matrix(kPk1)->GatherTo(indices, 9, w->flat(), 0, d);
      if (status == StatusCode::kOk)
        status = Matrix::MultiplyVector(b, true, w->flat(), w->point_force(), d);
      if (status == StatusCode::kOk) status = w->point_force()->Scale(q, d);
      if (status == StatusCode::kOk)
        status = g == 0 ? force->CopyFrom(w->point_force(), d)
                        : force->AddScaled(1, w->point_force(), d);
      if (status != StatusCode::kOk) return fail(status, g);
    }
    if (want_tangent) {
      auto* point_k = w->element_matrix(kExact);
      status = Matrix::Multiply(w->element_matrix(kPointJacobian), false,
          b, false, w->element_matrix(kProduct), d);
      if (status == StatusCode::kOk)
        status = Matrix::Multiply(b, true, w->element_matrix(kProduct), false,
                                  point_k, d);
      if (status == StatusCode::kOk) status = point_k->Scale(q, d);
      if (status == StatusCode::kOk && exact)
        status = g == 0 ? exact->CopyFrom(point_k, d)
                        : exact->AddScaled(1, point_k, d);
      if (status == StatusCode::kOk)
        status = RateCorrection(gradients[g], q, material, w, point_k, d,
                                use_committed_geometric_stress);
      if (status == StatusCode::kOk)
        status = g == 0 ? symmetric->CopyFrom(point_k, d)
                        : symmetric->AddScaled(1, point_k, d);
      if (status != StatusCode::kOk) return fail(status, g);
    }
    const int components[] = {0, 4, 8, 1, 2, 5};
    status = w->element_matrix(kStress)->GatherTo(components, 6, w->tensor(), 0, d);
    if (status == StatusCode::kOk) status = stress->SetRow(g, w->tensor(), d);
    if (status != StatusCode::kOk) return fail(status, g);
  }
  return StatusCode::kOk;
}
}  // namespace fesa
