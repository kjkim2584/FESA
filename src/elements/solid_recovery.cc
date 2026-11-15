#include "elements/solid_recovery.h"

#include <cmath>
#include <memory>

#include "elements/solid_access.h"
#include "fesa/math/dense_factorization.h"

namespace fesa {
namespace {
StatusCode InvalidRecovery(Diagnostic* d) {
  return ReportFailure(StatusCode::kInvalidArgument,
                       DiagnosticStage::kValidation,
                       "Solid tensor recovery", d);
}
// F11, five-point full Pyramid: square rational interpolation. One LU
// serves twelve tensor RHS. N at the natural nodes is identity, including
// the apex limit, so coefficients are already the element nodal results.
StatusCode RecoverPyramid(const SolidReferenceData* data,
                          const Matrix* ip_stress, const Matrix* ip_strain,
                          Matrix* nodal_stress, Matrix* nodal_strain,
                          Diagnostic* d) {
  Matrix *rhs = nullptr, *coefficients = nullptr;
  auto status = Matrix::Create(5, 12, &rhs, d);
  std::unique_ptr<Matrix> ro(rhs);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Create(5, 12, &coefficients, d);
  std::unique_ptr<Matrix> co(coefficients);
  if (status != StatusCode::kOk) return status;
  status = rhs->CopyBlockFrom(ip_stress, 0, 0, 5, 6, 0, 0, d);
  if (status != StatusCode::kOk) return status;
  status = rhs->CopyBlockFrom(ip_strain, 0, 0, 5, 6, 0, 6, d);
  if (status != StatusCode::kOk) return status;
  DenseLuFactorization* factor = nullptr;
  ObserveSolidEvaluation(SolidEvaluationEvent::kRecoveryFactor);
  status = DenseLuFactorization::Create(data->shape_values(), &factor, d);
  std::unique_ptr<DenseLuFactorization> fo(factor);
  if (status != StatusCode::kOk) return status;
  status = factor->Solve(rhs, false, coefficients, d);
  if (status != StatusCode::kOk) return status;
  status = nodal_stress->CopyBlockFrom(coefficients, 0, 0, 5, 6, 0, 0, d);
  if (status != StatusCode::kOk) return status;
  return nodal_strain->CopyBlockFrom(coefficients, 0, 6, 5, 6, 0, 0, d);
}
}  // namespace
StatusCode ValidateSolidTensorRecovery(const SolidReferenceData* data,
    const Matrix* ip_stress, const Matrix* ip_strain,
    const Matrix* nodal_stress, const Matrix* nodal_strain, Diagnostic* d) {
  if (!data || data->purpose() != IntegrationPurpose::kStiffness)
    return InvalidRecovery(d);
  const Matrix* fields[] = {ip_stress, ip_strain, nodal_stress, nodal_strain};
  for (int i = 0; i < 4; ++i) {
    const int rows = i < 2 ? data->point_count() : data->node_count();
    if (!fields[i] || fields[i]->Columns() != 6 ||
        fields[i]->Rows() != rows ||
        fields[i] == data->shape_values()) return InvalidRecovery(d);
    for (int j = 0; j < i; ++j)
      if (fields[i] == fields[j]) return InvalidRecovery(d);
  }
  return StatusCode::kOk;
}
StatusCode ExtrapolateSolidTensors(const SolidReferenceData* data,
    const Matrix* ip_stress, const Matrix* ip_strain,
    Matrix* nodal_stress, Matrix* nodal_strain, Diagnostic* d) {
  auto status = ValidateSolidTensorRecovery(data, ip_stress, ip_strain,
                                             nodal_stress, nodal_strain, d);
  if (status != StatusCode::kOk) return status;
  for (int i = 0; i < 6 * data->point_count(); ++i) {
    if (!std::isfinite(ip_stress->Data()[i]) ||
        !std::isfinite(ip_strain->Data()[i])) return InvalidRecovery(d);
  }
  ObserveSolidEvaluation(SolidEvaluationEvent::kRecovery);
  if (data->formulation() == SolidFormulation::kPyramid5Full)
    return RecoverPyramid(data, ip_stress, ip_strain, nodal_stress,
                          nodal_strain, d);
  if (data->point_count() == 1) {
    Matrix* constant = nullptr;
    status = Matrix::Create(data->node_count(), 1, &constant, d);
    std::unique_ptr<Matrix> constant_owner(constant);
    if (status != StatusCode::kOk) return status;
    for (int node = 0; node < data->node_count(); ++node) {
      status = constant->Set(node, 0, 1, d);
      if (status != StatusCode::kOk) return status;
    }
    status =
        Matrix::Multiply(constant, false, ip_stress, false, nodal_stress, d);
    if (status != StatusCode::kOk) return status;
    return Matrix::Multiply(constant, false, ip_strain, false, nodal_strain,
                            d);
  }
  // One call-local interpolation LU serves all twelve tensor RHS.
  const bool hexa = data->formulation() == SolidFormulation::kHexa8Full;
  const bool wedge = data->formulation() == SolidFormulation::kWedge6;
  const int samples = hexa ? 8 : wedge ? 2 : 4;
  Matrix *basis = nullptr, *rhs = nullptr, *coefficients = nullptr,
         *evaluation = nullptr, *result = nullptr;
  status = Matrix::Create(samples, samples, &basis, d);
  std::unique_ptr<Matrix> basis_owner(basis);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Create(samples, 12, &rhs, d);
  std::unique_ptr<Matrix> rhs_owner(rhs);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Create(samples, 12, &coefficients, d);
  std::unique_ptr<Matrix> coefficients_owner(coefficients);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Create(data->node_count(), samples, &evaluation, d);
  std::unique_ptr<Matrix> evaluation_owner(evaluation);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Create(data->node_count(), 12, &result, d);
  std::unique_ptr<Matrix> result_owner(result);
  if (status != StatusCode::kOk) return status;
  if (hexa) {
    status = basis->CopyFrom(data->shape_values(), d);
    if (status != StatusCode::kOk) return status;
  } else if (wedge) {
    for (int p = 0; p < 2; ++p) {
      status = basis->Set(p, 0, 1, d);
      if (status != StatusCode::kOk) return status;
      status =
          basis->Set(p, 1, data->natural_coordinates()->Data()[3 * p + 2], d);
      if (status != StatusCode::kOk) return status;
    }
  } else {
    for (int p = 0; p < 4; ++p) {
      const double* x = data->natural_coordinates()->Data() + 3 * p;
      const double l[] = {1 - x[0] - x[1] - x[2], x[0], x[1], x[2]};
      for (int i = 0; i < 4; ++i) {
        status = basis->Set(p, i, l[i], d);
        if (status != StatusCode::kOk) return status;
      }
    }
  }
  status = rhs->CopyBlockFrom(ip_stress, 0, 0, samples, 6, 0, 0, d);
  if (status != StatusCode::kOk) return status;
  status = rhs->CopyBlockFrom(ip_strain, 0, 0, samples, 6, 0, 6, d);
  if (status != StatusCode::kOk) return status;
  if (hexa) {
    for (int node = 0; node < 8; ++node) {
      status = evaluation->Set(node, node, 1, d);
      if (status != StatusCode::kOk) return status;
    }
  } else if (wedge) {
    for (int node = 0; node < 6; ++node) {
      status = evaluation->Set(node, 0, 1, d);
      if (status != StatusCode::kOk) return status;
      status = evaluation->Set(node, 1, node < 3 ? -1 : 1, d);
      if (status != StatusCode::kOk) return status;
    }
  } else {
    const int edges[6][2] = {{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}};
    for (int node = 0; node < data->node_count(); ++node)
      for (int i = 0; i < 4; ++i) {
        const double value =
            node < 4
                ? (node == i ? 1 : 0)
                : ((edges[node - 4][0] == i || edges[node - 4][1] == i) ? .5
                                                                        : 0);
        status = evaluation->Set(node, i, value, d);
        if (status != StatusCode::kOk) return status;
      }
  }
  DenseLuFactorization* factor = nullptr;
  ObserveSolidEvaluation(SolidEvaluationEvent::kRecoveryFactor);
  status = DenseLuFactorization::Create(basis, &factor, d);
  std::unique_ptr<DenseLuFactorization> factor_owner(factor);
  if (status != StatusCode::kOk) return status;
  status = factor->Solve(rhs, false, coefficients, d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Multiply(evaluation, false, coefficients, false, result, d);
  if (status != StatusCode::kOk) return status;
  status = nodal_stress->CopyBlockFrom(result, 0, 0, data->node_count(),
                                       6, 0, 0, d);
  if (status != StatusCode::kOk) return status;
  return nodal_strain->CopyBlockFrom(result, 0, 6, data->node_count(),
                                      6, 0, 0, d);
}
}  // namespace fesa
