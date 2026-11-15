#include "fesa/model/thermal_element.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <memory>
#include <new>

#include "fesa/model/domain.h"
#include "fesa/model/film_property_table.h"
#include "model/model_storage.h"

namespace fesa {
StatusCode CollectThermalNodeIds(const Domain* domain, std::vector<int>* ids,
                                 Diagnostic* d) {
  if (!domain || !ids)
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation,
                         "CollectThermalNodeIds null argument", d);
  try {
    std::vector<int> candidate(0);
    for (const auto entry : *domain->ThermalElements()) {
      const auto* element = entry.second;
      for (int local = 0; local < element->node_count(); ++local)
        candidate.push_back(element->GetNode(local)->Id());
    }
    std::sort(candidate.begin(), candidate.end());
    candidate.erase(std::unique(candidate.begin(), candidate.end()),
                    candidate.end());
    ids->swap(candidate);
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return ReportFailure(StatusCode::kAllocationFailure,
                         DiagnosticStage::kValidation,
                         "CollectThermalNodeIds allocation", d);
  }
}
ThermalElement::ThermalElement(int id) : id_(id) {}
int ThermalElement::Id() const { return id_; }
int ThermalElement::dof_count() const { return node_count() * dofs_per_node(); }
StatusCode ThermalElement::GetShellThickness(double*, Diagnostic* d) const {
  return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                       "Element has no shell section thickness", d);
}
StatusCode ThermalElement::ComputeSurfaceWeights(const char*, Vector*,
                                                 Diagnostic* d) const {
  return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                       "Thermal face weights unsupported", d);
}
namespace {
// Nodal radiation law and its full tangent, formulation8.6/8.6.1. Scalar
// constitutive coefficients are separate from Matrix::Diagonal's operation.
StatusCode EvaluateRadiation(const ThermalElement* element,
                             const ThermalLoadContext* load,
                             const Vector* theta, Vector* flux,
                             Matrix* jacobian, Diagnostic* d) {
  const int n = element->dof_count();
  if (!theta || theta->Size() != n || !std::isfinite(load->emissivity) ||
      load->emissivity < 0 || load->emissivity > 1 ||
      !std::isfinite(load->sigma) || load->sigma <= 0 ||
      !std::isfinite(load->absolute_zero) ||
      !std::isfinite(load->ambient_temperature) ||
      !std::isfinite(load->ambient_baseline) ||
      load->ambient_temperature < load->absolute_zero ||
      load->ambient_baseline < load->absolute_zero ||
      !std::isfinite(load->step_scale) || load->step_scale < 0 ||
      load->step_scale > 1 ||
      (load->step_scale < 1 && !load->step_start_temperature) ||
      (load->step_start_temperature &&
       load->step_start_temperature->Size() != n) ||
      (load->previous_outward_flux && load->previous_outward_flux->Size() != n))
    return model_storage::Invalid("Thermal radiation local arguments", d);
  auto status = element->ComputeSurfaceWeights(load->face_label, flux, d);
  Vector* raw = nullptr;
  if (status == StatusCode::kOk && jacobian)
    status = Vector::Create(n, &raw, d);
  std::unique_ptr<Vector> diagonal(raw);
  const double ambient = load->ambient_temperature - load->absolute_zero;
  const double baseline = load->ambient_baseline - load->absolute_zero;
  const double ambient_squared = ambient * ambient;
  const double baseline_squared = baseline * baseline;
  for (int i = 0; i < n && status == StatusCode::kOk; ++i) {
    const double area = flux->Data()[i];
    if (area == 0) continue;
    const double absolute = theta->Data()[i] - load->absolute_zero;
    if (!std::isfinite(absolute) ||
        (absolute < 0 && (!load->extrapolated_free_dofs ||
                          load->extrapolated_free_dofs[i] != 1))) {
      if (d) {
        d->element_id = element->Id();
        d->node_id = element->GetNode(i / element->dofs_per_node())->Id();
        d->condition_dof = 11 + i % element->dofs_per_node();
      }
      return ReportFailure(StatusCode::kLibraryFailure,
                           DiagnosticStage::kComputation,
                           "Invalid radiation absolute trial temperature", d);
    }
    // Supplied frozen state is strict even for STEP/carried loads (A=1).
    double start = 0;
    if (load->step_start_temperature) {
      start = load->step_start_temperature->Data()[i] - load->absolute_zero;
      if (!std::isfinite(start) || start < 0) {
        if (d) {
          d->element_id = element->Id();
          d->node_id = element->GetNode(i / element->dofs_per_node())->Id();
          d->condition_dof = 11 + i % element->dofs_per_node();
        }
        return ReportFailure(StatusCode::kLibraryFailure,
                             DiagnosticStage::kComputation,
                             "Invalid radiation absolute start temperature", d);
      }
    }
    const double squared = absolute * absolute;
    const double weight = area * load->emissivity * load->sigma;
    double outward =
        weight * (squared * squared - ambient_squared * ambient_squared);
    if (load->step_scale < 1) {
      const double start_squared = start * start;
      const double prior = load->previous_outward_flux
                               ? load->previous_outward_flux->Data()[i]
                               : 0;
      outward += (1 - load->step_scale) *
                 (prior - weight * (start_squared * start_squared -
                                    baseline_squared * baseline_squared));
    }
    const double tangent = 4 * weight * squared * absolute;
    if (!std::isfinite(outward) || !std::isfinite(tangent))
      return ReportFailure(StatusCode::kLibraryFailure,
                           DiagnosticStage::kComputation,
                           "Nonfinite thermal radiation law", d);
    status = flux->Set(i, -outward, d);
    if (status == StatusCode::kOk && diagonal)
      status = diagonal->Set(i, tangent, d);
  }
  if (status == StatusCode::kOk && jacobian)
    status = Matrix::Diagonal(diagonal.get(), jacobian, d);
  return status;
}
// Pure local constitutive coefficients; Matrix::Diagonal performs the matrix
// operation. Flux initially contains area weights, then receives signed Q.
StatusCode EvaluateFilm(const ThermalElement* element,
                        const ThermalLoadContext* load, const Vector* theta,
                        Vector* flux, Matrix* jacobian, Diagnostic* d) {
  const int n = element->dof_count();
  if (!theta || theta->Size() != n || !std::isfinite(load->coefficient) ||
      load->coefficient < 0 || !std::isfinite(load->sink_temperature) ||
      !std::isfinite(load->sink_baseline) || !std::isfinite(load->step_scale) ||
      load->step_scale < 0 || load->step_scale > 1 ||
      (load->step_scale < 1 && (!load->step_start_temperature ||
                                load->step_start_temperature->Size() != n)) ||
      (load->previous_outward_flux &&
       load->previous_outward_flux->Size() != n) ||
      (load->first_accepted_temperature &&
       load->first_accepted_temperature->Size() != n))
    return model_storage::Invalid("Thermal film local arguments", d);
  for (int i = 0; i < n; ++i)
    if (!std::isfinite(theta->Data()[i]))
      return ReportFailure(StatusCode::kLibraryFailure,
                           DiagnosticStage::kComputation,
                           "Nonfinite thermal film trial temperature", d);
  auto status = element->ComputeSurfaceWeights(load->face_label, flux, d);
  Vector* raw = nullptr;
  if (status == StatusCode::kOk && jacobian)
    status = Vector::Create(n, &raw, d);
  std::unique_ptr<Vector> diagonal(raw);
  for (int i = 0; i < n && status == StatusCode::kOk; ++i) {
    const double area = flux->Data()[i];
    if (area == 0) continue;
    const double temperature = theta->Data()[i];
    double h = load->coefficient, slope = 0;
    if (load->film_property)
      status = load->film_property->Evaluate(temperature, &h, &slope, d);
    if (status != StatusCode::kOk) return status;
    const double difference = temperature - load->sink_temperature;
    double outward = area * h * difference;
    if (load->step_scale < 1) {
      const double start = load->step_start_temperature->Data()[i];
      const double prior = load->previous_outward_flux
                               ? load->previous_outward_flux->Data()[i]
                               : 0;
      double correction_h = h, ignored_slope = 0;
      if (load->film_property && load->first_accepted_temperature) {
        status = load->film_property->Evaluate(
            load->first_accepted_temperature->Data()[i], &correction_h,
            &ignored_slope, d);
        if (status != StatusCode::kOk) return status;
      }
      outward += (1 - load->step_scale) *
                 (prior - area * correction_h * (start - load->sink_baseline));
    }
    const double tangent = area * (h + slope * difference);
    if (!std::isfinite(outward) || !std::isfinite(tangent))
      return ReportFailure(StatusCode::kLibraryFailure,
                           DiagnosticStage::kComputation,
                           "Nonfinite thermal film law", d);
    status = flux->Set(i, -outward, d);
    if (status == StatusCode::kOk && diagonal)
      status = diagonal->Set(i, tangent, d);
  }
  if (status == StatusCode::kOk && jacobian)
    status = Matrix::Diagonal(diagonal.get(), jacobian, d);
  return status;
}
}  // namespace
StatusCode ThermalElement::EvaluateLoadContribution(
    const ThermalLoadContext* load, const Vector* theta, Vector* flux,
    Matrix* jacobian, Diagnostic* d) const {
  if (!load || !load->face_label || !std::isfinite(load->magnitude) || !flux ||
      flux->Size() != dof_count() ||
      (jacobian &&
       (jacobian->Rows() != dof_count() || jacobian->Columns() != dof_count())))
    return model_storage::Invalid("Thermal flux contribution arguments", d);
  if (load->kind == ThermalConditionKind::kFilm)
    return EvaluateFilm(this, load, theta, flux, jacobian, d);
  if (load->kind == ThermalConditionKind::kRadiation)
    return EvaluateRadiation(this, load, theta, flux, jacobian, d);
  if (load->kind != ThermalConditionKind::kBodyFlux &&
      load->kind != ThermalConditionKind::kSurfaceFlux)
    return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                         "Thermal load contribution unsupported", d);
  if (load->kind == ThermalConditionKind::kBodyFlux &&
      std::strcmp(load->face_label, "BF") != 0)
    return model_storage::Invalid("Thermal body flux label", d);
  auto status = load->kind == ThermalConditionKind::kBodyFlux
                    ? ComputeCapacityWeights(flux, d)
                    : ComputeSurfaceWeights(load->face_label, flux, d);
  if (status == StatusCode::kOk) status = flux->Scale(load->magnitude, d);
  if (status == StatusCode::kOk && jacobian) {
    Matrix* raw = nullptr;
    status = Matrix::Create(dof_count(), dof_count(), &raw, d);
    std::unique_ptr<Matrix> zero(raw);
    if (status == StatusCode::kOk) status = jacobian->CopyFrom(zero.get(), d);
  }
  return status;
}
}  // namespace fesa
