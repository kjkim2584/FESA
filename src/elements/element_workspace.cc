#include "fesa/elements/element_workspace.h"

#include <memory>
#include <new>

#include "elements/beam3d_dynamic_metrics.h"
#include "elements/beam3d_kinematics.h"
#include "elements/beam3d_time_kinematics.h"
#include "elements/incremental_elastic_workspace.h"
#include "fesa/elements/solid_reference_data.h"

namespace fesa {

ElementWorkspace::~ElementWorkspace() {
  delete beam_rayleigh_;
  delete beam_kinematics_;
  delete beam_time_kinematics_;
  if (incremental_)
    for (int i = 0; i < reference_count_; ++i) delete incremental_[i];
  delete[] incremental_;
  delete local_k_;
  delete transform_;
  delete product_;
  delete recovery_;
  delete local_u_;
  delete local_force_;
  delete end_results_;
}

StatusCode ElementWorkspace::Create(int element_dofs, ElementWorkspace** out,
                                    Diagnostic* d) {
  if (!out)
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation,
                         "ElementWorkspace::Create output", d);
  *out = nullptr;
  if (element_dofs != 12)
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation,
                         "ElementWorkspace::Create requires 12 DOFs", d);
  std::unique_ptr<ElementWorkspace> candidate(new (std::nothrow)
                                                  ElementWorkspace);
  if (!candidate)
    return ReportFailure(StatusCode::kAllocationFailure,
                         DiagnosticStage::kAllocation,
                         "ElementWorkspace::Create", d);
  StatusCode status = Matrix::Create(12, 12, &candidate->local_k_, d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Create(12, 12, &candidate->transform_, d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Create(12, 12, &candidate->product_, d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Create(16, 12, &candidate->recovery_, d);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(12, &candidate->local_u_, d);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(12, &candidate->local_force_, d);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(16, &candidate->end_results_, d);
  if (status != StatusCode::kOk) return status;
  *out = candidate.release();
  return StatusCode::kOk;
}

StatusCode ElementWorkspace::CreateForSolid(const SolidReferenceData* stiffness,
                                            const SolidReferenceData* mass,
                                            const SolidReferenceData* geometry,
                                            ElementWorkspace** out,
                                            Diagnostic* d) {
  if (!out)
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation, "Solid context output",
                         d);
  *out = nullptr;
  if (!stiffness && !mass && !geometry)
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation,
                         "Solid context requires reference data", d);
  auto* context = new (std::nothrow) ElementWorkspace;
  if (!context)
    return ReportFailure(StatusCode::kAllocationFailure,
                         DiagnosticStage::kAllocation, "Solid call context", d);
  context->stiffness_reference_ = stiffness;
  context->mass_reference_ = mass;
  context->geometry_reference_ = geometry;
  *out = context;
  return StatusCode::kOk;
}
StatusCode ElementWorkspace::CreateForSolidCatalog(
    const SolidReferenceData* const* references, int count,
    ElementWorkspace** out, Diagnostic* d) {
  if (!out)
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation, "Solid catalog output",
                         d);
  *out = nullptr;
  if (!references || count <= 0)
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation, "Solid catalog entries",
                         d);
  for (int i = 0; i < count; ++i) {
    if (!references[i])
      return ReportFailure(StatusCode::kInvalidArgument,
                           DiagnosticStage::kValidation,
                           "Null solid catalog table", d);
    for (int j = 0; j < i; ++j)
      if (references[i]->formulation() == references[j]->formulation() &&
          references[i]->purpose() == references[j]->purpose())
        return ReportFailure(StatusCode::kInvalidArgument,
                             DiagnosticStage::kValidation,
                             "Duplicate solid catalog key", d);
  }
  auto* result = new (std::nothrow) ElementWorkspace;
  if (!result)
    return ReportFailure(StatusCode::kAllocationFailure,
                         DiagnosticStage::kAllocation, "Solid catalog context",
                         d);
  result->reference_catalog_ = references;
  result->reference_count_ = count;
  *out = result;
  return StatusCode::kOk;
}
StatusCode ElementWorkspace::CreateForIncrementalSolidCatalog(
    const SolidReferenceData* const* references, int count,
    ElementWorkspace** out, Diagnostic* d) {
  auto status = CreateForSolidCatalog(references, count, out, d);
  if (status != StatusCode::kOk) return status;
  (*out)->incremental_ =
      new (std::nothrow) IncrementalElasticWorkspace*[count]();
  if (!(*out)->incremental_)
    status = ReportFailure(StatusCode::kAllocationFailure,
        DiagnosticStage::kAllocation, "Incremental scratch catalog", d);
  for (int i = 0; i < count && status == StatusCode::kOk; ++i) {
    const auto* data = references[i];
    if (data->purpose() != IntegrationPurpose::kStiffness) continue;
    status = IncrementalElasticWorkspace::CreateForSolid(
        data->node_count(), data->point_count(), &(*out)->incremental_[i], d);
  }
  if (status != StatusCode::kOk) {
    delete *out;
    *out = nullptr;
  }
  return status;
}
IncrementalElasticWorkspace* ElementWorkspace::incremental_solid_workspace()
    const {
  return incremental_solid_workspace(SolidFormulation::kTetra4Full);
}
IncrementalElasticWorkspace* ElementWorkspace::incremental_solid_workspace(
    SolidFormulation formulation) const {
  if (!incremental_) return nullptr;
  for (int i = 0; i < reference_count_; ++i)
    if (reference_catalog_[i]->formulation() == formulation &&
        reference_catalog_[i]->purpose() == IntegrationPurpose::kStiffness)
      return incremental_[i];
  return nullptr;
}
const SolidReferenceData* ElementWorkspace::reference(
    SolidFormulation formulation, IntegrationPurpose purpose) const {
  if (reference_catalog_) {
    for (int i = 0; i < reference_count_; ++i) {
      const auto* item = reference_catalog_[i];
      if (item->formulation() == formulation && item->purpose() == purpose)
        return item;
    }
    return nullptr;
  }
  switch (purpose) {
    case IntegrationPurpose::kStiffness:
      return stiffness_reference_;
    case IntegrationPurpose::kMass:
      return mass_reference_;
    case IntegrationPurpose::kGeometry:
      return geometry_reference_;
    case IntegrationPurpose::kTrialGeometry:
      return nullptr;  // Trial samples require the exact-key catalog.
  }
  return nullptr;
}
const SolidReferenceData* ElementWorkspace::stiffness_reference() const {
  return stiffness_reference_;
}
const SolidReferenceData* ElementWorkspace::mass_reference() const {
  return mass_reference_;
}
const SolidReferenceData* ElementWorkspace::geometry_reference() const {
  return geometry_reference_;
}
}  // namespace fesa
