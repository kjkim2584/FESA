#include "fesa/analysis/heat_flux_residual.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <new>
#include <vector>

#include "analysis/heat_flux_measurement.h"
#include "analysis/thermal_assembler.h"
#include "fesa/analysis/thermal_history.h"
#include "fesa/io/control_settings.h"

namespace fesa {
namespace {
StatusCode Failure(StatusCode code, const char* message, Diagnostic* d) {
  const auto stage =
      code == StatusCode::kAllocationFailure ? DiagnosticStage::kAllocation
      : code == StatusCode::kLibraryFailure  ? DiagnosticStage::kComputation
                                             : DiagnosticStage::kValidation;
  return ReportFailure(code, stage, message, d);
}
}  // namespace
StatusCode MeasureThermalFluxOrigins(const Vector* element,
                                     const Vector* external,
                                     HeatFluxHistory* history,
                                     const FieldControls* c,
                                     HeatFluxResidual* out, Diagnostic* d) {
  if (history) history->DiscardCandidate();
  if (!history || !c || !out)
    return Failure(StatusCode::kInvalidArgument,
                   "Null thermal flux measurement argument", d);
  if (!std::isfinite(c->epsilon_l) || c->epsilon_l < 0 ||
      (c->qu_explicit && (!std::isfinite(c->qu) || c->qu <= 0)) ||
      (!c->qu_explicit && c->q0_explicit &&
       (!std::isfinite(c->q0) || c->q0 <= 0)))
    return Failure(StatusCode::kInputError,
                   "Invalid thermal flux measurement controls", d);
  const int ne = element ? element->Size() : 0,
            nx = external ? external->Size() : 0;
  if (nx > std::numeric_limits<int>::max() - ne)
    return Failure(StatusCode::kInvalidArgument,
                   "Thermal flux origin count overflow", d);
  int external_count = 0;
  for (int i = 0; i < ne; ++i)
    if (!std::isfinite(element->Data()[i]))
      return Failure(StatusCode::kLibraryFailure,
                     "Nonfinite unassembled thermal element flux", d);
  for (int i = 0; i < nx; ++i) {
    if (!std::isfinite(external->Data()[i]))
      return Failure(StatusCode::kLibraryFailure,
                     "Nonfinite external thermal flux", d);
    if (external->Data()[i] != 0) ++external_count;
  }
  double maximum = 0, value = 0;
  int index = 0;
  auto status = StatusCode::kOk;
  if (element) status = element->MaximumAbs(&maximum, &index, d);
  if (status != StatusCode::kOk) return status;
  if (external) status = external->MaximumAbs(&value, &index, d);
  if (status != StatusCode::kOk) return status;
  maximum = std::max(maximum, value);
  const double previous = history->PreviousMaximumMean();
  const double threshold = c->epsilon_l * previous;
  const bool all = maximum < .1 * previous || threshold == 0;
  if (!all && !std::isfinite(threshold))
    return Failure(StatusCode::kLibraryFailure,
                   "Thermal active flux threshold overflow", d);
  double mean = 0;
  int selected_count = 0;
  try {
    std::vector<int> selected_element(0), selected_external(0);
    for (int i = 0; i < ne; ++i)
      if (all || std::abs(element->Data()[i]) >= threshold)
        selected_element.push_back(i);
    for (int i = 0; i < nx; ++i)
      if (external->Data()[i] != 0 &&
          (all || std::abs(external->Data()[i]) >= threshold))
        selected_external.push_back(i);
    selected_count =
        static_cast<int>(selected_element.size() + selected_external.size());
    if (selected_count == 0 && ne + external_count > 0)
      return Failure(StatusCode::kLibraryFailure,
                     "Thermal active flux selection is empty", d);
    if (selected_count > 0) {
      Vector* raw = nullptr;
      status = Vector::Create(selected_count, &raw, d);
      std::unique_ptr<Vector> active(raw);
      if (status != StatusCode::kOk) return status;
      const int count = static_cast<int>(selected_element.size());
      if (count > 0)
        status =
            active->GatherFrom(element, selected_element.data(), count, 0, d);
      if (status == StatusCode::kOk && !selected_external.empty())
        status = active->GatherFrom(external, selected_external.data(),
                                    static_cast<int>(selected_external.size()),
                                    count, d);
      if (status == StatusCode::kOk)
        status = active->Scale(1. / selected_count, d);
      if (status == StatusCode::kOk) status = active->AbsoluteSum(&mean, d);
      if (status != StatusCode::kOk) return status;
    }
  } catch (const std::bad_alloc&) {
    return Failure(StatusCode::kAllocationFailure,
                   "Thermal flux selection allocation", d);
  }
  double used = 0;
  status = history->Preview(mean, maximum, &used, d);
  if (status != StatusCode::kOk) return status;
  out->contribution_count = ne + external_count;
  out->active_contribution_count = selected_count;
  out->spatial_mean = mean;
  out->maximum_flux = maximum;
  out->time_mean_used = used;
  return StatusCode::kOk;
}
StatusCode MeasureHeatFlux(const ThermalAssembly* a, HeatFluxHistory* h,
                           const FieldControls* controls, HeatFluxResidual* out,
                           Diagnostic* d) {
  if (h) h->DiscardCandidate();
  if (!a || !h || !controls || !out || !a->IsValid() ||
      h->StepNumber() != a->DofMap()->StepNumber())
    return Failure(StatusCode::kInvalidArgument,
                   "Invalid thermal residual measurement", d);
  const auto* map = a->DofMap();
  double maximum = 0, increment = 0;
  int index = 0, node = 0, label = 0;
  auto status = StatusCode::kOk;
  if (a->FreeResidual()) {
    status = a->FreeResidual()->MaximumAbs(&maximum, &index, d);
    if (status != StatusCode::kOk) return status;
    const int global = map->GlobalDofForEquation(index);
    for (int row = 0; row < map->NodeCount(); ++row) {
      const int id = map->NodeId(row), base = map->GlobalDof(id, 11);
      if (global >= base && global < base + map->NodeDofCount(id)) {
        node = id;
        label = 11 + global - base;
        break;
      }
    }
  }
  status = a->ParticipatingIncrement()->MaximumAbs(&increment, &index, d);
  if (status != StatusCode::kOk) return status;
  HeatFluxResidual measured;
  status = MeasureThermalFluxOrigins(a->ElementFluxes(), a->ExternalFluxes(), h,
                                     controls, &measured, d);
  if (status != StatusCode::kOk) return status;
  const double phi = maximum / measured.time_mean_used;
  if (!std::isfinite(phi)) {
    h->DiscardCandidate();
    return Failure(StatusCode::kLibraryFailure,
                   "Thermal normalized residual overflow", d);
  }
  out->dof_count = map->FreeCount();
  out->maximum_node = node;
  out->maximum_label = label;
  out->maximum = maximum;
  out->increment_max = increment;
  out->correction_max = 0;
  out->phi = phi;
  out->contribution_count = measured.contribution_count;
  out->active_contribution_count = measured.active_contribution_count;
  out->spatial_mean = measured.spatial_mean;
  out->maximum_flux = measured.maximum_flux;
  out->time_mean_used = measured.time_mean_used;
  out->branch = map->FreeCount() == 0 ? ThermalResidualBranch::kNoDofs
                                      : ThermalResidualBranch::kGeneral;
  out->passed = false;
  return StatusCode::kOk;
}
}  // namespace fesa
