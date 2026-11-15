#include "fesa/io/analysis_log.h"

#include <cmath>
#include <cstdio>
#include <new>
#include <string>

namespace fesa {
namespace {
const char* EventName(IncrementEvent event) {
  switch (event) {
    case IncrementEvent::kAttemptStarted:
      return "attempt-start";
    case IncrementEvent::kIteration:
      return "iteration";
    case IncrementEvent::kAttemptFailed:
      return "attempt-failed";
    case IncrementEvent::kCutback:
      return "cutback";
    case IncrementEvent::kCommitted:
      return "committed";
  }
  return nullptr;
}
const char* InspectionName(GeometryInspection result) {
  switch (result) {
    case GeometryInspection::kNotRequested:
      return "not-requested";
    case GeometryInspection::kUnavailableCandidate:
      return "candidate-unavailable";
    case GeometryInspection::kUnsupported:
      return "unsupported";
    case GeometryInspection::kPassed:
      return "passed";
    case GeometryInspection::kGeometry:
      return "geometry";
    case GeometryInspection::kFailed:
      return "failed";
  }
  return nullptr;
}
void OptionalNumber(bool available, double value, char* out) {
  if (available)
    std::snprintf(out, 40, "%.17g", value);
  else
    std::snprintf(out, 40, "NA");
}
StatusCode WriteField(const char* name, const NewtonIterationInfo* info,
                      const FieldResidual* field, int index, AnalysisLog* log,
                      Diagnostic* d) {
  char text[2048] = {};
  if (!info || !info->residual_valid) {
    std::snprintf(text, sizeof(text), "%s residual=NA\n", name);
  } else {
    char correction[40] = {}, estimate[40] = {}, prediction[40] = {};
    OptionalNumber(info->correction_valid, field->correction_max, correction);
    OptionalNumber(info->estimate_available[index], info->estimate_max[index],
                   estimate);
    OptionalNumber(info->prediction_available[index],
                   info->predicted_total[index], prediction);
    std::snprintf(
        text, sizeof(text),
        "%s targets=%d contributions=%d branch=%d maximum=%.17g "
        "node=%d dof=%d element=%d internal_end=%d "
        "l2=%.17g rms=%.17g spatial=%.17g raw=%.17g used=%.17g "
        "maximum_flux=%.17g converted=%d phi=%.17g increment=%.17g "
        "correction=%s estimate=%s prediction=%s residual_limit=%.17g "
        "correction_limit=%.17g passed=%d\n",
        name, field->dof_count, field->contribution_count,
        static_cast<int>(field->branch), field->maximum, field->maximum_node,
        field->maximum_dof, field->maximum_element, field->maximum_internal_end,
        field->l2, field->rms, field->spatial_mean,
        field->raw_mean, field->used_mean, field->maximum_flux,
        field->converted ? 1 : 0, field->phi, field->increment_max, correction,
        estimate, prediction, field->residual_limit, field->correction_limit,
        field->passed ? 1 : 0);
  }
  return log->Write(text, d);
}
}  // namespace

StatusCode AnalysisLog::WriteIncrementProgress(
    IncrementEvent event, const IncrementAttemptInfo* attempt,
    const NewtonIterationInfo* iteration, const FieldResidual* displacement,
    const FieldResidual* rotation, const Diagnostic* primary, Diagnostic* d) {
  const char* event_name = EventName(event);
  if (!event_name || !attempt || (primary && primary == d) ||
      !InspectionName(attempt->geometry_inspection) ||
      (attempt->geometry_diagnostic && attempt->geometry_diagnostic == d) ||
      attempt->step_number <= 0 || attempt->increment_number <= 0 ||
      attempt->attempt_number <= 0 || !std::isfinite(attempt->step_time) ||
      attempt->step_time < 0 || !std::isfinite(attempt->total_time) ||
      attempt->total_time < 0 || !std::isfinite(attempt->increment_size) ||
      attempt->increment_size <= 0 ||
      !std::isfinite(attempt->nominal_increment) ||
      attempt->nominal_increment <= 0 ||
      !std::isfinite(attempt->characteristic_length) ||
      attempt->characteristic_length < 0 ||
      (attempt->predictor_used && (!std::isfinite(attempt->predictor_ratio) ||
                                   attempt->predictor_ratio <= 0)) ||
      (attempt->proposed_available &&
       (!std::isfinite(attempt->proposed_increment) ||
        attempt->proposed_increment <= 0)) ||
      (event == IncrementEvent::kIteration && !iteration) ||
      (iteration && iteration->residual_valid && (!displacement || !rotation)))
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation,
                         "WriteIncrementProgress input", d);
  char text[1024] = {}, ratio[40] = {}, proposed[40] = {};
  OptionalNumber(attempt->predictor_used, attempt->predictor_ratio, ratio);
  OptionalNumber(attempt->proposed_available, attempt->proposed_increment,
                 proposed);
  std::snprintf(text, sizeof(text),
                "Increment event=%s step=%d increment=%d attempt=%d "
                "step_time=%.17g total_time=%.17g size=%.17g nominal=%.17g "
                "complete=%d Lc=%.17g predictor=%d ratio=%s proposed=%s\n",
                event_name, attempt->step_number, attempt->increment_number,
                attempt->attempt_number, attempt->step_time,
                attempt->total_time, attempt->increment_size,
                attempt->nominal_increment, attempt->completes_step ? 1 : 0,
                attempt->characteristic_length, attempt->predictor_used ? 1 : 0,
                ratio, proposed);
  auto status = Write(text, d);
  if (status != StatusCode::kOk) return status;
  if (iteration) {
    std::snprintf(
        text, sizeof(text),
        "Newton iteration=%d attempted_iteration=%d status=%d "
        "numeric_status=%d "
        "outcome=%d stage=%d solves=%d residual_evaluations=%lld "
        "residual_valid=%d correction_valid=%d logging_failed=%d "
        "residual_seconds=%.9g tangent_seconds=%.9g solve_seconds=%.9g "
        "update_seconds=%.9g criteria_seconds=%.9g\n",
        iteration->iteration, iteration->attempted_iteration,
        static_cast<int>(iteration->status),
        static_cast<int>(iteration->numerical_status),
        static_cast<int>(iteration->outcome),
        static_cast<int>(iteration->stage), iteration->linear_solver_calls,
        static_cast<long long>(iteration->residual_evaluations),
        iteration->residual_valid ? 1 : 0, iteration->correction_valid ? 1 : 0,
        iteration->logging_failed ? 1 : 0, iteration->residual_seconds,
        iteration->tangent_seconds, iteration->solve_seconds,
        iteration->update_seconds, iteration->criteria_seconds);
  } else {
    std::snprintf(text, sizeof(text), "Newton iteration=NA\n");
  }
  status = Write(text, d);
  if (status != StatusCode::kOk) return status;
  status = WriteField("displacement", iteration, displacement, 0, this, d);
  if (status != StatusCode::kOk) return status;
  status = WriteField("rotation", iteration, rotation, 1, this, d);
  if (status != StatusCode::kOk) return status;
  try {
    std::string diagnostic;
    if (primary && primary->code != StatusCode::kOk) {
      status = FormatDiagnostic(primary, &diagnostic);
      if (status != StatusCode::kOk)
        return ReportFailure(status, DiagnosticStage::kAllocation,
                             "Format increment diagnostic", d);
      status = Write(diagnostic.c_str(), d);
      if (status != StatusCode::kOk) return status;
    }
    if (attempt->geometry_inspection == GeometryInspection::kNotRequested)
      return StatusCode::kOk;
    std::snprintf(text, sizeof(text),
                  "geometry_inspection=%s original_outcome=%d\n",
                  InspectionName(attempt->geometry_inspection),
                  static_cast<int>(attempt->original_outcome));
    status = Write(text, d);
    if (status != StatusCode::kOk || !attempt->geometry_diagnostic)
      return status;
    std::snprintf(
        text, sizeof(text), "geometry_measurement_available=%d\n",
        attempt->geometry_diagnostic->geometry_value_available ? 1 : 0);
    status = Write(text, d);
    if (status != StatusCode::kOk) return status;
    status = FormatDiagnostic(attempt->geometry_diagnostic, &diagnostic);
    if (status != StatusCode::kOk)
      return ReportFailure(status, DiagnosticStage::kAllocation,
                           "Format geometry inspection diagnostic", d);
    return Write(diagnostic.c_str(), d);
  } catch (const std::bad_alloc&) {
    return ReportFailure(StatusCode::kAllocationFailure,
                         DiagnosticStage::kAllocation,
                         "Format increment diagnostic allocation", d);
  }
}
}  // namespace fesa
