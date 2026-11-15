#include <algorithm>
#include <chrono>
#include <cstring>
#include <cmath>
#include <cstdio>
#include <limits>
#include <memory>
#include <vector>

#include "adapters/abaqus/thermal_model_builder.h"
#include "fesa/analysis/heat_transfer_analysis.h"
#include "fesa/solver/fesa_solver.h"
#include "model/model_storage.h"
#include "solver/solver_access.h"

namespace fesa {
namespace {
bool SameThermalName(const char* a, const char* b) {
  while (*a || *b) {
    const int x = *a >= 'a' && *a <= 'z' ? *a - 'a' + 'A' : *a;
    const int y = *b >= 'a' && *b <= 'z' ? *b - 'a' + 'A' : *b;
    if (x != y) return false;
    if (!*a) return true;
    ++a;
    ++b;
  }
  return true;
}
double ThermalWallTime() {
  return std::chrono::duration<double>(
             std::chrono::steady_clock::now().time_since_epoch())
      .count();
}
const char* EventName(IncrementEvent event) {
  switch (event) {
    case IncrementEvent::kAttemptStarted:
      return "attempt";
    case IncrementEvent::kIteration:
      return "iteration";
    case IncrementEvent::kAttemptFailed:
      return "failed";
    case IncrementEvent::kCutback:
      return "cutback";
    case IncrementEvent::kCommitted:
      return "committed";
  }
  return "unknown";
}
}  // namespace
StatusCode FESASolver::PrepareThermal(Diagnostic* d) {
  for (const auto& item : *domain_->Steps())
    if (!item.second->IsHeatTransfer())
      return ReportFailure(StatusCode::kUnsupported,
                           DiagnosticStage::kValidation,
                           "Mixed mechanical/thermal steps", d);
  for (const auto* block : *deck_->Blocks()) {
    if (!block->compatibility_warning) continue;
    char text[256] = {};
    if (*block->keyword == "SURFACE") {
      bool referenced = false;
      for (const auto& item : *domain_->ThermalConditions()) {
        const auto* spec = item.second->specification();
        if (spec->target == ThermalConditionTarget::kSurface &&
            (std::strncmp(spec->scope, "PART.", 5) == 0) ==
                (block->scope_kind == InputScope::kPart) &&
            SameThermalName(spec->name, block->Option("NAME"))) {
          referenced = true;
          break;
        }
      }
      std::snprintf(text, sizeof(text), "Input line=%d: *SURFACE ",
                    block->line);
      for (const char* part :
           {static_cast<const char*>(text), block->Option("NAME"),
            referenced ? " validated; referenced by distributed load\n"
                       : " validated; inactive in M1\n"}) {
        const auto status = log_->Write(part, d);
        if (status != StatusCode::kOk) return status;
      }
      continue;
    }
    if (*block->keyword == "PHYSICALCONSTANTS") {
      bool radiation = false;
      for (const auto& row : *domain_->ThermalConditions())
        radiation = radiation || row.second->specification()->kind ==
                                     ThermalConditionKind::kRadiation;
      std::snprintf(text, sizeof(text),
                    radiation ? "Input line=%d: *PHYSICAL CONSTANTS stored; "
                                "used by radiation validation/evaluation\n"
                              : "Input line=%d: *PHYSICAL CONSTANTS stored; "
                                "unused by M1 solution\n",
                    block->line);
    } else
      std::snprintf(
          text, sizeof(text),
          "Input line=%d step=%d: request replaced by FESA fixed output.\n",
          block->line, block->step);
    const auto status = log_->Write(text, d);
    if (status != StatusCode::kOk) return status;
  }
  auto status = log_->Write(
      "Heat transfer fixed output: NT, RFL; HFL/GRADT/TEMP/CFL not written.\n",
      d);
  if (status != StatusCode::kOk) return status;
  std::vector<int> participating(0);
  status = CollectThermalNodeIds(domain_, &participating, d);
  if (status != StatusCode::kOk) return status;
  for (const auto* block : *deck_->Blocks()) {
    const bool initial = *block->keyword == "INITIALCONDITIONS";
    if (!initial && *block->keyword != "BOUNDARY" && *block->keyword != "CFLUX")
      continue;
    for (std::size_t row_index = 0; row_index < block->rows->size();) {
      const auto* row = block->rows->at(row_index);
      int total = 0, included = 0;
      if (initial) {
        std::vector<int> ids(0);
        std::vector<double> values(0);
        std::size_t next = row_index;
        int source_line = 0;
        status = ReadThermalInitialRow(deck_, domain_, block, row_index,
                                       &participating, &ids, &values, &next,
                                       &source_line, d);
        if (status != StatusCode::kOk) return status;
        row_index = next;
        total = static_cast<int>(ids.size());
        for (int id : ids)
          if (std::binary_search(participating.begin(), participating.end(),
                                 id))
            ++included;
      } else {
        status = CountThermalTargetParticipation(
            deck_, domain_, block, row, &participating, &total, &included, d);
        if (status != StatusCode::kOk) return status;
        ++row_index;
      }
      const int excluded = total - included;
      // Empty BC/CFLUX targets are errors in the shared condition preflight,
      // never warnings. Only an IC target can be an accepted empty no-op.
      if (excluded == 0 || (!initial && included == 0)) continue;
      char text[240] = {};
      std::snprintf(
          text, sizeof(text),
          "Input line=%d: thermal target total=%d included=%d "
          "excluded_nonparticipating=%d%s\n",
          row->line, total, included, excluded,
          initial && included == 0 ? "; empty initial temperature target" : "");
      status = log_->Write(text, d);
      if (status != StatusCode::kOk) return status;
    }
  }
  ControlSettings* raw_controls = nullptr;
  status = CloneControlSettings(controls_, &raw_controls, d);
  std::unique_ptr<ControlSettings> validation(raw_controls);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<ThermalConditions> previous(nullptr);
  double total_time = 0;
  std::int64_t expected = 1;
  NotifySolverStage(SolverStage::kBeforePreflight);
  for (const auto& entry : *domain_->Steps()) {
    if (entry.first != expected || expected > std::numeric_limits<int>::max())
      return model_storage::Invalid("Sequential thermal step numbering", d);
    ++expected;
    d->source_file = domain_->SourceFile();
    d->path_kind = DiagnosticPathKind::kLiteralFile;
    d->step = entry.first;
    total_time += entry.second->Period();
    if (!std::isfinite(total_time))
      return ReportFailure(StatusCode::kInputError,
                           DiagnosticStage::kValidation,
                           "Cumulative thermal step time overflow", d);
    status = ApplyStepControls(deck_, entry.first, validation.get(), d);
    if (status != StatusCode::kOk) return status;
    ThermalConditions* raw_conditions = nullptr;
    status = BuildThermalConditions(domain_, entry.first, previous.get(),
                                    &raw_conditions, d);
    std::unique_ptr<ThermalConditions> current(raw_conditions);
    if (status != StatusCode::kOk) return status;
    ThermalModel* raw_model = nullptr;
    status = ThermalModel::Create(&raw_model, d);
    std::unique_ptr<ThermalModel> model(raw_model);
    if (status != StatusCode::kOk) return status;
    status = model->BindStep(domain_, current.get(), d);
    if (status == StatusCode::kOk)
      status = HeatTransferAnalysis::ValidateStep(
          model.get(), model->DofMap(), current.get(),
          entry.second->HeatTransfer(), d);
    if (status != StatusCode::kOk) return status;
    model.reset();
    previous.swap(current);
  }
  status = ValidateResultModel(domain_, d);
  if (status != StatusCode::kOk) return status;
  return log_->Write(
      "Thermal input, controls, conditions, mapping and metadata validated.\n",
      d);
}
StatusCode FESASolver::ExecuteThermal(Diagnostic* d) {
  auto status = ThermalState::CreateForDomain(domain_, &thermal_state_, d);
  if (status == StatusCode::kOk)
    status = ThermalModel::Create(&thermal_model_, d);
  if (status != StatusCode::kOk) return status;
  analysis_report_ = new AnalysisReport;
  for (const auto& entry : *domain_->Steps()) {
    analysis_report_->Reset();
    logged_assembly_ = logged_solve_ = logged_recovery_ = logged_dispatch_ = 0;
    delete thermal_analysis_;
    thermal_analysis_ = nullptr;
    thermal_model_->Clear();
    delete trial_thermal_conditions_;
    trial_thermal_conditions_ = nullptr;
    d->source_file = domain_->SourceFile();
    d->source_line = d->first_source_line = d->node_id = d->element_id = 0;
    d->path_kind = DiagnosticPathKind::kLiteralFile;
    d->step = entry.first;
    status = ApplyStepControls(deck_, entry.first, controls_, d);
    if (status == StatusCode::kOk)
      status = log_->WriteControls(controls_, entry.first, d);
    if (status == StatusCode::kOk)
      status = BuildThermalConditions(domain_, entry.first, thermal_conditions_,
                                      &trial_thermal_conditions_, d);
    if (status == StatusCode::kOk)
      status = thermal_model_->BindStep(domain_, trial_thermal_conditions_, d);
    if (status == StatusCode::kOk)
      status = thermal_state_->ValidateMapping(thermal_model_->DofMap(), d);
    if (status == StatusCode::kOk)
      status = HeatTransferAnalysis::Create(entry.second->HeatTransfer(),
                                            thermal_conditions_,
                                            &thermal_analysis_, d);
    if (status != StatusCode::kOk) return status;
    char text[240] = {};
    const auto* parameters = entry.second->HeatTransfer();
    std::snprintf(text, sizeof(text),
                  "Thermal step=%d procedure=%s initial_dt=%.17g period=%.17g "
                  "DELTMX=%s\n",
                  entry.first, parameters->steady ? "steady" : "transient",
                  parameters->initial_increment, parameters->period,
                  parameters->has_deltmx ? "active" : "absent");
    status = log_->Write(text, d);
    if (status != StatusCode::kOk) return status;
    calculation_started_ = ThermalWallTime();
    status = thermal_analysis_->SolveStep(
        thermal_model_, thermal_model_->DofMap(), trial_thermal_conditions_,
        controls_, thermal_state_, writer_, next_frame_number_,
        completed_step_total_time_, analysis_report_, this, d);
    if (status != StatusCode::kOk) return status;
    if (!analysis_report_->step_completed)
      return model_storage::Invalid("Thermal analysis step did not complete",
                                    d);
  }
  return FinishOutput(d);
}
StatusCode FESASolver::LogThermalCalculation(const AnalysisReport* report,
                                             Diagnostic* d) {
  if (!report->residual_valid || !report->heat_flux_field)
    return model_storage::Invalid("Missing thermal residual report", d);
  NotifySolverStage(SolverStage::kSolved);
  const auto* heat = report->heat_flux_field;
  char text[800] = {};
  std::snprintf(
      text, sizeof(text),
      "heat flux residual: targets=%d contributions=%d branch=%d max=%.17g "
      "node=%d "
      "label=%d spatial=%.17g used=%.17g phi=%.17g increment=%.17g "
      "correction=%.17g passed=%d\n"
      "linear solver calls=%d; reduced equations=%d\n"
      "Calculation passed; elapsed_seconds=%.9g; output pending\n",
      heat->dof_count, heat->contribution_count, static_cast<int>(heat->branch),
      heat->maximum, heat->maximum_node, heat->maximum_label,
      heat->spatial_mean, heat->time_mean_used, heat->phi, heat->increment_max,
      heat->correction_max, heat->passed ? 1 : 0, report->linear_solver_calls,
      report->equation_count, ThermalWallTime() - calculation_started_);
  const auto status = log_->Write(text, d);
  if (status == StatusCode::kOk) output_started_ = ThermalWallTime();
  return status;
}
StatusCode FESASolver::OnThermalIncrementProgress(
    IncrementEvent event, const IncrementAttemptInfo* attempt,
    const ThermalIterationInfo* iteration, const HeatFluxResidual* heat,
    const Diagnostic* primary, Diagnostic* d) {
  if (!heat_transfer_ || !log_ || !thermal_state_ || !analysis_report_ ||
      !attempt || attempt->step_number != analysis_report_->step_number ||
      (primary && primary == d))
    return model_storage::Invalid("FESASolver thermal progress context", d);
  char text[900] = {};
  std::snprintf(text, sizeof(text),
                "thermal %s step=%d increment=%d attempt=%d dt=%.17g "
                "step_time=%.17g total_time=%.17g "
                "proposed_available=%d proposed_dt=%.17g status=%d\n",
                EventName(event), attempt->step_number,
                attempt->increment_number, attempt->attempt_number,
                attempt->increment_size, attempt->step_time,
                attempt->total_time, attempt->proposed_available ? 1 : 0,
                attempt->proposed_increment,
                primary ? static_cast<int>(primary->code) : 0);
  auto status = log_->Write(text, d);
  if (status != StatusCode::kOk) return status;
  if (attempt->thermal_change_available)
    std::snprintf(text, sizeof(text),
                  "thermal accuracy S_available=1 S=%.17g\n",
                  attempt->thermal_free_change);
  else
    std::snprintf(text, sizeof(text), "thermal accuracy S_available=0\n");
  status = log_->Write(text, d);
  if (status != StatusCode::kOk || !iteration) return status;
  std::snprintf(
      text, sizeof(text),
      "thermal iteration i=%d outcome=%d numerical_status=%d residual_valid=%d "
      "phi=%.17g branch=%d correction=%.17g increment=%.17g "
      "estimate_available=%d estimate=%.17g\n",
      iteration->iteration, static_cast<int>(iteration->outcome),
      static_cast<int>(iteration->numerical_status),
      iteration->residual_valid ? 1 : 0, iteration->phi,
      static_cast<int>(iteration->branch), iteration->correction_max,
      iteration->increment_max, iteration->estimate_available ? 1 : 0,
      iteration->estimate_max);
  status = log_->Write(text, d);
  if (status != StatusCode::kOk || !heat) return status;
  std::snprintf(text, sizeof(text),
                "heat flux spatial=%.17g time_mean=%.17g maximum=%.17g node=%d "
                "label=%d contributions=%d\n",
                heat->spatial_mean, heat->time_mean_used, heat->maximum,
                heat->maximum_node, heat->maximum_label,
                heat->contribution_count);
  return log_->Write(text, d);
}
}  // namespace fesa
