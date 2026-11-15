#include "fesa/solver/fesa_solver.h"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <memory>
#include <new>
#include <string>
#include <utility>

#include "adapters/abaqus/dynamic_settings_reader.h"
#include "adapters/abaqus/initial_velocity_reader.h"
#include "adapters/abaqus/input_syntax.h"
#include "adapters/abaqus/model_builder.h"
#include "adapters/openmp/thread_runtime.h"
#include "analysis/equilibrium_check.h"
#include "analysis/nonlinear_assembler.h"
#include "analysis/nonlinear_condition_validator.h"
#include "analysis/static_assembler.h"
#include "fesa/analysis/heat_transfer_analysis.h"
#include "fesa/analysis/linear_dynamic_analysis.h"
#include "fesa/analysis/nonlinear_dynamic_analysis.h"
#include "fesa/analysis/nonlinear_static_analysis.h"
#include "fesa/analysis/static_analysis.h"
#include "fesa/elements/beam3d.h"
#include "io/file_access.h"
#include "model/model_storage.h"
#include "solver/solver_access.h"
namespace fesa {
namespace {
thread_local void (*observer)(SolverStage) = nullptr;
}
void SetSolverObserver(void (*callback)(SolverStage)) { observer = callback; }
void NotifySolverStage(SolverStage stage) {
  if (observer) observer(stage);
}
namespace {
bool FiniteBeam(const Domain* domain) {
  if (!domain || domain->Elements()->empty() || domain->Steps()->empty() ||
      !domain->Steps()->begin()->second->HasFiniteGeometry()) return false;
  for (const auto& item : *domain->Elements())
    if (!dynamic_cast<const Beam3D*>(item.second)) return false;
  return true;
}
// Keep validation and execution on the same family-specific mass contract.
MassMode DynamicMassMode(MassSelection selected, bool nonlinear,
                         bool finite_beam) {
  if (selected == MassSelection::kConsistent) return MassMode::kConsistent;
  if (selected == MassSelection::kLumped) return MassMode::kRowSum;
  if (finite_beam) return MassMode::kConsistent;
  return nonlinear ? MassMode::kRowSum : MassMode::kDefault;
}
double WallSeconds() {
  return std::chrono::duration<double>(
             std::chrono::steady_clock::now().time_since_epoch())
      .count();
}
void Console(const Diagnostic* d) {
  try {
    std::string text("");
    if (FormatDiagnostic(d, &text) == StatusCode::kOk) {
      std::fputs(text.c_str(), stderr);
      return;
    }
  } catch (const std::bad_alloc&) {
  }
  std::fprintf(stderr, "FESA: %s (diagnostic formatting failed)\n",
               d && d->operation ? d->operation : "runtime failure");
}
bool SameDiagnosticText(const char* first, const char* second) {
  return first == second ||
         (first && second && std::strcmp(first, second) == 0);
}
// Exact copied-event metadata identity; floating metadata use bit identity,
// not numerical tolerance. Never allocates or modifies either diagnostic.
bool SameDiagnosticEvent(const Diagnostic* first, const Diagnostic* second) {
  if (!first || !second) return first == second;
  return first->code == second->code && first->stage == second->stage &&
         SameDiagnosticText(first->operation, second->operation) &&
         SameDiagnosticText(first->source_file, second->source_file) &&
         first->path_kind == second->path_kind &&
         first->source_line == second->source_line &&
         first->first_source_line == second->first_source_line &&
         first->step == second->step && first->node_id == second->node_id &&
         first->element_id == second->element_id &&
         first->condition_dof == second->condition_dof &&
         first->duplicate_boundary == second->duplicate_boundary &&
         first->boundary_dof == second->boundary_dof &&
         std::memcmp(&first->first_boundary_value,
                     &second->first_boundary_value, sizeof(double)) == 0 &&
         std::memcmp(&first->boundary_value, &second->boundary_value,
                     sizeof(double)) == 0 &&
         first->active_flux_failure == second->active_flux_failure &&
         first->flux_field == second->flux_field &&
         first->flux_count == second->flux_count &&
         first->active_flux_count == second->active_flux_count &&
         std::memcmp(&first->flux_epsilon_l, &second->flux_epsilon_l,
                     sizeof(double)) == 0 &&
         std::memcmp(&first->flux_history_max, &second->flux_history_max,
                     sizeof(double)) == 0 &&
         std::memcmp(&first->flux_threshold, &second->flux_threshold,
                     sizeof(double)) == 0 &&
         std::memcmp(&first->flux_current_max, &second->flux_current_max,
                     sizeof(double)) == 0 &&
         first->geometry_point_available == second->geometry_point_available &&
         std::memcmp(&first->natural_xi, &second->natural_xi, sizeof(double)) ==
             0 &&
         std::memcmp(&first->natural_eta, &second->natural_eta,
                     sizeof(double)) == 0 &&
         std::memcmp(&first->natural_zeta, &second->natural_zeta,
                     sizeof(double)) == 0 &&
         first->geometry_value_available == second->geometry_value_available &&
         first->geometry_quantity == second->geometry_quantity &&
         std::memcmp(&first->geometry_value, &second->geometry_value,
                     sizeof(double)) == 0 &&
         std::memcmp(&first->geometry_limit, &second->geometry_limit,
                     sizeof(double)) == 0 &&
         first->vendor_code == second->vendor_code;
}
// Derives the cfg path for valid decks and conservative error-log safety.
StatusCode SelectedConfigPath(const RunOptions* request, std::string* path,
                              Diagnostic* d) {
  d->source_file = request->config_path ? request->config_path->c_str()
                   : request->executable_directory
                       ? request->executable_directory->c_str()
                       : nullptr;
  d->path_kind = request->config_path
                     ? DiagnosticPathKind::kLiteralFile
                     : DiagnosticPathKind::kDefaultConfigDirectory;
  if (request->config_path)
    path->assign(*request->config_path);
  else if (request->executable_directory &&
           !request->executable_directory->empty())
    io::DefaultConfigPath(request->executable_directory->c_str(), path);
  else
    return model_storage::Invalid("Missing executable cfg context", d);
  return StatusCode::kOk;
}
StatusCode WriteThreadRequest(const RunOptions* options, int available_cpus,
                              AnalysisLog* log, Diagnostic* d) {
  char line[120] = {};
  std::snprintf(line, sizeof(line),
                "threads requested=%d; analysis runtime request\n",
                options->threads);
  auto status = log->Write(line, d);
  if (status != StatusCode::kOk) return status;
  std::snprintf(line, sizeof(line),
                "thread validation upper_bound=%d (caller validation context; "
                "0=unspecified)\n",
                available_cpus);
  status = log->Write(line, d);
  if (status != StatusCode::kOk) return status;
  for (const char* key : {"OMP_NUM_THREADS", "MKL_NUM_THREADS"}) {
    std::size_t size = 0;
    char* value = nullptr;
    const errno_t result = _dupenv_s(&value, &size, key);
    std::unique_ptr<char, decltype(&std::free)> owner(value, &std::free);
    if (result != 0)
      return ReportFailure(StatusCode::kAllocationFailure,
                           DiagnosticStage::kAllocation,
                           "Read thread environment", d);
    for (const char* text : {key, "=", value ? value : "<unset>", "\n"}) {
      status = log->Write(text, d);
      if (status != StatusCode::kOk) return status;
    }
  }
  return StatusCode::kOk;
}
StatusCode LogResidual(const FieldResidual* field, const char* name,
                       AnalysisLog* log, Diagnostic* d) {
  char line[1024] = {};
  std::snprintf(
      line, sizeof(line),
      "%s residual: targets=%d contributions=%d branch=%d max=%.17g "
      "node=%d dof=%d l2=%.17g rms=%.17g spatial=%.17g raw=%.17g "
      "used=%.17g converted=%d phi=%.17g increment=%.17g correction=%.17g "
      "residual_limit=%.17g correction_limit=%.17g passed=%d\n",
      name, field->dof_count, field->contribution_count,
      static_cast<int>(field->branch), field->maximum, field->maximum_node,
      field->maximum_dof, field->l2, field->rms, field->spatial_mean,
      field->raw_mean, field->used_mean, field->converted ? 1 : 0, field->phi,
      field->increment_max, field->correction_max, field->residual_limit,
      field->correction_limit, field->passed ? 1 : 0);
  return log->Write(line, d);
}
}  // namespace
StatusCode FESASolver::Create(FESASolver** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("FESASolver output", d);
  *out = new (std::nothrow) FESASolver;
  return *out ? StatusCode::kOk
              : model_storage::AllocationFailure("FESASolver allocation", d);
}
FESASolver::~FESASolver() { Reset(); }
void FESASolver::Reset() {
  heat_transfer_ = false;
  writer_finalized_ = false;
  writer_succeeded_ = false;
  log_finalized_ = false;
  last_completed_step_ = 0;
  next_frame_number_ = 1;
  completed_step_total_time_ = 0;
  calculation_started_ = output_started_ = output_seconds_ = 0;
  logged_assembly_ = logged_solve_ = logged_recovery_ = logged_dispatch_ = 0;
  delete writer_;
  writer_ = nullptr;
  if (analysis_report_) analysis_report_->Reset();
  delete analysis_report_;
  analysis_report_ = nullptr;
  delete analysis_;
  analysis_ = nullptr;
  delete dynamic_settings_;
  dynamic_settings_ = nullptr;
  delete initial_velocity_;
  initial_velocity_ = nullptr;
  delete model_;
  model_ = nullptr;
  delete dofs_;
  dofs_ = nullptr;
  delete conditions_;
  conditions_ = nullptr;
  delete trial_conditions_;
  trial_conditions_ = nullptr;
  delete state_;
  state_ = nullptr;
  delete thermal_analysis_;
  thermal_analysis_ = nullptr;
  delete thermal_model_;
  thermal_model_ = nullptr;
  delete trial_thermal_conditions_;
  trial_thermal_conditions_ = nullptr;
  delete thermal_conditions_;
  thermal_conditions_ = nullptr;
  delete thermal_state_;
  thermal_state_ = nullptr;
  delete domain_;
  domain_ = nullptr;
  delete controls_;
  controls_ = nullptr;
  delete deck_;
  deck_ = nullptr;
  delete log_;
  log_ = nullptr;
  delete output_;
  output_ = nullptr;
  delete input_;
  input_ = nullptr;
  delete request_;
  request_ = nullptr;
}
StatusCode FESASolver::Run(const RunOptions* options, Diagnostic* diagnostic) {
  Diagnostic local;
  Diagnostic* d = diagnostic ? diagnostic : &local;
  // A new Run ends previous diagnostic borrowing before destroying its owners.
  d->source_file = nullptr;
  d->source_line = d->first_source_line = d->node_id = d->element_id = 0;
  d->step = 0;
  d->duplicate_boundary = false;
  d->condition_dof = 0;
  d->path_kind = DiagnosticPathKind::kLiteralFile;
  Reset();
  const auto start = std::chrono::steady_clock::now();
  StatusCode status = StatusCode::kOk;
  try {
    status = Prepare(options, d);
    if (status == StatusCode::kOk) {
      char text[100] = {};
      const double seconds = std::chrono::duration<double>(
          std::chrono::steady_clock::now() - start).count();
      std::snprintf(text, sizeof(text),
                    "Preparation elapsed_seconds=%.9g\n", seconds);
      status = log_->Write(text, d);
    }
    if (status == StatusCode::kOk) status = Execute(d);
    if (status == StatusCode::kOk) {
      const double elapsed = std::chrono::duration<double>(
                                 std::chrono::steady_clock::now() - start)
                                 .count();
      char text[160] = {};
      std::snprintf(text, sizeof(text),
                    "Run elapsed_seconds=%.9g (includes HDF5 flush/close); "
                    "last_completed_step=%d\n",
                    elapsed, last_completed_step_);
      status = log_->Write(text, d);
    }
    if (status == StatusCode::kOk) {
      log_finalized_ = true;
      status = log_->Finish(d);
    }
  } catch (const std::bad_alloc&) {
    status = model_storage::AllocationFailure("FESASolver allocation", d);
  }
  return status == StatusCode::kOk ? status : Fail(status, d);
}
StatusCode FESASolver::Prepare(const RunOptions* options, Diagnostic* d) {
  if (!options || !options->input_path || options->input_path->empty() ||
      options->help || options->threads < 1 ||
      (options->available_cpus > 0 &&
       options->threads > options->available_cpus))
    return model_storage::Invalid("FESASolver RunOptions", d);
  request_ = new RunOptions;
  request_->input_path = new std::string(*options->input_path);
  request_->threads = options->threads;
  request_->b33_material_first = options->b33_material_first;
  request_->available_cpus = options->available_cpus;
  d->source_file = request_->input_path->c_str();
  if (options->config_path)
    request_->config_path = new std::string(*options->config_path);
  if (options->executable_directory)
    request_->executable_directory =
        new std::string(*options->executable_directory);
  input_ = new std::string("");
  auto status = io::AbsolutePath(request_->input_path->c_str(), input_, d);
  if (status != StatusCode::kOk) return status;
  d->source_file = input_->c_str();
  bool exists = false, regular = false;
  status = io::InspectPath(input_->c_str(), &exists, &regular, d);
  if (status != StatusCode::kOk) return status;
  if (!exists || !regular)
    return ReportFailure(StatusCode::kIoFailure, DiagnosticStage::kValidation,
                         "Input path must be an existing regular file", d);
  output_ = new std::string("");
  io::ResultPath(input_->c_str(), output_);
  std::string log_path("");
  io::AnalysisLogPath(input_->c_str(), &log_path);
  status = io::ValidateDistinctPaths(input_->c_str(), output_->c_str(), d);
  if (status != StatusCode::kOk) return status;
  status = io::ValidateDistinctPaths(input_->c_str(), log_path.c_str(), d);
  if (status != StatusCode::kOk) return status;
  status = ReadInputDeck(input_->c_str(), &deck_, d);
  if (status != StatusCode::kOk) {
    // Preserve the input error. A malformed deck cannot establish cfg bypass,
    // so conservatively protect the declared source before error-log
    // replacement.
    Diagnostic secondary;
    try {
      std::string cfg("");
      auto safe = SelectedConfigPath(request_, &cfg, &secondary);
      if (safe == StatusCode::kOk)
        safe = io::ValidateDistinctPaths(cfg.c_str(), log_path.c_str(),
                                         &secondary);
      if (safe == StatusCode::kOk)
        safe = OpenAnalysisLog(input_->c_str(), &log_, &secondary);
      if (safe != StatusCode::kOk) Console(&secondary);
    } catch (const std::bad_alloc&) {
      model_storage::AllocationFailure("Input error log safety allocation",
                                       &secondary);
      Console(&secondary);
    }
    return status;
  }
  {
    std::string cfg("");
    status = SelectedConfigPath(request_, &cfg, d);
    if (status != StatusCode::kOk) return status;
    status = io::ValidateDistinctPaths(cfg.c_str(), output_->c_str(), d);
    if (status != StatusCode::kOk) return status;
    status = io::ValidateDistinctPaths(cfg.c_str(), log_path.c_str(), d);
    if (status != StatusCode::kOk) return status;
  }
  status = OpenAnalysisLog(input_->c_str(), &log_, d);
  if (status != StatusCode::kOk) return status;
  for (const char* text : {"input=", input_->c_str(),
                           "\nPreparation; no HDF5 result was opened.\n"}) {
    status = log_->Write(text, d);
    if (status != StatusCode::kOk) return status;
  }
  status = WriteThreadRequest(request_, request_->available_cpus, log_, d);
  if (status != StatusCode::kOk) return status;
  status = SelectInitialControls(deck_, request_,
                                 request_->executable_directory
                                     ? request_->executable_directory->c_str()
                                     : nullptr,
                                 &controls_, d);
  if (status != StatusCode::kOk) return status;
  {
    // Missing cfg may now exist; repeat native identity checking before output.
    d->source_file = controls_->config_path->c_str();
    d->path_kind = DiagnosticPathKind::kLiteralFile;
    status = io::ValidateDistinctPaths(controls_->config_path->c_str(),
                                       output_->c_str(), d);
    if (status != StatusCode::kOk) return status;
    status = io::ValidateDistinctPaths(controls_->config_path->c_str(),
                                       log_path.c_str(), d);
    if (status != StatusCode::kOk) return status;
  }
  if (controls_->config_created) {
    status =
        log_->Write("Created missing cfg with automatic q0/qu defaults.\n", d);
    if (status != StatusCode::kOk) return status;
  }
  for (const auto* block : *deck_->Blocks())
    if (*block->keyword == "HEATTRANSFER") heat_transfer_ = true;
  if (heat_transfer_) controls_->temperature = new FieldControls;
  status = log_->WriteControls(controls_, 0, d);
  if (status != StatusCode::kOk) return status;
  status = ConfigureThreadRuntime(request_->threads, d);
  if (status != StatusCode::kOk) return status;
  char threads[100] = {};
  std::snprintf(threads, sizeof(threads), "runtime threads applied=%d\n",
                request_->threads);
  status = log_->Write(threads, d);
  if (status != StatusCode::kOk) return status;
  int cnr_ceiling = 0;
  status = PardisoCnrThreadCeiling(&cnr_ceiling, d);
  if (status != StatusCode::kOk) return status;
  char cnr[180] = {};
  std::snprintf(cnr, sizeof(cnr),
                "PARDISO CNR ceiling=%d; source=current process affinity "
                "snapshot; PARDISO actual workers=unobserved\n",
                cnr_ceiling);
  status = log_->Write(cnr, d);
  if (status != StatusCode::kOk) return status;
  if (!heat_transfer_) {
    status = ResolveSolidConstitutiveModel(
        deck_, controls_->solid_constitutive_model,
        &controls_->solid_constitutive_model, d);
    if (status != StatusCode::kOk) return status;
  }
  status = BuildDomain(deck_, &domain_, d, controls_->solid_constitutive_model);
  if (status != StatusCode::kOk) return status;
  if (request_->b33_material_first) {
    bool supported = FiniteBeam(domain_);
    for (const auto& entry : *domain_->Steps())
      supported = supported && !entry.second->IsDynamic() &&
                  entry.second->HasFiniteGeometry() &&
                  entry.second->NonlinearParameters();
    if (!supported)
      return ReportFailure(StatusCode::kUnsupported,
                           DiagnosticStage::kValidation,
                           "--b33-material-first requires pure nonlinear "
                           "static B33 steps",
                           d);
    status = log_->Write(
        "B33 tangent policy=material-first; later corrections=consistent; "
        "exact residual and convergence unchanged\n",
        d);
    if (status != StatusCode::kOk) return status;
  }
  if (heat_transfer_ || !domain_->ThermalElements()->empty()) {
    heat_transfer_ = true;
    return PrepareThermal(d);
  }
  for (const auto* block : *deck_->Blocks())
    if (*block->keyword == "INITIALCONDITIONS" &&
        !domain_->Steps()->begin()->second->IsDynamic())
      return abaqus::Error(deck_->SourceFile(), block->line, block->step,
                           "Initial velocity requires direct dynamic start", d,
                           StatusCode::kUnsupported);
  for (const auto& entry : *domain_->Elements()) {
    const auto* element = entry.second;
    if (controls_->solid_constitutive_model ==
            SolidConstitutiveModel::kIncrementalElastic &&
        !element->SupportsTrialEvaluation()) {
      if (d) d->element_id = entry.first;
      return ReportFailure(StatusCode::kUnsupported,
          DiagnosticStage::kValidation,
          "Incremental Solid trial implementation is not ready", d);
    }
    const bool default_mass = controls_->mass == MassSelection::kDefault;
    if (default_mass && controls_->solid_constitutive_model ==
                            SolidConstitutiveModel::kIncrementalElastic)
      continue;
    // Legacy linear Beam has no mass API. Preserve its default execution;
    // explicit unsupported requests must fail before any result is replaced.
    if (default_mass && element->family() == ElementFamily::kBeam &&
        !FiniteBeam(domain_)) continue;
    const auto mode = controls_->mass == MassSelection::kConsistent
        ? MassMode::kConsistent
        : controls_->mass == MassSelection::kLumped
            ? MassMode::kRowSum
            : element->SupportsMass(MassMode::kRowSum)
                ? MassMode::kRowSum : MassMode::kConsistent;
    // Finite B33 supplies the full 14-variable dynamic evaluation, not the
    // legacy physical 12-by-12 mass API used by linear Beam consumers.
    const bool supported = FiniteBeam(domain_)
        ? mode == MassMode::kConsistent &&
              element->SupportsDynamicEvaluation(
                  DynamicKinematics::kFiniteRotation)
        : element->SupportsMass(mode);
    if (!supported) {
      d->source_file = controls_->config_path->c_str();
      d->source_line = controls_->mass_line;
      d->path_kind = DiagnosticPathKind::kLiteralFile;
      d->element_id = entry.first;
      return ReportFailure(StatusCode::kUnsupported,
                           DiagnosticStage::kValidation,
                           "Unsupported requested mass mode", d);
    }
    char mass_text[120] = {};
    std::snprintf(mass_text, sizeof(mass_text),
                  "mass element=%d selected=%s (used by dynamic solve only)\n",
                  entry.first,
                  mode == MassMode::kConsistent ? "consistent" : "lumped");
    status = log_->Write(mass_text, d);
    if (status != StatusCode::kOk) return status;
  }
  for (const auto* block : *deck_->Blocks()) {
    if (!block->compatibility_warning) continue;
    char text[180] = {};
    std::snprintf(
        text, sizeof(text),
        "Input line=%d step=%d: request replaced by FESA fixed output.\n",
        block->line, block->step);
    status = log_->Write(text, d);
    if (status != StatusCode::kOk) return status;
  }
  ControlSettings* validation_raw = nullptr;
  status = CloneControlSettings(controls_, &validation_raw, d);
  std::unique_ptr<ControlSettings> validation(validation_raw);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<EffectiveConditions> previous;
  const bool nonlinear = !domain_->Steps()->empty() &&
                         domain_->Steps()->begin()->second->HasFiniteGeometry();
  std::unique_ptr<NonlinearConditionValidator> condition_validator;
  if (nonlinear) {
    NonlinearConditionValidator* raw = nullptr;
    status = NonlinearConditionValidator::Create(domain_, &raw, d);
    condition_validator.reset(raw);
    if (status != StatusCode::kOk) return status;
  }
  double total_time = 0;
  std::int64_t expected = 1;
  NotifySolverStage(SolverStage::kBeforePreflight);
  for (const auto& entry : *domain_->Steps()) {
    const int number = entry.first;
    d->step = number;
    d->source_file = domain_->SourceFile();
    d->path_kind = DiagnosticPathKind::kLiteralFile;
    if (number != expected || expected > std::numeric_limits<int>::max())
      return model_storage::Invalid("Sequential step numbering", d);
    ++expected;
    if (entry.second->HasFiniteGeometry() != nonlinear)
      return ReportFailure(StatusCode::kUnsupported,
                           DiagnosticStage::kValidation,
                           "Linear/nonlinear step transition unsupported", d);
    total_time += entry.second->Period();
    if (!std::isfinite(total_time))
      return ReportFailure(StatusCode::kInputError,
                           DiagnosticStage::kValidation,
                           "Cumulative step time overflow", d);
    status = ApplyStepControls(deck_, number, validation.get(), d);
    if (status != StatusCode::kOk) return status;
    EffectiveConditions* raw = nullptr;
    status = BuildEffectiveConditions(domain_, number, previous.get(), &raw, d);
    std::unique_ptr<EffectiveConditions> current(raw);
    if (status != StatusCode::kOk) return status;
    if (condition_validator) {
      status = condition_validator->ValidateStep(current.get(), d);
      if (status != StatusCode::kOk) return status;
    }
    DofManager* raw_dofs = nullptr;
    status = FiniteBeam(domain_)
        ? DofManager::CreateForNonlinearBeam(domain_, current.get(),
                                             &raw_dofs, d)
        : DofManager::Create(domain_, current.get(), &raw_dofs, d);
    std::unique_ptr<DofManager> map(raw_dofs);
    if (status != StatusCode::kOk) return status;
    AnalysisModel* raw_model = nullptr;
    status = AnalysisModel::Create(&raw_model, d);
    std::unique_ptr<AnalysisModel> view(raw_model);
    if (status != StatusCode::kOk) return status;
    status = view->BindStep(domain_, current.get(), d);
    if (status != StatusCode::kOk) return status;
    if (entry.second->IsDynamic()) {
      const auto mode = DynamicMassMode(validation->mass, nonlinear,
                                        FiniteBeam(domain_));
      status = nonlinear ? NonlinearDynamicAnalysis::ValidateStep(
                               view.get(), map.get(), mode, previous.get(), d)
                         : LinearDynamicAnalysis::ValidateStep(
                               view.get(), mode, previous.get(), d);
      if (status != StatusCode::kOk) return status;
      if (number == 1) {
        Vector* raw_velocity = nullptr;
        status = abaqus::ResolveInitialVelocity(deck_, domain_, map.get(),
                                                &raw_velocity, d);
        std::unique_ptr<Vector> velocity(raw_velocity);
        if (status != StatusCode::kOk) return status;
      }
    } else if (previous &&
               domain_->Steps()->at(previous->StepNumber())->IsDynamic()) {
      return ReportFailure(StatusCode::kUnsupported,
                           DiagnosticStage::kValidation,
                           "Dynamic to static transition unsupported", d);
    }
    status = FiniteBeam(domain_)
        ? ValidateNonlinearBeamAssemblyCounts(view.get(), map.get(), d)
        : ValidateAssemblyCounts(view.get(), map.get(), d);
    if (status != StatusCode::kOk) return status;
    status = ValidateEquilibriumCounts(view.get(), current.get(), d);
    if (status != StatusCode::kOk) return status;
    view.reset();
    map.reset();
    previous.swap(current);
  }
  status = ValidateResultModel(domain_, d);
  if (status != StatusCode::kOk) return status;
  return log_->Write(
      "Input, controls, conditions, mapping and metadata validated.\n", d);
}
StatusCode FESASolver::LogCalculation(const AnalysisReport* report,
                                      Diagnostic* d) {
  char text[320] = {};
  std::snprintf(text, sizeof(text),
                "Stages assembly_inclusive=%.9g solve=%.9g "
                "equilibrium_recovery_inclusive=%.9g element_dispatch=%.9g "
                "observed_element_team_max=%d\n",
                report->assembly_seconds - logged_assembly_,
                report->solve_seconds - logged_solve_,
                report->equilibrium_recovery_seconds - logged_recovery_,
                report->element_dispatch_seconds - logged_dispatch_,
                report->observed_element_team_max);
  auto status = log_->Write(text, d);
  if (status != StatusCode::kOk) return status;
  logged_assembly_ = report->assembly_seconds;
  logged_solve_ = report->solve_seconds;
  logged_recovery_ = report->equilibrium_recovery_seconds;
  logged_dispatch_ = report->element_dispatch_seconds;
  if (report->calculation_status != StatusCode::kOk) return StatusCode::kOk;
  if (heat_transfer_) return LogThermalCalculation(report, d);
  if (dynamic_settings_ && report->increment_number == 0) {
    status = log_->Write(
        "Dynamic initialization and recovery passed; output pending.\n", d);
    if (status == StatusCode::kOk) output_started_ = WallSeconds();
    return status;
  }
  if (!report->residual_valid || !report->displacement_field ||
      !report->rotation_field)
    return model_storage::Invalid("Missing completed residual report", d);
  NotifySolverStage(SolverStage::kSolved);
  std::snprintf(text, sizeof(text),
                "linear solver calls=%d; reduced equations=%d\n",
                report->linear_solver_calls, report->equation_count);
  status = log_->Write(text, d);
  if (status != StatusCode::kOk) return status;
  for (int field = 0; field < 2; ++field) {
    status = LogResidual(field == 0 ? report->displacement_field
                                    : report->rotation_field,
                         field == 0 ? "displacement" : "rotation", log_, d);
    if (status != StatusCode::kOk) return status;
  }
  std::snprintf(text, sizeof(text),
                "Calculation passed; characteristic_length=%.17g "
                "elapsed_seconds=%.9g; output pending\n",
                report->characteristic_length,
                WallSeconds() - calculation_started_);
  status = log_->Write(text, d);
  if (status == StatusCode::kOk) output_started_ = WallSeconds();
  return status;
}
StatusCode FESASolver::OnProgress(AnalysisEvent event,
                                  const AnalysisReport* report, Diagnostic* d) {
  if (!report || report != analysis_report_ || !log_ ||
      (!state_ && !thermal_state_))
    return model_storage::Invalid("FESASolver progress context", d);
  if (event == AnalysisEvent::kCalculationFinished)
    return LogCalculation(report, d);
  if (event != AnalysisEvent::kIncrementCommitted)
    return model_storage::Invalid("FESASolver progress event", d);
  if (report->step_completed) {
    if (heat_transfer_) {
      delete thermal_conditions_;
      thermal_conditions_ = trial_thermal_conditions_;
      trial_thermal_conditions_ = nullptr;
    } else {
      std::swap(conditions_, trial_conditions_);
    }
    // The completed Analysis still borrows the previous conditions. Retain
    // them until that Analysis is destroyed at the next step or Reset.
    last_completed_step_ = report->step_number;
    completed_step_total_time_ = report->last_committed_total_time;
  }
  next_frame_number_ = report->last_committed_frame + 1;
  output_seconds_ += WallSeconds() - output_started_;
  NotifySolverStage(SolverStage::kCommitted);
  char text[160] = {};
  std::snprintf(text, sizeof(text),
                "Step %d increment %d completed and committed.\n",
                report->step_number, report->increment_number);
  return log_->Write(text, d);
}
StatusCode FESASolver::OnIncrementProgress(
    IncrementEvent event, const IncrementAttemptInfo* attempt,
    const NewtonIterationInfo* iteration, const FieldResidual* force,
    const FieldResidual* moment, const Diagnostic* primary, Diagnostic* d) {
  if (!log_ || !state_ || !analysis_report_ || !attempt ||
      attempt->step_number != analysis_report_->step_number)
    return model_storage::Invalid("FESASolver nonlinear progress context", d);
  return log_->WriteIncrementProgress(event, attempt, iteration, force, moment,
                                      primary, d);
}
StatusCode FESASolver::Execute(Diagnostic* d) {
  // All-step preflight context does not describe pre-step runtime work.
  d->source_file = domain_->SourceFile();
  d->path_kind = DiagnosticPathKind::kLiteralFile;
  d->source_line = d->first_source_line = d->node_id = d->element_id = 0;
  d->duplicate_boundary = false;
  d->condition_dof = 0;
  d->step = 0;
  auto status =
      OpenResultWriter(output_->c_str(), input_->c_str(), &writer_, d);
  if (status != StatusCode::kOk) return status;
  status = writer_->WriteModel(domain_, d);
  if (status != StatusCode::kOk) return status;
  NotifySolverStage(SolverStage::kModelWritten);
  if (heat_transfer_) return ExecuteThermal(d);
  if (domain_->Steps()->begin()->second->HasFiniteGeometry())
    status = AnalysisState::CreateNonlinearForDomain(domain_, &state_, d);
  else
    status = AnalysisState::CreateForDomain(domain_, &state_, d);
  if (status != StatusCode::kOk) return status;
  analysis_report_ = new AnalysisReport;
  for (const auto& entry : *domain_->Steps()) {
    const int number = entry.first;
    const auto* step = entry.second;
    analysis_report_->Reset();
    logged_assembly_ = logged_solve_ = logged_recovery_ = logged_dispatch_ = 0;
    delete analysis_;
    analysis_ = nullptr;
    delete dynamic_settings_;
    dynamic_settings_ = nullptr;
    delete initial_velocity_;
    initial_velocity_ = nullptr;
    delete trial_conditions_;
    trial_conditions_ = nullptr;
    delete model_;
    model_ = nullptr;
    delete dofs_;
    dofs_ = nullptr;
    d->source_file = domain_->SourceFile();
    d->path_kind = DiagnosticPathKind::kLiteralFile;
    d->source_line = d->first_source_line = d->node_id = d->element_id = 0;
    d->duplicate_boundary = false;
    d->condition_dof = 0;
    d->step = number;
    status = ApplyStepControls(deck_, number, controls_, d);
    if (status != StatusCode::kOk) return status;
    status = log_->WriteControls(controls_, number, d);
    if (status != StatusCode::kOk) return status;
    status = BuildEffectiveConditions(domain_, number, conditions_,
                                      &trial_conditions_, d);
    if (status != StatusCode::kOk) return status;
    status = FiniteBeam(domain_)
        ? DofManager::CreateForNonlinearBeam(domain_, trial_conditions_,
                                             &dofs_, d)
        : DofManager::Create(domain_, trial_conditions_, &dofs_, d);
    if (status != StatusCode::kOk) return status;
    status = state_->ValidateMapping(dofs_, d);
    if (status != StatusCode::kOk) return status;
    status = AnalysisModel::Create(&model_, d);
    if (status != StatusCode::kOk) return status;
    status = model_->BindStep(domain_, trial_conditions_, d);
    if (status != StatusCode::kOk) return status;
    if (step->IsDynamic()) {
      status = abaqus::ReadDynamicSettings(deck_, number, &dynamic_settings_, d);
      if (status != StatusCode::kOk) return status;
      bool explicit_automatic = false;
      for (const auto* block : *deck_->Blocks())
        if (block->step == number && *block->keyword == "DYNAMIC")
          explicit_automatic = block->HasOption("INCREMENTATION");
      const auto* time = dynamic_settings_->Parameters();
      const bool automatic = dynamic_settings_->IsAutomatic();
      char time_text[700] = {};
      std::snprintf(
          time_text, sizeof(time_text),
          "HHT control=%s selection=%s source_line=%d alpha=%.17g "
          "initial_input=%.17g minimum_input=%.17g maximum_input=%.17g "
          "minimum=%.17g maximum=%.17g nominal=%.17g INC=%d\n",
          automatic ? "AGGRESSIVE" : "DIRECT",
          automatic ? (explicit_automatic ? "INP" : "FESA-default") : "INP",
          dynamic_settings_->SourceLine(), time->alpha, time->initial_increment,
          time->minimum, time->maximum, dynamic_settings_->MinimumIncrement(),
          dynamic_settings_->MaximumIncrement(),
          dynamic_settings_->NominalIncrement(), time->increment_limit);
      status = log_->Write(time_text, d);
      if (status != StatusCode::kOk) return status;
      if (automatic) {
        const auto* increments = controls_->increments;
        std::snprintf(
            time_text, sizeof(time_text),
            "AGGRESSIVE convergence history only; no half-increment accuracy "
            "control. FESA default differs from Abaqus TF default. "
            "effective_IL=%d effective_IG=%d\n",
            increments->slow_success_iterations == 0
                ? 10 : increments->slow_success_iterations,
            increments->quick_iterations == 0 ? 4 : increments->quick_iterations);
        status = log_->Write(time_text, d);
        if (status != StatusCode::kOk) return status;
      }
      if (number == 1) {
        status = abaqus::ResolveInitialVelocity(deck_, domain_, dofs_,
                                                &initial_velocity_, d);
        if (status != StatusCode::kOk) return status;
      }
      const auto mode = DynamicMassMode(controls_->mass,
                                        step->HasFiniteGeometry(),
                                        FiniteBeam(domain_));
      if (step->HasFiniteGeometry()) {
        NonlinearDynamicAnalysis* created = nullptr;
        status = NonlinearDynamicAnalysis::Create(
            dynamic_settings_, mode, conditions_, initial_velocity_, &created,
            d);
        analysis_ = created;
      } else {
        LinearDynamicAnalysis* created = nullptr;
        status =
            LinearDynamicAnalysis::Create(dynamic_settings_, mode, conditions_,
                                          initial_velocity_, &created, d);
        analysis_ = created;
      }
    } else if (step->NonlinearParameters()) {
      NonlinearStaticAnalysis* created = nullptr;
      status = NonlinearStaticAnalysis::Create(&created, d,
                                               request_->b33_material_first);
      analysis_ = created;
    } else {
      StaticAnalysis* created = nullptr;
      status = StaticAnalysis::Create(&created, d);
      analysis_ = created;
    }
    if (status != StatusCode::kOk) return status;
    char start_text[160] = {};
    std::snprintf(
        start_text, sizeof(start_text),
        "Step %d start; increment=%d step_time=%.17g total_time=%.17g\n",
        number, step->IsDynamic() ? 0 : 1,
        step->IsDynamic() ? 0 : step->Period(),
        completed_step_total_time_ + (step->IsDynamic() ? 0 : step->Period()));
    status = log_->Write(start_text, d);
    if (status != StatusCode::kOk) return status;
    calculation_started_ = WallSeconds();
    status = analysis_->SolveStep(model_, dofs_, trial_conditions_, controls_,
                                  state_, writer_, next_frame_number_,
                                  completed_step_total_time_, analysis_report_,
                                  this, d);
    if (status != StatusCode::kOk) return status;
    if (!analysis_report_->step_completed)
      return model_storage::Invalid("Analysis step did not complete", d);
  }
  return FinishOutput(d);
}
StatusCode FESASolver::FinishOutput(Diagnostic* d) {
  StatusCode status = StatusCode::kOk;
  writer_finalized_ = true;
  const auto finish_start = std::chrono::steady_clock::now();
  status = writer_->Finish(true, last_completed_step_, d);
  if (status != StatusCode::kOk) return status;
  writer_succeeded_ = true;
  output_seconds_ += std::chrono::duration<double>(
                         std::chrono::steady_clock::now() - finish_start)
                         .count();
  char summary[240] = {};
  std::snprintf(summary, sizeof(summary),
                "Output completed; elapsed_seconds=%.9g "
                "(payload/completion/final flush and close)\n",
                output_seconds_);
  return log_->Write(summary, d);
}
StatusCode FESASolver::Fail(StatusCode status, Diagnostic* d) {
  if (thermal_state_) thermal_state_->DiscardCandidate();
  if (thermal_model_) thermal_model_->Clear();
  delete trial_thermal_conditions_;
  trial_thermal_conditions_ = nullptr;
  if (state_) state_->DiscardCandidate();
  Console(d);
  if (heat_transfer_ && analysis_report_ &&
      analysis_report_->logging_diagnostic &&
      !SameDiagnosticEvent(d, analysis_report_->logging_diagnostic)) {
    std::fputs("First logging failure: ", stderr);
    Console(analysis_report_->logging_diagnostic);
  }
  Diagnostic secondary;
  if (writer_ && (!writer_finalized_ || writer_succeeded_) &&
      writer_->Finish(false, last_completed_step_, &secondary) !=
          StatusCode::kOk) {
    Console(&secondary);
  }
  if (writer_) {
    // Any failed run with opened output requires reanalysis, including frame
    // failure with a successful failed-file close and primary Finish failure.
    constexpr char kGuidance[] =
        "Output completeness is not guaranteed; rerun analysis.\n";
    std::fputs(kGuidance, stderr);
    if (log_ && !log_finalized_ && log_->Write(kGuidance, &secondary) !=
                                      StatusCode::kOk)
      Console(&secondary);
  }
  writer_succeeded_ = false;
  if (log_ && !log_finalized_) {
    char completed[80] = {};
    std::snprintf(completed, sizeof(completed), "last_completed_step=%d\n",
                  last_completed_step_);
    if (log_->Write(completed, &secondary) != StatusCode::kOk)
      Console(&secondary);
    if (!heat_transfer_ && analysis_report_ &&
        analysis_report_->residual_valid &&
        status == StatusCode::kEquilibriumFailure) {
      if (LogResidual(analysis_report_->displacement_field,
                      "displacement", log_, &secondary) != StatusCode::kOk)
        Console(&secondary);
      if (LogResidual(analysis_report_->rotation_field, "rotation", log_,
                      &secondary) != StatusCode::kOk)
        Console(&secondary);
    }
    try {
      std::string text("");
      const auto formatted = FormatDiagnostic(d, &text);
      if (formatted != StatusCode::kOk)
        Console(d);
      else if (log_->Write(text.c_str(), &secondary) != StatusCode::kOk) {
        Console(d);
        Console(&secondary);
      }
    } catch (const std::bad_alloc&) {
      Console(d);
    }
    if (log_->Finish(&secondary) != StatusCode::kOk) Console(&secondary);
  }
  return status;
}
}  // namespace fesa
