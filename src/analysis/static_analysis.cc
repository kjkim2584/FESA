#include "fesa/analysis/static_analysis.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <memory>
#include <new>
#include <vector>

#include "adapters/mkl/pardiso_solver.h"
#include "adapters/openmp/thread_runtime.h"
#include "analysis/analysis_access.h"
#include "analysis/equilibrium_check.h"
#include "analysis/linear_result_recovery.h"
#include "analysis/static_assembler.h"
#include "fesa/elements/solid_element.h"
#include "fesa/elements/solid_reference_data.h"
#include "model/model_storage.h"
namespace fesa {
namespace {
thread_local void (*analysis_observer)(AnalysisStage) = nullptr;
void NotifyAnalysisStage(AnalysisStage stage) {
  if (analysis_observer) analysis_observer(stage);
}
}  // namespace
void SetAnalysisObserver(void (*callback)(AnalysisStage)) {
  analysis_observer = callback;
}
namespace {
StatusCode NotifyProgress(AnalysisProgress* progress, AnalysisEvent event,
                          const AnalysisReport* report, Diagnostic* d) {
  try {
    return progress->OnProgress(event, report, d);
  } catch (const std::bad_alloc&) {
    return model_storage::AllocationFailure("Analysis progress allocation", d);
  }
}
StatusCode PublishLinearStep(const FrameInfo* frame, const Domain* domain,
                             AnalysisState* state, ResultWriter* writer,
                             AnalysisReport* report, AnalysisProgress* progress,
                             Diagnostic* d) {
  const auto started = std::chrono::steady_clock::now();
  auto status = writer->WriteFrame(frame, domain, state, d);
  if (status == StatusCode::kOk) {
    state->Commit(report->step_number, report->increment_number);
    report->step_completed = true;
    report->last_committed_step = report->step_number;
    report->last_committed_increment = report->increment_number;
    report->last_committed_frame = frame->frame_number;
    report->last_committed_total_time = frame->total_time;
  }
  report->output_seconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - started)
          .count();
  if (status == StatusCode::kOk)
    status =
        NotifyProgress(progress, AnalysisEvent::kIncrementCommitted, report, d);
  return status;
}
}  // namespace
StatusCode StaticAnalysis::Create(StaticAnalysis** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("StaticAnalysis output", d);
  *out = nullptr;
  try {
    std::unique_ptr<StaticAnalysis> result(new StaticAnalysis);
    result->solver_ = new PardisoSolver;
    const auto status = LinearAlgorithm::Create(&result->algorithm_, d);
    if (status != StatusCode::kOk) return status;
    *out = result.release();
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return model_storage::AllocationFailure("StaticAnalysis allocation", d);
  }
}
StaticAnalysis::~StaticAnalysis() {
  ClearScratch();
  delete algorithm_;
  delete solver_;
}
void StaticAnalysis::ClearScratch() {
  ReleaseWorkers();
  ReleaseAttemptResources();
  delete report_;
  report_ = nullptr;
}
void StaticAnalysis::ReleaseWorkers() noexcept {
  for (int i = 0; i < worker_count_; ++i) delete workers_[i];
  delete[] workers_;
  workers_ = nullptr;
  worker_count_ = 0;
  for (int i = 0; i < solid_reference_count_; ++i) delete solid_references_[i];
  delete[] solid_references_;
  solid_references_ = nullptr;
  solid_reference_count_ = 0;
}
void StaticAnalysis::ReleaseEquationSystem() noexcept {
  delete equations_;
  equations_ = nullptr;
}
void StaticAnalysis::ReleaseSolution() noexcept {
  if (workspace_) workspace_->ReleaseReducedBuffers();
}
void StaticAnalysis::ReleaseAttemptResources() noexcept {
  ReleaseEquationSystem();
  delete workspace_;
  workspace_ = nullptr;
}
StatusCode StaticAnalysis::SolveStep(
    const AnalysisModel* model, const DofManager* dofs,
    const EffectiveConditions* conditions, const ControlSettings* controls,
    AnalysisState* state, ResultWriter* writer, std::int64_t first_frame_number,
    double total_time_before_step, AnalysisReport* report,
    AnalysisProgress* progress, Diagnostic* d) {
  if (report) report->Reset();
  ClearScratch();
  StatusCode status = StatusCode::kOk;
  assembly_seconds_ = solve_seconds_ = recovery_seconds_ = 0;
  try {
    if (!model || !dofs || !conditions || !controls || !state || !writer ||
        !report || !progress || first_frame_number <= 0 ||
        !std::isfinite(total_time_before_step) || total_time_before_step < 0)
      status = model_storage::Invalid("StaticAnalysis input", d);
    else if (state->StepNumber() == std::numeric_limits<int>::max() ||
             conditions->StepNumber() != state->StepNumber() + 1)
      status =
          ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                        "StaticAnalysis requires sequential committed step", d);
    if (status == StatusCode::kOk) {
      const auto* domain = conditions->SourceDomain();
      const auto found = domain->Steps()->find(conditions->StepNumber());
      if (found != domain->Steps()->end() && found->second->IsDynamic()) {
        status = ReportFailure(StatusCode::kUnsupported,
                               DiagnosticStage::kValidation,
                               "StaticAnalysis cannot execute dynamic step", d);
      } else if (found == domain->Steps()->end() ||
                 !std::isfinite(total_time_before_step +
                                found->second->Period())) {
        status = model_storage::Invalid("StaticAnalysis step/time", d);
      } else {
        const auto* step = found->second;
        report->step_number = conditions->StepNumber();
        report->increment_number = 1;
        report->attempt_count = 1;
        report->last_committed_step = state->StepNumber();
        report->last_committed_increment = state->IncrementNumber();
        report->last_committed_frame = first_frame_number - 1;
        report->last_committed_total_time = total_time_before_step;
        status = CalculateWithReport(model, dofs, conditions, controls, state,
                                     report, d);
        Diagnostic secondary;
        const auto notified =
            NotifyProgress(progress, AnalysisEvent::kCalculationFinished,
                           report, status == StatusCode::kOk ? d : &secondary);
        if (status == StatusCode::kOk) status = notified;
        if (status == StatusCode::kOk) {
          FrameInfo frame;
          frame.frame_number = first_frame_number;
          frame.step_number = conditions->StepNumber();
          frame.increment_number = 1;
          frame.step_name = step->Name();
          frame.step_time = step->Period();
          frame.total_time = total_time_before_step + step->Period();
          status = PublishLinearStep(&frame, domain, state, writer, report,
                                     progress, d);
        }
      }
    }
  } catch (const std::bad_alloc&) {
    status = model_storage::AllocationFailure("StaticAnalysis allocation", d);
  }
  if (status != StatusCode::kOk) {
    if (state) state->DiscardCandidate();
    if (status == StatusCode::kEquilibriumFailure) {
      auto* retained = report_;
      report_ = nullptr;
      ClearScratch();
      report_ = retained;
    } else {
      if (report) {
        report->residual = nullptr;
        report->displacement_field = report->rotation_field = nullptr;
        report->residual_valid = false;
      }
      ClearScratch();
    }
  }
  if (report) report->status = status;
  return status;
}
StatusCode StaticAnalysis::CalculateWithReport(
    const AnalysisModel* model, const DofManager* dofs,
    const EffectiveConditions* conditions, const ControlSettings* controls,
    AnalysisState* state, AnalysisReport* report, Diagnostic* d) {
  const auto started = std::chrono::steady_clock::now();
  const double element_started = ElementDispatchSeconds();
  const auto solves_started = PardisoSolveCallCount();
  StatusCode status = StatusCode::kOk;
  try {
    status = Calculate(model, dofs, conditions, controls, state, d);
  } catch (const std::bad_alloc&) {
    status = model_storage::AllocationFailure("StaticAnalysis allocation", d);
  }
  // Official recovery has finished (or failed). No later consumer borrows
  // the worker context/catalog or the global residual scratch.
  ReleaseAttemptResources();
  ReleaseWorkers();
  model_ = nullptr;
  dofs_ = nullptr;
  conditions_ = nullptr;
  controls_ = nullptr;
  state_ = nullptr;
  report->calculation_seconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - started)
          .count();
  report->assembly_seconds = assembly_seconds_;
  report->solve_seconds = solve_seconds_;
  report->equilibrium_recovery_seconds = recovery_seconds_;
  report->element_dispatch_seconds = ElementDispatchSeconds() - element_started;
  report->observed_element_team_max = ObservedElementTeamMaximum();
  // A direct linear step enters the reduced solver at most once.
  report->linear_solver_calls =
      static_cast<int>(PardisoSolveCallCount() - solves_started);
  report->calculation_finished = true;
  report->equation_count = dofs->FreeCount();
  report->status = report->calculation_status = status;
  if (status == StatusCode::kOk || status == StatusCode::kEquilibriumFailure) {
    report->residual = report_;
    report->residual_valid = report_ != nullptr;
    if (report_) {
      report->displacement_field = report_->Displacement();
      report->rotation_field = report_->Rotation();
      report->characteristic_length = report_->CharacteristicLength();
    }
  }
  return status;
}
StatusCode StaticAnalysis::Calculate(const AnalysisModel* model,
                                     const DofManager* dofs,
                                     const EffectiveConditions* conditions,
                                     const ControlSettings* controls,
                                     AnalysisState* state, Diagnostic* d) {
  auto status = StaticWorkspace::Create(dofs, &workspace_, d);
  if (status != StatusCode::kOk) return status;
  const int count = ConfiguredThreadCount();
  std::vector<SolidFormulation> formulations(0);
  for (const auto& entry : *model->Elements()) {
    const auto* solid = dynamic_cast<const SolidElement*>(entry.second);
    if (solid && std::find(formulations.begin(), formulations.end(),
                           solid->formulation()) == formulations.end())
      formulations.push_back(solid->formulation());
  }
  if (!formulations.empty()) {
    const int stiffness_count = static_cast<int>(formulations.size());
    const bool uniform_hexa =
        std::find(formulations.begin(), formulations.end(),
                  SolidFormulation::kHexa8Reduced) != formulations.end();
    const int reference_count = stiffness_count + (uniform_hexa ? 1 : 0);
    solid_references_ = new SolidReferenceData* [reference_count] {};
    solid_reference_count_ = reference_count;
    for (int i = 0; i < stiffness_count; ++i) {
      status = SolidReferenceData::Create(formulations[i],
                                          IntegrationPurpose::kStiffness,
                                          &solid_references_[i], d);
      if (status != StatusCode::kOk) return status;
    }
    if (uniform_hexa) {
      // Uniform strain needs averaged geometry, independently of its single
      // physical constitutive point. Static analysis never requests mass.
      status = SolidReferenceData::Create(
          SolidFormulation::kHexa8Reduced, IntegrationPurpose::kGeometry,
          &solid_references_[stiffness_count], d);
      if (status != StatusCode::kOk) return status;
    }
    workers_ = new ElementWorkspace* [count] {};
    worker_count_ = count;
    for (int i = 0; i < count; ++i) {
      status = ElementWorkspace::CreateForSolidCatalog(
          solid_references_, solid_reference_count_, &workers_[i], d);
      if (status != StatusCode::kOk) return status;
    }
  } else {
    status = PrepareWorkerWorkspaces(count, 12, &workers_, d);
    if (status != StatusCode::kOk) return status;
    worker_count_ = count;
  }
  status = ResidualReport::Create(dofs->GlobalCount(), &report_, d);
  if (status != StatusCode::kOk) return status;
  status = state->PrepareCandidate(d);
  if (status != StatusCode::kOk) return status;
  for (const auto* boundary : *conditions->Boundaries()) {
    const int index =
        dofs->GlobalDof(boundary->node_id, boundary->definition->Dof());
    status = state->CandidateDisplacement()->Set(
        index, boundary->definition->Value(), d);
    if (status != StatusCode::kOk) return status;
  }
  model_ = model;
  dofs_ = dofs;
  conditions_ = conditions;
  controls_ = controls;
  state_ = state;
  status = algorithm_->Solve(this, solver_, d);
  solve_seconds_ = algorithm_->SolveSeconds();
  model_ = nullptr;
  dofs_ = nullptr;
  conditions_ = nullptr;
  controls_ = nullptr;
  state_ = nullptr;
  if (status != StatusCode::kOk) return status;
  const auto recovery_start = std::chrono::steady_clock::now();
  status = state->StageFluxHistory(report_->Displacement()->raw_mean,
                                   report_->Rotation()->raw_mean,
                                   report_->Displacement()->maximum_flux,
                                   report_->Rotation()->maximum_flux, d);
  if (status != StatusCode::kOk) return status;
  status = state->CandidateReaction()->CopyFrom(report_->Reaction(), d);
  if (status != StatusCode::kOk) return status;
  status = RecoverLinearResults(model, dofs, state, state->StepNumber() + 1, 1,
                                workers_, worker_count_, d, NotifyAnalysisStage);
  recovery_seconds_ += std::chrono::duration<double>(
      std::chrono::steady_clock::now() - recovery_start).count();
  return status;
}
StatusCode StaticAnalysis::AssembleSystem(
    const EquationSystem** system, Vector** solution, bool* solve_required,
    Diagnostic* d) {
  *system = nullptr;
  *solution = nullptr;
  *solve_required = false;
  const auto started = std::chrono::steady_clock::now();
  const auto status = AssembleReducedSystem(model_, dofs_, conditions_,
      workspace_, &equations_, d, workers_, worker_count_);
  assembly_seconds_ = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - started).count();
  if (status != StatusCode::kOk) return status;
  *solve_required = dofs_->FreeCount() > 0;
  *system = equations_;
  *solution = workspace_->ReducedSolution();
  return StatusCode::kOk;
}
StatusCode StaticAnalysis::ApplySolution(const Vector* solution, Diagnostic* d) {
  std::unique_ptr<int[]> indices(new int[dofs_->FreeCount()]);
  for (int i = 0; i < dofs_->FreeCount(); ++i)
    indices[i] = dofs_->GlobalDofForEquation(i);
  Vector* old_raw = nullptr;
  auto status = Vector::Create(dofs_->FreeCount(), &old_raw, d);
  std::unique_ptr<Vector> old(old_raw);
  if (status != StatusCode::kOk) return status;
  status = old->GatherFrom(state_->CandidateDisplacement(), indices.get(),
                           dofs_->FreeCount(), 0, d);
  if (status != StatusCode::kOk) return status;
  status = old->Scale(-1, d);
  if (status != StatusCode::kOk) return status;
  status = state_->CandidateDisplacement()->ScatterAdd(
      old.get(), indices.get(), dofs_->FreeCount(), d);
  if (status != StatusCode::kOk) return status;
  return state_->CandidateDisplacement()->ScatterAdd(
      solution, indices.get(), dofs_->FreeCount(), d);
}
StatusCode StaticAnalysis::EvaluateResidual(Diagnostic* d) {
  const auto started = std::chrono::steady_clock::now();
  const auto status = CheckEquilibrium(model_, dofs_, conditions_, controls_,
      state_->CandidateDisplacement(), state_, workspace_, report_, d,
      workers_, worker_count_);
  recovery_seconds_ = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - started).count();
  return status;
}
}  // namespace fesa
