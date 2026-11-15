#include "fesa/analysis/dynamic_analysis.h"

#include <chrono>
#include <cmath>
#include <limits>
#include <memory>
#include <new>
#include <utility>

#include "analysis/direct_increment.h"
#include "base/diagnostic_copy.h"
#include "fesa/analysis/dynamic_increment_controller.h"
#include "fesa/analysis/dynamic_state.h"

namespace fesa {
namespace {
double Seconds() {
  return std::chrono::duration<double>(
             std::chrono::steady_clock::now().time_since_epoch())
      .count();
}
}  // namespace
DynamicAnalysis::~DynamicAnalysis() {
  delete integrator_;
  delete logging_error_;
}
StatusCode DynamicAnalysis::InitializeIntegration(
    const DynamicStepSettings* settings, Diagnostic* d) {
  if (!settings || settings_ || integrator_)
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation, "Dynamic settings", d);
  try {
    TimeIntegrator* raw = nullptr;
    const auto status =
        CreateHhtIntegrator(settings->IntegratorParameters(), &raw, d);
    std::unique_ptr<TimeIntegrator> integrator(raw);
    if (status != StatusCode::kOk) return status;
    std::unique_ptr<Diagnostic> error(new Diagnostic);
    integrator_ = integrator.release();
    logging_error_ = error.release();
    settings_ = settings;
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return ReportFailure(StatusCode::kAllocationFailure,
                         DiagnosticStage::kAllocation, "Dynamic integration",
                         d);
  }
}
const DynamicStepSettings* DynamicAnalysis::Settings() const {
  return settings_;
}
const TimeIntegrator* DynamicAnalysis::Integrator() const {
  return integrator_;
}
StatusCode DynamicAnalysis::Notify(const DynamicStepContext* context,
                                   AnalysisEvent event, StatusCode primary,
                                   Diagnostic* d) {
  StatusCode status;
  try {
    status =
        context->progress->OnProgress(event, context->report, logging_error_);
  } catch (const std::bad_alloc&) {
    status = ReportFailure(StatusCode::kAllocationFailure,
                           DiagnosticStage::kAllocation, "Dynamic progress",
                           logging_error_);
  }
  if (status != StatusCode::kOk) {
    if (!context->report->logging_diagnostic)
      context->report->logging_diagnostic = logging_error_;
    if (primary == StatusCode::kOk) {
      CopyDiagnosticMetadata(logging_error_, d);
      return status;
    }
  }
  return primary;
}
StatusCode DynamicAnalysis::CalculateFrame(const DynamicStepContext* context,
                                           bool initial, double increment_size,
                                           double step_time, double total_time,
                                           bool* numerical_failure,
                                           Diagnostic* d) {
  if (numerical_failure) *numerical_failure = false;
  auto* report = context->report;
  report->calculation_finished = false;
  report->residual_valid = false;
  report->residual = nullptr;
  report->displacement_field = report->rotation_field = nullptr;
  const double started = Seconds();
  StatusCode status;
  try {
    status = initial ? StatusCode::kOk
                     : NotifyIncrement(context, IncrementEvent::kAttemptStarted,
                                        StatusCode::kOk, d);
    if (status == StatusCode::kOk) {
      status = initial ? InitializeCandidate(context, d)
                       : AdvanceCandidate(context, increment_size, step_time,
                                          total_time, d);
      if (numerical_failure)
        *numerical_failure = !initial && status != StatusCode::kOk;
    }
    if (status == StatusCode::kOk) status = RecoverCandidate(context, d);
  } catch (const std::bad_alloc&) {
    status =
        ReportFailure(StatusCode::kAllocationFailure,
                      DiagnosticStage::kAllocation, "Dynamic calculation", d);
  }
  report->calculation_seconds += Seconds() - started;
  report->calculation_status = status;
  report->calculation_finished = true;
  return Notify(context, AnalysisEvent::kCalculationFinished, status, d);
}
StatusCode DynamicAnalysis::NotifyIncrement(const DynamicStepContext* context,
                                            IncrementEvent event,
                                            StatusCode primary,
                                            Diagnostic* d) {
  // Linear receivers are not required to implement the detailed Newton log.
  if (!context->state->FieldHistory()) return primary;
  StatusCode status;
  try {
    status = context->progress->OnIncrementProgress(
        event, context->attempt, context->newton,
        context->report->displacement_field, context->report->rotation_field,
        primary == StatusCode::kOk ? nullptr : d, logging_error_);
  } catch (const std::bad_alloc&) {
    status = ReportFailure(StatusCode::kAllocationFailure,
                           DiagnosticStage::kAllocation,
                           "Dynamic increment progress", logging_error_);
  }
  if (status != StatusCode::kOk) {
    context->report->logging_diagnostic = logging_error_;
    if (primary == StatusCode::kOk) {
      CopyDiagnosticMetadata(logging_error_, d);
      return status;
    }
  }
  return primary;
}
StatusCode DynamicAnalysis::PublishFrame(const DynamicStepContext* context,
                                         const FrameInfo* frame,
                                         Diagnostic* d) {
  auto* state = context->state;
  const auto* motion = state->Dynamics();
  if (!motion || !motion->CandidateVelocity() ||
      motion->CandidateStepTime() != frame->step_time ||
      motion->CandidateTotalTime() != frame->total_time ||
      state->Results() != ResultState::kUnconfirmed ||
      state->ResultStepNumber() != frame->step_number ||
      state->ResultIncrementNumber() != frame->increment_number)
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation,
                         "Dynamic frame candidate", d);
  const double started = Seconds();
  StatusCode status;
  try {
    status = context->writer->WriteFrame(
        frame, context->conditions->SourceDomain(), state, d);
  } catch (const std::bad_alloc&) {
    status = ReportFailure(StatusCode::kAllocationFailure,
                           DiagnosticStage::kAllocation, "Dynamic output", d);
  }
  if (status == StatusCode::kOk) {
    auto* report = context->report;
    state->Commit(report->step_number, report->increment_number);
    report->last_committed_step = report->step_number;
    report->last_committed_increment = report->increment_number;
    report->last_committed_frame = frame->frame_number;
    report->last_committed_total_time = frame->total_time;
    report->step_completed = frame->completes_step;
  }
  context->report->output_seconds += Seconds() - started;
  return status;
}
StatusCode DynamicAnalysis::SolveStep(
    const AnalysisModel* model, const DofManager* dofs,
    const EffectiveConditions* conditions, const ControlSettings* controls,
    AnalysisState* state, ResultWriter* writer, std::int64_t first_frame_number,
    double total_time_before_step, AnalysisReport* report,
    AnalysisProgress* progress, Diagnostic* d) {
  if (report) report->Reset();
  if (logging_error_) {
    Diagnostic empty;
    CopyDiagnosticMetadata(&empty, logging_error_);
  }
  StatusCode status = StatusCode::kOk;
  DynamicStepContext context;
  context.model = model;
  context.dofs = dofs;
  context.conditions = conditions;
  context.controls = controls;
  context.state = state;
  context.writer = writer;
  context.report = report;
  context.progress = progress;
  context.total_time_before_step = total_time_before_step;
  const StepDefinition* step = nullptr;
  if (!settings_ || !model || !dofs || !conditions || !controls || !state ||
      !writer || !report || !progress ||
      (settings_ && settings_->IsAutomatic() &&
       (!controls || !controls->increments)) || first_frame_number <= 0 ||
      !std::isfinite(total_time_before_step) || total_time_before_step < 0)
    status =
        ReportFailure(StatusCode::kInvalidArgument,
                      DiagnosticStage::kValidation, "Dynamic step inputs", d);
  if (status == StatusCode::kOk) {
    status = conditions->SourceDomain()->FindStepDefinition(
        conditions->StepNumber(), &step, d);
    if (status == StatusCode::kOk &&
        (!step->IsDynamic() || step->Id() != settings_->Step() ||
         step->Period() != settings_->Parameters()->period ||
         state->StepNumber() == std::numeric_limits<int>::max() ||
         state->StepNumber() + 1 != step->Id() ||
         !std::isfinite(total_time_before_step + step->Period())))
      status = ReportFailure(StatusCode::kInvalidArgument,
                             DiagnosticStage::kValidation,
                             "Dynamic step identity/time", d);
  }
  IncrementClock first, second;
  IncrementClock* clock = &first;
  IncrementClock* next = &second;
  DynamicIncrementController automatic(settings_);
  if (status == StatusCode::kOk) {
    const auto* parameters = settings_->Parameters();
    status = settings_->IsAutomatic()
                 ? automatic.Initialize(clock, d)
                 : InitializeDirectIncrementClock(
                       parameters->initial_increment, parameters->period,
                       parameters->minimum, parameters->maximum, clock, d);
    report->step_number = step->Id();
    report->last_committed_step = state->StepNumber();
    report->last_committed_increment = state->IncrementNumber();
    report->last_committed_frame = first_frame_number - 1;
    report->last_committed_total_time = total_time_before_step;
    report->equation_count = dofs->EquationCount();
    const double started = Seconds();
    if (status == StatusCode::kOk) {
      try {
        status = BeginStep(&context, d);
      } catch (const std::bad_alloc&) {
        status = ReportFailure(StatusCode::kAllocationFailure,
                               DiagnosticStage::kAllocation,
                               "Dynamic step preparation", d);
      }
    }
    report->calculation_seconds += Seconds() - started;
    if (status != StatusCode::kOk) {
      report->calculation_status = status;
      report->calculation_finished = true;
      status = Notify(&context, AnalysisEvent::kCalculationFinished, status, d);
    }
  }
  FrameInfo frame;
  if (status == StatusCode::kOk) {
    frame.frame_number = first_frame_number;
    frame.step_number = step->Id();
    frame.step_name = step->Name();
    frame.increment_number = 0;
    frame.step_time = 0;
    frame.total_time = total_time_before_step;
    frame.procedure = FrameProcedure::kDynamic;
    frame.geometry = state->FieldHistory() ? FrameGeometry::kNonlinear
                                           : FrameGeometry::kLinear;
    frame.kind = FrameKind::kInitial;
    frame.completes_step = false;
    status = CalculateFrame(&context, true, 0, 0, total_time_before_step,
                             nullptr, d);
    if (status == StatusCode::kOk) status = PublishFrame(&context, &frame, d);
    if (status == StatusCode::kOk)
      status = Notify(&context, AnalysisEvent::kIncrementCommitted, status, d);
  }
  while (status == StatusCode::kOk && clock->time < step->Period()) {
    double nominal = clock->nominal;
    int attempt_number = 1;
    while (status == StatusCode::kOk) {
      IncrementTarget target;
      status = settings_->IsAutomatic()
                   ? automatic.Propose(controls->increments, clock, nominal,
                                        &target, d)
                   : ProposeDirectIncrement(
                         step->Period(), settings_->Parameters()->increment_limit,
                         clock, &target, d);
      if (status != StatusCode::kOk) break;
      const double total_time = total_time_before_step + target.time;
      if (report->last_committed_frame ==
              std::numeric_limits<std::int64_t>::max() ||
          !std::isfinite(total_time) ||
          total_time <= report->last_committed_total_time) {
        status = ReportFailure(StatusCode::kLibraryFailure,
                               DiagnosticStage::kComputation,
                               "Dynamic frame/time advancement", d);
        break;
      }
      frame.frame_number = report->last_committed_frame + 1;
      frame.increment_number = clock->successful + 1;
      frame.step_time = target.time;
      frame.total_time = total_time;
      frame.kind = FrameKind::kIncrementEnd;
      frame.completes_step = target.completes_step;
      report->increment_number = clock->successful + 1;
      report->attempt_count = attempt_number;
      report->iteration_count = 0;
      IncrementAttemptInfo attempt;
      attempt.step_number = step->Id();
      attempt.increment_number = report->increment_number;
      attempt.attempt_number = attempt_number;
      attempt.step_time = target.time;
      attempt.total_time = total_time;
      attempt.increment_size = target.size;
      attempt.nominal_increment = target.nominal;
      attempt.completes_step = target.completes_step;
      NewtonIterationInfo outcome;
      context.attempt = &attempt;
      context.newton = &outcome;
      bool numerical_failure = false;
      status = CalculateFrame(&context, false, target.size, target.time,
                               total_time, &numerical_failure, d);
      if (status != StatusCode::kOk) {
        if (!numerical_failure ||
            outcome.outcome == NewtonOutcome::kNotStarted) {
          outcome.outcome = NewtonOutcome::kFatal;
          outcome.status = outcome.numerical_status = status;
        }
        if (!report->logging_diagnostic)
          status = NotifyIncrement(&context, IncrementEvent::kAttemptFailed,
                                    status, d);
        if (!settings_->IsAutomatic() || !numerical_failure ||
            report->logging_diagnostic)
          break;
        state->DiscardCandidate();
        double retry = 0;
        const auto cutback = automatic.Cutback(
            controls->increments, clock, &target, &outcome, attempt_number,
            &retry, d);
        if (cutback != StatusCode::kOk) {
          status = cutback;
          break;
        }
        attempt.proposed_available = true;
        attempt.proposed_increment = retry;
        status = NotifyIncrement(&context, IncrementEvent::kCutback,
                                  StatusCode::kOk, d);
        if (status != StatusCode::kOk) break;
        nominal = retry;
        ++attempt_number;
        continue;
      }
      if (settings_->IsAutomatic()) {
        const auto* actual = outcome.outcome == NewtonOutcome::kNotStarted
                                 ? nullptr : &outcome;
        status = automatic.Accept(controls->increments, clock, &target, actual,
                                   next, d);
      } else {
        status = AcceptDirectIncrement(
            step->Period(), settings_->Parameters()->increment_limit, clock,
            &target, status, next, d);
      }
      if (status == StatusCode::kOk) status = PublishFrame(&context, &frame, d);
      if (status == StatusCode::kOk) {
        // Same no-fail commit gate as PublishFrame's physical/history commit.
        std::swap(clock, next);
        status = NotifyIncrement(&context, IncrementEvent::kCommitted,
                                  StatusCode::kOk, d);
        if (status == StatusCode::kOk)
          status = Notify(&context, AnalysisEvent::kIncrementCommitted, status, d);
      }
      break;
    }
  }
  if (status != StatusCode::kOk && state) state->DiscardCandidate();
  ReleaseStepResources();
  if (report) report->status = status;
  return status;
}
}  // namespace fesa
