#include "fesa/analysis/newton_raphson.h"

#include <algorithm>
#include <chrono>
#include <limits>
#include <memory>
#include <new>

#include "analysis/nonlinear_convergence.h"
#include "base/diagnostic_copy.h"
#include "model/model_storage.h"
namespace fesa {
namespace {
class AttemptCleanup final {
 public:
  explicit AttemptCleanup(NonlinearProblem* problem) : problem_(problem) {}
  ~AttemptCleanup() {
    if (problem_) problem_->ReleaseAttemptResources();
  }
  AttemptCleanup(const AttemptCleanup&) = delete;
  AttemptCleanup& operator=(const AttemptCleanup&) = delete;
  AttemptCleanup(AttemptCleanup&&) = delete;
  AttemptCleanup& operator=(AttemptCleanup&&) = delete;

 private:
  NonlinearProblem* problem_;  // Borrowed for the enclosing Solve call.
};
double Elapsed(const std::chrono::steady_clock::time_point* start) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                       *start)
      .count();
}
StatusCode Notify(NewtonProgress* progress, const NewtonReport* report,
                  Diagnostic* d) {
  try {
    return progress->OnIteration(report, d);
  } catch (const std::bad_alloc&) {
    return model_storage::AllocationFailure("Newton progress allocation", d);
  }
}
StatusCode Emit(NewtonProgress* progress, NewtonReport* report,
                NewtonIterationInfo* info, Diagnostic* d) {
  if (!progress) return info->status;
  Diagnostic secondary;
  const bool numeric_failure = info->numerical_status != StatusCode::kOk;
  if (!numeric_failure) CopyDiagnosticMetadata(d, &secondary);
  const auto log_status = Notify(progress, report, &secondary);
  if (log_status != StatusCode::kOk) {
    report->RecordLoggingFailure(log_status, &secondary);
    if (!numeric_failure) {
      CopyDiagnosticMetadata(&secondary, d);
      info->status = log_status;
      info->outcome = NewtonOutcome::kFatal;
    }
  }
  return info->status;
}
NewtonOutcome FailureOutcome(StatusCode status, NonlinearFailure failure) {
  if (status == StatusCode::kEquilibriumFailure) {
    if (failure == NonlinearFailure::kGeometry) return NewtonOutcome::kGeometry;
    if (failure == NonlinearFailure::kRotationBranch)
      return NewtonOutcome::kRotationBranch;
  }
  return NewtonOutcome::kFatal;
}
}  // namespace
StatusCode NewtonRaphson::Create(NewtonRaphson** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("NewtonRaphson output", d);
  *out = new (std::nothrow) NewtonRaphson;
  return *out ? StatusCode::kOk
              : model_storage::AllocationFailure("NewtonRaphson", d);
}
StatusCode NewtonRaphson::Solve(NonlinearProblem* problem, LinearSolver* solver,
                                const ControlSettings* controls,
                                const NewtonControls* counts,
                                NewtonReport* report, NewtonProgress* progress,
                                Diagnostic* d, NewtonExecution execution) {
  AttemptCleanup cleanup(problem);
  if (report) report->Reset();
  if (!problem || !solver || !controls || !controls->displacement ||
      !controls->rotation || !counts || counts->divergence_start < 3 ||
      counts->slow_start <= 0 || counts->nonquadratic_start <= 0 ||
      counts->prediction_limit <= 0 || counts->iteration_limit <= 0 ||
      !report) {
    const auto status =
        model_storage::Invalid("NewtonRaphson bindings/controls", d);
    if (report) {
      report->info_->status = report->info_->numerical_status = status;
      report->info_->outcome = NewtonOutcome::kFatal;
    }
    return status;
  }
  NewtonControls effective;
  if (execution == NewtonExecution::kDirect) {
    effective.divergence_start = counts->divergence_start;
    effective.slow_start = std::numeric_limits<int>::max();
    effective.nonquadratic_start = counts->nonquadratic_start;
    effective.prediction_limit = counts->prediction_limit;
    effective.iteration_limit = counts->prediction_limit;
    counts = &effective;
  }
  auto* info = report->info_;
  info->outcome = NewtonOutcome::kInProgress;
  auto fail = [&](StatusCode status, NonlinearFailure typed) {
    info->status = info->numerical_status = status;
    info->outcome = FailureOutcome(status, typed);
    return Emit(progress, report, info, d);
  };
  NewtonResidualHistory history;
  auto started = std::chrono::steady_clock::now();
  try {
    int i = 0;
    while (true) {
      info->iteration = i;
      info->attempted_iteration = i;
      info->correction_valid = i > 0;
      info->residual_valid = false;
      info->stage = NewtonStage::kResidual;
      for (int field = 0; field < 2; ++field) {
        info->estimate_available[field] = false;
        info->prediction_available[field] = false;
        info->estimate_max[field] = info->predicted_total[field] = 0;
      }
      ++info->residual_evaluations;
      double length = 0;
      NonlinearFailure typed = NonlinearFailure::kNone;
      started = std::chrono::steady_clock::now();
      auto status = problem->EvaluateResidual(
          i, report->force_, report->moment_, &length, &typed, d);
      info->residual_seconds += Elapsed(&started);
      problem->ReleaseCorrectionAssessment();
      if (status != StatusCode::kOk) return fail(status, typed);
      if (typed != NonlinearFailure::kNone)
        return fail(
            model_storage::Invalid("Residual success with failure flag", d),
            NonlinearFailure::kNone);
      info->stage = NewtonStage::kCriteria;
      NewtonConvergenceResult convergence;
      started = std::chrono::steady_clock::now();
      status = EvaluateNewtonConvergence(controls, counts, length, i, i > 0,
                                         &history, report->force_,
                                         report->moment_, &convergence, d);
      info->criteria_seconds += Elapsed(&started);
      if (status != StatusCode::kOk)
        return fail(status, NonlinearFailure::kNone);
      info->residual_valid = true;
      for (int field = 0; field < 2; ++field) {
        info->estimate_available[field] = convergence.estimate_available[field];
        info->estimate_max[field] = convergence.estimate_max[field];
        info->prediction_available[field] =
            convergence.prediction_available[field];
        info->predicted_total[field] = convergence.predicted_total[field];
      }
      // A dynamic kinematic predictor has not solved equilibrium. Keep public
      // scalar results, but honor the Analysis execution contract before
      // commit.
      if (convergence.decision == NewtonDecision::kConverged && i == 0 &&
          problem->RequiresInitialCorrection() &&
          (report->force_->dof_count > 0 || report->moment_->dof_count > 0))
        convergence.decision = NewtonDecision::kContinue;
      if (convergence.decision == NewtonDecision::kConverged) {
        info->outcome = NewtonOutcome::kConverged;
        info->stage = NewtonStage::kFinished;
        return Emit(progress, report, info, d);
      }
      if (convergence.decision != NewtonDecision::kContinue) {
        const char* operation = nullptr;
        if (convergence.decision == NewtonDecision::kDivergence) {
          info->outcome = NewtonOutcome::kDivergence;
          operation = "Newton diverging residual";
        } else if (convergence.decision == NewtonDecision::kSlowConvergence) {
          info->outcome = NewtonOutcome::kSlowConvergence;
          operation = "Newton slow residual convergence";
        } else {
          info->outcome = NewtonOutcome::kIterationLimit;
          operation = "Newton iteration limit";
        }
        info->status = info->numerical_status =
            ReportFailure(StatusCode::kEquilibriumFailure,
                          DiagnosticStage::kComputation, operation, d);
        return Emit(progress, report, info, d);
      }
      status = Emit(progress, report, info, d);
      if (status != StatusCode::kOk) return status;
      const FieldResidual* fields[2] = {report->force_, report->moment_};
      for (int field = 0; field < 2; ++field) {
        history.maximum[field][1] = history.maximum[field][0];
        history.maximum[field][0] = fields[field]->maximum;
      }
      history.count = std::min(history.count + 1, 2);
      info->attempted_iteration = i + 1;
      info->stage = NewtonStage::kTangent;
      const EquationSystem* system = nullptr;
      Vector* correction = nullptr;
      typed = NonlinearFailure::kNone;
      started = std::chrono::steady_clock::now();
      status = problem->AssembleTangent(&system, &correction, &typed, d);
      info->tangent_seconds += Elapsed(&started);
      if (status != StatusCode::kOk) return fail(status, typed);
      if (typed != NonlinearFailure::kNone || !system || !correction ||
          correction == system->Rhs() ||
          correction->Size() != system->Rhs()->Size())
        return fail(model_storage::Invalid("Newton correction system/flag", d),
                    NonlinearFailure::kNone);
      info->stage = NewtonStage::kLinearSolve;
      ++info->linear_solver_calls;
      started = std::chrono::steady_clock::now();
      status = solver->Solve(system, correction, d);
      info->solve_seconds += Elapsed(&started);
      problem->ReleaseEquationSystem();
      if (status != StatusCode::kOk)
        return fail(status, NonlinearFailure::kNone);
      info->stage = NewtonStage::kUpdate;
      typed = NonlinearFailure::kNone;
      started = std::chrono::steady_clock::now();
      status = problem->ApplyCorrection(correction, &typed, d);
      info->update_seconds += Elapsed(&started);
      problem->ReleaseFreeCorrection();
      if (status != StatusCode::kOk) return fail(status, typed);
      if (typed != NonlinearFailure::kNone)
        return fail(
            model_storage::Invalid("Update success with failure flag", d),
            NonlinearFailure::kNone);
      ++i;
    }
  } catch (const std::bad_alloc&) {
    const double failed_seconds = Elapsed(&started);
    switch (info->stage) {
      case NewtonStage::kResidual:
        info->residual_seconds += failed_seconds;
        break;
      case NewtonStage::kCriteria:
        info->criteria_seconds += failed_seconds;
        break;
      case NewtonStage::kTangent:
        info->tangent_seconds += failed_seconds;
        break;
      case NewtonStage::kLinearSolve:
        info->solve_seconds += failed_seconds;
        break;
      case NewtonStage::kUpdate:
        info->update_seconds += failed_seconds;
        break;
      default:
        break;
    }
    return fail(
        model_storage::AllocationFailure("Newton attempt allocation", d),
        NonlinearFailure::kNone);
  }
}
}  // namespace fesa
