#include "fesa/analysis/newton_report.h"

#include <memory>
#include <new>

#include "base/diagnostic_copy.h"
#include "model/model_storage.h"
namespace fesa {
namespace {
void ResetField(FieldResidual* f) {
  f->dof_count = f->contribution_count = f->maximum_node = f->maximum_dof = 0;
  f->maximum_element = 0;
  f->maximum_internal_end = -1;
  f->maximum = f->l2 = f->rms = f->spatial_mean = f->raw_mean = f->used_mean =
      0;
  f->maximum_flux = f->phi = f->increment_max = f->correction_max = 0;
  f->residual_limit = f->correction_limit = 0;
  f->converted = false;
  f->passed = true;
  f->branch = ResidualBranch::kNoDofs;
}
}  // namespace
StatusCode NewtonReport::Create(NewtonReport** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("NewtonReport output", d);
  *out = nullptr;
  std::unique_ptr<NewtonReport> report(new (std::nothrow) NewtonReport);
  if (!report) return model_storage::AllocationFailure("NewtonReport", d);
  report->info_ = new (std::nothrow) NewtonIterationInfo;
  report->force_ = new (std::nothrow) FieldResidual;
  report->moment_ = new (std::nothrow) FieldResidual;
  report->logging_diagnostic_ = new (std::nothrow) Diagnostic;
  if (!report->info_ || !report->force_ || !report->moment_ ||
      !report->logging_diagnostic_)
    return model_storage::AllocationFailure("NewtonReport fields", d);
  *out = report.release();
  return StatusCode::kOk;
}
NewtonReport::~NewtonReport() {
  delete info_;
  delete force_;
  delete moment_;
  delete logging_diagnostic_;
}
void NewtonReport::Reset() noexcept {
  Diagnostic empty;
  CopyDiagnosticMetadata(&empty, logging_diagnostic_);
  auto* i = info_;
  i->status = i->numerical_status = StatusCode::kOk;
  i->outcome = NewtonOutcome::kNotStarted;
  i->stage = NewtonStage::kIdle;
  i->iteration = i->attempted_iteration = i->linear_solver_calls = 0;
  i->residual_evaluations = 0;
  i->residual_valid = i->correction_valid = i->logging_failed = false;
  i->residual_seconds = i->tangent_seconds = i->solve_seconds = 0;
  i->update_seconds = i->criteria_seconds = 0;
  for (int field = 0; field < 2; ++field) {
    i->estimate_available[field] = i->prediction_available[field] = false;
    i->estimate_max[field] = i->predicted_total[field] = 0;
  }
  ResetField(force_);
  ResetField(moment_);
}
const NewtonIterationInfo* NewtonReport::Info() const { return info_; }
const FieldResidual* NewtonReport::Displacement() const { return force_; }
const FieldResidual* NewtonReport::Rotation() const { return moment_; }
const Diagnostic* NewtonReport::LoggingDiagnostic() const {
  return info_->logging_failed ? logging_diagnostic_ : nullptr;
}
void NewtonReport::RecordLoggingFailure(StatusCode status,
                                        const Diagnostic* diagnostic) noexcept {
  if (status == StatusCode::kOk || info_->logging_failed) return;
  CopyDiagnosticMetadata(diagnostic, logging_diagnostic_);
  logging_diagnostic_->code = status;
  info_->logging_failed = true;
}
}  // namespace fesa
