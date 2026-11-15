#include "fesa/io/analysis_log.h"

#include <cerrno>
#include <cstdio>
#include <memory>
#include <new>
#include <string>

#include "fesa/io/control_settings.h"
#include "io/file_access.h"
#include "io/log_stdio.h"

namespace fesa {
namespace io {
namespace {
thread_local LogFault fault = LogFault::kNone;
thread_local int remaining = 0;
thread_local unsigned long long write_calls = 0;
}  // namespace
void FailLogAfter(LogFault point, int skip) {
  fault = point;
  remaining = skip;
}
unsigned long long LogWriteCallCount() { return write_calls; }
bool FailLogBoundary(LogFault point) {
  if (point == LogFault::kWrite) ++write_calls;
  if (fault != point) return false;
  if (remaining > 0) {
    --remaining;
    return false;
  }
  fault = LogFault::kNone;
  errno = EIO;
  return true;
}
}  // namespace io
namespace {
StatusCode LogFailure(const char* input, const char* operation, Diagnostic* d,
                      StatusCode code = StatusCode::kIoFailure,
                      int native = 0) {
  if (d) {
    d->source_file = input;
    d->duplicate_boundary = false;
    d->condition_dof = 0;
    d->path_kind = DiagnosticPathKind::kAnalysisInputFile;
    d->source_line = 0;
    d->first_source_line = 0;
    d->step = 0;
    d->node_id = d->element_id = 0;
    d->Record(code,
              code == StatusCode::kAllocationFailure
                  ? DiagnosticStage::kAllocation
                  : DiagnosticStage::kValidation,
              operation, native);
  }
  return code;
}
StatusCode FieldText(const FieldControls* field, const char* name,
                     AnalysisLog* log, Diagnostic* d) {
  char q0[100] = {}, qu[100] = {}, line[1024] = {};
  if (field->q0_explicit) {
    std::snprintf(q0, sizeof(q0), "%.17g (explicit)", field->q0);
  } else {
    std::snprintf(q0, sizeof(q0),
                  "automatic (Step1=.01; then previous committed raw mean)");
  }
  if (field->qu_explicit) {
    std::snprintf(qu, sizeof(qu), "%.17g (explicit)", field->qu);
  } else {
    std::snprintf(qu, sizeof(qu), "automatic");
  }
  std::snprintf(
      line, sizeof(line),
      "%s: Rn=%.17g Cn=%.17g q0=%s qu=%s Rp=%.17g "
      "epsilon=%.17g Cepsilon=%.17g Rl=%.17g Cf=%.17g epsilon_l=%.17g",
      name, field->rn, field->cn, q0, qu, field->rp, field->epsilon,
      field->cepsilon, field->rl, field->cf, field->epsilon_l);
  auto status = log->Write(line, d);
  if (status != StatusCode::kOk) return status;
  if (name[0] == 'd') {
    std::snprintf(line, sizeof(line), " epsilon_d=%.17g", field->epsilon_d);
    status = log->Write(line, d);
    if (status != StatusCode::kOk) return status;
  }
  return log->Write("\n", d);
}
}  // namespace
AnalysisLog::AnalysisLog(std::string* input_path)
    : stream_(nullptr), input_path_(input_path) {}
AnalysisLog::~AnalysisLog() {
  if (stream_) std::fclose(static_cast<FILE*>(stream_));
  delete input_path_;
}
StatusCode AnalysisLog::Write(const char* text, Diagnostic* d) {
  if (!text || !stream_)
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation,
                         "AnalysisLog.Write input", d);
  FILE* file = static_cast<FILE*>(stream_);
  if (io::FailLogBoundary(io::LogFault::kWrite) ||
      std::fputs(text, file) == EOF ||
      io::FailLogBoundary(io::LogFault::kFlush) || std::fflush(file) != 0)
    return LogFailure(input_path_->c_str(), "Write analysis log", d,
                      StatusCode::kIoFailure, errno);
  return StatusCode::kOk;
}
StatusCode AnalysisLog::WriteControls(const ControlSettings* settings, int step,
                                      Diagnostic* d) {
  if (!settings || !settings->source_path || !settings->displacement ||
      !settings->rotation)
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation, "WriteControls input",
                         d);
  if (settings->temperature) {
    char thermal[360] = {};
    std::snprintf(thermal, sizeof(thermal),
                  "thermal controls step=%d source=Abaqus defaults/INP; "
                  "mechanical cfg convergence unused\n",
                  step);
    auto status = Write(thermal, d);
    if (status != StatusCode::kOk) return status;
    if (settings->config_path) {
      for (const char* text :
           {"config validated=", settings->config_path->c_str(), "\n"}) {
        status = Write(text, d);
        if (status != StatusCode::kOk) return status;
      }
    }
    status = FieldText(settings->temperature, "temperature", this, d);
    if (status != StatusCode::kOk) return status;
    if (settings->newton) {
      const auto* n = settings->newton;
      std::snprintf(thermal, sizeof(thermal),
                    "thermal Newton I0=%d IR=%d IP=%d IC=%d ISc=%d; "
                    "dedicated heat increment controller, fixed coefficients\n",
                    n->divergence_start, n->slow_start, n->nonquadratic_start,
                    n->prediction_limit, n->iteration_limit);
      return Write(thermal, d);
    }
    return StatusCode::kOk;
  }
  char prefix[100] = {};
  std::snprintf(prefix, sizeof(prefix),
                "controls step=%d source=%s path=", step,
                settings->source == ControlSource::kAbaqus ? "Abaqus" : "cfg");
  for (const char* text : {static_cast<const char*>(prefix),
                           settings->source_path->c_str(), "\n"}) {
    const auto status = Write(text, d);
    if (status != StatusCode::kOk) return status;
  }
  if (settings->config_path) {
    const char* mass = settings->mass == MassSelection::kDefault ? "default"
                       : settings->mass == MassSelection::kConsistent
                           ? "consistent" : "lumped";
    for (const char* text : {"config validated=", settings->config_path->c_str(),
                             " mass requested=", mass, "\n"}) {
      const auto status = Write(text, d);
      if (status != StatusCode::kOk) return status;
    }
  }
  auto status = FieldText(settings->displacement, "displacement", this, d);
  if (status != StatusCode::kOk) return status;
  status = FieldText(settings->rotation, "rotation", this, d);
  if (status != StatusCode::kOk) return status;
  if (settings->newton && settings->increments) {
    const auto* n = settings->newton;
    const auto* i = settings->increments;
    char line[640] = {};
    std::snprintf(line, sizeof(line),
                  "Newton I0=%d IR=%d IP=%d IC=%d ISc=%d\n"
                  "Increment IA=%d IL=%d IG=%d Df=%.17g DC=%.17g DH=%.17g "
                  "DB=%.17g DD=%.17g DM=%.17g DE=%.17g DT=%.17g\n",
                  n->divergence_start, n->slow_start, n->nonquadratic_start,
                  n->prediction_limit, n->iteration_limit, i->attempt_limit,
                  i->slow_success_iterations, i->quick_iterations,
                  i->divergence_cutback, i->slow_cutback, i->geometry_cutback,
                  i->success_scale, i->growth, i->maximum_growth,
                  i->extrapolation_cutoff, i->endpoint_stretch);
    status = Write(line, d);
    if (status != StatusCode::kOk) return status;
    return Write("Rp applies at Newton iteration>=IP; unused by direct linear "
                 "solve.\n", d);
  }
  return Write("Rp is parsed but unused by the single direct solve.\n", d);
}
StatusCode AnalysisLog::Finish(Diagnostic* d) {
  if (!stream_) return StatusCode::kOk;
  FILE* file = static_cast<FILE*>(stream_);
  stream_ = nullptr;
  const int closed = std::fclose(file);
  if (io::FailLogBoundary(io::LogFault::kClose) || closed != 0)
    return LogFailure(input_path_->c_str(), "Close analysis log", d,
                      StatusCode::kIoFailure, errno);
  return StatusCode::kOk;
}
StatusCode OpenAnalysisLog(const char* input_path, AnalysisLog** out,
                           Diagnostic* d) {
  if (!out)
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation, "OpenAnalysisLog output",
                         d);
  *out = nullptr;
  if (!input_path || !*input_path)
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation, "OpenAnalysisLog input",
                         d);
  try {
    std::string path("");
    io::AnalysisLogPath(input_path, &path);
    std::unique_ptr<std::string> copied(new std::string(input_path));
    std::unique_ptr<AnalysisLog> log(new AnalysisLog(copied.get()));
    copied.release();
    FILE* file = nullptr;
    Diagnostic native;
    const auto status =
        io::OpenFile(path.c_str(), io::FileMode::kReplace, &file, &native);
    if (status != StatusCode::kOk)
      return LogFailure(input_path, "Open analysis log", d, status,
                        native.vendor_code);
    log->stream_ = file;
    *out = log.release();
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return LogFailure(input_path, "OpenAnalysisLog allocation", d,
                      StatusCode::kAllocationFailure);
  }
}
StatusCode FormatDiagnostic(const Diagnostic* d, std::string* text) {
  if (!d || !text) return StatusCode::kInvalidArgument;
  try {
    text->assign("FESA: ");
    text->append(d->operation ? d->operation : "runtime failure");
    if (d->source_file) {
      std::string path("");
      if (d->path_kind == DiagnosticPathKind::kDefaultConfigDirectory) {
        io::DefaultConfigPath(d->source_file, &path);
      } else if (d->path_kind == DiagnosticPathKind::kAnalysisInputFile) {
        io::AnalysisLogPath(d->source_file, &path);
      } else {
        path.assign(d->source_file);
      }
      text->append(" file=");
      text->append(path);
    }
    char context[256] = {};
    std::snprintf(context, sizeof(context),
                  " line=%d first_line=%d step=%d node=%d element=%d stage=%d "
                  "native=%d condition_dof=%d\n",
                  d->source_line, d->first_source_line, d->step, d->node_id,
                  d->element_id, static_cast<int>(d->stage), d->vendor_code,
                  d->condition_dof);
    text->append(context);
    if (d->geometry_point_available) {
      std::snprintf(context, sizeof(context),
                    "natural xi=%.17g eta=%.17g zeta=%.17g\n", d->natural_xi,
                    d->natural_eta, d->natural_zeta);
      text->append(context);
    }
    if (d->geometry_value_available) {
      const char* quantity =
          d->geometry_quantity == GeometryQuantity::kNormalizedCurrentJacobian
              ? "normalized_current_jacobian"
              : (d->geometry_quantity == GeometryQuantity::kDeformationJacobian
                     ? "det_f"
                     : "none");
      std::snprintf(context, sizeof(context),
                    "geometry quantity=%s value=%.17g limit=%.17g\n", quantity,
                    d->geometry_value, d->geometry_limit);
      text->append(context);
    }
    if (d->duplicate_boundary) {
      std::snprintf(context, sizeof(context),
                    "boundary dof=%d first_value=%.17g value=%.17g\n",
                    d->boundary_dof, d->first_boundary_value,
                    d->boundary_value);
      text->append(context);
    }
    if (d->active_flux_failure) {
      std::snprintf(context, sizeof(context),
                    "active flux field=%s epsilon_l=%.17g history_max=%.17g "
                    "threshold=%.17g current_max=%.17g count=%d active=%d\n",
                    d->flux_field == 0 ? "force" : "moment", d->flux_epsilon_l,
                    d->flux_history_max, d->flux_threshold, d->flux_current_max,
                    d->flux_count, d->active_flux_count);
      text->append(context);
    }
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return StatusCode::kAllocationFailure;
  }
}
}  // namespace fesa
