#include "fesa/io/control_settings.h"

#include <cstdio>
#include <memory>
#include <new>
#include <string>

#include "adapters/abaqus/control_reader.h"
#include "adapters/abaqus/input_syntax.h"
#include "io/file_access.h"

namespace fesa {
ControlSettings::~ControlSettings() {
  delete newton;
  delete increments;
  delete displacement;
  delete rotation;
  delete temperature;
  delete source_path;
  delete config_path;
}
namespace {
// q0/qu are deliberately omitted: writing a numeric q0 would replace the
// previous-committed-mean policy with a constant after the first step.
constexpr char kDefaultConfig[] =
    "# FESA convergence controls; omitted keys use Abaqus defaults.\n"
    "# Nonlinear Solid defaults to incremental_elastic; StVK is deprecated.\n"
    "# Ratios, Cf and epsilon thresholds are dimensionless.\n"
    "# Rn: general residual; Cn: correction (Cn=1 disables that check).\n"
    "# Rp: parsed/logged only; unused by one direct solve.\n"
    "# epsilon: zero-flux residual/classification; Cepsilon: correction.\n"
    "# Rl: linear branch; Cf: field conversion (0 disables conversion).\n"
    "# epsilon_l: active contribution threshold (0 includes all).\n"
    "# epsilon_d: displacement increment / characteristic length.\n"
    "# q0: automatic .01 in Step 1, previous committed raw mean thereafter.\n"
    "# qu: automatic average. Numeric q0/qu must be positive; no auto syntax.\n"
    "# global.q0/global.qu apply the same number in force AND moment units.\n"
    "# Use displacement.q0/qu or rotation.q0/qu for field-specific units.\n"
    "mass = default\n"
    "global.Rn = .005\nglobal.Cn = .01\nglobal.Rp = .02\n"
    "global.epsilon = 1e-5\nglobal.Cepsilon = .001\nglobal.Rl = 1e-8\n"
    "global.Cf = 1\nglobal.epsilon_l = 1e-5\nglobal.epsilon_d = 1e-8\n";
StatusCode ConfigError(StatusCode code, const Diagnostic* context, int line,
                       const char* operation, Diagnostic* d, int first = 0) {
  if (d) {
    d->source_file = context->source_file;
    d->duplicate_boundary = false;
    d->condition_dof = 0;
    d->path_kind = context->path_kind;
    d->source_line = line;
    d->first_source_line = first;
    d->step = 0;
    d->Record(code,
              code == StatusCode::kAllocationFailure
                  ? DiagnosticStage::kAllocation
                  : DiagnosticStage::kValidation,
              operation, context->vendor_code);
  }
  return code;
}
StatusCode ParseConfig(FILE* file, ControlSettings* settings,
                       Diagnostic* context, Diagnostic* d) {
  double values[3][abaqus::kControlCount] = {};
  int lines[3][abaqus::kControlCount] = {};
  std::string raw(""), text(""), key(""), number("");
  int line = 0;
  while (io::ReadLine(file, &raw)) {
    ++line;
    context->source_line = line;
    if (raw.find('\0') != std::string::npos)
      return ConfigError(StatusCode::kInputError, context, line,
                         "Embedded NUL in cfg", d);
    if (line == 1 && raw.compare(0, 3, "\xEF\xBB\xBF") == 0) raw.erase(0, 3);
    abaqus::Trim(raw.c_str(), &text);
    if (text.empty() || text[0] == '#') continue;
    const auto equal = text.find('=');
    if (equal == std::string::npos)
      return ConfigError(StatusCode::kInputError, context, line,
                         "cfg requires key = value", d);
    abaqus::Trim(text.substr(0, equal).c_str(), &key);
    abaqus::Trim(text.substr(equal + 1).c_str(), &number);
    if (key == "solid.constitutive_model") {
      if (settings->solid_constitutive_model_line)
        return ConfigError(StatusCode::kInputError, context, line,
            "Duplicate cfg key", d, settings->solid_constitutive_model_line);
      if (number == "stvk")
        return ConfigError(StatusCode::kUnsupported, context, line,
            "StVK nonlinear Solid is deprecated; use incremental_elastic", d);
      else if (number == "incremental_elastic")
        settings->solid_constitutive_model =
            SolidConstitutiveModel::kIncrementalElastic;
      else
        return ConfigError(StatusCode::kInputError, context, line,
                           "Invalid cfg solid constitutive model", d);
      settings->solid_constitutive_model_line = line;
      continue;
    }
    if (key == "mass") {
      if (settings->mass_line)
        return ConfigError(StatusCode::kInputError, context, line,
                           "Duplicate cfg key", d, settings->mass_line);
      if (number == "default") settings->mass = MassSelection::kDefault;
      else if (number == "consistent") settings->mass = MassSelection::kConsistent;
      else if (number == "lumped") settings->mass = MassSelection::kLumped;
      else
        return ConfigError(StatusCode::kInputError, context, line,
                           "Invalid cfg mass value", d);
      settings->mass_line = line;
      continue;
    }
    const auto dot = key.find('.');
    const std::string field(key.substr(0, dot));
    const int group = field == "global"         ? 0
                      : field == "displacement" ? 1
                      : field == "rotation"     ? 2
                                                : -1;
    int index = 0;
    for (; index < abaqus::kControlCount; ++index)
      if (dot != std::string::npos &&
          key.compare(dot + 1, std::string::npos,
                      abaqus::ControlSymbol(index)) == 0)
        break;
    if (group < 0 || index == abaqus::kControlCount ||
        (group == 2 && index == 10))
      return ConfigError(StatusCode::kInputError, context, line,
                         "Unknown cfg key", d);
    if (lines[group][index])
      return ConfigError(StatusCode::kInputError, context, line,
                         "Duplicate cfg key", d, lines[group][index]);
    double value = 0;
    if (!abaqus::Real(number.c_str(), &value) || value < 0 ||
        ((index == 2 || index == 3) && value == 0))
      return ConfigError(StatusCode::kInputError, context, line,
                         "Invalid cfg value", d);
    lines[group][index] = line;
    values[group][index] = value;
  }
  if (std::ferror(file))
    return ConfigError(StatusCode::kIoFailure, context, line, "Read cfg", d);
  // All cfg entries have been validated, including values overridden by INP.
  if (settings->source == ControlSource::kAbaqus) return StatusCode::kOk;
  // Apply common settings before field overrides, independent of source lines.
  for (int group = 0; group < 3; ++group) {
    for (int index = 0; index < abaqus::kControlCount; ++index) {
      if (!lines[group][index]) continue;
      for (int field = 1; field <= 2; ++field) {
        if ((group != 0 && group != field) || (field == 2 && index == 10))
          continue;
        FieldControls* target =
            field == 1 ? settings->displacement : settings->rotation;
        *abaqus::ControlValue(target, index) = values[group][index];
        if (index == 2) target->q0_explicit = true;
        if (index == 3) target->qu_explicit = true;
      }
    }
  }
  return StatusCode::kOk;
}
StatusCode LoadConfig(const char* path, ControlSettings* settings,
                      Diagnostic* context, Diagnostic* d) {
  bool exists = false, regular = false;
  const auto inspected = io::InspectPath(path, &exists, &regular, context);
  if (inspected != StatusCode::kOk)
    return ConfigError(inspected, context, 0, "Inspect cfg", d);
  FILE* raw = nullptr;
  const auto mode = exists ? io::FileMode::kRead : io::FileMode::kCreateNew;
  StatusCode status = io::OpenFile(path, mode, &raw, context);
  if (status != StatusCode::kOk)
    return ConfigError(status, context, 0, exists ? "Open cfg" : "Create cfg",
                       d);
  std::unique_ptr<FILE, decltype(&std::fclose)> file(raw, &std::fclose);
  if (exists) {
    status = ParseConfig(file.get(), settings, context, d);
  } else {
    settings->config_created = true;
    if (std::fputs(kDefaultConfig, file.get()) == EOF ||
        std::fflush(file.get()) != 0)
      status = ConfigError(StatusCode::kIoFailure, context, 0, "Write cfg", d);
  }
  const int closed = std::fclose(file.release());
  if (status == StatusCode::kOk && closed != 0)
    return ConfigError(StatusCode::kIoFailure, context, 0, "Close cfg", d);
  return status;
}
}  // namespace
StatusCode SelectInitialControls(const InputDeck* deck,
                                 const RunOptions* options,
                                 const char* executable_directory,
                                 ControlSettings** out, Diagnostic* d) {
  if (!out)
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation,
                         "SelectInitialControls output", d);
  *out = nullptr;
  if (!deck || !options ||
      (!options->config_path &&
       (!executable_directory || !*executable_directory)))
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation,
                         "SelectInitialControls input", d);
  Diagnostic context;
  const bool has_controls = deck->HasControls();
  context.source_file = options->config_path ? options->config_path->c_str()
                                              : executable_directory;
  context.path_kind = options->config_path
                          ? DiagnosticPathKind::kLiteralFile
                          : DiagnosticPathKind::kDefaultConfigDirectory;
  try {
    std::unique_ptr<ControlSettings> settings(new ControlSettings);
    settings->displacement = new FieldControls;
    settings->rotation = new FieldControls;
    settings->newton = new NewtonControls;
    settings->increments = new IncrementControls;
    settings->source = has_controls ? ControlSource::kAbaqus
                                    : ControlSource::kConfig;
    settings->config_path = new std::string("");
    if (options->config_path)
      settings->config_path->assign(*options->config_path);
    else
      io::DefaultConfigPath(executable_directory, settings->config_path);
    settings->source_path = new std::string(
        has_controls ? deck->SourceFile() : settings->config_path->c_str());
    const auto status = LoadConfig(settings->config_path->c_str(),
                                   settings.get(), &context, d);
    if (status != StatusCode::kOk) return status;
    *out = settings.release();
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return ConfigError(StatusCode::kAllocationFailure, &context,
                       context.source_line, "SelectInitialControls allocation",
                       d);
  }
}
StatusCode ApplyStepControls(const InputDeck* deck, int step_number,
                             ControlSettings* settings, Diagnostic* d) {
  if (!deck || !settings || !settings->displacement || !settings->rotation ||
      step_number < 1)
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation,
                         "ApplyStepControls input", d);
  if (settings->temperature)
    return abaqus::ReadThermalFieldControls(deck, step_number,
                                            settings->temperature, d);
  if (settings->source == ControlSource::kConfig) return StatusCode::kOk;
  try {
    FieldControls displacement, rotation;
    NewtonControls counts;
    IncrementControls increments;
    if (settings->newton) abaqus::CopyNewtonControls(settings->newton, &counts);
    if (settings->increments)
      abaqus::CopyIncrementControls(settings->increments, &increments);
    abaqus::CopyFieldControls(settings->displacement, &displacement);
    abaqus::CopyFieldControls(settings->rotation, &rotation);
    auto status = abaqus::ReadStepControls(deck, step_number, &displacement,
                                           &rotation, d);
    if (status != StatusCode::kOk) return status;
    status = abaqus::ReadStepTimeControls(
        deck, step_number, settings->newton ? &counts : nullptr,
        settings->increments ? &increments : nullptr, d);
    if (status != StatusCode::kOk) return status;
    if (settings->newton) abaqus::CopyNewtonControls(&counts, settings->newton);
    if (settings->increments)
      abaqus::CopyIncrementControls(&increments, settings->increments);
    abaqus::CopyFieldControls(&displacement, settings->displacement);
    abaqus::CopyFieldControls(&rotation, settings->rotation);
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return abaqus::Error(deck->SourceFile(), 0, step_number,
                         "ApplyStepControls allocation", d,
                         StatusCode::kAllocationFailure);
  }
}
StatusCode CloneControlSettings(const ControlSettings* source,
                                ControlSettings** out, Diagnostic* d) {
  if (!out)
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation,
                         "CloneControlSettings output", d);
  *out = nullptr;
  if (!source || !source->displacement || !source->rotation ||
      !source->source_path)
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation,
                         "CloneControlSettings input", d);
  try {
    std::unique_ptr<ControlSettings> copy(new ControlSettings);
    copy->displacement = new FieldControls;
    copy->rotation = new FieldControls;
    if (source->temperature) {
      copy->temperature = new FieldControls;
      abaqus::CopyFieldControls(source->temperature, copy->temperature);
    }
    if (source->newton) {
      copy->newton = new NewtonControls;
      abaqus::CopyNewtonControls(source->newton, copy->newton);
    }
    if (source->increments) {
      copy->increments = new IncrementControls;
      abaqus::CopyIncrementControls(source->increments, copy->increments);
    }
    copy->source_path = new std::string(*source->source_path);
    abaqus::CopyFieldControls(source->displacement, copy->displacement);
    abaqus::CopyFieldControls(source->rotation, copy->rotation);
    copy->source = source->source;
    copy->config_created = source->config_created;
    if (source->config_path)
      copy->config_path = new std::string(*source->config_path);
    copy->mass = source->mass;
    copy->mass_line = source->mass_line;
    copy->solid_constitutive_model = source->solid_constitutive_model;
    copy->solid_constitutive_model_line = source->solid_constitutive_model_line;
    *out = copy.release();
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return ReportFailure(StatusCode::kAllocationFailure,
                         DiagnosticStage::kAllocation,
                         "CloneControlSettings allocation", d);
  }
}
}  // namespace fesa
