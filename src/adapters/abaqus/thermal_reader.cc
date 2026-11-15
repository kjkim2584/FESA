#include "adapters/abaqus/thermal_reader.h"

#include <limits>
#include <map>
#include <string>
#include <vector>

#include "adapters/abaqus/amplitude_reader.h"
#include "adapters/abaqus/input_syntax.h"

namespace fesa {
namespace abaqus {
const char* ThermalAllowedOptions(const char* k) {
  if (Equal(k, "HEADING") || Equal(k, "NODE") || Equal(k, "ENDPART") ||
      Equal(k, "ENDASSEMBLY") || Equal(k, "ENDINSTANCE") ||
      Equal(k, "ENDSTEP") || Equal(k, "DENSITY") || Equal(k, "SPECIFICHEAT") ||
      Equal(k, "NODEOUTPUT"))
    return "";
  if (Equal(k, "PREPRINT")) return "ECHO,MODEL,HISTORY,CONTACT";
  if (Equal(k, "PART") || Equal(k, "ASSEMBLY") || Equal(k, "MATERIAL"))
    return "NAME";
  if (Equal(k, "INSTANCE")) return "NAME,PART";
  if (Equal(k, "ELEMENT")) return "TYPE,ELSET";
  if (Equal(k, "NSET")) return "NSET,GENERATE,INSTANCE";
  if (Equal(k, "ELSET")) return "ELSET,GENERATE,INSTANCE";
  if (Equal(k, "SOLIDSECTION") || Equal(k, "SHELLSECTION"))
    return "ELSET,MATERIAL";
  if (Equal(k, "CONDUCTIVITY")) return "TYPE";
  if (Equal(k, "INITIALCONDITIONS")) return "TYPE,SECTIONSPECIFICATION";
  if (Equal(k, "STEP")) return "NAME,INC,NLGEOM,AMPLITUDE";
  if (Equal(k, "HEATTRANSFER")) return "STEADYSTATE,DELTMX,END";
  if (Equal(k, "BOUNDARY") || Equal(k, "CFLUX") || Equal(k, "DFLUX") ||
      Equal(k, "DSFLUX"))
    return "OP,AMPLITUDE";
  if (Equal(k, "FILM") || Equal(k, "SFILM"))
    return "OP,FILMAMPLITUDE,SINKAMPLITUDE";
  if (Equal(k, "RADIATE") || Equal(k, "SRADIATE")) return "OP,AMBIENTAMPLITUDE";
  if (Equal(k, "FILMPROPERTY")) return "NAME,DEPENDENCIES";
  if (Equal(k, "AMPLITUDE")) return "NAME,DEFINITION,TIME,VALUE";
  if (Equal(k, "CONTROLS")) return "RESET,PARAMETERS,FIELD,CONVERGENCECHECK";
  if (Equal(k, "PHYSICALCONSTANTS")) return "ABSOLUTEZERO,STEFANBOLTZMANN";
  if (Equal(k, "SURFACE")) return "NAME,TYPE";
  if (Equal(k, "OUTPUT")) return "FIELD,HISTORY,VARIABLE";
  if (Equal(k, "ELEMENTOUTPUT")) return "DIRECTIONS";
  if (Equal(k, "RESTART")) return "WRITE,FREQUENCY";
  return nullptr;
}
namespace {
bool Named(const InputBlock* b, const char* key) {
  return b->HasOption(key) && Name(b->Option(key));
}
bool Option(const InputBlock* b, const char* key, const char* expected) {
  return !b->HasOption(key) || Equal(b->Option(key), expected);
}
bool OneOf(const char* value, const char* list) {
  std::vector<std::string> choices(0);
  Fields(list, &choices);
  for (const auto& item : choices)
    if (Equal(value, item.c_str())) return true;
  return false;
}
StatusCode Rows(const char* source, InputBlock* b, Diagnostic* d) {
  const auto& k = *b->keyword;
  const auto count = b->rows->size();
  if (k == "HEADING") return StatusCode::kOk;
  if (k == "NODE" || k == "ELEMENT" || k == "NSET" || k == "ELSET")
    return ValidateInputMeshRows(source, b, d);
  if (k == "AMPLITUDE") {
    int points = 0;
    return ValidateAmplitudeRows(source, b, &points, d);
  }
  const bool data =
      k == "CONDUCTIVITY" || k == "SPECIFICHEAT" || k == "DENSITY" ||
      k == "SOLIDSECTION" || k == "SHELLSECTION" || k == "INITIALCONDITIONS" ||
      k == "BOUNDARY" || k == "CFLUX" || k == "DFLUX" || k == "DSFLUX" ||
      k == "FILM" || k == "SFILM" || k == "RADIATE" || k == "SRADIATE" ||
      k == "FILMPROPERTY" || k == "HEATTRANSFER" || k == "CONTROLS" ||
      k == "SURFACE" || k == "NODEOUTPUT" || k == "ELEMENTOUTPUT";
  if (!data && count)
    return Error(source, b->rows->front()->line, b->step,
                 "Keyword accepts no data", d);
  if (k == "DENSITY" && count != 1)
    return Error(source, b->line, b->step,
                 "M1 requires one constant property row", d,
                 StatusCode::kUnsupported);
  if ((k == "BOUNDARY" || k == "CFLUX" || k == "DFLUX" || k == "DSFLUX" ||
       k == "FILM" || k == "SFILM" || k == "RADIATE" || k == "SRADIATE") &&
      count == 0 && !Equal(b->Option("OP"), "NEW"))
    return Error(source, b->line, b->step, "Empty condition requires NEW", d);
  if ((k == "INITIALCONDITIONS" || k == "SURFACE" || k == "NODEOUTPUT" ||
       k == "ELEMENTOUTPUT") &&
      count == 0)
    return Error(source, b->line, b->step, "Required data missing", d);
  if (k == "SHELLSECTION" && count != 1)
    return Error(source, b->line, b->step,
                 "Shell section requires one layer row", d);
  if (k == "SOLIDSECTION" && count > 1)
    return Error(source, b->line, b->step, "Invalid solid section data", d);
  if (k == "CONTROLS") {
    if (b->HasOption("RESET")) {
      if (b->options->size() != 1 || count)
        return Error(source, b->line, b->step, "Invalid CONTROLS RESET", d);
      return StatusCode::kOk;
    }
    if (Equal(b->Option("PARAMETERS"), "TIMEINCREMENTATION")) {
      if (count || b->HasOption("FIELD"))
        return Error(source, b->line, b->step,
                     "Thermal time controls use defaults only", d,
                     StatusCode::kUnsupported);
      return StatusCode::kOk;
    }
    if (!Equal(b->Option("PARAMETERS"), "FIELD") || count < 1 || count > 2)
      return Error(source, b->line, b->step, "Invalid field controls", d);
  }
  if ((k == "CONDUCTIVITY" || k == "SPECIFICHEAT" || k == "FILMPROPERTY") &&
      count == 0)
    return Error(source, b->line, b->step, "Thermal property data missing", d);
  double previous_temperature = 0;
  bool has_previous_temperature = false;
  for (std::size_t r = 0; r < count; ++r) {
    const auto* row = b->rows->at(r);
    std::vector<std::string> f(0);
    Fields(row->text->c_str(), &f);
    double v = 0;
    int first = 11, last = 11;
    bool valid = true;
    StatusCode code = StatusCode::kInputError;
    if (k == "DENSITY") {
      valid = f.size() == 1 && Real(f[0].c_str(), &v) && v > 0;
      if (f.size() > 1) code = StatusCode::kUnsupported;
    } else if (k == "CONDUCTIVITY" || k == "SPECIFICHEAT" ||
               k == "FILMPROPERTY") {
      double temperature = 0;
      valid = !f.empty() && f.size() <= 2 && Real(f[0].c_str(), &v) &&
              (k == "FILMPROPERTY" ? v >= 0 : v > 0) &&
              (f.size() == 2 ? Real(f[1].c_str(), &temperature) : count == 1);
      if (valid && has_previous_temperature)
        valid = temperature > previous_temperature;
      if (f.size() > 2) code = StatusCode::kUnsupported;
      previous_temperature = temperature;
      has_previous_temperature = true;
    } else if (k == "SOLIDSECTION") {
      valid = f.size() <= 2 &&
              (f[0].empty() || (Real(f[0].c_str(), &v) && v == 1)) &&
              (f.size() < 2 || f[1].empty());
    } else if (k == "SHELLSECTION") {
      int ns = 5;
      valid = f.size() >= 1 && f.size() <= 2 && Real(f[0].c_str(), &v) &&
              v > 0 &&
              (f.size() < 2 || f[1].empty() ||
               (Integer(f[1].c_str(), 1, 19, &ns) && ns % 2 == 1));
    } else if (k == "INITIALCONDITIONS") {
      const bool uniform = Equal(b->Option("SECTIONSPECIFICATION"), "UNIFORM");
      const bool numeric_first = Real(f[0].c_str(), &v);
      valid = uniform ? f.size() == 2 && Name(f[0].c_str(), true)
                      : !f.empty() && f.size() <= 8 &&
                            (numeric_first ||
                             (f.size() >= 2 && Name(f[0].c_str(), true)));
      for (std::size_t i = 1; valid && i < f.size(); ++i)
        valid = Real(f[i].c_str(), &v);
    } else if (k == "BOUNDARY" || k == "CFLUX") {
      valid = f.size() >= 2 && f.size() <= 4 && Name(f[0].c_str(), true);
      if (k == "CFLUX") {
        valid = valid && f.size() == 3 &&
                (f[1].empty() || Equal(f[1].c_str(), "0") ||
                 Integer(f[1].c_str(), 11, 30, &first)) &&
                Real(f[2].c_str(), &v);
      } else {
        valid = valid && Integer(f[1].c_str(), 11, 30, &first) &&
                (f.size() < 3 || f[2].empty() ||
                 Integer(f[2].c_str(), first, 30, &last)) &&
                (f.size() < 4 || f[3].empty() || Real(f[3].c_str(), &v));
        if (valid && b->step == 0 && v != 0)
          return Error(source, row->line, 0,
                       "Model BOUNDARY magnitude must be zero", d);
      }
    } else if (k == "RADIATE" || k == "SRADIATE") {
      double emissivity = 0;
      valid = f.size() == 4 && Name(f[0].c_str(), true) &&
              Real(f[2].c_str(), &v) && Real(f[3].c_str(), &emissivity) &&
              emissivity >= 0 && emissivity <= 1 &&
              OneOf(f[1].c_str(),
                    k == "RADIATE" ? "R1,R2,R3,R4,R5,R6,RNEG,RPOS" : "R");
    } else if (k == "FILM" || k == "SFILM") {
      double coefficient = 0;
      valid = f.size() == 4 && Name(f[0].c_str(), true) &&
              Real(f[2].c_str(), &v) &&
              (Real(f[3].c_str(), &coefficient) ? coefficient >= 0
                                                : Name(f[3].c_str())) &&
              OneOf(f[1].c_str(),
                    k == "FILM" ? "F1,F2,F3,F4,F5,F6,FNEG,FPOS" : "F");
    } else if (k == "DFLUX" || k == "DSFLUX") {
      valid = f.size() == 3 && Name(f[0].c_str(), true) &&
              Real(f[2].c_str(), &v) &&
              OneOf(f[1].c_str(),
                    k == "DFLUX" ? "BF,S1,S2,S3,S4,S5,S6,SNEG,SPOS" : "S");
      code = StatusCode::kUnsupported;
    } else if (k == "SURFACE") {
      valid = f.size() == 2 && Name(f[0].c_str(), true) && !f[1].empty();
    } else if (k == "NODEOUTPUT" || k == "ELEMENTOUTPUT") {
      for (const auto& item : f)
        valid =
            valid && OneOf(item.c_str(),
                           k == "NODEOUTPUT" ? "NT,RFL,CFL" : "HFL,TEMP,GRADT");
      code = StatusCode::kUnsupported;
    } else if (k == "CONTROLS") {
      valid = f.size() <= (r == 0 ? 8u : 2u);
      for (std::size_t i = 0; valid && i < f.size(); ++i)
        valid = f[i].empty() || (Real(f[i].c_str(), &v) && v >= 0 &&
                                 !(r == 0 && (i == 2 || i == 3) && v <= 0));
      if (f.size() > (r == 0 ? 8u : 2u)) code = StatusCode::kUnsupported;
    }
    if (!valid)
      return Error(source, row->line, b->step, "Invalid thermal keyword data",
                   d, code);
  }
  return StatusCode::kOk;
}
}  // namespace
StatusCode ReadThermalStep(const char* source, const InputBlock* step,
                           const InputBlock* heat,
                           HeatTransferStepParameters* out, Diagnostic* d) {
  if (!step || !heat || !out)
    return Error(source, 0, 0, "Thermal step reader input", d);
  HeatTransferStepParameters p;
  p.steady = heat->HasOption("STEADYSTATE");
  if (heat->rows->size() != 1)
    return Error(source, heat->line, heat->step, "Heat step needs one time row",
                 d);
  std::vector<std::string> f(0);
  Fields(heat->rows->front()->text->c_str(), &f);
  bool valid = f.size() >= 2 && f.size() <= 4 &&
               Real(f[0].c_str(), &p.initial_increment) &&
               Real(f[1].c_str(), &p.period);
  p.has_minimum = f.size() > 2 && !f[2].empty();
  p.has_maximum = f.size() > 3 && !f[3].empty();
  if (p.has_minimum) valid = valid && Real(f[2].c_str(), &p.minimum);
  if (p.has_maximum) valid = valid && Real(f[3].c_str(), &p.maximum);
  p.has_deltmx = heat->HasOption("DELTMX");
  if (p.has_deltmx) valid = valid && Real(heat->Option("DELTMX"), &p.deltmx);
  if (step->HasOption("INC"))
    valid =
        valid && Integer(step->Option("INC"), 1,
                         std::numeric_limits<int>::max(), &p.max_increments);
  p.has_default_amplitude = step->HasOption("AMPLITUDE");
  p.default_amplitude = Equal(step->Option("AMPLITUDE"), "RAMP")
                            ? ThermalStepAmplitude::kRamp
                            : ThermalStepAmplitude::kStep;
  if (p.steady && p.has_deltmx)
    return Error(source, heat->line, heat->step, "Steady DELTMX unsupported", d,
                 StatusCode::kUnsupported);
  if (!valid ||
      ValidateHeatTransferStepParameters(&p, nullptr) != StatusCode::kOk)
    return Error(source, heat->rows->front()->line, heat->step,
                 "Invalid thermal time settings", d);
  out->steady = p.steady;
  out->initial_increment = p.initial_increment;
  out->period = p.period;
  out->has_minimum = p.has_minimum;
  out->minimum = p.minimum;
  out->has_maximum = p.has_maximum;
  out->maximum = p.maximum;
  out->has_deltmx = p.has_deltmx;
  out->deltmx = p.deltmx;
  out->has_default_amplitude = p.has_default_amplitude;
  out->default_amplitude = p.default_amplitude;
  out->max_increments = p.max_increments;
  return StatusCode::kOk;
}
StatusCode ValidateThermalDeck(const char* source,
                               std::vector<InputBlock*>* blocks,
                               Diagnostic* d) {
  bool part = false, assembly = false, instance = false, active = false;
  bool had_part = false, had_assembly = false, had_instance = false;
  bool field = false, material = false;
  int step_number = 0, procedures = 0;
  std::string part_name(""), assembly_name("");
  const InputBlock* step = nullptr;
  std::map<std::string, int> amplitude_names;
  for (auto* b : *blocks) {
    const auto& k = *b->keyword;
    b->step = active ? step_number : 0;
    b->scope_kind = part       ? InputScope::kPart
                    : assembly ? InputScope::kAssembly
                    : active   ? InputScope::kStep
                               : InputScope::kModel;
    *b->scope = part ? part_name : assembly ? assembly_name : "";
    bool scope = true, options = true;
    if (k == "PART") {
      scope = !part && !assembly && !active && !had_part && step_number == 0;
      options = Named(b, "NAME");
      part = had_part = true;
      part_name = b->Option("NAME");
      b->scope_kind = InputScope::kPart;
      *b->scope = part_name;
    } else if (k == "ENDPART") {
      scope = part;
      part = false;
    } else if (k == "ASSEMBLY") {
      scope = had_part && !part && !assembly && !active && !had_assembly;
      options = Named(b, "NAME");
      assembly = had_assembly = true;
      assembly_name = b->Option("NAME");
      b->scope_kind = InputScope::kAssembly;
      *b->scope = assembly_name;
    } else if (k == "ENDASSEMBLY") {
      scope = assembly && !instance;
      assembly = false;
    } else if (k == "INSTANCE") {
      scope = assembly && !instance && !had_instance;
      options = Named(b, "NAME") && Named(b, "PART");
      instance = had_instance = true;
    } else if (k == "ENDINSTANCE") {
      scope = instance;
      instance = false;
    } else if (k == "NODE" || k == "ELEMENT" || k == "SOLIDSECTION" ||
               k == "SHELLSECTION") {
      scope = part;
      if (k == "ELEMENT")
        options = OneOf(b->Option("TYPE"), "DC3D4,DC3D8,DS4") &&
                  (!b->HasOption("ELSET") || Named(b, "ELSET"));
      if (k == "SOLIDSECTION" || k == "SHELLSECTION")
        options = Named(b, "ELSET") && Named(b, "MATERIAL");
    } else if (k == "NSET" || k == "ELSET") {
      scope = part || (assembly && !instance);
      options = Named(b, k.c_str()) && (!b->HasOption("INSTANCE") ||
                                        (assembly && Named(b, "INSTANCE")));
    } else if (k == "SURFACE") {
      scope = part || (assembly && !instance);
      options = Named(b, "NAME") && Equal(b->Option("TYPE"), "ELEMENT");
    } else if (k == "MATERIAL") {
      scope = !part && !assembly && !active && step_number == 0;
      options = Named(b, "NAME");
      material = true;
    } else if (k == "CONDUCTIVITY" || k == "SPECIFICHEAT" || k == "DENSITY") {
      scope = material;
      options = k != "CONDUCTIVITY" || Option(b, "TYPE", "ISOTROPIC");
    } else if (k == "AMPLITUDE" || k == "INITIALCONDITIONS" ||
               k == "PHYSICALCONSTANTS") {
      scope = !part && !assembly && !active && step_number == 0;
      if (k == "AMPLITUDE") {
        options = Named(b, "NAME") && Option(b, "DEFINITION", "TABULAR") &&
                  Option(b, "VALUE", "RELATIVE") &&
                  (!b->HasOption("TIME") ||
                   OneOf(b->Option("TIME"), "STEPTIME,TOTALTIME"));
        std::string name("");
        Canonical(b->Option("NAME"), &name);
        if (!amplitude_names.emplace(name, b->line).second)
          return Error(source, b->line, 0, "Duplicate amplitude name", d);
      } else if (k == "INITIALCONDITIONS") {
        options = Equal(b->Option("TYPE"), "TEMPERATURE") &&
                  Option(b, "SECTIONSPECIFICATION", "UNIFORM");
      }
    } else if (k == "FILMPROPERTY") {
      scope = !part && !assembly && !active && step_number == 0;
      options = Named(b, "NAME") && Option(b, "DEPENDENCIES", "0");
    } else if (k == "STEP") {
      scope = had_instance && !part && !assembly && !active;
      options = Option(b, "NLGEOM", "NO") &&
                (!b->HasOption("NAME") || Named(b, "NAME")) &&
                (!b->HasOption("AMPLITUDE") ||
                 OneOf(b->Option("AMPLITUDE"), "STEP,RAMP"));
      active = true;
      b->step = ++step_number;
      step = b;
      b->scope_kind = InputScope::kStep;
      procedures = 0;
      field = false;
    } else if (k == "ENDSTEP") {
      scope = active && procedures == 1;
      active = false;
    } else if (k == "BOUNDARY" || k == "CFLUX" || k == "DFLUX" ||
               k == "DSFLUX") {
      scope = !part && !assembly &&
              (active || (k == "BOUNDARY" && step_number == 0));
      options = (!b->HasOption("OP") || OneOf(b->Option("OP"), "MOD,NEW")) &&
                (!b->HasOption("AMPLITUDE") || Named(b, "AMPLITUDE")) &&
                (active || b->options->empty());
    } else if (k == "RADIATE" || k == "SRADIATE") {
      scope = active && !part && !assembly;
      options =
          (!b->HasOption("OP") || OneOf(b->Option("OP"), "MOD,NEW")) &&
          (!b->HasOption("AMBIENTAMPLITUDE") || Named(b, "AMBIENTAMPLITUDE"));
    } else if (k == "FILM" || k == "SFILM") {
      scope = active && !part && !assembly;
      options = (!b->HasOption("OP") || OneOf(b->Option("OP"), "MOD,NEW")) &&
                (!b->HasOption("FILMAMPLITUDE") || Named(b, "FILMAMPLITUDE")) &&
                (!b->HasOption("SINKAMPLITUDE") || Named(b, "SINKAMPLITUDE"));
    } else if (k == "PREPRINT") {
      scope = !part && !assembly && !active && step_number == 0;
      for (const auto& option : *b->options)
        options = options && Equal(option.second.c_str(), "NO");
    } else if (k == "HEADING") {
      scope = b == blocks->front();
    } else {
      scope = active;
      if (k == "HEATTRANSFER") {
        ++procedures;
        options = Option(b, "END", "PERIOD");
        HeatTransferStepParameters settings;
        if (scope && options) {
          auto status = ReadThermalStep(source, step, b, &settings, d);
          if (status != StatusCode::kOk) return status;
        }
      } else if (k == "CONTROLS") {
        options = Option(b, "CONVERGENCECHECK", "STRICT") &&
                  (!b->HasOption("FIELD") ||
                   OneOf(b->Option("FIELD"), "GLOBAL,TEMPERATURE"));
      } else if (k == "RESTART") {
        options = b->HasOption("WRITE") && Equal(b->Option("FREQUENCY"), "0");
      } else if (k == "OUTPUT") {
        options = b->HasOption("FIELD") != b->HasOption("HISTORY");
        options = options && Option(b, "VARIABLE", "PRESELECT");
        field = b->HasOption("FIELD");
      } else if (k == "NODEOUTPUT" || k == "ELEMENTOUTPUT") {
        scope = scope && field;
        options = Option(b, "DIRECTIONS", "YES");
      }
    }
    if (k != "MATERIAL" && k != "CONDUCTIVITY" && k != "DENSITY" &&
        k != "SPECIFICHEAT")
      material = false;
    if (!scope)
      return Error(source, b->line, b->step,
                   "Invalid thermal keyword scope/order", d);
    if (!options)
      return Error(source, b->line, b->step, "Unsupported thermal option", d,
                   StatusCode::kUnsupported);
    b->compatibility_warning = k == "SURFACE" || k == "PHYSICALCONSTANTS" ||
                               k == "OUTPUT" || k == "NODEOUTPUT" ||
                               k == "ELEMENTOUTPUT" || k == "RESTART";
    const auto status = Rows(source, b, d);
    if (status != StatusCode::kOk) return status;
  }
  if (part || assembly || instance || active || !had_instance ||
      step_number == 0)
    return Error(source, blocks->empty() ? 1 : blocks->back()->line,
                 step_number, "Incomplete thermal Part/Assembly/Instance/Step",
                 d);
  return StatusCode::kOk;
}
}  // namespace abaqus
}  // namespace fesa
