#include "adapters/abaqus/inp_reader.h"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <limits>
#ifdef _WIN32
#define NOMINMAX
#include <Windows.h>
#endif
#include <memory>
#include <new>

#include "adapters/abaqus/amplitude_reader.h"
#include "adapters/abaqus/control_reader.h"
#include "adapters/abaqus/input_syntax.h"
#include "adapters/abaqus/thermal_reader.h"

namespace fesa {
InputLine::~InputLine() { delete text; }
InputBlock::~InputBlock() {
  if (elements)
    for (auto* element : *elements) delete element;
  delete elements;
  if (rows)
    for (auto* row : *rows) delete row;
  delete rows;
  delete options;
  delete keyword;
  delete original;
  delete scope;
}
const char* InputBlock::Option(const char* key) const {
  const auto it = options->find(key);
  return it == options->end() ? "" : it->second.c_str();
}
bool InputBlock::HasOption(const char* key) const {
  return options->find(key) != options->end();
}
InputDeck::InputDeck(std::string* source, std::vector<InputBlock*>* blocks)
    : source_(source), blocks_(blocks) {}
InputDeck::~InputDeck() {
  for (auto* block : *blocks_) delete block;
  delete blocks_;
  delete source_;
}
bool InputDeck::HasControls() const {
  for (const auto* block : *blocks_)
    if (*block->keyword == "CONTROLS") return true;
  return false;
}
const char* InputDeck::SourceFile() const { return source_->c_str(); }
const std::vector<InputBlock*>* InputDeck::Blocks() const { return blocks_; }

namespace {
using abaqus::Canonical;
using abaqus::Equal;
using abaqus::Error;
using abaqus::Fields;
using abaqus::Integer;
using abaqus::Name;
using abaqus::Real;
using abaqus::Trim;

// Every accepted option is enumerated here. Unknown syntax is never discarded.
const char* AllowedOptions(const char* keyword) {
  if (Equal(keyword, "HEADING") || Equal(keyword, "NODE") ||
      Equal(keyword, "ENDPART") || Equal(keyword, "ENDINSTANCE") ||
      Equal(keyword, "ENDASSEMBLY") || Equal(keyword, "STATIC") ||
      Equal(keyword, "ENDSTEP") || Equal(keyword, "NODEOUTPUT") ||
      Equal(keyword, "DENSITY"))
    return "";
  if (Equal(keyword, "INITIALCONDITIONS")) return "TYPE";
  if (Equal(keyword, "AMPLITUDE")) return "NAME,DEFINITION,TIME,VALUE";
  if (Equal(keyword, "DYNAMIC"))
    return "DIRECT,APPLICATION,TIMEINTEGRATOR,INCREMENTATION,ALPHA,NOHAF";
  if (Equal(keyword, "PREPRINT")) return "ECHO,MODEL,HISTORY,CONTACT";
  if (Equal(keyword, "PART") || Equal(keyword, "ASSEMBLY") ||
      Equal(keyword, "MATERIAL"))
    return "NAME";
  if (Equal(keyword, "INSTANCE")) return "NAME,PART";
  if (Equal(keyword, "ELEMENT")) return "TYPE,ELSET";
  if (Equal(keyword, "ELASTIC")) return "TYPE";
  if (Equal(keyword, "DAMPING")) return "ALPHA,BETA";
  if (Equal(keyword, "SOLIDSECTION")) return "ELSET,MATERIAL";
  if (Equal(keyword, "SHELLSECTION"))
    return "ELSET,MATERIAL,SECTIONINTEGRATION";
  if (Equal(keyword, "NSET")) return "NSET,GENERATE,INSTANCE";
  if (Equal(keyword, "ELSET")) return "ELSET,GENERATE,INSTANCE";
  if (Equal(keyword, "BEAMGENERALSECTION")) return "ELSET,MATERIAL,SECTION";
  if (Equal(keyword, "BOUNDARY")) return "OP,TYPE";
  if (Equal(keyword, "CLOAD")) return "OP,AMPLITUDE";
  if (Equal(keyword, "DLOAD")) return "OP";
  if (Equal(keyword, "STEP"))
    return "NAME,NLGEOM,INC,EXTRAPOLATION,AMPLITUDE,CONVERTSDI";
  if (Equal(keyword, "CONTROLS"))
    return "RESET,PARAMETERS,FIELD,CONVERGENCECHECK";
  if (Equal(keyword, "RESTART")) return "WRITE,FREQUENCY";
  if (Equal(keyword, "OUTPUT")) return "FIELD,HISTORY,VARIABLE";
  if (Equal(keyword, "ELEMENTOUTPUT")) return "DIRECTIONS";
  if (Equal(keyword, "CONTACTOUTPUT")) return "VARIABLE";
  return nullptr;
}
bool Listed(const char* value, const char* list) {
  std::vector<std::string> fields(0);
  Fields(list, &fields);
  for (const auto& field : fields)
    if (Equal(value, field.c_str())) return true;
  return false;
}
int OptionLine(const char* text, std::size_t index, int line) {
  std::size_t field = 0;
  for (const char* p = text; *p; ++p) {
    if (*p == '\n') ++line;
    if (*p == ',') {
      ++field;
      continue;
    }
    if (field == index && *p != ' ' && *p != '\t' && *p != '\n') return line;
  }
  return line;
}
StatusCode Keyword(const char* source, const char* text, int line,
                   int dynamic_step, bool thermal,
                   InputBlock** out, Diagnostic* d) {
  std::unique_ptr<InputBlock> block(new InputBlock);
  block->line = line;
  block->original = new std::string(text);
  block->keyword = new std::string("");
  block->scope = new std::string("");
  block->options = new std::map<std::string, std::string>;
  block->rows = new std::vector<InputLine*>(0);
  std::vector<std::string> fields(0);
  std::map<std::string, int> option_lines;
  Fields(text + 1, &fields);
  Canonical(fields.front().c_str(), block->keyword);
  const bool step_keyword =
      *block->keyword == "DYNAMIC" || *block->keyword == "CLOAD";
  const int step = step_keyword ? dynamic_step : 0;
  const char* allowed = thermal
      ? abaqus::ThermalAllowedOptions(block->keyword->c_str())
      : AllowedOptions(block->keyword->c_str());
  if (!allowed)
    return Error(source, line, 0, "Unsupported keyword", d,
                 StatusCode::kUnsupported);
  for (std::size_t i = 1; i < fields.size(); ++i) {
    const auto equal = fields[i].find('=');
    std::string key(""), value("");
    Canonical(fields[i].substr(0, equal).c_str(), &key);
    if (equal != std::string::npos)
      Trim(fields[i].substr(equal + 1).c_str(), &value);
    const int option_line = OptionLine(text, i, line);
    if (!Listed(key.c_str(), allowed))
      return Error(source, step_keyword ? option_line : line,
                   step, "Unsupported keyword option", d,
                   StatusCode::kUnsupported);
    const auto inserted = option_lines.emplace(key, option_line);
    if (!inserted.second)
      return Error(source, option_line, step, "Duplicate keyword option", d,
                   StatusCode::kInputError, inserted.first->second);
    block->options->emplace(key, value);
    const bool flag =
        key == "GENERATE" || key == "RESET" || key == "WRITE" ||
        (thermal && key == "STEADYSTATE" &&
         *block->keyword == "HEATTRANSFER") ||
        ((key == "DIRECT" || key == "NOHAF") &&
         *block->keyword == "DYNAMIC") ||
        ((key == "FIELD" || key == "HISTORY") && *block->keyword == "OUTPUT");
    if ((flag && equal != std::string::npos) ||
        (!flag && ((equal == std::string::npos &&
                    !(key == "NLGEOM" && *block->keyword == "STEP")) ||
                   (equal != std::string::npos && value.empty()))))
      return Error(source, step_keyword ? option_line : line,
                   step, "Invalid option value or flag", d,
                   *block->keyword == "DYNAMIC" && key == "DIRECT" &&
                           !value.empty()
                       ? StatusCode::kUnsupported : StatusCode::kInputError);
  }
  *out = block.release();
  return StatusCode::kOk;
}
bool OptionIs(const InputBlock* b, const char* key, const char* value) {
  return !b->HasOption(key) || Equal(b->Option(key), value);
}
bool Named(const InputBlock* b, const char* key) {
  return b->HasOption(key) && Name(b->Option(key));
}
bool SetEndpoint(const char* text, int* id) {
  const std::string target(text);
  const auto dot = target.find('.');
  if (dot == std::string::npos) return Integer(text, 1, 999999999, id);
  return Name(target.substr(0, dot).c_str()) &&
         Integer(target.substr(dot + 1).c_str(), 1, 999999999, id);
}
StatusCode ParseElements(const char* source, InputBlock* block, Diagnostic* d) {
  block->elements = new std::vector<InputElement*>(0);
  const int nodes = Equal(block->Option("TYPE"), "B33") ? 2
                    : (Equal(block->Option("TYPE"), "C3D4") ||
                       Equal(block->Option("TYPE"), "DC3D4") ||
                       Equal(block->Option("TYPE"), "DS4") ||
                       Equal(block->Option("TYPE"), "S4"))
                        ? 4
                    : (Equal(block->Option("TYPE"), "C3D8") ||
                       Equal(block->Option("TYPE"), "DC3D8") ||
                       Equal(block->Option("TYPE"), "C3D8R"))
                        ? 8
                    : Equal(block->Option("TYPE"), "C3D6") ? 6
                    : (Equal(block->Option("TYPE"), "C3D5") ||
                       Equal(block->Option("TYPE"), "FESA_C3D5R"))
                        ? 5
                        : 10;
  for (std::size_t r = 0; r < block->rows->size();) {
    const auto* first = block->rows->at(r);
    std::unique_ptr<InputElement> record(new (std::nothrow) InputElement);
    if (!record)
      return Error(source, first->line, block->step, "Input element allocation",
                   d, StatusCode::kAllocationFailure);
    record->line = first->line;
    auto fail = [&](const InputLine* row, const char* message,
                    StatusCode code) {
      if (d) d->element_id = record->id;
      return Error(source, row->line, block->step, message, d, code);
    };
    bool first_row = true;
    while (true) {
      const auto* row = block->rows->at(r);
      std::vector<std::string> fields(0);
      Fields(row->text->c_str(), &fields);
      std::size_t start = 0;
      if (first_row) {
        if (fields.empty() ||
            !Integer(fields[0].c_str(), 1, 999999999, &record->id))
          return fail(row, "Invalid element ID", StatusCode::kInputError);
        start = 1;
        first_row = false;
      }
      if (fields.size() - start >
          static_cast<std::size_t>(nodes - record->node_count))
        return fail(row, "Too many element nodes", StatusCode::kUnsupported);
      for (std::size_t i = start; i < fields.size(); ++i) {
        int node = 0;
        if (!Integer(fields[i].c_str(), 1, 999999999, &node))
          return fail(row, "Invalid element node", StatusCode::kInputError);
        const int index = record->node_count++;
        record->node_ids[index] = node;
        record->node_lines[index] = row->line;
      }
      record->last_line = row->line;
      ++r;
      if (record->node_count == nodes) break;
      const char* begin = row->text->c_str();
      const char* end = begin + row->text->size();
      while (end > begin && (end[-1] == ' ' || end[-1] == '\t')) --end;
      if (nodes != 10 || end == begin || end[-1] != ',' ||
          r == block->rows->size())
        return fail(row, "Incomplete element connectivity",
                    StatusCode::kInputError);
    }
    block->elements->push_back(record.get());
    record.release();
  }
  return StatusCode::kOk;
}
StatusCode ValidateDynamic(const char* source, const InputBlock* b,
                           Diagnostic* d) {
  if (!OptionIs(b, "APPLICATION", "TRANSIENT FIDELITY") ||
      !OptionIs(b, "TIMEINTEGRATOR", "HHT-TF") ||
      !OptionIs(b, "INCREMENTATION", "AGGRESSIVE"))
    return Error(source, b->line, b->step, "Unsupported DYNAMIC option", d,
                 StatusCode::kUnsupported);
  double alpha = -.05;
  if (b->HasOption("ALPHA") &&
      (!Real(b->Option("ALPHA"), &alpha) || alpha < -1. / 3. || alpha > 0))
    return Error(source, b->line, b->step, "Invalid HHT alpha", d);
  if (b->rows->size() != 1)
    return Error(source, b->rows->size() > 1 ? b->rows->at(1)->line : b->line,
                 b->step, "DYNAMIC requires one data row", d);
  const auto* row = b->rows->front();
  std::vector<std::string> fields;
  Fields(row->text->c_str(), &fields);
  if (fields.size() < 2 || fields.size() > 4)
    return Error(source, row->line, b->step, "Invalid DYNAMIC data fields", d);
  for (std::size_t i = 0; i < fields.size(); ++i) {
    if (i >= 2 && fields[i].empty()) continue;
    double value = 0;
    if (!Real(fields[i].c_str(), &value) ||
        (i < 2 ? value <= 0 : value < 0))
      return Error(source, row->line, b->step, "Invalid DYNAMIC time value", d);
  }
  return StatusCode::kOk;
}
StatusCode ValidateRows(const char* source, InputBlock* b, bool nonlinear,
                        bool dynamic,
                        Diagnostic* d) {
  const std::string* k = b->keyword;
  const std::size_t count = b->rows->size();
  if (*k == "HEADING") return StatusCode::kOk;
  if (*k == "DYNAMIC") return ValidateDynamic(source, b, d);
  if (*k == "AMPLITUDE") {
    int point_count = 0;
    return abaqus::ValidateAmplitudeRows(source, b, &point_count, d);
  }
  const bool data =
      *k == "NODE" || *k == "ELEMENT" || *k == "NSET" || *k == "ELSET" ||
      *k == "BEAMGENERALSECTION" || *k == "ELASTIC" || *k == "DENSITY" ||
      *k == "SOLIDSECTION" || *k == "SHELLSECTION" || *k == "BOUNDARY" ||
      *k == "CLOAD" || *k == "DLOAD" || *k == "STATIC" || *k == "CONTROLS" ||
      *k == "NODEOUTPUT" || *k == "ELEMENTOUTPUT" || *k == "INITIALCONDITIONS";
  if (!data && count)
    return Error(
        source, b->rows->front()->line, b->step, "Keyword does not accept data",
        d,
        *k == "INSTANCE" ? StatusCode::kUnsupported : StatusCode::kInputError);
  if (*k == "INITIALCONDITIONS" && count == 0)
    return Error(source, b->line, b->step, "Initial velocity requires data", d);
  if (*k == "SHELLSECTION" && count != 1)
    return Error(source, b->line, b->step,
                 "Shell section requires one homogeneous thickness row", d,
                 StatusCode::kUnsupported);
  if (((*k == "ELASTIC" || *k == "DENSITY") && count != 1) ||
      (*k == "SOLIDSECTION" && count > 1) ||
      (*k == "BEAMGENERALSECTION" && (count < 1 || count > 2)) ||
      (*k == "STATIC" && count > 1) ||
      ((*k == "NODEOUTPUT" || *k == "ELEMENTOUTPUT") && count == 0) ||
      ((*k == "BOUNDARY" || *k == "CLOAD" || *k == "DLOAD") && count == 0 &&
       !Equal(b->Option("OP"), "NEW")))
    return Error(source, b->line, b->step, "Invalid data row count", d);
  if (*k == "CONTROLS") {
    if (b->HasOption("RESET")) {
      if (b->options->size() != 1 || count)
        return Error(source, b->line, b->step,
                     "RESET accepts no other options or data", d);
    } else if (Equal(b->Option("PARAMETERS"), "TIMEINCREMENTATION")) {
      if ((!nonlinear && !dynamic) || b->HasOption("FIELD"))
        return Error(source, b->line, b->step,
                     "TIME controls require nonlinear step without FIELD", d,
                     StatusCode::kUnsupported);
      NewtonControls n;
      IncrementControls c;
      return abaqus::ReadTimeControlsBlock(source, b, &n, &c, d, dynamic);
    } else if (!Equal(b->Option("PARAMETERS"), "FIELD") || count < 1 ||
               count > 2) {
      return Error(source, b->line, b->step, "CONTROLS requires FIELD data", d);
    }
  }
  if (*k == "STATIC") {
    double period = 0;
    NonlinearStaticParameters parameters;
    if (!(nonlinear ? abaqus::NonlinearStaticData(b, &parameters)
                    : abaqus::StaticPeriod(b, &period)))
      return Error(source, count ? b->rows->front()->line : b->line, b->step,
                   "Invalid STATIC data for selected step procedure", d);
  }
  if (*k == "ELEMENT") return ParseElements(source, b, d);
  for (std::size_t r = 0; r < count; ++r) {
    const InputLine* row = b->rows->at(r);
    std::vector<std::string> f(0);
    Fields(row->text->c_str(), &f);
    bool valid = true;
    StatusCode code = StatusCode::kInputError;
    int id = 0;
    double value = 0;
    if (*k == "NODE") {
      valid = f.size() <= 4 && Integer(f[0].c_str(), 1, 999999999, &id);
      for (std::size_t i = 1; valid && i < f.size(); ++i)
        valid = f[i].empty() || Real(f[i].c_str(), &value);
    } else if (*k == "NSET" || *k == "ELSET") {
      if (b->HasOption("GENERATE")) {
        int start = 0, end = 0, increment = 1;
        valid = f.size() >= 2 && f.size() <= 3 &&
                SetEndpoint(f[0].c_str(), &start) &&
                SetEndpoint(f[1].c_str(), &end) &&
                (f.size() == 2 || f[2].empty() ||
                 Integer(f[2].c_str(), 1, 999999999, &increment)) &&
                start <= end && (end - start) % increment == 0;
      } else {
        valid = f.size() <= 16;
        for (const auto& item : f) {
          const auto dot = item.find('.');
          valid = valid && (dot == std::string::npos
                                ? Integer(item.c_str(), 1, 999999999, &id)
                                : Name(item.substr(0, dot).c_str()) &&
                                      Integer(item.substr(dot + 1).c_str(), 1,
                                              999999999, &id));
        }
      }
    } else if (*k == "DENSITY") {
      valid = f.size() == 1 && Real(f[0].c_str(), &value) && value > 0;
    } else if (*k == "SOLIDSECTION") {
      valid = f.size() <= 2 &&
              std::all_of(f.begin(), f.end(),
                          [](const std::string& item) { return item.empty(); });
      code = StatusCode::kUnsupported;
    } else if (*k == "SHELLSECTION") {
      // Explicit two-point Gauss integration is required. In particular the
      // Abaqus default five-point Simpson row must never be reinterpreted.
      int points = 0;
      if (f.size() != 2 || !Integer(f[1].c_str(), 2, 2, &points))
        return Error(source, row->line, b->step,
                     "Shell section supports exactly two Gauss points", d,
                     StatusCode::kUnsupported);
      valid = Real(f[0].c_str(), &value) && value > 0.0;
    } else if (*k == "ELASTIC" || *k == "BEAMGENERALSECTION") {
      const std::size_t expected = *k == "ELASTIC" ? 2 : r == 0 ? 5 : 3;
      valid = f.size() == expected;
      const bool blank_direction =
          *k == "BEAMGENERALSECTION" && r == 1 &&
          std::all_of(f.begin(), f.end(),
                      [](const std::string& v) { return v.empty(); });
      if (blank_direction)
        valid = f.size() <= 3;
      else
        for (const auto& item : f) valid = valid && Real(item.c_str(), &value);
    } else if (*k == "DLOAD") {
      valid = f.size() == 3 && Name(f[0].c_str(), true) &&
              Listed(f[1].c_str(), "PX,PY,PZ") && Real(f[2].c_str(), &value);
      if (f.size() >= 2 && !Listed(f[1].c_str(), "PX,PY,PZ"))
        code = StatusCode::kUnsupported;
    } else if (*k == "INITIALCONDITIONS") {
      valid = f.size() == 3 && Name(f[0].c_str(), true) &&
              Integer(f[1].c_str(), 1, 6, &id) && Real(f[2].c_str(), &value);
    } else if (*k == "BOUNDARY" || *k == "CLOAD") {
      int first = 0, last = 0;
      valid = f.size() >= 2 && f.size() <= 4 && Name(f[0].c_str(), true) &&
              Integer(f[1].c_str(), 1, 6, &first);
      if (*k == "CLOAD")
        valid = valid && f.size() == 3 && Real(f[2].c_str(), &value);
      else {
        last = first;
        valid = valid &&
                (f.size() < 3 || f[2].empty() ||
                 Integer(f[2].c_str(), first, 6, &last)) &&
                (f.size() < 4 || f[3].empty() || Real(f[3].c_str(), &value));
        if (valid && b->step == 0 && value != 0) {
          valid = false;
          code = StatusCode::kUnsupported;
        }
      }
    } else if (*k == "CONTROLS") {
      valid = f.size() <= (r == 0 ? 8u : 3u);
      for (std::size_t i = 0; valid && i < f.size(); ++i) {
        if (f[i].empty()) continue;
        valid = Real(f[i].c_str(), &value) && value >= 0;
        if (r == 0 && (i == 2 || i == 3)) valid = valid && value > 0;
        if (r == 1 && i == 2 && Equal(b->Option("FIELD"), "ROTATION"))
          valid = false;
      }
    } else if (*k == "NODEOUTPUT" || *k == "ELEMENTOUTPUT") {
      for (const auto& item : f)
        valid = valid && Listed(item.c_str(),
                                *k == "NODEOUTPUT"
                                    ? "CF,RF,TF,U"
                                    : "LE,NFORC,NFORCSO,PE,PEEQ,PEMAG,S,SF");
      code = StatusCode::kUnsupported;
    }
    if (!valid)
      return Error(source, row->line, b->step, "Invalid keyword data", d, code);
  }
  return StatusCode::kOk;
}

StatusCode Validate(const char* source, std::vector<InputBlock*>* blocks,
                    Diagnostic* d) {
  bool part = false, assembly = false, instance = false, active_step = false;
  bool had_part = false, had_assembly = false, had_instance = false;
  bool field_output = false;
  bool nonlinear = false, dynamic = false;
  int step = 0, procedures = 0;
  std::map<std::string, int> amplitude_lines;
  std::string part_name(""), assembly_name(""), previous("");
  for (std::size_t index = 0; index < blocks->size(); ++index) {
    auto* b = blocks->at(index);
    const std::string* k = b->keyword;
    bool scope_ok = true, options_ok = true;
    b->step = active_step ? step : 0;
    b->scope_kind = part          ? InputScope::kPart
                    : assembly    ? InputScope::kAssembly
                    : active_step ? InputScope::kStep
                                  : InputScope::kModel;
    *b->scope = part ? part_name : assembly ? assembly_name : "";
    if (*k == "HEADING")
      scope_ok = b == blocks->front();
    else if (*k == "PART") {
      scope_ok = !part && !assembly && step == 0 && !had_part;
      options_ok = Named(b, "NAME");
      part = had_part = true;
      part_name = b->Option("NAME");
      *b->scope = part_name;
      b->scope_kind = InputScope::kPart;
    } else if (*k == "ENDPART") {
      scope_ok = part;
      part = false;
    } else if (*k == "ASSEMBLY") {
      scope_ok = !part && !assembly && step == 0 && !had_assembly;
      options_ok = Named(b, "NAME");
      assembly = had_assembly = true;
      assembly_name = b->Option("NAME");
      *b->scope = assembly_name;
      b->scope_kind = InputScope::kAssembly;
    } else if (*k == "ENDASSEMBLY") {
      scope_ok = assembly && !instance;
      assembly = false;
    } else if (*k == "INSTANCE") {
      scope_ok = assembly && !instance && !had_instance;
      options_ok = Named(b, "NAME") && Named(b, "PART");
      instance = had_instance = true;
    } else if (*k == "ENDINSTANCE") {
      scope_ok = instance;
      instance = false;
    } else if (*k == "NODE" || *k == "ELEMENT" || *k == "BEAMGENERALSECTION" ||
               *k == "SOLIDSECTION" || *k == "SHELLSECTION") {
      scope_ok = part;
      if (*k == "ELEMENT" &&
          (Equal(b->Option("TYPE"), "FESA_C3D10R") ||
           Equal(b->Option("TYPE"), "C3D10R")))
        return Error(source, b->line, b->step,
                     "FESA_C3D10R retired; use C3D10", d,
                     StatusCode::kUnsupported);
      if (*k == "ELEMENT")
        options_ok = (Equal(b->Option("TYPE"), "B33") ||
                      Equal(b->Option("TYPE"), "C3D4") ||
                      Equal(b->Option("TYPE"), "S4") ||
                      Equal(b->Option("TYPE"), "C3D8") ||
                      Equal(b->Option("TYPE"), "C3D8R") ||
                      Equal(b->Option("TYPE"), "C3D6") ||
                      Equal(b->Option("TYPE"), "C3D5") ||
                      Equal(b->Option("TYPE"), "FESA_C3D5R") ||
                      Equal(b->Option("TYPE"), "C3D10")) &&
                     (!b->HasOption("ELSET") || Named(b, "ELSET"));
      if (*k == "SOLIDSECTION")
        options_ok = Named(b, "ELSET") && Named(b, "MATERIAL");
      if (*k == "SHELLSECTION")
        options_ok = Named(b, "ELSET") && Named(b, "MATERIAL") &&
                     Equal(b->Option("SECTIONINTEGRATION"), "GAUSS");
      if (*k == "BEAMGENERALSECTION")
        options_ok = Named(b, "ELSET") && Named(b, "MATERIAL") &&
                     OptionIs(b, "SECTION", "GENERAL");
    } else if (*k == "NSET" || *k == "ELSET") {
      scope_ok = part || (assembly && !instance);
      options_ok = Named(b, k->c_str()) && (!b->HasOption("INSTANCE") ||
                                            (!part && Named(b, "INSTANCE")));
    } else if (*k == "INITIALCONDITIONS") {
      scope_ok = !part && !assembly && step == 0;
      options_ok = Equal(b->Option("TYPE"), "VELOCITY");
    } else if (*k == "AMPLITUDE") {
      scope_ok = !part && !assembly && step == 0;
      options_ok = Named(b, "NAME") && OptionIs(b, "DEFINITION", "TABULAR") &&
                   OptionIs(b, "VALUE", "RELATIVE") &&
                   (OptionIs(b, "TIME", "STEP TIME") ||
                    Equal(b->Option("TIME"), "TOTAL TIME"));
      if (scope_ok && options_ok) {
        std::string name("");
        Canonical(b->Option("NAME"), &name);
        const auto inserted = amplitude_lines.emplace(name, b->line);
        if (!inserted.second)
          return Error(source, b->line, 0, "Duplicate amplitude name", d,
                       StatusCode::kInputError, inserted.first->second);
      }
    } else if (*k == "MATERIAL") {
      scope_ok = !part && !assembly && step == 0;
      options_ok = Named(b, "NAME");
    } else if (*k == "ELASTIC" || *k == "DENSITY" || *k == "DAMPING") {
      scope_ok = previous == "MATERIAL" || previous == "ELASTIC" ||
                 previous == "DENSITY" || previous == "DAMPING";
      if (*k == "ELASTIC") options_ok = OptionIs(b, "TYPE", "ISOTROPIC");
      if (*k == "DAMPING") {
        for (const char* option : {"ALPHA", "BETA"}) {
          double coefficient = 0;
          if (b->HasOption(option) &&
              (!Real(b->Option(option), &coefficient) || coefficient < 0))
            return Error(source, b->line, b->step,
                         "Damping requires finite nonnegative constants", d);
        }
      }
    } else if (*k == "PREPRINT") {
      scope_ok = !part && !assembly && step == 0;
      for (const auto& option : *b->options)
        options_ok = options_ok && Equal(option.second.c_str(), "NO");
    } else if (*k == "STEP") {
      scope_ok = !part && !assembly && !active_step;
      dynamic = false;
      // STEP options depend on its procedure, which occurs later in the file.
      for (std::size_t next = index + 1; next < blocks->size(); ++next) {
        const auto* keyword = blocks->at(next)->keyword;
        if (*keyword == "STEP" || *keyword == "ENDSTEP") break;
        if (*keyword == "DYNAMIC") dynamic = true;
      }
      const bool activate =
          b->HasOption("NLGEOM") &&
          (!*b->Option("NLGEOM") || Equal(b->Option("NLGEOM"), "YES"));
      options_ok = (!b->HasOption("NLGEOM") || activate ||
                    Equal(b->Option("NLGEOM"), "NO")) &&
                   (!b->HasOption("NAME") || Named(b, "NAME"));
      if (activate && step > 0 && !nonlinear)
        return Error(source, b->line, step + 1,
                     "Linear to nonlinear transition unsupported", d,
                     StatusCode::kUnsupported);
      nonlinear = nonlinear || activate;
      if (dynamic) {
        int limit = 0;
        options_ok = options_ok &&
                     (!b->HasOption("INC") ||
                      Integer(b->Option("INC"), 1,
                              std::numeric_limits<int>::max(), &limit)) &&
                     !b->HasOption("EXTRAPOLATION") &&
                     !b->HasOption("AMPLITUDE") && !b->HasOption("CONVERTSDI");
      } else if (nonlinear) {
        int limit = 0;
        options_ok = options_ok &&
                     (!b->HasOption("INC") ||
                      Integer(b->Option("INC"), 1,
                              std::numeric_limits<int>::max(), &limit)) &&
                     (!b->HasOption("EXTRAPOLATION") ||
                      Listed(b->Option("EXTRAPOLATION"), "LINEAR,NO")) &&
                     OptionIs(b, "AMPLITUDE", "RAMP") &&
                     OptionIs(b, "CONVERTSDI", "YES");
        bool has_shell = false;
        bool has_solid = false;
        bool has_beam = false;
        for (const auto* definition : *blocks) {
          if (*definition->keyword != "ELEMENT") continue;
          const auto* type = definition->Option("TYPE");
          if (Equal(type, "B33")) {
            has_beam = true;
          } else if (Equal(type, "S4")) {
            has_shell = true;
          } else if (Listed(type, "C3D4,C3D10,C3D8,C3D6,C3D5")) {
            has_solid = true;
          } else {
            return Error(source, b->line, step + 1,
                         "Nonlinear step requires B33, S4 or full Solid elements", d,
                         StatusCode::kUnsupported);
          }
        }
        if ((has_shell && has_solid) ||
            (has_beam && (has_shell || has_solid)))
          return Error(source, b->line, step + 1,
                       "Nonlinear element-family mixture unsupported", d,
                       StatusCode::kUnsupported);
      } else {
        options_ok = options_ok && !b->HasOption("INC") &&
                     !b->HasOption("EXTRAPOLATION") &&
                     !b->HasOption("AMPLITUDE") && !b->HasOption("CONVERTSDI");
      }
      active_step = true;
      b->step = ++step;
      procedures = 0;
      field_output = false;
      b->scope_kind = InputScope::kStep;
    } else if (*k == "ENDSTEP") {
      scope_ok = active_step && procedures == 1;
      active_step = false;
    } else if (*k == "BOUNDARY" || *k == "CLOAD" || *k == "DLOAD") {
      scope_ok = !part && !assembly &&
                 (active_step || (*k == "BOUNDARY" && step == 0));
      options_ok =
          (!b->HasOption("OP") || Listed(b->Option("OP"), "MOD,NEW")) &&
          OptionIs(b, "TYPE", "DISPLACEMENT") &&
          (active_step || b->options->empty());
      if (*k == "CLOAD" && b->HasOption("AMPLITUDE")) {
        if (!dynamic)
          return Error(source, b->line, b->step,
                       "Static CLOAD amplitude unsupported", d,
                       StatusCode::kUnsupported);
        if (!Name(b->Option("AMPLITUDE")))
          return Error(source, b->line, b->step,
                       "Invalid CLOAD amplitude name", d);
        std::string amplitude_name;
        Canonical(b->Option("AMPLITUDE"), &amplitude_name);
        if (amplitude_lines.find(amplitude_name) == amplitude_lines.end())
          return Error(source, b->line, b->step,
                       "Undefined CLOAD amplitude", d);
      }
    } else {
      scope_ok = active_step;
      if (*k == "STATIC" || *k == "DYNAMIC")
        ++procedures;
      else if (*k == "CONTROLS")
        options_ok =
            OptionIs(b, "CONVERGENCECHECK", "STRICT") &&
            (!b->HasOption("FIELD") ||
             Listed(b->Option("FIELD"), "GLOBAL,DISPLACEMENT,ROTATION"));
      else if (*k == "RESTART")
        options_ok =
            b->HasOption("WRITE") && Equal(b->Option("FREQUENCY"), "0");
      else if (*k == "OUTPUT") {
        options_ok = b->HasOption("FIELD") != b->HasOption("HISTORY");
        options_ok = options_ok && OptionIs(b, "VARIABLE", "PRESELECT");
        field_output = b->HasOption("FIELD");
      } else if (*k == "NODEOUTPUT" || *k == "ELEMENTOUTPUT" ||
                 *k == "CONTACTOUTPUT") {
        scope_ok = scope_ok && field_output;
        options_ok = *k == "CONTACTOUTPUT"
                         ? Equal(b->Option("VARIABLE"), "PRESELECT")
                         : OptionIs(b, "DIRECTIONS", "YES");
      }
    }
    if (!scope_ok)
      return Error(source, b->line, b->step, "Invalid keyword scope/order", d);
    if (((*k == "PART" || *k == "ASSEMBLY" || *k == "MATERIAL" ||
          *k == "INSTANCE" || *k == "AMPLITUDE") &&
         !b->HasOption("NAME")) ||
        (*k == "INSTANCE" && !b->HasOption("PART")) ||
        (*k == "ELEMENT" && !b->HasOption("TYPE")) ||
        ((*k == "NSET" || *k == "ELSET") && !b->HasOption(k->c_str())) ||
        ((*k == "BEAMGENERALSECTION" || *k == "SOLIDSECTION" ||
          *k == "SHELLSECTION") &&
         (!b->HasOption("ELSET") || !b->HasOption("MATERIAL"))))
      return Error(source, b->line, b->step, "Missing required option", d);
    if (!options_ok)
      return Error(source, b->line, b->step,
                   "Unsupported or missing option value", d,
                   StatusCode::kUnsupported);
    b->compatibility_warning = *k == "RESTART" || *k == "OUTPUT" ||
                               *k == "NODEOUTPUT" || *k == "ELEMENTOUTPUT" ||
                               *k == "CONTACTOUTPUT";
    const StatusCode status = ValidateRows(source, b, nonlinear, dynamic, d);
    if (status != StatusCode::kOk) return status;
    previous = *k;
  }
  if (part || assembly || instance || active_step || !had_part ||
      !had_assembly || !had_instance || step == 0)
    return Error(source, blocks->empty() ? 1 : blocks->back()->line, step,
                 "Incomplete Part/Assembly/Instance/Step", d);
  return StatusCode::kOk;
}
}  // namespace

namespace abaqus {
StatusCode ValidateInputMeshRows(const char* source, InputBlock* block,
                                  Diagnostic* d) {
  return ValidateRows(source, block, false, false, d);
}
}  // namespace abaqus

namespace {
FILE* OpenInput(const char* path) {
#ifdef _WIN32
  const int length =
      MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, path, -1, nullptr, 0);
  if (length == 0) return nullptr;
  std::unique_ptr<wchar_t[]> wide(new wchar_t[length]);
  if (!MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, path, -1, wide.get(),
                           length))
    return nullptr;
  FILE* file = nullptr;
  if (_wfopen_s(&file, wide.get(), L"rb") != 0) return nullptr;
  return file;
#else
  return std::fopen(path, "rb");
#endif
}
bool ReadLine(FILE* file, std::string* line) {
  line->clear();
  int c = 0;
  while ((c = std::fgetc(file)) != EOF) {
    if (c == '\n') return true;
    line->push_back(static_cast<char>(c));
  }
  return !line->empty();
}
}  // namespace
StatusCode ReadInputDeck(const char* input_path, InputDeck** out,
                         Diagnostic* d) {
  if (!out)
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation, "ReadInputDeck output",
                         d);
  *out = nullptr;
  if (!input_path)
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation, "ReadInputDeck path", d);
  int line_number = 0;
  try {
    std::unique_ptr<std::string> source(new std::string(input_path));
    std::unique_ptr<std::vector<InputBlock*>> blocks(
        new std::vector<InputBlock*>(0));
    std::unique_ptr<InputDeck> deck(new InputDeck(source.get(), blocks.get()));
    source.release();
    blocks.release();
    std::unique_ptr<FILE, decltype(&std::fclose)> file(OpenInput(input_path),
                                                       &std::fclose);
    if (!file)
      return Error(input_path, 0, 0, "Cannot open input file", d,
                   StatusCode::kIoFailure);
    // Select one whitelist from actual keyword lines; comments never dispatch.
    bool thermal = false;
    std::string scan(""), scan_trimmed(""), scan_keyword("");
    while (ReadLine(file.get(), &scan)) {
      Trim(scan.c_str(), &scan_trimmed);
      if (scan_trimmed.empty() || scan_trimmed[0] != '*' ||
          scan_trimmed.compare(0, 2, "**") == 0) continue;
      Canonical(scan_trimmed.substr(1, scan_trimmed.find(',') - 1).c_str(),
                &scan_keyword);
      if (scan_keyword == "HEATTRANSFER") thermal = true;
    }
    std::rewind(file.get());
    std::string raw(""), pending("");
    int pending_line = 0, lexical_step = 0;
    bool lexical_active_step = false;
    InputBlock* current = nullptr;
    while (ReadLine(file.get(), &raw)) {
      ++line_number;
      if (raw.find('\0') != std::string::npos)
        return Error(input_path, line_number, 0, "Embedded NUL in input", d);
      if (line_number == 1 && raw.compare(0, 3, "\xef\xbb\xbf") == 0)
        raw.erase(0, 3);
      if (!raw.empty() && raw.back() == '\r') raw.pop_back();
      std::string trimmed("");
      Trim(raw.c_str(), &trimmed);
      if (trimmed.empty() || trimmed.compare(0, 2, "**") == 0) {
        // OptionLine counts these physical rows even though their text is not
        // part of the continued keyword's option syntax.
        if (!pending.empty()) pending.push_back('\n');
        continue;
      }
      if (!pending.empty() && trimmed[0] == '*')
        return Error(input_path, line_number, 0,
                     "Keyword continuation interrupted by keyword", d);
      if (trimmed[0] == '*' || !pending.empty()) {
        if (pending.empty()) {
          pending = trimmed;
          pending_line = line_number;
        } else {
          pending.push_back('\n');
          pending.append(trimmed);
        }
        if (trimmed.back() == ',') continue;
        InputBlock* parsed = nullptr;
        const StatusCode status =
            Keyword(input_path, pending.c_str(), pending_line,
                    lexical_active_step ? lexical_step : 0, thermal, &parsed, d);
        if (status != StatusCode::kOk) return status;
        if (*parsed->keyword == "STEP") {
          ++lexical_step;
          lexical_active_step = true;
        } else if (*parsed->keyword == "ENDSTEP") {
          lexical_active_step = false;
        }
        std::unique_ptr<InputBlock> guard(parsed);
        deck->blocks_->push_back(parsed);
        current = guard.release();
        pending.clear();
      } else {
        if (!current)
          return Error(input_path, line_number, 0, "Data before keyword", d);
        std::unique_ptr<InputLine> row(new InputLine);
        row->line = line_number;
        row->text = new std::string(raw);
        current->rows->push_back(row.get());
        row.release();
      }
    }
    if (std::ferror(file.get()))
      return Error(input_path, line_number, 0, "Input read failure", d,
                   StatusCode::kIoFailure);
    if (!pending.empty())
      return Error(input_path, pending_line, 0,
                   "Incomplete keyword continuation", d);
    const StatusCode status = thermal
        ? abaqus::ValidateThermalDeck(input_path, deck->blocks_, d)
        : Validate(input_path, deck->blocks_, d);
    if (status != StatusCode::kOk) return status;
    *out = deck.release();
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return Error(input_path, line_number, 0, "ReadInputDeck allocation", d,
                 StatusCode::kAllocationFailure);
  }
}
}  // namespace fesa
