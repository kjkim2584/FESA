#include "adapters/abaqus/control_reader.h"

#include <new>
#include <string>
#include <vector>

#include "adapters/abaqus/input_syntax.h"

namespace fesa {
namespace abaqus {
StatusCode ReadThermalFieldControls(const InputDeck* deck, int step,
                                    FieldControls* out, Diagnostic* d) {
  if (!deck || !out || step < 1)
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation,
                         "Thermal field control input", d);
  try {
    FieldControls candidate;
    CopyFieldControls(out, &candidate);
    for (const auto* block : *deck->Blocks()) {
      if (block->step != step || *block->keyword != "CONTROLS") continue;
      FieldControls field;
      if (block->HasOption("RESET")) {
        CopyFieldControls(&field, &candidate);
        continue;
      }
      if (!Equal(block->Option("PARAMETERS"), "FIELD")) continue;
      if (block->HasOption("FIELD") &&
          !Equal(block->Option("FIELD"), "GLOBAL") &&
          !Equal(block->Option("FIELD"), "TEMPERATURE"))
        return Error(deck->SourceFile(), block->line, step,
                     "Thermal controls require temperature field", d,
                     StatusCode::kUnsupported);
      for (std::size_t row = 0; row < block->rows->size(); ++row) {
        const auto* data = block->rows->at(row);
        std::vector<std::string> values(0);
        Fields(data->text->c_str(), &values);
        for (std::size_t i = 0; i < values.size(); ++i) {
          if (values[i].empty()) continue;
          const int index = static_cast<int>(i) + (row == 0 ? 0 : 8);
          double value = 0;
          if (index >= 10 || !Real(values[i].c_str(), &value))
            return Error(deck->SourceFile(), data->line, step,
                         "Invalid thermal FIELD data", d);
          *ControlValue(&field, index) = value;
          if (index == 2) field.q0_explicit = true;
          if (index == 3) field.qu_explicit = true;
        }
      }
      CopyFieldControls(&field, &candidate);
    }
    CopyFieldControls(&candidate, out);
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return Error(deck->SourceFile(), 0, step,
                 "Thermal field controls allocation", d,
                 StatusCode::kAllocationFailure);
  }
}
const char* ControlSymbol(int index) {
  static const char* const symbols[] = {"Rn", "Cn",        "q0",       "qu",
                                        "Rp", "epsilon",   "Cepsilon", "Rl",
                                        "Cf", "epsilon_l", "epsilon_d"};
  return symbols[index];
}
double* ControlValue(FieldControls* f, int index) {
  double* values[] = {&f->rn, &f->cn,        &f->q0,       &f->qu,
                      &f->rp, &f->epsilon,   &f->cepsilon, &f->rl,
                      &f->cf, &f->epsilon_l, &f->epsilon_d};
  return values[index];
}
void CopyFieldControls(const FieldControls* from, FieldControls* to) {
  to->rn = from->rn;
  to->cn = from->cn;
  to->q0 = from->q0;
  to->qu = from->qu;
  to->rp = from->rp;
  to->epsilon = from->epsilon;
  to->cepsilon = from->cepsilon;
  to->rl = from->rl;
  to->cf = from->cf;
  to->epsilon_l = from->epsilon_l;
  to->epsilon_d = from->epsilon_d;
  to->q0_explicit = from->q0_explicit;
  to->qu_explicit = from->qu_explicit;
}
StatusCode ReadStepControls(const InputDeck* deck, int step,
                            FieldControls* displacement,
                            FieldControls* rotation, Diagnostic* d) {
  for (const auto* block : *deck->Blocks()) {
    if (block->step != step || *block->keyword != "CONTROLS") continue;
    FieldControls field;
    if (block->HasOption("RESET")) {
      CopyFieldControls(&field, displacement);
      CopyFieldControls(&field, rotation);
      continue;
    }
    if (!Equal(block->Option("PARAMETERS"), "FIELD")) continue;
    for (std::size_t row = 0; row < block->rows->size(); ++row) {
      const auto* data = block->rows->at(row);
      std::vector<std::string> values(0);
      Fields(data->text->c_str(), &values);
      for (std::size_t i = 0; i < values.size(); ++i) {
        if (values[i].empty()) continue;
        const int index = static_cast<int>(i) + (row == 0 ? 0 : 8);
        double value = 0;
        if (index >= kControlCount || !Real(values[i].c_str(), &value))
          return Error(deck->SourceFile(), data->line, step,
                       "Invalid validated CONTROLS data", d);
        *ControlValue(&field, index) = value;
        if (index == 2) field.q0_explicit = true;
        if (index == 3) field.qu_explicit = true;
      }
    }
    const char* name = block->Option("FIELD");
    if (!*name || Equal(name, "GLOBAL") || Equal(name, "DISPLACEMENT"))
      CopyFieldControls(&field, displacement);
    if (!*name || Equal(name, "GLOBAL") || Equal(name, "ROTATION"))
      CopyFieldControls(&field, rotation);
  }
  return StatusCode::kOk;
}
}  // namespace abaqus
}  // namespace fesa
