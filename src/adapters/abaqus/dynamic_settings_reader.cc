#include "adapters/abaqus/dynamic_settings_reader.h"

#include <climits>
#include <new>
#include <string>
#include <vector>

#include "adapters/abaqus/input_syntax.h"

namespace fesa {
namespace abaqus {
namespace {
StatusCode ReadSettings(const InputDeck* deck, int step, bool direct_only,
                        DynamicStepSettings** out, Diagnostic* d) {
  if (out) *out = nullptr;
  if (!deck || !out || step <= 0)
    return Error(nullptr, 0, step, "DIRECT settings arguments", d,
                 StatusCode::kInvalidArgument);
  const InputBlock* definition = nullptr;
  const InputBlock* dynamic = nullptr;
  for (const auto* block : *deck->Blocks()) {
    if (block->step != step) continue;
    if (*block->keyword == "STEP") definition = block;
    if (*block->keyword == "DYNAMIC") dynamic = block;
  }
  if (!definition)
    return Error(deck->SourceFile(), 0, step, "Undefined dynamic step", d);
  try {
    if (!dynamic || (direct_only && !dynamic->HasOption("DIRECT")))
      return Error(deck->SourceFile(), definition->line, step,
                   "DIRECT procedure required", d, StatusCode::kUnsupported);
    DynamicStepParameters parameters;
    parameters.alpha_specified = dynamic->HasOption("ALPHA");
    if (parameters.alpha_specified)
      Real(dynamic->Option("ALPHA"), &parameters.alpha);
    parameters.increment_limit_specified = definition->HasOption("INC");
    if (parameters.increment_limit_specified)
      Integer(definition->Option("INC"), 1, INT_MAX,
              &parameters.increment_limit);
    std::vector<std::string> fields(0);
    const auto* row = dynamic->rows->front();
    Fields(row->text->c_str(), &fields);
    Real(fields[0].c_str(), &parameters.initial_increment);
    Real(fields[1].c_str(), &parameters.period);
    parameters.minimum_specified = fields.size() > 2 && !fields[2].empty();
    parameters.maximum_specified = fields.size() > 3 && !fields[3].empty();
    if (parameters.minimum_specified)
      Real(fields[2].c_str(), &parameters.minimum);
    if (parameters.maximum_specified)
      Real(fields[3].c_str(), &parameters.maximum);
    const auto status = dynamic->HasOption("DIRECT")
                            ? DynamicStepSettings::CreateDirect(
                                  &parameters, deck->SourceFile(), dynamic->line,
                                  step, out, d)
                            : DynamicStepSettings::CreateAutomatic(
                                  &parameters, deck->SourceFile(), dynamic->line,
                                  step, out, d);
    if (status == StatusCode::kOk) return status;
    return Error(
        deck->SourceFile(), row->line, step, "Invalid dynamic settings", d,
        status == StatusCode::kInvalidArgument ? StatusCode::kInputError
                                               : status);
  } catch (const std::bad_alloc&) {
    return Error(deck->SourceFile(), dynamic->line, step,
                 "DIRECT settings allocation", d,
                 StatusCode::kAllocationFailure);
  }
}
}  // namespace
StatusCode ReadDirectDynamicSettings(const InputDeck* deck, int step,
                                     DynamicStepSettings** out, Diagnostic* d) {
  return ReadSettings(deck, step, true, out, d);
}
StatusCode ReadDynamicSettings(const InputDeck* deck, int step,
                               DynamicStepSettings** out, Diagnostic* d) {
  return ReadSettings(deck, step, false, out, d);
}
}  // namespace abaqus
}  // namespace fesa
