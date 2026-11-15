#include <limits>
#include <string>
#include <vector>

#include "adapters/abaqus/control_reader.h"
#include "adapters/abaqus/input_syntax.h"

namespace fesa {
namespace abaqus {
void CopyNewtonControls(const NewtonControls* from, NewtonControls* to) {
  to->divergence_start = from->divergence_start;
  to->slow_start = from->slow_start;
  to->nonquadratic_start = from->nonquadratic_start;
  to->prediction_limit = from->prediction_limit;
  to->iteration_limit = from->iteration_limit;
}
void CopyIncrementControls(const IncrementControls* from,
                           IncrementControls* to) {
  to->attempt_limit = from->attempt_limit;
  to->slow_success_iterations = from->slow_success_iterations;
  to->quick_iterations = from->quick_iterations;
  to->divergence_cutback = from->divergence_cutback;
  to->slow_cutback = from->slow_cutback;
  to->geometry_cutback = from->geometry_cutback;
  to->success_scale = from->success_scale;
  to->growth = from->growth;
  to->maximum_growth = from->maximum_growth;
  to->extrapolation_cutoff = from->extrapolation_cutoff;
  to->endpoint_stretch = from->endpoint_stretch;
}
namespace {
int* Count(int column, NewtonControls* n, IncrementControls* c) {
  switch (column) {
    case 0:
      return &n->divergence_start;
    case 1:
      return &n->slow_start;
    case 2:
      return &n->nonquadratic_start;
    case 3:
      return &n->prediction_limit;
    case 4:
      return &c->slow_success_iterations;
    case 5:
      return &c->quick_iterations;
    case 7:
      return &c->attempt_limit;
    case 10:
      return &n->iteration_limit;
    default:
      return nullptr;
  }
}
double* Factor(int row, int column, IncrementControls* c) {
  if (row == 1) {
    switch (column) {
      case 0:
        return &c->divergence_cutback;
      case 1:
        return &c->slow_cutback;
      case 2:
        return &c->success_scale;
      case 5:
        return &c->geometry_cutback;
      case 6:
        return &c->growth;
      default:
        return nullptr;
    }
  }
  if (row == 2 && column == 1) return &c->maximum_growth;
  if (row == 2 && column == 5) return &c->extrapolation_cutoff;
  if (row == 3 && column == 0) return &c->endpoint_stretch;
  return nullptr;
}
bool FactorDomain(int row, int column, double value) {
  if (row == 1 && column != 6)
    return value > 0 && (column == 2 ? value <= 1 : value < 1);
  if (row == 2 && column == 5) return value >= 0;
  return value >= 1;
}
}  // namespace
StatusCode ReadTimeControlsBlock(const char* source, const InputBlock* block,
                                 NewtonControls* counts,
                                 IncrementControls* increments, Diagnostic* d,
                                 bool dynamic) {
  if (!block || !counts || !increments)
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation, "TIME controls inputs",
                         d);
  if (block->rows->empty() || block->rows->size() > 4)
    return Error(source, block->line, block->step, "TIME controls row count",
                 d);
  NewtonControls n;
  IncrementControls c;
  const int widths[] = {13, 8, 8, 1};
  for (std::size_t row = 0; row < block->rows->size(); ++row) {
    const auto* data = block->rows->at(row);
    std::vector<std::string> fields(0);
    Fields(data->text->c_str(), &fields);
    if (fields.size() > static_cast<std::size_t>(widths[row]))
      return Error(source, data->line, block->step, "TIME controls field count",
                   d);
    for (std::size_t column = 0; column < fields.size(); ++column) {
      if (fields[column].empty()) continue;
      const int index = static_cast<int>(column);
      if (dynamic && ((row == 2 && column == 5) || row == 3))
        return Error(source, data->line, block->step,
                     "Dynamic DE/DT controls unsupported", d,
                     StatusCode::kUnsupported);
      if (row == 0) {
        int* target = Count(index, &n, &c);
        if (!target)
          return Error(source, data->line, block->step,
                       "Unsupported TIME iteration control", d,
                       StatusCode::kUnsupported);
        const int minimum = index == 0 ? 3 : (index == 4 || index == 5) ? 0 : 1;
        if (!Integer(fields[column].c_str(), minimum,
                     std::numeric_limits<int>::max(), target))
          return Error(source, data->line, block->step,
                       "Invalid TIME iteration count", d);
      } else {
        double* target = Factor(static_cast<int>(row), index, &c);
        if (!target)
          return Error(source, data->line, block->step,
                       "Unsupported TIME factor", d, StatusCode::kUnsupported);
        if (!Real(fields[column].c_str(), target) ||
            !FactorDomain(static_cast<int>(row), index, *target))
          return Error(source, data->line, block->step, "Invalid TIME factor",
                       d);
      }
    }
  }
  CopyNewtonControls(&n, counts);
  CopyIncrementControls(&c, increments);
  return StatusCode::kOk;
}
StatusCode ReadStepTimeControls(const InputDeck* deck, int step,
                                NewtonControls* counts,
                                IncrementControls* increments, Diagnostic* d) {
  bool dynamic = false;
  for (const auto* block : *deck->Blocks())
    if (block->step == step && *block->keyword == "DYNAMIC") dynamic = true;
  for (const auto* block : *deck->Blocks()) {
    if (block->step != step || *block->keyword != "CONTROLS") continue;
    if (block->HasOption("RESET")) {
      NewtonControls n;
      IncrementControls c;
      if (counts) CopyNewtonControls(&n, counts);
      if (increments) CopyIncrementControls(&c, increments);
    } else if (Equal(block->Option("PARAMETERS"), "TIMEINCREMENTATION")) {
      const auto status = ReadTimeControlsBlock(deck->SourceFile(), block,
                                                counts, increments, d, dynamic);
      if (status != StatusCode::kOk) return status;
    }
  }
  return StatusCode::kOk;
}
}  // namespace abaqus
}  // namespace fesa
