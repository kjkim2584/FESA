#include "adapters/abaqus/amplitude_reader.h"

#include <limits>
#include <memory>
#include <new>

#include "adapters/abaqus/input_syntax.h"

namespace fesa {
namespace abaqus {
namespace {
StatusCode ReadRows(const char* source, const InputBlock* block, Vector* times,
                     Vector* values, int* count, Diagnostic* d) {
  int points = 0;
  double previous = 0;
  bool single_pair = false;
  std::vector<std::string> fields(0);
  for (std::size_t row_index = 0; row_index < block->rows->size(); ++row_index) {
    const auto* row = block->rows->at(row_index);
    Fields(row->text->c_str(), &fields);
    if (row_index == 0) single_pair = fields.size() == 2;
    const bool last = row_index + 1 == block->rows->size();
    if (fields.size() < 2 || fields.size() > 8 || fields.size() % 2 != 0 ||
        (single_pair && fields.size() != 2) ||
        (!single_pair && !last && fields.size() != 8))
      return Error(source, row->line, block->step, "Amplitude data row width", d);
    for (std::size_t i = 0; i < fields.size(); i += 2) {
      double time = 0, value = 0;
      if (!Real(fields[i].c_str(), &time) ||
          !Real(fields[i + 1].c_str(), &value) ||
          (points && time <= previous))
        return Error(source, row->line, block->step,
                     "Amplitude requires finite increasing time/value pairs", d);
      if (points == std::numeric_limits<int>::max())
        return Error(source, row->line, block->step, "Amplitude point limit", d);
      if (times) {
        auto status = times->Set(points, time, d);
        if (status == StatusCode::kOk) status = values->Set(points, value, d);
        if (status != StatusCode::kOk) return status;
      }
      previous = time;
      ++points;
    }
  }
  if (!points)
    return Error(source, block->line, block->step, "Empty amplitude table", d);
  *count = points;
  return StatusCode::kOk;
}
}  // namespace
StatusCode ValidateAmplitudeRows(const char* source, const InputBlock* block,
                                  int* count, Diagnostic* d) {
  if (!source || !block || !block->rows || !count)
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation, "Amplitude rows", d);
  try {
    return ReadRows(source, block, nullptr, nullptr, count, d);
  } catch (const std::bad_alloc&) {
    return Error(source, block->line, block->step, "Amplitude row allocation",
                 d, StatusCode::kAllocationFailure);
  }
}
StatusCode BuildAmplitude(const char* source, const InputBlock* block, int id,
                           Amplitude** out, Diagnostic* d) {
  if (!out)
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation, "Amplitude output", d);
  *out = nullptr;
  int count = 0;
  auto status = ValidateAmplitudeRows(source, block, &count, d);
  if (status != StatusCode::kOk) return status;
  try {
    Vector* raw_times = nullptr;
    status = Vector::Create(count, &raw_times, d);
    std::unique_ptr<Vector> times(raw_times);
    Vector* raw_values = nullptr;
    if (status == StatusCode::kOk)
      status = Vector::Create(count, &raw_values, d);
    std::unique_ptr<Vector> values(raw_values);
    if (status == StatusCode::kOk)
      status = ReadRows(source, block, times.get(), values.get(), &count, d);
    if (status == StatusCode::kOk)
      status = Amplitude::Create(
          id, block->Option("NAME"),
          Equal(block->Option("TIME"), "TOTAL TIME")
              ? AmplitudeTimeBasis::kTotalTime : AmplitudeTimeBasis::kStepTime,
          times.get(), values.get(), out, d);
  } catch (const std::bad_alloc&) {
    return Error(source, block->line, block->step, "Amplitude build allocation",
                 d, StatusCode::kAllocationFailure);
  }
  if (status != StatusCode::kOk && d) {
    d->source_file = source;
    d->source_line = block->line;
    d->step = block->step;
    d->path_kind = DiagnosticPathKind::kLiteralFile;
  }
  return status;
}
}  // namespace abaqus
}  // namespace fesa
