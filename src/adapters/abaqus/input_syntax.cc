#include "adapters/abaqus/input_syntax.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <cstring>

namespace fesa {
namespace abaqus {
namespace {
bool Space(char c) { return c == ' ' || c == '\t' || c == '\r' || c == '\n'; }
bool Digit(char c) { return c >= '0' && c <= '9'; }
}  // namespace
void Canonical(const char* text, std::string* out) {
  out->clear();
  for (; *text; ++text) {
    if (!Space(*text))
      out->push_back(*text >= 'a' && *text <= 'z'
                         ? static_cast<char>(*text - 'a' + 'A')
                         : *text);
  }
}
void Trim(const char* text, std::string* out) {
  while (Space(*text)) ++text;
  const char* end = text + std::strlen(text);
  while (end > text && Space(end[-1])) --end;
  out->assign(text, end);
}
void Fields(const char* text, std::vector<std::string>* out) {
  out->clear();
  // Reserve before inserting strings: MSVC Debug's noexcept string move can
  // allocate an iterator proxy while vector grows, defeating status recovery.
  out->reserve(static_cast<std::size_t>(
                   std::count(text, text + std::strlen(text), ',')) +
               1);
  const char* start = text;
  for (const char* p = text;; ++p) {
    if (*p == ',' || *p == '\0') {
      std::string raw(start, p), field("");
      Trim(raw.c_str(), &field);
      out->push_back(field);
      if (*p == '\0') break;
      start = p + 1;
    }
  }
  if (out->size() > 1 && out->back().empty()) out->pop_back();
}
bool Integer(const char* text, int minimum, int maximum, int* out) {
  const char* p = text;
  if (*p == '+' || *p == '-') ++p;
  if (!Digit(*p)) return false;
  while (Digit(*p)) ++p;
  if (*p) return false;
  errno = 0;
  const long long value = std::strtoll(text, nullptr, 10);
  if (errno == ERANGE || value < minimum || value > maximum) return false;
  *out = static_cast<int>(value);
  return true;
}
bool Real(const char* text, double* out) {
  const char* p = text;
  if (*p == '+' || *p == '-') ++p;
  bool digits = false;
  while (Digit(*p)) {
    digits = true;
    ++p;
  }
  if (*p == '.') {
    ++p;
    while (Digit(*p)) {
      digits = true;
      ++p;
    }
  }
  if (!digits) return false;
  if (*p == 'e' || *p == 'E' || *p == 'd' || *p == 'D') {
    ++p;
    if (*p == '+' || *p == '-') ++p;
    if (!Digit(*p)) return false;
    while (Digit(*p)) ++p;
  }
  if (*p) return false;
  std::string normalized(text);
  for (char& c : normalized)
    if (c == 'd' || c == 'D') c = 'E';
  errno = 0;
  char* end = nullptr;
  const double value = std::strtod(normalized.c_str(), &end);
  if (errno == ERANGE || *end || !std::isfinite(value)) return false;
  *out = value;
  return true;
}
bool Name(const char* text, bool qualified) {
  bool present = false;
  for (; *text; ++text) {
    const unsigned char c = static_cast<unsigned char>(*text);
    if (c >= 127 || c < 32 || c == '"' || c == '\'' || c == ',' || c == '=' ||
        (!qualified && c == '.'))
      return false;
    if (!Space(*text)) present = true;
  }
  return present;
}
bool Equal(const char* left, const char* right) {
  std::string a(""), b("");
  Canonical(left, &a);
  Canonical(right, &b);
  return a == b;
}
StatusCode Error(const char* source, int line, int step, const char* operation,
                 Diagnostic* d, StatusCode code, int first_line) {
  if (d) {
    d->source_file = source;
    d->duplicate_boundary = false;
    d->condition_dof = 0;
    d->path_kind = DiagnosticPathKind::kLiteralFile;
    d->source_line = line;
    d->first_source_line = first_line;
    d->step = step;
  }
  return ReportFailure(code,
                       code == StatusCode::kAllocationFailure
                           ? DiagnosticStage::kAllocation
                           : DiagnosticStage::kValidation,
                       operation, d);
}
bool NonlinearStaticData(const InputBlock* block,
                         NonlinearStaticParameters* out) {
  if (!block || !out || block->rows->size() > 1) return false;
  double v[4] = {};
  bool maximum = false;
  if (!block->rows->empty()) {
    std::vector<std::string> fields(0);
    Fields(block->rows->front()->text->c_str(), &fields);
    if (fields.size() > 4) return false;
    for (std::size_t i = 0; i < fields.size(); ++i) {
      if (!fields[i].empty() && !Real(fields[i].c_str(), &v[i])) return false;
      if (i == 3 && !fields[i].empty()) maximum = true;
    }
  }
  if (v[1] == 0) v[1] = 1;
  if (v[0] == 0) v[0] = v[1];
  if (v[2] == 0) v[2] = std::min(v[0], 1e-5 * v[1]);
  NonlinearStaticParameters parsed;
  parsed.period = v[1];
  parsed.initial_increment = v[0];
  parsed.minimum_increment = v[2];
  parsed.has_maximum = maximum;
  parsed.maximum_increment = v[3];
  if (ValidateNonlinearStaticParameters(&parsed, nullptr) != StatusCode::kOk)
    return false;
  out->period = parsed.period;
  out->initial_increment = parsed.initial_increment;
  out->minimum_increment = parsed.minimum_increment;
  out->has_maximum = parsed.has_maximum;
  out->maximum_increment = parsed.maximum_increment;
  return true;
}
bool StaticPeriod(const InputBlock* block, double* period) {
  NonlinearStaticParameters parsed;
  if (!period || !NonlinearStaticData(block, &parsed) ||
      parsed.initial_increment != parsed.period)
    return false;
  *period = parsed.period;
  return true;
}

}  // namespace abaqus
}  // namespace fesa
