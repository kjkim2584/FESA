#include "fesa/loads/beam_line_load.h"

#include <algorithm>
#include <cerrno>
#include <charconv>
#include <cmath>
#include <cstdlib>
#include <map>
#include <memory>
#include <new>
#include <string>
#include <vector>

#include "fesa/elements/beam3d.h"

namespace fesa {
namespace {

// Owns one immutable source definition. Program owns all rows; snapshots may
// borrow these rows, but never another snapshot's buffers.
struct LineRow {
  LineRow() = default;
  ~LineRow() {
    delete key;
    delete token;
    delete label;
    delete ids;
    delete density;
  }
  LineRow(const LineRow&) = delete;
  LineRow& operator=(const LineRow&) = delete;
  LineRow(LineRow&&) = delete;
  LineRow& operator=(LineRow&&) = delete;
  std::string* key = nullptr;
  std::string* token = nullptr;
  std::string* label = nullptr;
  std::vector<int>* ids = nullptr;
  Vector* density = nullptr;
  int step = 0;
  int line = 0;
  int component = 0;
};

// Independent current origin. Definition pointers borrow the immutable program.
struct LineOrigin {
  LineOrigin() = default;
  ~LineOrigin() { delete density; }
  LineOrigin(const LineOrigin&) = delete;
  LineOrigin& operator=(const LineOrigin&) = delete;
  LineOrigin(LineOrigin&&) = delete;
  LineOrigin& operator=(LineOrigin&&) = delete;
  const LineRow* original = nullptr;
  const LineRow* first = nullptr;
  const LineRow* current = nullptr;
  Vector* density = nullptr;
  int multiplicity = 0;
};

using OriginMap = std::map<std::string, LineOrigin*>;

void ClearOrigins(OriginMap* origins) {
  if (!origins) return;
  for (auto entry : *origins) delete entry.second;
  origins->clear();
}

StatusCode Fail(StatusCode code, const char* operation, Diagnostic* d) {
  return ReportFailure(code, DiagnosticStage::kValidation, operation, d);
}

void Source(const char* path, int step, int line, Diagnostic* d) {
  if (!d) return;
  d->source_file = path;
  d->path_kind = DiagnosticPathKind::kLiteralFile;
  d->step = step;
  d->source_line = line;
  d->first_source_line = 0;
  d->node_id = 0;
  d->element_id = 0;
  d->condition_dof = 0;
  d->duplicate_boundary = false;
  d->boundary_dof = 0;
}

void Canonical(std::string* value) {
  std::size_t write = 0;
  for (char c : *value) {
    if (c == ' ' || c == '\t' || c == '\r' || c == '\n') continue;
    if (c >= 'a' && c <= 'z') c = static_cast<char>(c - 'a' + 'A');
    (*value)[write++] = c;
  }
  value->resize(write);
}

void Trim(std::string* value) {
  const auto first = value->find_first_not_of(" \t\r\n");
  if (first == std::string::npos) {
    value->clear();
    return;
  }
  const auto last = value->find_last_not_of(" \t\r\n");
  *value = value->substr(first, last - first + 1);
}

void Split(const std::string* text, std::vector<std::string>* fields) {
  // Avoid noexcept Debug string moves when a populated vector reallocates.
  fields->reserve(static_cast<std::size_t>(
                      std::count(text->begin(), text->end(), ',')) + 1);
  std::size_t start = 0;
  while (true) {
    const auto end = text->find(',', start);
    std::string field = text->substr(start, end - start);
    Trim(&field);
    fields->push_back(field);
    if (end == std::string::npos) return;
    start = end + 1;
  }
}

StatusCode CheckElement(const Domain* domain, int id, Diagnostic* d) {
  const Element* element = nullptr;
  if (domain->FindElement(id, &element, nullptr) != StatusCode::kOk) {
    if (d) d->element_id = id;
    return Fail(StatusCode::kInputError, "DLOAD element target missing", d);
  }
  const auto* beam = dynamic_cast<const Beam3D*>(element);
  const Matrix* axes = beam ? beam->LocalAxes() : nullptr;
  if (!element->GetNode(0) || !element->GetNode(1) || !axes ||
      axes->Rows() != 3 || axes->Columns() != 3) {
    if (d) d->element_id = id;
    return Fail(StatusCode::kUnsupported, "DLOAD requires two-node beam", d);
  }
  return StatusCode::kOk;
}

StatusCode Target(const Domain* domain, const std::string* token, int component,
                  LineRow* row, Diagnostic* d) {
  std::string name = *token;
  Canonical(&name);
  bool qualified = false;
  const auto dot = name.find('.');
  if (dot != std::string::npos) {
    std::string instance = domain->InstanceName();
    Canonical(&instance);
    if (name.substr(0, dot) != instance ||
        name.find('.', dot + 1) != std::string::npos) {
      return Fail(StatusCode::kInputError, "DLOAD invalid instance scope", d);
    }
    qualified = true;
    name.erase(0, dot + 1);
  }
  if (name.empty()) {
    return Fail(StatusCode::kInputError, "DLOAD empty target", d);
  }
  row->ids = new std::vector<int>(0);
  row->key = new std::string("");
  if ((name.front() >= '0' && name.front() <= '9') || name.front() == '+' ||
      name.front() == '-') {
    if (name.front() == '+') name.erase(0, 1);
    int id = 0;
    const auto result =
        std::from_chars(name.data(), name.data() + name.size(), id);
    if (result.ec != std::errc() || result.ptr != name.data() + name.size() ||
        id <= 0 || id > 999999999) {
      return Fail(StatusCode::kInputError, "DLOAD invalid element ID", d);
    }
    StatusCode code = CheckElement(domain, id, d);
    if (code != StatusCode::kOk) return code;
    row->key->assign("ELEMENT:");
    row->key->append(std::to_string(id));
    row->ids->push_back(id);
  } else {
    std::string scope = qualified ? "PART." : "ASSEMBLY.";
    scope += qualified ? domain->PartName() : domain->AssemblyName();
    Canonical(&scope);
    const IndexSet* set = nullptr;
    const StatusCode lookup = domain->FindIndexSetByName(
        scope.c_str(), name.c_str(), true, &set, d);
    if (lookup == StatusCode::kAllocationFailure) return lookup;
    if (lookup != StatusCode::kOk) {
      return Fail(StatusCode::kInputError, "DLOAD scoped ELSET missing", d);
    }
    for (int i = 0; i < set->Size(); ++i) {
      StatusCode code = CheckElement(domain, set->Ids()[i], d);
      if (code != StatusCode::kOk) return code;
      row->ids->push_back(set->Ids()[i]);
    }
    row->key->assign(scope);
    row->key->append(":ELSET:");
    row->key->append(name);
  }
  row->key->push_back(':');
  row->key->append(std::to_string(component));
  return StatusCode::kOk;
}

StatusCode MakeOrigin(const LineRow* row, const LineOrigin* previous,
                      LineOrigin** out, Diagnostic* d) {
  std::unique_ptr<LineOrigin> origin(new LineOrigin);
  origin->original = previous ? previous->original : row;
  origin->first = row;
  origin->current = row;
  origin->multiplicity = 1;
  StatusCode code = Vector::Create(3, &origin->density, d);
  if (code != StatusCode::kOk) return code;
  code = origin->density->CopyFrom(row->density, d);
  if (code != StatusCode::kOk) return code;
  *out = origin.release();
  return StatusCode::kOk;
}

StatusCode CloneOrigins(const OriginMap* from, OriginMap* to, Diagnostic* d) {
  for (auto entry : *from) {
    std::unique_ptr<LineOrigin> origin(new LineOrigin);
    origin->original = entry.second->original;
    origin->first = entry.second->first;
    origin->current = entry.second->current;
    origin->multiplicity = entry.second->multiplicity;
    StatusCode code = Vector::Create(3, &origin->density, d);
    if (code != StatusCode::kOk) return code;
    code = origin->density->CopyFrom(entry.second->density, d);
    if (code != StatusCode::kOk) return code;
    auto inserted = to->emplace(entry.first, origin.get());
    if (!inserted.second) {
      return Fail(StatusCode::kInvalidArgument, "DLOAD duplicate snapshot key",
                  d);
    }
    origin.release();
  }
  return StatusCode::kOk;
}

}  // namespace

class BeamLineLoadProgramData {
 public:
  BeamLineLoadProgramData() = default;
  ~BeamLineLoadProgramData() {
    if (rows_)
      for (LineRow* row : *rows_) delete row;
    delete rows_;
    delete operations_;
    delete source_;
    delete source_text_;
  }
  BeamLineLoadProgramData(const BeamLineLoadProgramData&) = delete;
  BeamLineLoadProgramData& operator=(const BeamLineLoadProgramData&) = delete;
  BeamLineLoadProgramData(BeamLineLoadProgramData&&) = delete;
  BeamLineLoadProgramData& operator=(BeamLineLoadProgramData&&) = delete;

 private:
  friend class BeamLineLoadProgram;
  const Domain* domain_ =
      nullptr;  // Borrowed throughout program/state lifetime.
  std::string* source_ = nullptr;          // Owned diagnostic path.
  std::vector<LineRow*>* rows_ = nullptr;  // Owned list and rows.
  std::vector<std::string>* source_text_ =
      nullptr;                                 // Owned original DLOAD text.
  std::map<int, bool>* operations_ = nullptr;  // Owned; true is NEW.
};

class BeamLineLoadStateData {
 public:
  BeamLineLoadStateData() = default;
  ~BeamLineLoadStateData() {
    ClearOrigins(origins_);
    delete origins_;
    if (sums_)
      for (auto entry : *sums_) delete entry.second;
    delete sums_;
    delete zero_;
  }
  BeamLineLoadStateData(const BeamLineLoadStateData&) = delete;
  BeamLineLoadStateData& operator=(const BeamLineLoadStateData&) = delete;
  BeamLineLoadStateData(BeamLineLoadStateData&&) = delete;
  BeamLineLoadStateData& operator=(BeamLineLoadStateData&&) = delete;

 private:
  friend class BeamLineLoadProgram;
  friend class BeamLineLoadState;
  const BeamLineLoadProgram* program_ =
      nullptr;                      // Borrowed owner of source rows.
  const Domain* domain_ = nullptr;  // Borrowed immutable model.
  OriginMap* origins_ = nullptr;    // Owned list and independent entries.
  std::map<int, Vector*>* sums_ = nullptr;  // Owned vectors, global F/L.
  Vector* zero_ = nullptr;                  // Owned absent-density source.
  int step_ = 0;
};

class BeamLineLoadWorkspaceData {
 public:
  BeamLineLoadWorkspaceData() = default;
  ~BeamLineLoadWorkspaceData() {
    delete map_;
    delete transform_;
    delete first_;
    delete second_;
  }
  BeamLineLoadWorkspaceData(const BeamLineLoadWorkspaceData&) = delete;
  BeamLineLoadWorkspaceData& operator=(const BeamLineLoadWorkspaceData&) =
      delete;
  BeamLineLoadWorkspaceData(BeamLineLoadWorkspaceData&&) = delete;
  BeamLineLoadWorkspaceData& operator=(BeamLineLoadWorkspaceData&&) = delete;

 private:
  friend class BeamLineLoadWorkspace;
  friend StatusCode ComputeBeamLineLoad(const Element*, const Vector*,
                                        BeamLineLoadWorkspace*, Vector*,
                                        Vector*, Diagnostic*);
  Matrix* map_ = nullptr;        // Owned L1 global 12x3 coefficient map.
  Matrix* transform_ = nullptr;  // Owned block-diagonal 12x12 T.
  Vector* first_ = nullptr;      // Owned coordinate/difference scratch, size 3.
  Vector* second_ = nullptr;     // Owned endpoint coordinate scratch, size 3.
};

BeamLineLoadProgram::BeamLineLoadProgram(BeamLineLoadProgramData* data)
    : data_(data) {}
BeamLineLoadProgram::~BeamLineLoadProgram() { delete data_; }
const Domain* BeamLineLoadProgram::SourceDomain() const {
  return data_->domain_;
}

StatusCode BeamLineLoadProgram::Create(const InputDeck* deck,
                                       const Domain* domain,
                                       BeamLineLoadProgram** out,
                                       Diagnostic* d) {
  if (out) *out = nullptr;
  if (!deck || !domain || !out) {
    return Fail(StatusCode::kInvalidArgument, "BeamLineLoadProgram::Create", d);
  }
  Source(deck->SourceFile(), 0, 0, d);
  try {
    std::unique_ptr<BeamLineLoadProgramData> data(new BeamLineLoadProgramData);
    data->domain_ = domain;
    data->source_ = new std::string(deck->SourceFile());
    data->rows_ = new std::vector<LineRow*>(0);
    data->source_text_ = new std::vector<std::string>(0);
    data->operations_ = new std::map<int, bool>;
    std::size_t source_count = deck->Blocks()->size();
    for (const InputBlock* block : *deck->Blocks())
      source_count += block->rows->size();
    // Copy source strings into reserved storage; Debug string move is noexcept
    // despite proxy allocation, so vector growth must not relocate strings.
    data->source_text_->reserve(source_count);
    for (const InputBlock* block : *deck->Blocks()) {
      if (*block->keyword != "DLOAD") continue;
      data->source_text_->push_back(*block->original);
      Source(deck->SourceFile(), block->step, block->line, d);
      const StepDefinition* step = nullptr;
      if (block->scope_kind != InputScope::kStep || block->step <= 0 ||
          domain->FindStepDefinition(block->step, &step, nullptr) !=
              StatusCode::kOk) {
        return Fail(StatusCode::kInputError, "DLOAD requires valid step scope",
                    d);
      }
      for (auto option : *block->options) {
        if (option.first != "OP") {
          return Fail(StatusCode::kUnsupported, "DLOAD unsupported option", d);
        }
      }
      std::string operation =
          block->HasOption("OP") ? block->Option("OP") : "MOD";
      Canonical(&operation);
      if (operation != "MOD" && operation != "NEW") {
        return Fail(StatusCode::kUnsupported, "DLOAD unsupported OP", d);
      }
      auto inserted =
          data->operations_->emplace(block->step, operation == "NEW");
      if (!inserted.second && inserted.first->second != (operation == "NEW")) {
        return Fail(StatusCode::kInputError, "DLOAD mixed step operations", d);
      }
      for (const InputLine* input : *block->rows) {
        Source(deck->SourceFile(), block->step, input->line, d);
        data->source_text_->push_back(*input->text);
        std::vector<std::string> fields(0);
        Split(input->text, &fields);
        if (fields.size() == 4 && fields.back().empty()) fields.pop_back();
        if (fields.size() != 3 || fields[0].empty() || fields[2].empty()) {
          return Fail(StatusCode::kInputError, "DLOAD requires target,label,q",
                      d);
        }
        std::string original_label = fields[1];
        Canonical(&fields[1]);
        int component = -1;
        if (fields[1] == "PX") component = 0;
        if (fields[1] == "PY") component = 1;
        if (fields[1] == "PZ") component = 2;
        if (component < 0) {
          return Fail(StatusCode::kUnsupported, "DLOAD unsupported label", d);
        }
        if (d) d->condition_dof = component + 1;
        // Approved decimal syntax includes Fortran D/d exponents.
        for (char& c : fields[2])
          if (c == 'D' || c == 'd') c = 'E';
        errno = 0;
        char* end = nullptr;
        const double magnitude = std::strtod(fields[2].c_str(), &end);
        if (end != fields[2].c_str() + fields[2].size() ||
            !std::isfinite(magnitude) || errno == ERANGE) {
          return Fail(StatusCode::kInputError, "DLOAD invalid finite density",
                      d);
        }
        std::unique_ptr<LineRow> row(new LineRow);
        row->token = new std::string(fields[0]);
        row->label = new std::string(original_label);
        row->step = block->step;
        row->line = input->line;
        row->component = component;
        StatusCode code = Target(domain, &fields[0], component, row.get(), d);
        if (code != StatusCode::kOk) return code;
        code = Vector::Create(3, &row->density, d);
        if (code != StatusCode::kOk) return code;
        code = row->density->Set(component, magnitude, d);
        if (code != StatusCode::kOk) return code;
        data->rows_->push_back(row.get());
        row.release();
      }
    }
    std::unique_ptr<BeamLineLoadProgram> program(
        new BeamLineLoadProgram(data.get()));
    data.release();
    *out = program.release();
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return ReportFailure(StatusCode::kAllocationFailure,
                         DiagnosticStage::kAllocation,
                         "BeamLineLoadProgram allocation", d);
  }
}

BeamLineLoadState::BeamLineLoadState(BeamLineLoadStateData* data)
    : data_(data) {}
BeamLineLoadState::~BeamLineLoadState() { delete data_; }
int BeamLineLoadState::StepNumber() const { return data_->step_; }
const Domain* BeamLineLoadState::SourceDomain() const { return data_->domain_; }

StatusCode BeamLineLoadProgram::ResolveStep(int step_number,
                                            const BeamLineLoadState* previous,
                                            BeamLineLoadState** out,
                                            Diagnostic* d) const {
  if (out) *out = nullptr;
  Source(data_->source_->c_str(), step_number, 0, d);
  if (!out)
    return Fail(StatusCode::kInvalidArgument, "DLOAD null state output", d);
  const StepDefinition* step = nullptr;
  if (step_number <= 0 || data_->domain_->FindStepDefinition(
                              step_number, &step, nullptr) != StatusCode::kOk) {
    return Fail(StatusCode::kInputError, "DLOAD unknown step", d);
  }
  if ((step_number == 1 && previous) ||
      (step_number > 1 && (!previous || previous->data_->program_ != this ||
                           previous->StepNumber() != step_number - 1))) {
    return Fail(StatusCode::kInvalidArgument,
                "DLOAD preceding snapshot required", d);
  }
  try {
    std::unique_ptr<BeamLineLoadStateData> state(new BeamLineLoadStateData);
    state->program_ = this;
    state->domain_ = data_->domain_;
    state->step_ = step_number;
    state->origins_ = new OriginMap;
    state->sums_ = new std::map<int, Vector*>;
    StatusCode code = Vector::Create(3, &state->zero_, d);
    if (code != StatusCode::kOk) return code;
    const auto operation = data_->operations_->find(step_number);
    const bool clear =
        operation != data_->operations_->end() && operation->second;
    if (previous && !clear) {
      code = CloneOrigins(previous->data_->origins_, state->origins_, d);
      if (code != StatusCode::kOk) return code;
    }
    // Current-step markers live separately from previous snapshot storage.
    std::map<std::string, LineOrigin*> changes;
    for (const LineRow* row : *data_->rows_) {
      if (row->step != step_number) continue;
      Source(data_->source_->c_str(), step_number, row->line, d);
      if (d) d->condition_dof = row->component + 1;
      auto changed = changes.find(*row->key);
      if (changed != changes.end()) {
        // DLOAD source multiplicity includes zero-valued definition rows.
        // A zero adds no density, but remains ambiguous for a later MOD.
        code = changed->second->density->AddScaled(1, row->density, d);
        if (code != StatusCode::kOk) return code;
        ++changed->second->multiplicity;
        changed->second->current = row;
        continue;
      }
      auto old = state->origins_->find(*row->key);
      if (old != state->origins_->end() && old->second->multiplicity > 1) {
        if (d) d->first_source_line = old->second->first->line;
        return Fail(StatusCode::kInputError, "DLOAD MOD repeated original key",
                    d);
      }
      LineOrigin* raw = nullptr;
      code = MakeOrigin(
          row, old == state->origins_->end() ? nullptr : old->second, &raw, d);
      if (code != StatusCode::kOk) return code;
      std::unique_ptr<LineOrigin> origin(raw);
      if (old != state->origins_->end()) {
        delete old->second;
        old->second = origin.get();
      } else {
        auto inserted = state->origins_->emplace(*row->key, origin.get());
        if (!inserted.second) {
          return Fail(StatusCode::kInvalidArgument, "DLOAD snapshot insertion",
                      d);
        }
      }
      origin.release();  // The independent snapshot now owns the replacement.
      auto marked = changes.emplace(*row->key, raw);
      if (!marked.second) {
        return Fail(StatusCode::kInvalidArgument, "DLOAD current key insertion",
                    d);
      }
    }
    for (auto entry : *state->origins_) {
      const LineOrigin* origin = entry.second;
      for (int id : *origin->current->ids) {
        auto found = state->sums_->find(id);
        if (found == state->sums_->end()) {
          Vector* raw = nullptr;
          code = Vector::Create(3, &raw, d);
          if (code != StatusCode::kOk) return code;
          std::unique_ptr<Vector> sum(raw);
          auto inserted = state->sums_->emplace(id, raw);
          if (!inserted.second) {
            return Fail(StatusCode::kInvalidArgument, "DLOAD density insertion",
                        d);
          }
          sum.release();
          found = inserted.first;
        }
        Source(data_->source_->c_str(), step_number, origin->current->line, d);
        if (d) d->element_id = id;
        code = found->second->AddScaled(1, origin->density, d);
        if (code != StatusCode::kOk) return code;
      }
    }
    std::unique_ptr<BeamLineLoadState> result(
        new BeamLineLoadState(state.get()));
    state.release();
    *out = result.release();
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return ReportFailure(StatusCode::kAllocationFailure,
                         DiagnosticStage::kAllocation,
                         "BeamLineLoadState allocation", d);
  }
}

StatusCode BeamLineLoadState::CopyDensity(int element_id, Vector* out,
                                          Diagnostic* d) const {
  if (!out || out->Size() != 3) {
    return Fail(StatusCode::kInvalidArgument, "DLOAD density output dimension",
                d);
  }
  const Element* element = nullptr;
  if (data_->domain_->FindElement(element_id, &element, nullptr) !=
      StatusCode::kOk) {
    if (d) d->element_id = element_id;
    return Fail(StatusCode::kInputError, "DLOAD density unknown element", d);
  }
  const auto found = data_->sums_->find(element_id);
  return out->CopyFrom(
      found == data_->sums_->end() ? data_->zero_ : found->second, d);
}

BeamLineLoadWorkspace::BeamLineLoadWorkspace(BeamLineLoadWorkspaceData* data)
    : data_(data) {}
BeamLineLoadWorkspace::~BeamLineLoadWorkspace() { delete data_; }

StatusCode BeamLineLoadWorkspace::Create(BeamLineLoadWorkspace** out,
                                         Diagnostic* d) {
  if (out) *out = nullptr;
  if (!out)
    return Fail(StatusCode::kInvalidArgument, "DLOAD null workspace", d);
  try {
    std::unique_ptr<BeamLineLoadWorkspaceData> data(
        new BeamLineLoadWorkspaceData);
    StatusCode code = Matrix::Create(12, 3, &data->map_, d);
    if (code != StatusCode::kOk) return code;
    code = Matrix::Create(12, 12, &data->transform_, d);
    if (code != StatusCode::kOk) return code;
    code = Vector::Create(3, &data->first_, d);
    if (code != StatusCode::kOk) return code;
    code = Vector::Create(3, &data->second_, d);
    if (code != StatusCode::kOk) return code;
    std::unique_ptr<BeamLineLoadWorkspace> result(
        new BeamLineLoadWorkspace(data.get()));
    data.release();
    *out = result.release();
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return ReportFailure(StatusCode::kAllocationFailure,
                         DiagnosticStage::kAllocation,
                         "BeamLineLoadWorkspace allocation", d);
  }
}

StatusCode ComputeBeamLineLoad(const Element* element,
                               const Vector* global_density,
                               BeamLineLoadWorkspace* workspace,
                               Vector* global_equivalent,
                               Vector* local_equivalent, Diagnostic* d) {
  if (!element || !global_density || !workspace || !global_equivalent ||
      !local_equivalent || global_density->Size() != 3 ||
      global_equivalent->Size() != 12 || local_equivalent->Size() != 12 ||
      global_equivalent == local_equivalent ||
      global_density == global_equivalent ||
      global_density == local_equivalent) {
    return Fail(StatusCode::kInvalidArgument, "ComputeBeamLineLoad arguments",
                d);
  }
  const auto* beam = dynamic_cast<const Beam3D*>(element);
  const Matrix* axes = beam ? beam->LocalAxes() : nullptr;
  if (!axes || axes->Rows() != 3 || axes->Columns() != 3 ||
      !element->GetNode(0) || !element->GetNode(1)) {
    return Fail(StatusCode::kUnsupported, "ComputeBeamLineLoad beam axes", d);
  }
  if (global_equivalent == beam->InternalForce() ||
      local_equivalent == beam->InternalForce() ||
      global_equivalent == beam->EndResults() ||
      local_equivalent == beam->EndResults())
    return Fail(StatusCode::kInvalidArgument,
                "ComputeBeamLineLoad result alias", d);
  BeamLineLoadWorkspaceData* data = workspace->data_;
  for (int i = 0; i < 3; ++i) {
    double a = 0, b = 0;
    StatusCode code = element->GetNode(0)->Coordinate(i, &a, d);
    if (code != StatusCode::kOk) return code;
    code = element->GetNode(1)->Coordinate(i, &b, d);
    if (code != StatusCode::kOk) return code;
    code = data->first_->Set(i, a, d);
    if (code != StatusCode::kOk) return code;
    code = data->second_->Set(i, b, d);
    if (code != StatusCode::kOk) return code;
  }
  StatusCode code = data->second_->AddScaled(-1, data->first_, d);
  if (code != StatusCode::kOk) return code;
  double length = 0;
  code = data->second_->Norm2(&length, d);
  if (code != StatusCode::kOk) return code;
  const double half_length = length / 2;
  const double moment = (length / 12) * length;
  if (!std::isfinite(length) || length <= 0 || !std::isfinite(moment) ||
      moment == 0) {
    return ReportFailure(StatusCode::kLibraryFailure,
                         DiagnosticStage::kComputation,
                         "DLOAD unrepresentable L1 coefficients", d);
  }
  // L1 map: translation qg*L/2 and moment +(L^2/12)[t]cross*qg
  // at end a, with the opposite moment at end b. These are scalar
  // coefficient definitions; the actual mapping is a common BLAS operation.
  const double* t = axes->Data();
  for (int row = 0; row < 12; ++row) {
    for (int col = 0; col < 3; ++col) {
      const int local_row = row % 6;
      double coefficient = 0;
      if (local_row < 3) {
        coefficient = local_row == col ? half_length : 0;
      } else if (local_row - 3 != col) {
        const int rotation = local_row - 3;
        const int tangent = 3 - rotation - col;
        const double sign = col == (rotation + 1) % 3 ? -1 : 1;
        coefficient = (row < 6 ? moment : -moment) * sign * t[tangent];
      }
      code = data->map_->Set(row, col, coefficient, d);
      if (code != StatusCode::kOk) return code;
    }
  }
  // Fixed zero off-diagonal blocks remain zero throughout workspace lifetime.
  for (int block = 0; block < 4; ++block) {
    for (int row = 0; row < 3; ++row) {
      for (int col = 0; col < 3; ++col) {
        code = data->transform_->Set(block * 3 + row, block * 3 + col,
                                     axes->Data()[row * 3 + col], d);
        if (code != StatusCode::kOk) return code;
      }
    }
  }
  code = Matrix::MultiplyVector(data->map_, false, global_density,
                                global_equivalent, d);
  if (code != StatusCode::kOk) return code;
  return Matrix::MultiplyVector(data->transform_, false, global_equivalent,
                                local_equivalent, d);
}

}  // namespace fesa
