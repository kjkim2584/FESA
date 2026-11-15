#include "adapters/abaqus/thermal_model_builder.h"

#include <algorithm>
#include <limits>
#include <map>
#include <memory>
#include <new>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "adapters/abaqus/amplitude_reader.h"
#include "adapters/abaqus/input_syntax.h"
#include "adapters/abaqus/model_builder.h"
#include "adapters/abaqus/thermal_reader.h"
#include "fesa/elements/thermal_hexa8.h"
#include "fesa/elements/thermal_shell4.h"
#include "fesa/elements/thermal_tetra4.h"
#include "fesa/model/film_property_table.h"

namespace fesa {
namespace {
using abaqus::Canonical;
using abaqus::Equal;
using abaqus::Error;
using abaqus::Fields;
using abaqus::Integer;
using abaqus::Real;
StatusCode At(StatusCode status, const InputDeck* deck, const InputBlock* b,
              int line, Diagnostic* d) {
  if (status == StatusCode::kOk) return status;
  if (status == StatusCode::kInvalidArgument) status = StatusCode::kInputError;
  if (d) {
    d->code = status;
    d->source_file = deck->SourceFile();
    d->source_line = line;
    d->step = b->step;
    d->path_kind = DiagnosticPathKind::kLiteralFile;
  }
  return status;
}
StatusCode PropertyTable(const InputBlock* property, TemperatureTable** out,
                         Diagnostic* d) {
  *out = nullptr;
  if (!property) return StatusCode::kOk;
  const int count = static_cast<int>(property->rows->size());
  Vector *v = nullptr, *t = nullptr;
  auto status = Vector::Create(count, &v, d);
  std::unique_ptr<Vector> values(v);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(count, &t, d);
  std::unique_ptr<Vector> temperatures(t);
  for (int i = 0; i < count && status == StatusCode::kOk; ++i) {
    std::vector<std::string> fields(0);
    Fields(property->rows->at(i)->text->c_str(), &fields);
    double value = 0, temperature = 0;
    Real(fields[0].c_str(), &value);
    if (fields.size() == 2) Real(fields[1].c_str(), &temperature);
    status = v->Set(i, value, d);
    if (status == StatusCode::kOk) status = t->Set(i, temperature, d);
  }
  if (status == StatusCode::kOk)
    status = TemperatureTable::Create(v, t, out, d);
  return status;
}
StatusCode FilmProperties(const InputDeck* deck, Domain* domain,
                          Diagnostic* d) {
  int id = 0;
  for (const auto* block : *deck->Blocks()) {
    if (*block->keyword != "FILMPROPERTY") continue;
    TemperatureTable* raw_table = nullptr;
    auto status = PropertyTable(block, &raw_table, d);
    std::unique_ptr<TemperatureTable> table(raw_table);
    FilmPropertyTable* raw = nullptr;
    if (status == StatusCode::kOk)
      status = FilmPropertyTable::Create(++id, block->Option("NAME"),
                                         table.get(), &raw, d);
    std::unique_ptr<FilmPropertyTable> property(raw);
    if (status == StatusCode::kOk)
      status = domain->RegisterFilmPropertyTable(raw, d);
    if (status != StatusCode::kOk)
      return At(status, deck, block, block->line, d);
    property.release();
  }
  return StatusCode::kOk;
}
StatusCode Materials(const InputDeck* deck, Domain* domain, Diagnostic* d) {
  int id = 0;
  for (std::size_t i = 0; i < deck->Blocks()->size(); ++i) {
    const auto* b = deck->Blocks()->at(i);
    if (*b->keyword != "MATERIAL") continue;
    const InputBlock *k = nullptr, *rho = nullptr, *c = nullptr;
    for (std::size_t j = i + 1; j < deck->Blocks()->size(); ++j) {
      const auto* p = deck->Blocks()->at(j);
      const InputBlock** slot = *p->keyword == "CONDUCTIVITY"   ? &k
                                : *p->keyword == "DENSITY"      ? &rho
                                : *p->keyword == "SPECIFICHEAT" ? &c
                                                                : nullptr;
      if (!slot) break;
      if (*slot)
        return Error(deck->SourceFile(), p->line, 0,
                     "Duplicate thermal property", d);
      *slot = p;
    }
    if (!k)
      return Error(deck->SourceFile(), b->line, 0,
                   "Thermal material needs conductivity", d);
    TemperatureTable *raw_k = nullptr, *raw_c = nullptr;
    auto status = PropertyTable(k, &raw_k, d);
    std::unique_ptr<TemperatureTable> conductivity(raw_k);
    if (status != StatusCode::kOk) return At(status, deck, k, k->line, d);
    status = PropertyTable(c, &raw_c, d);
    std::unique_ptr<TemperatureTable> capacity(raw_c);
    if (status != StatusCode::kOk) return At(status, deck, c, c->line, d);
    double density = 0;
    if (rho) {
      std::vector<std::string> f(0);
      Fields(rho->rows->front()->text->c_str(), &f);
      Real(f[0].c_str(), &density);
    }
    ThermalMaterial* raw = nullptr;
    status = ThermalMaterial::Create(++id, b->Option("NAME"), raw_k,
                                     rho != nullptr, density, raw_c, &raw, d);
    std::unique_ptr<ThermalMaterial> owner(raw);
    if (status == StatusCode::kOk) status = domain->RegisterMaterial(raw, d);
    if (status != StatusCode::kOk) return At(status, deck, b, b->line, d);
    owner.release();
  }
  return StatusCode::kOk;
}
// Resolves original qualification without changing the stored row spelling.
StatusCode Target(const InputDeck* deck, const Domain* domain,
                  const InputBlock* b, int line, const char* spelling,
                  bool element, std::vector<int>* ids, std::string* scope,
                  std::string* name, int* direct, Diagnostic* d) {
  ids->clear();
  *direct = 0;
  *name = spelling;
  *scope = b->scope_kind == InputScope::kPart ? "PART." : "ASSEMBLY.";
  scope->append(b->scope_kind == InputScope::kPart ? domain->PartName()
                                                   : domain->AssemblyName());
  const auto dot = name->find('.');
  if (dot != std::string::npos) {
    if (b->scope_kind == InputScope::kPart ||
        !Equal(name->substr(0, dot).c_str(), domain->InstanceName()))
      return Error(deck->SourceFile(), line, b->step,
                   "Undefined target Instance", d);
    *scope = "PART.";
    scope->append(domain->PartName());
    *name = name->substr(dot + 1);
  }
  if (Integer(name->c_str(), 1, 999999999, direct)) {
    const Node* node = nullptr;
    const auto status =
        element ? StatusCode::kOk : domain->FindNode(*direct, &node, nullptr);
    if (status != StatusCode::kOk)
      return Error(deck->SourceFile(), line, b->step, "Undefined target node",
                   d);
    if (element) {
      const ThermalElement* found = nullptr;
      if (domain->FindThermalElement(*direct, &found, nullptr) !=
          StatusCode::kOk)
        return Error(deck->SourceFile(), line, b->step,
                     "Undefined target element", d);
    }
    ids->push_back(*direct);
    name->clear();
    return StatusCode::kOk;
  }
  const IndexSet* set = nullptr;
  const auto status = domain->FindIndexSetByName(scope->c_str(), name->c_str(),
                                                 element, &set, d);
  if (status != StatusCode::kOk) return At(status, deck, b, line, d);
  if (set->Size() == 0)
    return Error(deck->SourceFile(), line, b->step, "Empty target set", d);
  ids->assign(set->Ids(), set->Ids() + set->Size());
  return StatusCode::kOk;
}
// Resolves section and node layout before creating immutable Node objects.
// Domain owns every adopted definition; the caller discards the entire
// candidate on any error. No active view exists during this construction.
StatusCode Sections(const InputDeck* deck, Domain* domain, bool transient,
                    const std::map<int, const InputBlock*>* element_blocks,
                    const std::map<int, const InputElement*>* records,
                    const std::map<int, int>* node_lines, Diagnostic* d) {
  std::map<int, ThermalSolidSection*> solids;
  std::map<int, ThermalShellSection*> shells;
  std::map<int, std::pair<int, bool>> layout;
  int section_id = 0;
  for (const auto* b : *deck->Blocks()) {
    const bool shell = *b->keyword == "SHELLSECTION";
    if (!shell && *b->keyword != "SOLIDSECTION") continue;
    const Material* base = nullptr;
    auto status = domain->FindMaterialByName(b->Option("MATERIAL"), &base, d);
    if (status != StatusCode::kOk) return At(status, deck, b, b->line, d);
    const auto* material = base->AsThermal();
    if (!material || (transient && !material->HasCapacity()))
      return Error(deck->SourceFile(), b->line, 0,
                   "Transient requires density and specific heat", d);
    int ns = 1;
    ThermalSolidSection* solid = nullptr;
    ThermalShellSection* shell_section = nullptr;
    if (shell) {
      std::vector<std::string> f(0);
      Fields(b->rows->front()->text->c_str(), &f);
      double h = 0;
      Real(f[0].c_str(), &h);
      ns = 5;
      if (f.size() > 1 && !f[1].empty()) Integer(f[1].c_str(), 1, 19, &ns);
      status = ThermalShellSection::Create(++section_id, h, ns, material,
                                           &shell_section, d);
      std::unique_ptr<ThermalShellSection> owner(shell_section);
      if (status == StatusCode::kOk)
        status = domain->RegisterThermalShellSection(shell_section, d);
      if (status != StatusCode::kOk) return At(status, deck, b, b->line, d);
      owner.release();
    } else {
      status = ThermalSolidSection::Create(++section_id, material, &solid, d);
      std::unique_ptr<ThermalSolidSection> owner(solid);
      if (status == StatusCode::kOk)
        status = domain->RegisterThermalSolidSection(solid, d);
      if (status != StatusCode::kOk) return At(status, deck, b, b->line, d);
      owner.release();
    }
    std::string scope("PART.");
    scope.append(*b->scope);
    const IndexSet* set = nullptr;
    status = domain->FindIndexSetByName(scope.c_str(), b->Option("ELSET"), true,
                                        &set, d);
    if (status == StatusCode::kAllocationFailure)
      return At(status, deck, b, b->line, d);
    if (status != StatusCode::kOk)
      return Error(deck->SourceFile(), b->line, 0,
                   "Undefined section element set", d);
    if (set->Size() == 0)
      return Error(deck->SourceFile(), b->line, 0, "Empty section element set",
                   d);
    for (int member = 0; member < set->Size(); ++member) {
      const int id = set->Ids()[member];
      if (solids.find(id) != solids.end() || shells.find(id) != shells.end())
        return Error(deck->SourceFile(), b->line, 0,
                     "Duplicate element section", d);
      const auto* record = records->find(id)->second;
      const auto* source = element_blocks->find(id)->second;
      if (shell != Equal(source->Option("TYPE"), "DS4"))
        return Error(deck->SourceFile(), b->line, 0,
                     "Thermal element section type mismatch", d);
      if ((shell && !solids.empty()) || (!shell && !shells.empty()))
        return Error(deck->SourceFile(), b->line, 0,
                     "Mixed solid/shell thermal model", d);
      if (shell)
        shells.emplace(id, shell_section);
      else
        solids.emplace(id, solid);
      for (int j = 0; j < record->node_count; ++j) {
        const int node = record->node_ids[j];
        if (node_lines->find(node) == node_lines->end())
          return Error(deck->SourceFile(), record->node_lines[j], 0,
                       "Undefined thermal connectivity node", d);
        const auto found = layout.find(node);
        if (found != layout.end() &&
            (found->second.first != ns || found->second.second != shell))
          return Error(deck->SourceFile(), record->node_lines[j], 0,
                       "Shared thermal node has incompatible section layout",
                       d);
        if (found == layout.end())
          layout.emplace(node, std::make_pair(ns, shell));
      }
    }
  }
  if (solids.size() + shells.size() != records->size())
    return Error(deck->SourceFile(), 1, 0,
                 "Every thermal element needs a section", d);
  // Nodes are created exactly once after the section layout is fully known.
  for (const auto* b : *deck->Blocks())
    if (*b->keyword == "NODE")
      for (const auto* row : *b->rows) {
        std::vector<std::string> f(0);
        Fields(row->text->c_str(), &f);
        int id = 0;
        Integer(f[0].c_str(), 1, 999999999, &id);
        double xyz[3] = {};
        for (std::size_t i = 1; i < f.size(); ++i)
          if (!f[i].empty()) Real(f[i].c_str(), &xyz[i - 1]);
        const auto found = layout.find(id);
        const int ns = found == layout.end() ? 1 : found->second.first;
        Node* node = nullptr;
        auto status =
            Node::CreateThermal(id, xyz[0], xyz[1], xyz[2], ns, &node, d);
        std::unique_ptr<Node> owner(node);
        if (status == StatusCode::kOk) status = domain->RegisterNode(node, d);
        if (status != StatusCode::kOk) return At(status, deck, b, row->line, d);
        owner.release();
      }
  for (const auto& entry : *records) {
    const int id = entry.first;
    const auto* record = entry.second;
    const auto* source = element_blocks->find(id)->second;
    const Node* nodes[8] = {};
    auto status = StatusCode::kOk;
    for (int j = 0; j < record->node_count; ++j) {
      status = domain->FindNode(record->node_ids[j], &nodes[j], d);
      if (status != StatusCode::kOk)
        return At(status, deck, source, record->node_lines[j], d);
    }
    ThermalElement* element = nullptr;
    if (Equal(source->Option("TYPE"), "DS4")) {
      ThermalShell4* raw = nullptr;
      status = ThermalShell4::Create(id, nodes, record->node_count,
                                     shells.find(id)->second, &raw, d);
      element = raw;
    } else if (Equal(source->Option("TYPE"), "DC3D8")) {
      ThermalHexa8* raw = nullptr;
      status = ThermalHexa8::Create(id, nodes, record->node_count,
                                    solids.find(id)->second, &raw, d);
      element = raw;
    } else {
      ThermalTetra4* raw = nullptr;
      status = ThermalTetra4::Create(id, nodes, record->node_count,
                                     solids.find(id)->second, &raw, d);
      element = raw;
    }
    std::unique_ptr<ThermalElement> owner(element);
    if (status == StatusCode::kOk)
      status = domain->RegisterThermalElement(element, d);
    if (status != StatusCode::kOk)
      return At(status, deck, source, record->line, d);
    owner.release();
  }
  return StatusCode::kOk;
}
StatusCode Metadata(const InputDeck* deck, Domain* domain, Diagnostic* d) {
  std::map<std::string, int> surfaces;
  for (const auto* b : *deck->Blocks()) {
    if (*b->keyword == "PHYSICALCONSTANTS") {
      if (domain->ThermalConstants())
        return Error(deck->SourceFile(), b->line, 0,
                     "Multiple physical-constants blocks unsupported", d,
                     StatusCode::kUnsupported);
      ThermalPhysicalConstants constants;
      constants.source_line = b->line;
      constants.has_absolute_zero = b->HasOption("ABSOLUTEZERO");
      constants.has_stefan_boltzmann = b->HasOption("STEFANBOLTZMANN");
      if ((constants.has_absolute_zero &&
           !Real(b->Option("ABSOLUTEZERO"), &constants.absolute_zero)) ||
          (constants.has_stefan_boltzmann &&
           !Real(b->Option("STEFANBOLTZMANN"), &constants.stefan_boltzmann)))
        return Error(deck->SourceFile(), b->line, 0,
                     "Invalid physical constants", d);
      const auto status = domain->SetThermalPhysicalConstants(&constants, d);
      if (status != StatusCode::kOk) return At(status, deck, b, b->line, d);
    } else if (*b->keyword == "SURFACE") {
      std::string key("");
      const std::string scope =
          b->scope_kind == InputScope::kPart ? "PART." : "ASSEMBLY.";
      std::string qualified(scope);
      qualified.append(*b->scope);
      qualified.append(":");
      qualified.append(b->Option("NAME"));
      Canonical(qualified.c_str(), &key);
      if (!surfaces.emplace(key, b->line).second)
        return Error(deck->SourceFile(), b->line, 0, "Duplicate surface name",
                     d);
      for (const auto* row : *b->rows) {
        std::vector<std::string> f(0);
        Fields(row->text->c_str(), &f);
        std::vector<int> ids(0);
        std::string target_scope(""), name("");
        int direct = 0;
        auto status = Target(deck, domain, b, row->line, f[0].c_str(), true,
                             &ids, &target_scope, &name, &direct, d);
        if (status != StatusCode::kOk) return status;
        for (int id : ids) {
          const ThermalElement* element = nullptr;
          status = domain->FindThermalElement(id, &element, d);
          if (status != StatusCode::kOk)
            return At(status, deck, b, row->line, d);
          int face = 0;
          const bool valid_face =
              element->type() == ThermalElementType::kDs4
                  ? Equal(f[1].c_str(), "SNEG") || Equal(f[1].c_str(), "SPOS")
                  : f[1].size() == 2 && (f[1][0] == 'S' || f[1][0] == 's') &&
                        Integer(f[1].substr(1).c_str(), 1,
                                element->type() == ThermalElementType::kDc3d8
                                    ? 6
                                    : 4,
                                &face);
          if (!valid_face)
            return Error(deck->SourceFile(), row->line, 0,
                         "Invalid surface face", d);
        }
      }
    }
  }
  return StatusCode::kOk;
}
// Expansion is copied into each original condition row, never retained from
// InputDeck. Original target spelling remains the temporal identity key.
int FaceCode(const char* label) {
  if (Equal(label, "BF")) return 0;
  if (Equal(label, "SNEG") || Equal(label, "FNEG") || Equal(label, "RNEG"))
    return -1;
  if (Equal(label, "SPOS") || Equal(label, "FPOS") || Equal(label, "RPOS"))
    return -2;
  return label && label[0] && label[1] && !label[2] ? label[1] - '0' : 99;
}
StatusCode DistributedRow(const InputDeck* deck, Domain* domain,
                          const InputBlock* block, const InputLine* row,
                          const Amplitude* amplitude, StepDefinition* step,
                          int id, Diagnostic* d) {
  std::vector<std::string> fields(0);
  Fields(row->text->c_str(), &fields);
  const bool surface = *block->keyword == "DSFLUX" ||
                       *block->keyword == "SFILM" ||
                       *block->keyword == "SRADIATE";
  const bool film = *block->keyword == "FILM" || *block->keyword == "SFILM";
  const bool radiation =
      *block->keyword == "RADIATE" || *block->keyword == "SRADIATE";
  std::vector<int> elements(0), faces(0), ids(0);
  std::string scope(""), name("");
  int direct = 0;
  StatusCode status = StatusCode::kOk;
  if (surface) {
    name = fields[0];
    const auto dot = name.find('.');
    const bool part = dot != std::string::npos;
    if (part) {
      if (!Equal(name.substr(0, dot).c_str(), domain->InstanceName()))
        return Error(deck->SourceFile(), row->line, block->step,
                     "Undefined surface Instance", d);
      name = name.substr(dot + 1);
    }
    scope = part ? "PART." : "ASSEMBLY.";
    scope.append(part ? domain->PartName() : domain->AssemblyName());
    const InputBlock* found = nullptr;
    for (const auto* candidate : *deck->Blocks()) {
      if (*candidate->keyword == "SURFACE" &&
          (candidate->scope_kind == InputScope::kPart) == part &&
          Equal(candidate->Option("NAME"), name.c_str())) {
        found = candidate;
        break;
      }
    }
    if (!found)
      return Error(deck->SourceFile(), row->line, block->step,
                   "Undefined thermal surface", d);
    std::set<std::pair<int, int>> unique;
    for (const auto* face_row : *found->rows) {
      std::vector<std::string> f(0);
      Fields(face_row->text->c_str(), &f);
      std::string target_scope(""), target_name("");
      int target = 0;
      status = Target(deck, domain, found, face_row->line, f[0].c_str(), true,
                      &ids, &target_scope, &target_name, &target, d);
      if (status != StatusCode::kOk) return status;
      const int face = FaceCode(f[1].c_str());
      for (int element : ids)
        if (unique.emplace(element, face).second) {
          elements.push_back(element);
          faces.push_back(face);
        }
    }
  } else {
    status = Target(deck, domain, block, row->line, fields[0].c_str(), true,
                    &elements, &scope, &name, &direct, d);
    if (status != StatusCode::kOk) return status;
    const int face = FaceCode(fields[1].c_str());
    for (int element_id : elements) {
      const ThermalElement* element = nullptr;
      status = domain->FindThermalElement(element_id, &element, d);
      if (status != StatusCode::kOk) return status;
      const bool valid =
          face == 0 ||
          (element->type() == ThermalElementType::kDs4
               ? face == -1 || face == -2
               : face >= 1 &&
                     face <= (element->type() == ThermalElementType::kDc3d8
                                  ? 6
                                  : 4));
      if (!valid)
        return Error(deck->SourceFile(), row->line, block->step,
                     "Flux label incompatible with thermal element", d,
                     StatusCode::kUnsupported);
      faces.push_back(face);
    }
  }
  if (elements.empty() ||
      elements.size() >
          static_cast<std::size_t>(std::numeric_limits<int>::max()))
    return Error(deck->SourceFile(), row->line, block->step,
                 "Empty or oversized thermal load expansion", d);
  ThermalConditionSpec spec;
  spec.id = id;
  spec.kind = radiation ? ThermalConditionKind::kRadiation
              : film    ? ThermalConditionKind::kFilm
              : !surface && Equal(fields[1].c_str(), "BF")
                  ? ThermalConditionKind::kBodyFlux
                  : ThermalConditionKind::kSurfaceFlux;
  spec.family = surface ? (radiation ? ThermalLoadFamily::kSurfaceRadiation
                           : film    ? ThermalLoadFamily::kSurfaceFilm
                                     : ThermalLoadFamily::kSurfaceFlux)
                        : (radiation ? ThermalLoadFamily::kElementRadiation
                           : film    ? ThermalLoadFamily::kElementFilm
                                     : ThermalLoadFamily::kElementFlux);
  spec.target = surface  ? ThermalConditionTarget::kSurface
                : direct ? ThermalConditionTarget::kElement
                         : ThermalConditionTarget::kElementSet;
  spec.target_id = direct;
  spec.scope = scope.c_str();
  spec.name = name.c_str();
  spec.face_label = fields[1].c_str();
  spec.expanded_count = static_cast<int>(elements.size());
  spec.expanded_element_ids = elements.data();
  spec.expanded_face_codes = faces.data();
  if (radiation) {
    Real(fields[2].c_str(), &spec.ambient_temperature);
    Real(fields[3].c_str(), &spec.emissivity);
  } else if (film) {
    Real(fields[2].c_str(), &spec.sink_temperature);
    if (!Real(fields[3].c_str(), &spec.coefficient)) {
      status = domain->FindFilmPropertyTableByName(fields[3].c_str(),
                                                   &spec.film_property, d);
      if (status != StatusCode::kOk)
        return At(status, deck, block, row->line, d);
    }
    if (block->HasOption("FILMAMPLITUDE")) {
      status = domain->FindAmplitudeByName(block->Option("FILMAMPLITUDE"),
                                           &spec.film_amplitude, d);
      if (status != StatusCode::kOk)
        return At(status, deck, block, row->line, d);
    }
  } else {
    Real(fields[2].c_str(), &spec.magnitude);
  }
  spec.amplitude = amplitude;
  spec.operation = Equal(block->Option("OP"), "NEW") ? ConditionOperation::kNew
                                                     : ConditionOperation::kMod;
  spec.step = block->step;
  spec.source_line = row->line;
  spec.source_file = deck->SourceFile();
  ThermalConditionDefinition* raw = nullptr;
  status = ThermalConditionDefinition::Create(&spec, &raw, d);
  std::unique_ptr<ThermalConditionDefinition> owner(raw);
  if (status == StatusCode::kOk)
    status = domain->RegisterThermalCondition(raw, d);
  if (status != StatusCode::kOk) return At(status, deck, block, row->line, d);
  owner.release();
  return step->AddThermalCondition(raw, d);
}
StatusCode InitialAndConditions(const InputDeck* deck, Domain* domain,
                                const std::map<int, StepDefinition*>* steps,
                                Diagnostic* d) {
  std::vector<int> participating(0);
  const auto participation_status =
      CollectThermalNodeIds(domain, &participating, d);
  if (participation_status != StatusCode::kOk) return participation_status;
  int condition_id = 0;
  std::map<std::tuple<int, int, int>, int> boundary_lines;
  for (const auto* b : *deck->Blocks()) {
    const bool initial = *b->keyword == "INITIALCONDITIONS";
    const bool flux = *b->keyword == "CFLUX";
    const bool film = *b->keyword == "FILM" || *b->keyword == "SFILM";
    const bool radiation =
        *b->keyword == "RADIATE" || *b->keyword == "SRADIATE";
    const bool distributed =
        radiation || film || *b->keyword == "DFLUX" || *b->keyword == "DSFLUX";
    if (!initial && !flux && !distributed && *b->keyword != "BOUNDARY")
      continue;
    const auto operation = Equal(b->Option("OP"), "NEW")
                               ? ConditionOperation::kNew
                               : ConditionOperation::kMod;
    StepDefinition* step =
        b->step == 0 ? nullptr : steps->find(b->step)->second;
    if (!initial && step) {
      const auto status =
          distributed
              ? step->RecordThermalLoadBlock(
                    radiation ? (*b->keyword == "SRADIATE"
                                     ? ThermalLoadFamily::kSurfaceRadiation
                                     : ThermalLoadFamily::kElementRadiation)
                    : film    ? (*b->keyword == "SFILM"
                                     ? ThermalLoadFamily::kSurfaceFilm
                                     : ThermalLoadFamily::kElementFilm)
                              : (*b->keyword == "DSFLUX"
                                     ? ThermalLoadFamily::kSurfaceFlux
                                     : ThermalLoadFamily::kElementFlux),
                    operation, b->line, d)
              : step->RecordConditionBlock(flux, operation, b->line, d);
      if (status != StatusCode::kOk) return At(status, deck, b, b->line, d);
    }
    const Amplitude* amplitude = nullptr;
    const char* amplitude_key = radiation ? "AMBIENTAMPLITUDE"
                                : film    ? "SINKAMPLITUDE"
                                          : "AMPLITUDE";
    if (b->HasOption(amplitude_key)) {
      const auto status =
          domain->FindAmplitudeByName(b->Option(amplitude_key), &amplitude, d);
      if (status != StatusCode::kOk) return At(status, deck, b, b->line, d);
    }
    if (distributed) {
      for (const auto* row : *b->rows) {
        const auto status = DistributedRow(deck, domain, b, row, amplitude,
                                           step, ++condition_id, d);
        if (status != StatusCode::kOk) return status;
      }
      continue;
    }
    for (std::size_t row_index = 0; row_index < b->rows->size();) {
      const auto* row = b->rows->at(row_index);
      if (initial) {
        std::vector<int> initial_ids(0);
        std::vector<double> initial_values(0);
        std::size_t next = row_index;
        int source_line = 0;
        auto status = ReadThermalInitialRow(
            deck, domain, b, row_index, &participating, &initial_ids,
            &initial_values, &next, &source_line, d);
        if (status != StatusCode::kOk) return status;
        row_index = next;
        const bool uniform =
            Equal(b->Option("SECTIONSPECIFICATION"), "UNIFORM");
        for (int id : initial_ids) {
          if (!std::binary_search(participating.begin(), participating.end(),
                                  id))
            continue;
          const Node* node = nullptr;
          status = domain->FindNode(id, &node, d);
          if (status != StatusCode::kOk)
            return At(status, deck, b, source_line, d);
          Vector* raw_values = nullptr;
          status =
              Vector::Create(node->temperature_dof_count(), &raw_values, d);
          std::unique_ptr<Vector> values(raw_values);
          if (status != StatusCode::kOk)
            return At(status, deck, b, source_line, d);
          for (int j = 0; j < node->temperature_dof_count(); ++j) {
            status = values->Set(j, initial_values[uniform ? 0 : j], d);
            if (status != StatusCode::kOk)
              return At(status, deck, b, source_line, d);
          }
          ThermalInitialTemperature* raw = nullptr;
          status = ThermalInitialTemperature::Create(
              id, values.get(), deck->SourceFile(), source_line, &raw, d);
          std::unique_ptr<ThermalInitialTemperature> owner(raw);
          if (status == StatusCode::kOk)
            status = domain->RegisterThermalInitialTemperature(raw, d);
          if (status != StatusCode::kOk)
            return At(status, deck, b, source_line, d);
          owner.release();
        }
        continue;
      }
      ++row_index;
      std::vector<std::string> f(0);
      Fields(row->text->c_str(), &f);
      std::vector<int> ids(0);
      std::string scope(""), name("");
      int direct = 0;
      auto status = Target(deck, domain, b, row->line, f[0].c_str(), false,
                           &ids, &scope, &name, &direct, d);
      if (status != StatusCode::kOk) return status;
      double value = 0;
      if (flux)
        Real(f[2].c_str(), &value);
      else if (f.size() > 3 && !f[3].empty())
        Real(f[3].c_str(), &value);
      int first = 11, last = 11;
      if (!f[1].empty() && !Equal(f[1].c_str(), "0"))
        Integer(f[1].c_str(), 11, 30, &first);
      last = first;
      if (!flux && f.size() > 2 && !f[2].empty())
        Integer(f[2].c_str(), first, 30, &last);
      if (!flux)
        for (int id : ids) {
          if (!std::binary_search(participating.begin(), participating.end(),
                                  id))
            continue;
          for (int label = first; label <= last; ++label)
            if (!boundary_lines
                     .emplace(std::make_tuple(b->step, id, label), row->line)
                     .second)
              return Error(deck->SourceFile(), row->line, b->step,
                           "Duplicate temperature boundary", d);
        }
      ThermalConditionSpec spec;
      spec.id = ++condition_id;
      spec.kind = flux ? ThermalConditionKind::kConcentratedFlux
                       : ThermalConditionKind::kTemperature;
      spec.target = direct ? ThermalConditionTarget::kNode
                           : ThermalConditionTarget::kNodeSet;
      spec.target_id = direct;
      spec.scope = scope.c_str();
      spec.name = name.c_str();
      spec.magnitude = value;
      spec.first = first;
      spec.last = last;
      spec.amplitude = amplitude;
      spec.source_file = deck->SourceFile();
      spec.source_line = row->line;
      spec.step = b->step;
      spec.operation = operation;
      ThermalConditionDefinition* raw = nullptr;
      status = ThermalConditionDefinition::Create(&spec, &raw, d);
      std::unique_ptr<ThermalConditionDefinition> owner(raw);
      if (status == StatusCode::kOk)
        status = domain->RegisterThermalCondition(raw, d);
      if (status != StatusCode::kOk) return At(status, deck, b, row->line, d);
      owner.release();
      if (step) {
        status = step->AddThermalCondition(raw, d);
        if (status != StatusCode::kOk) return At(status, deck, b, row->line, d);
      }
    }
  }
  return StatusCode::kOk;
}
}  // namespace
StatusCode ReadThermalInitialRow(const InputDeck* deck, const Domain* domain,
                                 const InputBlock* block, std::size_t first_row,
                                 const std::vector<int>* participating,
                                 std::vector<int>* ids,
                                 std::vector<double>* values,
                                 std::size_t* next_row, int* source_line,
                                 Diagnostic* d) {
  if (!deck || !domain || !block || !participating || !ids || !values ||
      !next_row || !source_line || *block->keyword != "INITIALCONDITIONS" ||
      first_row >= block->rows->size())
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation,
                         "Thermal initial row inputs", d);
  const auto* first = block->rows->at(first_row);
  try {
    std::vector<std::string> fields(0);
    Fields(first->text->c_str(), &fields);
    if (fields.size() < 2 || fields.size() > 8)
      return Error(
          deck->SourceFile(), first->line, 0,
          "Initial temperature target row requires one to seven values", d);
    std::vector<int> parsed_ids(0);
    std::string scope(""), name("");
    int direct = 0;
    auto status = Target(deck, domain, block, first->line, fields[0].c_str(),
                         false, &parsed_ids, &scope, &name, &direct, d);
    if (status != StatusCode::kOk) return status;
    if (direct && !std::binary_search(participating->begin(),
                                      participating->end(), direct))
      return Error(deck->SourceFile(), first->line, 0,
                   "Initial temperature node is not participating", d);
    const bool uniform =
        Equal(block->Option("SECTIONSPECIFICATION"), "UNIFORM");
    int ns = 0;
    for (int id : parsed_ids) {
      if (!std::binary_search(participating->begin(), participating->end(), id))
        continue;
      const Node* node = nullptr;
      status = domain->FindNode(id, &node, d);
      if (status != StatusCode::kOk)
        return At(status, deck, block, first->line, d);
      const int count = node->temperature_dof_count();
      if (ns != 0 && ns != count && !uniform)
        return Error(deck->SourceFile(), first->line, 0,
                     "Nonuniform initial target has different section counts",
                     d);
      ns = count;
    }
    std::vector<double> parsed_values(0);
    parsed_values.reserve(19);
    for (std::size_t i = 1; i < fields.size(); ++i) {
      double value = 0;
      if (!Real(fields[i].c_str(), &value))
        return Error(deck->SourceFile(), first->line, 0,
                     "Invalid initial temperature", d);
      parsed_values.push_back(value);
    }
    if (!uniform && ns == 0 && parsed_values.size() == 7)
      return Error(deck->SourceFile(), first->line, 0,
                   "Initial-temperature continuation cannot be resolved for a "
                   "target without participating nodes",
                   d, StatusCode::kUnsupported);
    const int first_count = uniform ? 1 : std::min(ns, 7);
    if ((uniform || ns > 0) &&
        parsed_values.size() != static_cast<std::size_t>(first_count))
      return Error(deck->SourceFile(), first->line, 0,
                   "Initial temperature first row count mismatch", d);
    std::size_t next = first_row + 1;
    while (!uniform && ns > static_cast<int>(parsed_values.size())) {
      if (next >= block->rows->size())
        return Error(deck->SourceFile(), first->line, 0,
                     "Missing initial temperature continuation", d);
      const auto* continuation = block->rows->at(next);
      Fields(continuation->text->c_str(), &fields);
      const int count =
          std::min(8, ns - static_cast<int>(parsed_values.size()));
      if (fields.size() != static_cast<std::size_t>(count))
        return Error(deck->SourceFile(), continuation->line, 0,
                     "Initial temperature continuation count mismatch", d);
      for (const auto& field : fields) {
        double value = 0;
        if (!Real(field.c_str(), &value))
          return Error(deck->SourceFile(), continuation->line, 0,
                       "Initial temperature continuation requires values only",
                       d);
        parsed_values.push_back(value);
      }
      ++next;
    }
    ids->swap(parsed_ids);
    values->swap(parsed_values);
    *next_row = next;
    *source_line = first->line;
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return Error(deck->SourceFile(), first->line, 0,
                 "Thermal initial row allocation", d,
                 StatusCode::kAllocationFailure);
  }
}
StatusCode CountThermalTargetParticipation(
    const InputDeck* deck, const Domain* domain, const InputBlock* block,
    const InputLine* row, const std::vector<int>* participating, int* total,
    int* included, Diagnostic* d) {
  if (!deck || !domain || !block || !row || !participating || !total ||
      !included || total == included)
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation,
                         "Thermal participation query inputs", d);
  try {
    std::vector<std::string> fields(0);
    Fields(row->text->c_str(), &fields);
    if (fields.empty())
      return Error(deck->SourceFile(), row->line, block->step,
                   "Empty thermal target row", d);
    std::vector<int> ids(0);
    std::string scope(""), name("");
    int direct = 0;
    const auto status =
        Target(deck, domain, block, row->line, fields[0].c_str(), false, &ids,
               &scope, &name, &direct, d);
    if (status != StatusCode::kOk) return status;
    int active = 0;
    for (int id : ids)
      if (std::binary_search(participating->begin(), participating->end(), id))
        ++active;
    if (direct && active == 0)
      return Error(deck->SourceFile(), row->line, block->step,
                   "Direct thermal target is not participating", d);
    *total = static_cast<int>(ids.size());
    *included = active;
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return Error(deck->SourceFile(), row->line, block->step,
                 "Thermal participation query allocation", d,
                 StatusCode::kAllocationFailure);
  }
}
StatusCode BuildThermalDomain(const InputDeck* deck, Domain** out,
                              Diagnostic* d) {
  if (!out)
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation, "Thermal Domain output",
                         d);
  *out = nullptr;
  if (!deck)
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation, "Thermal Domain input",
                         d);
  try {
    Domain* raw = nullptr;
    auto status = Domain::Create(&raw, d);
    std::unique_ptr<Domain> domain(raw);
    if (status != StatusCode::kOk) return status;
    const char *part = "", *assembly = "", *instance = "";
    for (const auto* b : *deck->Blocks()) {
      if (*b->keyword == "PART") part = b->Option("NAME");
      if (*b->keyword == "ASSEMBLY") assembly = b->Option("NAME");
      if (*b->keyword == "INSTANCE") {
        if (!Equal(b->Option("PART"), part))
          return Error(deck->SourceFile(), b->line, 0,
                       "Undefined Instance Part", d);
        instance = b->Option("NAME");
      }
    }
    status = domain->SetInputIdentity(deck->SourceFile(), part, assembly,
                                      instance, d);
    if (status != StatusCode::kOk) return status;
    std::map<int, int> node_lines, element_lines;
    std::map<int, const InputBlock*> element_blocks;
    std::map<int, const InputElement*> records;
    for (const auto* b : *deck->Blocks()) {
      if (*b->keyword == "NODE") {
        for (const auto* row : *b->rows) {
          std::vector<std::string> f(0);
          Fields(row->text->c_str(), &f);
          int id = 0;
          Integer(f[0].c_str(), 1, 999999999, &id);
          if (!node_lines.emplace(id, row->line).second)
            return Error(deck->SourceFile(), row->line, 0, "Duplicate node ID",
                         d);
        }
      } else if (*b->keyword == "ELEMENT") {
        for (const auto* record : *b->elements) {
          if (!element_lines.emplace(record->id, record->line).second)
            return Error(deck->SourceFile(), record->line, 0,
                         "Duplicate element ID", d);
          element_blocks.emplace(record->id, b);
          records.emplace(record->id, record);
        }
      }
    }
    if (node_lines.empty() || records.empty())
      return Error(deck->SourceFile(), 1, 0,
                   "Thermal model needs nodes and elements", d);
    status = FilmProperties(deck, domain.get(), d);
    if (status != StatusCode::kOk) return status;
    status = Materials(deck, domain.get(), d);
    if (status != StatusCode::kOk) return status;
    status = abaqus::BuildInputIndexSets(deck, domain.get(), &node_lines,
                                         &element_lines, d);
    if (status != StatusCode::kOk) return status;
    int amplitude_id = 0;
    std::map<int, StepDefinition*> steps;
    bool transient = false;
    for (const auto* b : *deck->Blocks()) {
      if (*b->keyword == "AMPLITUDE") {
        Amplitude* raw_amplitude = nullptr;
        status = abaqus::BuildAmplitude(deck->SourceFile(), b, ++amplitude_id,
                                        &raw_amplitude, d);
        std::unique_ptr<Amplitude> owner(raw_amplitude);
        if (status == StatusCode::kOk)
          status = domain->RegisterAmplitude(raw_amplitude, d);
        if (status != StatusCode::kOk) return At(status, deck, b, b->line, d);
        owner.release();
      } else if (*b->keyword == "STEP") {
        const InputBlock* heat = nullptr;
        for (const auto* candidate : *deck->Blocks())
          if (candidate->step == b->step &&
              *candidate->keyword == "HEATTRANSFER")
            heat = candidate;
        HeatTransferStepParameters parameters;
        status = abaqus::ReadThermalStep(deck->SourceFile(), b, heat,
                                         &parameters, d);
        if (status != StatusCode::kOk) return status;
        transient = transient || !parameters.steady;
        StepDefinition* raw_step = nullptr;
        status = StepDefinition::CreateHeatTransfer(b->step, b->Option("NAME"),
                                                    &parameters, &raw_step, d);
        std::unique_ptr<StepDefinition> owner(raw_step);
        if (status == StatusCode::kOk)
          status = domain->RegisterStepDefinition(raw_step, d);
        if (status != StatusCode::kOk) return At(status, deck, b, b->line, d);
        // Domain must own alone before this auxiliary lookup can allocate.
        owner.release();
        steps.emplace(b->step, raw_step);
      }
    }
    status = Sections(deck, domain.get(), transient, &element_blocks, &records,
                      &node_lines, d);
    if (status == StatusCode::kOk) status = Metadata(deck, domain.get(), d);
    if (status == StatusCode::kOk)
      status = InitialAndConditions(deck, domain.get(), &steps, d);
    if (status != StatusCode::kOk) return status;
    *out = domain.release();
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return Error(deck->SourceFile(), 0, 0, "Thermal model allocation", d,
                 StatusCode::kAllocationFailure);
  }
}
}  // namespace fesa
