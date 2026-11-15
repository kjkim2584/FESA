#include "adapters/abaqus/model_builder.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <new>
#include <string>
#include <vector>

#include "adapters/abaqus/amplitude_reader.h"
#include "adapters/abaqus/dynamic_settings_reader.h"
#include "adapters/abaqus/input_syntax.h"
#include "adapters/abaqus/thermal_model_builder.h"
#include "fesa/elements/beam3d.h"
#include "fesa/elements/hexa8.h"
#include "fesa/elements/incremental_tetra4.h"
#include "fesa/elements/incremental_tetra10.h"
#include "fesa/elements/incremental_wedge6.h"
#include "fesa/elements/incremental_pyramid5.h"
#include "fesa/elements/incremental_hexa8.h"
#include "fesa/elements/mitc4_shell.h"
#include "fesa/elements/pyramid5.h"
#include "fesa/elements/tetra10.h"
#include "fesa/elements/tetra4.h"
#include "fesa/elements/wedge6.h"
#include "fesa/loads/beam_line_load.h"
#include "fesa/model/elastic_material.h"

namespace fesa {
namespace {
using abaqus::Canonical;
using abaqus::Equal;
using abaqus::Error;
using abaqus::Fields;
using abaqus::Integer;
using abaqus::Real;

// Adapter-only accumulation, never a solver model value member.
struct SetRows {
  SetRows() = default;
  ~SetRows() { delete ids; }
  SetRows(const SetRows&) = delete;
  SetRows& operator=(const SetRows&) = delete;
  SetRows(SetRows&&) = delete;
  SetRows& operator=(SetRows&&) = delete;
  const InputBlock* first = nullptr;
  std::vector<int>* ids = nullptr;
};
void Scope(const InputBlock* block, std::string* out) {
  *out = block->scope_kind == InputScope::kPart ? "PART." : "ASSEMBLY.";
  *out += *block->scope;
}
void SetKey(const InputBlock* block, std::string* out) {
  std::string scope("");
  Scope(block, &scope);
  scope.push_back(':');
  const char* kind =
      *block->keyword == "ELEMENT" ? "ELSET" : block->keyword->c_str();
  scope.append(kind);
  scope.push_back(':');
  scope.append(block->Option(kind));
  Canonical(scope.c_str(), out);
}
StatusCode Context(StatusCode status, const InputDeck* deck,
                   const InputBlock* block, int line, Diagnostic* d,
                   int node = 0, int element = 0) {
  if (status == StatusCode::kOk) return status;
  if (status == StatusCode::kInvalidArgument) {
    status = StatusCode::kInputError;
    if (d) d->code = status;
  }
  if (d) {
    d->source_file = deck->SourceFile();
    d->duplicate_boundary = false;
    d->condition_dof = 0;
    d->path_kind = DiagnosticPathKind::kLiteralFile;
    d->source_line = line;
    d->first_source_line = 0;
    d->step = block->step;
    d->node_id = node;
    d->element_id = element;
  }
  return status;
}
StatusCode Duplicate(const InputDeck* deck, const InputBlock* b, int line,
                     int first, const char* what, Diagnostic* d, int node = 0,
                     int element = 0) {
  if (d) {
    d->node_id = node;
    d->element_id = element;
  }
  return Error(deck->SourceFile(), line, b->step, what, d,
               StatusCode::kInputError, first);
}
StatusCode Missing(const InputDeck* deck, const InputBlock* b, int line,
                   const char* what, Diagnostic* d, int node = 0,
                   int element = 0) {
  if (d) {
    d->node_id = node;
    d->element_id = element;
  }
  return Error(deck->SourceFile(), line, b->step, what, d);
}
bool MemberId(const char* text, bool assembly, const char* explicit_instance,
              const char* instance, int* id) {
  std::string target(text);
  const auto dot = target.find('.');
  if (assembly && !*explicit_instance) {
    return dot != std::string::npos &&
           Equal(target.substr(0, dot).c_str(), instance) &&
           Integer(target.substr(dot + 1).c_str(), 1, 999999999, id);
  }
  if (dot != std::string::npos) return false;
  return Integer(text, 1, 999999999, id);
}

// Select one physical family before any Node is constructed. All family
// and distributed-load errors precede numerical work and output opening.
StatusCode ModelDofs(const InputDeck* deck, int* out, Diagnostic* d) {
  const char* selected = nullptr;
  for (const auto* b : *deck->Blocks()) {
    if (*b->keyword != "ELEMENT") continue;
    const char* type = b->Option("TYPE");
    if (Equal(type, "FESA_C3D10R") || Equal(type, "C3D10R"))
      return Error(deck->SourceFile(), b->line, b->step,
                   "FESA_C3D10R retired; use C3D10", d,
                   StatusCode::kUnsupported);
    const char* current = Equal(type, "B33")  ? "B33"
                          : Equal(type, "S4") ? "S4"
                                              : "SOLID";
    if (selected && !Equal(selected, current))
      return Error(deck->SourceFile(), b->line, 0,
                   "Mixed element families are unsupported", d,
                   StatusCode::kUnsupported);
    selected = current;
  }
  const bool solid = selected && Equal(selected, "SOLID");
  const bool shell = selected && Equal(selected, "S4");
  *out = solid ? 3 : 6;
  for (const auto* b : *deck->Blocks()) {
    if ((solid && (*b->keyword == "BEAMGENERALSECTION" ||
                   *b->keyword == "SHELLSECTION" || *b->keyword == "DLOAD")) ||
        (shell && (*b->keyword == "BEAMGENERALSECTION" ||
                   *b->keyword == "SOLIDSECTION" || *b->keyword == "DLOAD")) ||
        (!solid && !shell &&
         (*b->keyword == "SOLIDSECTION" || *b->keyword == "SHELLSECTION")))
      return Error(deck->SourceFile(), b->line, b->step,
                   "Unsupported section/load for element family", d,
                   StatusCode::kUnsupported);
  }
  return StatusCode::kOk;
}

StatusCode RegisterNodesAndMaterials(const InputDeck* deck, Domain* domain,
                                     int node_dofs,
                                     std::map<int, int>* node_lines,
                                     std::map<int, int>* element_lines,
                                     Diagnostic* d) {
  std::map<std::string, int> material_lines;
  const auto* blocks = deck->Blocks();
  for (std::size_t bi = 0; bi < blocks->size(); ++bi) {
    const InputBlock* b = blocks->at(bi);
    if (*b->keyword == "ELEMENT") {
      for (const auto* row : *b->elements) {
        const auto inserted = element_lines->emplace(row->id, row->line);
        if (!inserted.second)
          return Duplicate(deck, b, row->line, inserted.first->second,
                           "Duplicate element", d, 0, row->id);
      }
    } else if (*b->keyword == "NODE") {
      for (const auto* row : *b->rows) {
        std::vector<std::string> f(0);
        Fields(row->text->c_str(), &f);
        int id = 0;
        Integer(f[0].c_str(), 1, 999999999, &id);
        const auto inserted = node_lines->emplace(id, row->line);
        if (!inserted.second)
          return Duplicate(deck, b, row->line, inserted.first->second,
                           "Duplicate node", d, id, 0);
        double xyz[3] = {};
        for (std::size_t i = 1; i < f.size(); ++i)
          if (!f[i].empty()) Real(f[i].c_str(), &xyz[i - 1]);
        Node* object = nullptr;
        auto status =
            Node::Create(id, xyz[0], xyz[1], xyz[2], node_dofs, &object, d);
        std::unique_ptr<Node> owner(object);
        if (status == StatusCode::kOk) status = domain->RegisterNode(object, d);
        if (status != StatusCode::kOk)
          return Context(status, deck, b, row->line, d, id);
        owner.release();
      }
    } else if (*b->keyword == "MATERIAL") {
      std::string key("");
      Canonical(b->Option("NAME"), &key);
      const auto inserted = material_lines.emplace(key, b->line);
      if (!inserted.second)
        return Duplicate(deck, b, b->line, inserted.first->second,
                         "Duplicate material", d);
      const InputBlock *elastic = nullptr, *density = nullptr,
                       *damping = nullptr;
      for (std::size_t property = bi + 1; property < blocks->size();
           ++property) {
        const auto* current = blocks->at(property);
        if (*current->keyword != "ELASTIC" && *current->keyword != "DENSITY" &&
            *current->keyword != "DAMPING")
          break;
        const InputBlock** slot =
            *current->keyword == "ELASTIC" ? &elastic
            : *current->keyword == "DENSITY" ? &density : &damping;
        if (*slot)
          return Duplicate(deck, current, current->line, (*slot)->line,
                           "Duplicate material property", d);
        *slot = current;
      }
      if (!elastic)
        return Missing(deck, b, b->line, "Material requires Elastic", d);
      std::vector<std::string> f(0);
      Fields(elastic->rows->front()->text->c_str(), &f);
      double young = 0, poisson = 0;
      Real(f[0].c_str(), &young);
      Real(f[1].c_str(), &poisson);
      double rho = 0;
      if (density) {
        Fields(density->rows->front()->text->c_str(), &f);
        Real(f[0].c_str(), &rho);
      }
      double rayleigh_alpha = 0, rayleigh_beta = 0;
      if (damping) {
        if (damping->HasOption("ALPHA"))
          Real(damping->Option("ALPHA"), &rayleigh_alpha);
        if (damping->HasOption("BETA"))
          Real(damping->Option("BETA"), &rayleigh_beta);
      }
      ElasticMaterial* object = nullptr;
      StatusCode status = ElasticMaterial::Create(
          static_cast<int>(material_lines.size()), b->Option("NAME"), young,
          poisson, density != nullptr, rho, rayleigh_alpha, rayleigh_beta,
          &object, d);
      std::unique_ptr<Material> owner(object);
      if (status == StatusCode::kOk)
        status = domain->RegisterMaterial(object, d);
      if (status != StatusCode::kOk)
        return Context(status, deck, elastic, elastic->rows->front()->line, d);
      owner.release();
    }
  }
  if (node_lines->empty())
    return Error(deck->SourceFile(), 1, 0, "Model has no nodes", d);
  return StatusCode::kOk;
}

StatusCode RegisterSets(const InputDeck* deck, Domain* domain,
                        const std::map<int, int>* node_lines,
                        const std::map<int, int>* element_lines,
                        Diagnostic* d) {
  const auto* blocks = deck->Blocks();
  std::map<std::string, std::unique_ptr<SetRows>> sets;
  for (const auto* b : *blocks) {
    const bool inline_elements =
        *b->keyword == "ELEMENT" && b->HasOption("ELSET");
    if (*b->keyword != "NSET" && *b->keyword != "ELSET" && !inline_elements)
      continue;
    const bool assembly = b->scope_kind == InputScope::kAssembly;
    if (assembly && b->HasOption("INSTANCE") &&
        !Equal(b->Option("INSTANCE"), domain->InstanceName()))
      return Missing(deck, b, b->line, "Undefined set Instance", d);
    std::string key("");
    SetKey(b, &key);
    auto found = sets.find(key);
    if (found == sets.end()) {
      std::unique_ptr<SetRows> rows(new SetRows);
      rows->first = b;
      rows->ids = new std::vector<int>(0);
      found = sets.emplace(key, std::move(rows)).first;
    }
    const auto* members = *b->keyword != "NSET" ? element_lines : node_lines;
    if (inline_elements) {
      for (const auto* row : *b->elements)
        found->second->ids->push_back(row->id);
      continue;
    }
    for (const auto* row : *b->rows) {
      std::vector<std::string> f(0);
      Fields(row->text->c_str(), &f);
      if (b->HasOption("GENERATE")) {
        int first = 0, last = 0, increment = 1;
        if (!MemberId(f[0].c_str(), assembly, b->Option("INSTANCE"),
                      domain->InstanceName(), &first) ||
            !MemberId(f[1].c_str(), assembly, b->Option("INSTANCE"),
                      domain->InstanceName(), &last))
          return Missing(deck, b, row->line,
                         "Incorrectly scoped GENERATE endpoints", d);
        if (f.size() > 2 && !f[2].empty())
          Integer(f[2].c_str(), 1, 999999999, &increment);
        for (int id = first;; id += increment) {
          if (members->find(id) == members->end())
            return Missing(deck, b, row->line, "Undefined generated set member",
                           d, *b->keyword == "NSET" ? id : 0,
                           *b->keyword == "ELSET" ? id : 0);
          found->second->ids->push_back(id);
          if (id == last) break;
        }
      } else {
        for (const auto& field : f) {
          int id = 0;
          if (!MemberId(field.c_str(), assembly, b->Option("INSTANCE"),
                        domain->InstanceName(), &id) ||
              members->find(id) == members->end())
            return Missing(deck, b, row->line,
                           "Undefined or incorrectly scoped set member", d,
                           *b->keyword == "NSET" ? id : 0,
                           *b->keyword == "ELSET" ? id : 0);
          found->second->ids->push_back(id);
        }
      }
    }
  }
  int set_id = 0;
  for (const auto& entry : sets) {
    const auto* b = entry.second->first;
    const auto* ids = entry.second->ids;
    std::string scope("");
    Scope(b, &scope);
    IndexSet* object = nullptr;
    StatusCode status = IndexSet::Create(
        ++set_id, scope.c_str(),
        b->Option(*b->keyword == "ELEMENT" ? "ELSET" : b->keyword->c_str()),
        *b->keyword != "NSET", static_cast<int>(ids->size()), ids->data(),
        &object, d);
    std::unique_ptr<IndexSet> owner(object);
    if (status == StatusCode::kOk) status = domain->RegisterIndexSet(object, d);
    if (status != StatusCode::kOk) return Context(status, deck, b, b->line, d);
    owner.release();
  }
  return StatusCode::kOk;
}

StatusCode RegisterSolidSectionsAndElements(const InputDeck* deck,
                                            Domain* domain, SolidConstitutiveModel model, Diagnostic* d) {
  std::map<int, const SolidSection*> assigned;
  std::map<int, int> section_lines;
  int section_id = 0;
  for (const auto* b : *deck->Blocks()) {
    if (*b->keyword != "SOLIDSECTION") continue;
    const Material* material = nullptr;
    auto status =
        domain->FindMaterialByName(b->Option("MATERIAL"), &material, d);
    if (status == StatusCode::kAllocationFailure)
      return Context(status, deck, b, b->line, d);
    if (status != StatusCode::kOk || !material->AsElastic())
      return Missing(deck, b, b->line, "Undefined solid elastic material", d);
    std::string scope("");
    Scope(b, &scope);
    const IndexSet* set = nullptr;
    status = domain->FindIndexSetByName(scope.c_str(), b->Option("ELSET"), true,
                                        &set, d);
    if (status == StatusCode::kAllocationFailure)
      return Context(status, deck, b, b->line, d);
    if (status != StatusCode::kOk)
      return Missing(deck, b, b->line, "Undefined solid section element set",
                     d);
    SolidSection* section = nullptr;
    status =
        SolidSection::Create(++section_id, material->AsElastic(), &section, d);
    std::unique_ptr<SolidSection> owner(section);
    if (status == StatusCode::kOk)
      status = domain->RegisterSolidSection(section, d);
    if (status != StatusCode::kOk) return Context(status, deck, b, b->line, d);
    owner.release();
    for (int i = 0; i < set->Size(); ++i) {
      const int id = set->Ids()[i];
      const auto inserted = section_lines.emplace(id, b->line);
      if (!inserted.second)
        return Duplicate(deck, b, b->line, inserted.first->second,
                         "Duplicate solid element section", d, 0, id);
      assigned.emplace(id, section);
    }
  }
  for (const auto* b : *deck->Blocks()) {
    if (*b->keyword != "ELEMENT") continue;
    for (const auto* row : *b->elements) {
      const int id = row->id;
      const Node* nodes[10] = {};
      for (int i = 0; i < row->node_count; ++i) {
        const int node_id = row->node_ids[i];
        if (domain->FindNode(node_id, &nodes[i], nullptr) != StatusCode::kOk)
          return Missing(deck, b, row->node_lines[i], "Undefined solid node", d,
                         node_id, id);
      }
      const auto section = assigned.find(id);
      if (section == assigned.end())
        return Missing(deck, b, row->line, "Solid element lacks section", d, 0,
                       id);
      Element* element = nullptr;
      auto status = StatusCode::kOk;
      if (row->node_count == 4) {
        if (model == SolidConstitutiveModel::kIncrementalElastic) {
          IncrementalTetra4* raw = nullptr;
          status = IncrementalTetra4::Create(id, nodes, 4, section->second,
              SolidFormulation::kTetra4Full, &raw, d);
          element = raw;
        } else {
          Tetra4* raw = nullptr;
          status = Tetra4::Create(id, nodes, 4, section->second,
              SolidFormulation::kTetra4Full, &raw, d);
          element = raw;
        }
      } else if (Equal(b->Option("TYPE"), "C3D8") ||
                 Equal(b->Option("TYPE"), "C3D8R")) {
        const auto formulation = Equal(b->Option("TYPE"), "C3D8R")
                                     ? SolidFormulation::kHexa8Reduced
                                     : SolidFormulation::kHexa8Full;
        if (model == SolidConstitutiveModel::kIncrementalElastic) {
          IncrementalHexa8* raw = nullptr;
          status = IncrementalHexa8::Create(id,nodes,8,section->second,
                                             formulation,&raw,d);
          element = raw;
        } else {
          Hexa8* raw = nullptr;
          status = Hexa8::Create(id,nodes,8,section->second,formulation,&raw,d);
          element = raw;
        }
      } else if (Equal(b->Option("TYPE"), "C3D6")) {
        if (model == SolidConstitutiveModel::kIncrementalElastic) {
          IncrementalWedge6* raw = nullptr;
          status = IncrementalWedge6::Create(id, nodes, 6, section->second,
              SolidFormulation::kWedge6, &raw, d);
          element = raw;
        } else {
          Wedge6* raw = nullptr;
          status = Wedge6::Create(id, nodes, 6, section->second,
                                  SolidFormulation::kWedge6, &raw, d);
          element = raw;
        }
      } else if (Equal(b->Option("TYPE"), "C3D5") ||
                 Equal(b->Option("TYPE"), "FESA_C3D5R")) {
        const auto formulation = Equal(b->Option("TYPE"), "FESA_C3D5R")
                                     ? SolidFormulation::kPyramid5Reduced
                                     : SolidFormulation::kPyramid5Full;
        if (model == SolidConstitutiveModel::kIncrementalElastic) {
          IncrementalPyramid5* raw = nullptr;
          status = IncrementalPyramid5::Create(id, nodes, 5, section->second,
              formulation, &raw, d);
          element = raw;
        } else {
          Pyramid5* raw = nullptr;
          status = Pyramid5::Create(id, nodes, 5, section->second, formulation,
                                    &raw, d);
          element = raw;
        }
      } else {
        if (model == SolidConstitutiveModel::kIncrementalElastic) {
          IncrementalTetra10* raw = nullptr;
          status = IncrementalTetra10::Create(id, nodes, 10, section->second,
              SolidFormulation::kTetra10Full, &raw, d);
          element = raw;
        } else {
          Tetra10* raw = nullptr;
          status = Tetra10::Create(id, nodes, 10, section->second,
                                   SolidFormulation::kTetra10Full, &raw, d);
          element = raw;
        }
      }
      std::unique_ptr<Element> owner(element);
      if (status == StatusCode::kOk)
        status = domain->RegisterElement(element, d);
      if (status != StatusCode::kOk)
        return Context(status, deck, b, row->line, d, d ? d->node_id : 0, id);
      owner.release();
    }
  }
  return StatusCode::kOk;
}

StatusCode RegisterShellSectionsAndElements(const InputDeck* deck,
                                            Domain* domain, Diagnostic* d) {
  std::map<int, const ShellSection*> assigned;
  std::map<int, int> section_lines;
  int section_id = 0;
  for (const auto* b : *deck->Blocks()) {
    if (*b->keyword != "SHELLSECTION") continue;
    const Material* material = nullptr;
    auto status =
        domain->FindMaterialByName(b->Option("MATERIAL"), &material, d);
    if (status == StatusCode::kAllocationFailure)
      return Context(status, deck, b, b->line, d);
    if (status != StatusCode::kOk || !material->AsElastic())
      return Missing(deck, b, b->line, "Undefined shell elastic material", d);
    std::string scope("");
    Scope(b, &scope);
    const IndexSet* set = nullptr;
    status = domain->FindIndexSetByName(scope.c_str(), b->Option("ELSET"), true,
                                        &set, d);
    if (status == StatusCode::kAllocationFailure)
      return Context(status, deck, b, b->line, d);
    if (status != StatusCode::kOk)
      return Missing(deck, b, b->line, "Undefined shell section element set",
                     d);
    std::vector<std::string> fields(0);
    Fields(b->rows->front()->text->c_str(), &fields);
    double thickness = 0.0;
    Real(fields[0].c_str(), &thickness);
    ShellSection* section = nullptr;
    status = ShellSection::Create(++section_id, thickness,
                                  material->AsElastic(), &section, d);
    std::unique_ptr<ShellSection> owner(section);
    if (status == StatusCode::kOk)
      status = domain->RegisterShellSection(section, d);
    if (status != StatusCode::kOk)
      return Context(status, deck, b, b->rows->front()->line, d);
    owner.release();
    for (int i = 0; i < set->Size(); ++i) {
      const int id = set->Ids()[i];
      const auto inserted = section_lines.emplace(id, b->line);
      if (!inserted.second)
        return Duplicate(deck, b, b->line, inserted.first->second,
                         "Duplicate shell element section", d, 0, id);
      assigned.emplace(id, section);
    }
  }
  for (const auto* b : *deck->Blocks()) {
    if (*b->keyword != "ELEMENT") continue;
    for (const auto* row : *b->rows) {
      std::vector<std::string> fields(0);
      Fields(row->text->c_str(), &fields);
      int id = 0;
      Integer(fields[0].c_str(), 1, 999999999, &id);
      const Node* nodes[4] = {};
      for (int i = 0; i < 4; ++i) {
        int node_id = 0;
        Integer(fields[i + 1].c_str(), 1, 999999999, &node_id);
        if (domain->FindNode(node_id, &nodes[i], nullptr) != StatusCode::kOk)
          return Missing(deck, b, row->line, "Undefined shell node", d, node_id,
                         id);
      }
      const auto section = assigned.find(id);
      if (section == assigned.end())
        return Missing(deck, b, row->line, "Shell element lacks section", d, 0,
                       id);
      Mitc4Shell* element = nullptr;
      auto status =
          Mitc4Shell::Create(id, nodes, 4, section->second, &element, d);
      std::unique_ptr<Mitc4Shell> owner(element);
      if (status == StatusCode::kOk)
        status = domain->RegisterElement(element, d);
      if (status != StatusCode::kOk)
        return Context(status, deck, b, row->line, d, 0, id);
      owner.release();
    }
  }
  return StatusCode::kOk;
}

StatusCode RegisterSectionsAndElements(const InputDeck* deck, Domain* domain,
                                       SolidConstitutiveModel model, Diagnostic* d) {
  const auto* blocks = deck->Blocks();
  if (domain->Nodes()->begin()->second->dof_count() == 3)
    return RegisterSolidSectionsAndElements(deck, domain, model, d);
  for (const auto* b : *blocks)
    if (*b->keyword == "ELEMENT" && Equal(b->Option("TYPE"), "S4"))
      return RegisterShellSectionsAndElements(deck, domain, d);
  std::map<int, const BeamSection*> assigned;
  std::map<int, int> section_lines;
  int section_id = 0;
  for (const auto* b : *blocks) {
    if (*b->keyword != "BEAMGENERALSECTION") continue;
    const Material* material = nullptr;
    const StatusCode material_status =
        domain->FindMaterialByName(b->Option("MATERIAL"), &material, d);
    if (material_status == StatusCode::kAllocationFailure)
      return Context(material_status, deck, b, b->line, d);
    if (material_status != StatusCode::kOk)
      return Missing(deck, b, b->line, "Undefined section material", d);
    if (!material->AsElastic())
      return Missing(deck, b, b->line, "Section requires elastic material", d);
    std::string scope("");
    Scope(b, &scope);
    const IndexSet* set = nullptr;
    const StatusCode set_status = domain->FindIndexSetByName(
        scope.c_str(), b->Option("ELSET"), true, &set, d);
    if (set_status == StatusCode::kAllocationFailure)
      return Context(set_status, deck, b, b->line, d);
    if (set_status != StatusCode::kOk)
      return Missing(deck, b, b->line, "Undefined section element set", d);
    std::vector<std::string> f(0);
    Fields(b->rows->front()->text->c_str(), &f);
    double c[5] = {0, 0, 0, 0, 0};
    for (int i = 0; i < 5; ++i) Real(f[i].c_str(), &c[i]);
    if (c[2] != 0)
      return Error(deck->SourceFile(), b->rows->front()->line, 0,
                   "Nonzero I12 is unsupported", d, StatusCode::kUnsupported);
    Vector* direction = nullptr;
    StatusCode status = Vector::Create(3, &direction, d);
    std::unique_ptr<Vector> direction_owner(direction);
    if (status != StatusCode::kOk) return Context(status, deck, b, b->line, d);
    double n[3] = {0, 0, -1};
    if (b->rows->size() == 2) {
      Fields(b->rows->at(1)->text->c_str(), &f);
      if (!f[0].empty())
        for (int i = 0; i < 3; ++i) Real(f[i].c_str(), &n[i]);
    }
    for (int i = 0; i < 3; ++i) {
      status = direction->Set(i, n[i], d);
      if (status != StatusCode::kOk)
        return Context(status, deck, b, b->line, d);
    }
    BeamSection* section = nullptr;
    status =
        BeamSection::Create(++section_id, material->AsElastic(), c[0], c[1],
                            c[2], c[3], c[4], direction, &section, d);
    std::unique_ptr<BeamSection> owner(section);
    if (status == StatusCode::kOk)
      status = domain->RegisterBeamSection(section, d);
    if (status != StatusCode::kOk)
      return Context(status, deck, b, b->rows->front()->line, d);
    owner.release();
    for (int i = 0; i < set->Size(); ++i) {
      const int id = set->Ids()[i];
      const auto inserted = section_lines.emplace(id, b->line);
      if (!inserted.second)
        return Duplicate(deck, b, b->line, inserted.first->second,
                         "Duplicate element section", d, 0, id);
      assigned.emplace(id, section);
    }
  }
  for (const auto* b : *blocks) {
    if (*b->keyword != "ELEMENT") continue;
    for (const auto* row : *b->elements) {
      const int id = row->id, a = row->node_ids[0], z = row->node_ids[1];
      const Node *first = nullptr, *second = nullptr;
      if (domain->FindNode(a, &first, nullptr) != StatusCode::kOk)
        return Missing(deck, b, row->line, "Undefined first endpoint", d, a,
                       id);
      if (domain->FindNode(z, &second, nullptr) != StatusCode::kOk)
        return Missing(deck, b, row->line, "Undefined second endpoint", d, z,
                       id);
      const auto section = assigned.find(id);
      if (section == assigned.end())
        return Missing(deck, b, row->line, "Element lacks section", d, 0, id);
      Beam3D* object = nullptr;
      StatusCode status =
          Beam3D::Create(id, first, second, section->second->GetMaterial(),
                         section->second, &object, d);
      std::unique_ptr<Beam3D> owner(object);
      if (status == StatusCode::kOk)
        status = domain->RegisterElement(object, d);
      if (status != StatusCode::kOk)
        return Context(status, deck, b, row->line, d, a, id);
      owner.release();
    }
  }
  return StatusCode::kOk;
}

StatusCode Conditions(const InputDeck* deck, Domain* domain, Diagnostic* d) {
  std::map<int, StepDefinition*> steps;
  std::map<std::string, int> names;
  int load_id = 0, boundary_id = 0;
  bool nonlinear = false;
  for (const auto* b : *deck->Blocks()) {
    if (*b->keyword == "STEP") {
      if (*b->Option("NAME")) {
        std::string key("");
        Canonical(b->Option("NAME"), &key);
        const auto inserted = names.emplace(key, b->line);
        if (!inserted.second)
          return Duplicate(deck, b, b->line, inserted.first->second,
                           "Duplicate step name", d);
      }
      nonlinear =
          nonlinear ||
          (b->HasOption("NLGEOM") &&
           (!*b->Option("NLGEOM") || Equal(b->Option("NLGEOM"), "YES")));
      double period = 1;
      bool dynamic = false;
      NonlinearStaticParameters parameters;
      for (const auto* other : *deck->Blocks()) {
        if (other->step == b->step && *other->keyword == "DYNAMIC") {
          dynamic = true;
          DynamicStepSettings* raw = nullptr;
          const auto resolved =
              abaqus::ReadDynamicSettings(deck, b->step, &raw, d);
          std::unique_ptr<DynamicStepSettings> settings(raw);
          if (resolved != StatusCode::kOk) return resolved;
          period = settings->Parameters()->period;
        }
        if (other->step != b->step || *other->keyword != "STATIC") continue;
        const bool valid = nonlinear
                               ? abaqus::NonlinearStaticData(other, &parameters)
                               : abaqus::StaticPeriod(other, &period);
        if (!valid)
          return Missing(deck, other, other->line, "Invalid STATIC settings",
                         d);
      }
      StepDefinition* step = nullptr;
      StatusCode status;
      if (dynamic) {
        status = nonlinear ? StepDefinition::CreateNonlinearDynamic(
                                 b->step, b->Option("NAME"), period, &step, d)
                           : StepDefinition::CreateDynamic(
                                 b->step, b->Option("NAME"), period, &step, d);
      } else if (nonlinear) {
        if (b->HasOption("INC") &&
            !Integer(b->Option("INC"), 1, std::numeric_limits<int>::max(),
                     &parameters.increment_limit))
          return Missing(deck, b, b->line, "Invalid STEP INC", d);
        parameters.extrapolate = !Equal(b->Option("EXTRAPOLATION"), "NO");
        status = StepDefinition::CreateNonlinearStatic(
            b->step, b->Option("NAME"), &parameters, &step, d);
      } else {
        status = StepDefinition::Create(b->step, b->Option("NAME"), period,
                                        &step, d);
      }
      std::unique_ptr<StepDefinition> owner(step);
      if (status == StatusCode::kOk)
        status = domain->RegisterStepDefinition(step, d);
      if (status != StatusCode::kOk)
        return Context(status, deck, b, b->line, d);
      owner.release();
      steps.emplace(b->step, step);
    }
    const bool load = *b->keyword == "CLOAD";
    if (!load && *b->keyword != "BOUNDARY") continue;
    const ConditionOperation operation = Equal(b->Option("OP"), "NEW")
                                             ? ConditionOperation::kNew
                                             : ConditionOperation::kMod;
    StepDefinition* step = b->step == 0 ? nullptr : steps.find(b->step)->second;
    if (step) {
      const int first_line =
          load ? step->LoadBlockLine() : step->BoundaryBlockLine();
      if (step->RecordConditionBlock(load, operation, b->line, nullptr) !=
          StatusCode::kOk)
        return Duplicate(deck, b, b->line, first_line,
                         "Mixed condition operations", d);
    }
    for (const auto* row : *b->rows) {
      std::vector<std::string> f(0);
      Fields(row->text->c_str(), &f);
      std::string scope("ASSEMBLY.");
      scope.append(domain->AssemblyName());
      std::string target = f[0];
      const auto dot = target.find('.');
      if (dot != std::string::npos) {
        if (!Equal(target.substr(0, dot).c_str(), domain->InstanceName()))
          return Missing(deck, b, row->line, "Undefined target Instance", d);
        scope.assign("PART.");
        scope.append(domain->PartName());
        target = target.substr(dot + 1);
      }
      int node_id = 0;
      const bool node_target = Integer(target.c_str(), 1, 999999999, &node_id);
      if (node_target) {
        const Node* node = nullptr;
        if (domain->FindNode(node_id, &node, nullptr) != StatusCode::kOk)
          return Missing(deck, b, row->line, "Undefined condition node", d,
                         node_id);
        // All direct node targets resolve to the unique Instance; keep their
        // original parsed qualification in OriginalTarget for diagnostics.
        // Fields trims exterior whitespace; InputDeck rows remain unchanged.
      } else {
        const IndexSet* set = nullptr;
        const StatusCode set_status = domain->FindIndexSetByName(
            scope.c_str(), target.c_str(), false, &set, d);
        if (set_status == StatusCode::kAllocationFailure)
          return Context(set_status, deck, b, row->line, d);
        if (set_status != StatusCode::kOk)
          return Missing(deck, b, row->line, "Undefined condition node set", d);
      }
      int first = 0, last = 0;
      double value = 0;
      Integer(f[1].c_str(), 1, 6, &first);
      last = first;
      if (load)
        Real(f[2].c_str(), &value);
      else {
        if (f.size() > 2 && !f[2].empty())
          Integer(f[2].c_str(), first, 6, &last);
        if (f.size() > 3 && !f[3].empty()) Real(f[3].c_str(), &value);
      }
      if (last > domain->Nodes()->begin()->second->dof_count()) {
        auto status = Error(deck->SourceFile(), row->line, b->step,
                            "Solid condition requires translation DOF 1-3", d,
                            StatusCode::kUnsupported);
        if (d) {
          d->condition_dof = std::max(first, 4);
          d->node_id = node_id;
        }
        return status;
      }
      for (int dof = first; dof <= last; ++dof) {
        StatusCode status = StatusCode::kOk;
        if (load) {
          LoadDefinition* object = nullptr;
          status = LoadDefinition::Create(
              ++load_id, scope.c_str(),
              node_target ? ConditionTarget::kNode : ConditionTarget::kNodeSet,
              node_id, node_target ? "" : target.c_str(), dof, value,
              deck->SourceFile(), row->line, b->step, operation, &object, d);
          std::unique_ptr<LoadDefinition> owner(object);
          if (status == StatusCode::kOk)
            status = object->SetOriginalTarget(f[0].c_str(), d);
          if (status == StatusCode::kOk && b->HasOption("AMPLITUDE")) {
            const Amplitude* amplitude = nullptr;
            status = domain->FindAmplitudeByName(b->Option("AMPLITUDE"),
                                                 &amplitude, d);
            if (status == StatusCode::kOk)
              status = object->BindDynamicAmplitude(amplitude, d);
          }
          if (status == StatusCode::kOk)
            status = domain->RegisterLoadDefinition(object, d);
          if (status == StatusCode::kOk) {
            owner.release();
            status = step->AddLoad(object, d);
          }
        } else {
          BoundaryDefinition* object = nullptr;
          status = BoundaryDefinition::Create(
              ++boundary_id, scope.c_str(),
              node_target ? ConditionTarget::kNode : ConditionTarget::kNodeSet,
              node_id, node_target ? "" : target.c_str(), dof, value,
              deck->SourceFile(), row->line, b->step, operation, &object, d);
          std::unique_ptr<BoundaryDefinition> owner(object);
          if (status == StatusCode::kOk)
            status = object->SetOriginalTarget(f[0].c_str(), d);
          if (status == StatusCode::kOk)
            status = domain->RegisterBoundaryDefinition(object, d);
          if (status == StatusCode::kOk) {
            owner.release();
            if (step) status = step->AddBoundary(object, d);
          }
        }
        if (status != StatusCode::kOk)
          return Context(status, deck, b, row->line, d, node_id);
      }
    }
  }
  return StatusCode::kOk;
}
}  // namespace

namespace abaqus {
StatusCode BuildInputIndexSets(const InputDeck* deck, Domain* domain,
    const std::map<int, int>* node_lines,
    const std::map<int, int>* element_lines, Diagnostic* d) {
  return RegisterSets(deck, domain, node_lines, element_lines, d);
}
}  // namespace abaqus

StatusCode ResolveSolidConstitutiveModel(const InputDeck* deck,
    SolidConstitutiveModel requested, SolidConstitutiveModel* out,
    Diagnostic* d) {
  if (!deck || !out)
    return ReportFailure(StatusCode::kInvalidArgument,
        DiagnosticStage::kValidation, "Solid constitutive resolution", d);
  if (requested != SolidConstitutiveModel::kStvk &&
      requested != SolidConstitutiveModel::kIncrementalElastic &&
      requested != SolidConstitutiveModel::kAutomatic)
    return ReportFailure(StatusCode::kInvalidArgument,
        DiagnosticStage::kValidation, "Solid constitutive selection", d);
  try {
    bool has_solid = false;
    const InputBlock* nonlinear_step = nullptr;
    for (const auto* block : *deck->Blocks()) {
      if (*block->keyword == "SOLIDSECTION") has_solid = true;
      if (*block->keyword == "STEP" && Equal(block->Option("NLGEOM"), "YES"))
        nonlinear_step = block;
    }
    if (has_solid && nonlinear_step &&
        requested == SolidConstitutiveModel::kStvk)
      return Error(deck->SourceFile(), nonlinear_step->line, nonlinear_step->step,
          "StVK nonlinear Solid is deprecated; use incremental_elastic", d,
          StatusCode::kUnsupported);
    *out = requested == SolidConstitutiveModel::kAutomatic
        ? (has_solid && nonlinear_step
               ? SolidConstitutiveModel::kIncrementalElastic
               : SolidConstitutiveModel::kStvk)
        : requested;
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return ReportFailure(StatusCode::kAllocationFailure,
        DiagnosticStage::kAllocation, "Solid constitutive resolution", d);
  }
}

StatusCode BuildDomain(const InputDeck* deck, Domain** out, Diagnostic* d,
                       SolidConstitutiveModel model) {
  if (!out)
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation, "BuildDomain output", d);
  *out = nullptr;
  if (!deck)
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation, "BuildDomain deck", d);
  for (const auto* block : *deck->Blocks())
    if (*block->keyword == "HEATTRANSFER") {
      if (model != SolidConstitutiveModel::kAutomatic)
        return Error(deck->SourceFile(), block->line, block->step,
                     "Mechanical constitutive option on thermal model", d,
                     StatusCode::kUnsupported);
      return BuildThermalDomain(deck, out, d);
    }
  try {
    const auto resolution = ResolveSolidConstitutiveModel(deck, model, &model, d);
    if (resolution != StatusCode::kOk) return resolution;
    if (model == SolidConstitutiveModel::kIncrementalElastic) {
      for (const auto* block : *deck->Blocks()) {
        if ((*block->keyword == "ELEMENT" &&
             !Equal(block->Option("TYPE"), "C3D4") &&
             !Equal(block->Option("TYPE"), "C3D10") &&
             !Equal(block->Option("TYPE"), "C3D6") &&
             !Equal(block->Option("TYPE"), "C3D5") &&
             !Equal(block->Option("TYPE"), "C3D8")) ||
            (*block->keyword == "STEP" &&
             !Equal(block->Option("NLGEOM"), "YES")) ||
            *block->keyword == "DYNAMIC")
          return Error(deck->SourceFile(), block->line, block->step,
              "Incremental elastic requires implemented full Solid nonlinear static",
              d, StatusCode::kUnsupported);
      }
    }
    bool has_dynamic = false;
    for (const auto* block : *deck->Blocks())
      if (*block->keyword == "DYNAMIC") has_dynamic = true;
    for (const auto* block : *deck->Blocks())
      if (*block->keyword == "INITIALCONDITIONS" && !has_dynamic)
        return Error(deck->SourceFile(), block->line, block->step,
                     "Initial velocity execution not yet connected", d,
                     StatusCode::kUnsupported);
    Domain* candidate = nullptr;
    StatusCode status = Domain::Create(&candidate, d);
    std::unique_ptr<Domain> owner(candidate);
    if (status != StatusCode::kOk) return status;
    const char *part = "", *assembly = "", *instance = "";
    const InputBlock* instance_block = nullptr;
    for (const auto* b : *deck->Blocks()) {
      if (*b->keyword == "PART") part = b->Option("NAME");
      if (*b->keyword == "ASSEMBLY") assembly = b->Option("NAME");
      if (*b->keyword == "INSTANCE") {
        instance = b->Option("NAME");
        instance_block = b;
      }
    }
    if (!instance_block || !Equal(instance_block->Option("PART"), part))
      return Error(deck->SourceFile(),
                   instance_block ? instance_block->line : 1, 0,
                   "Undefined Instance Part", d);
    status = candidate->SetInputIdentity(deck->SourceFile(), part, assembly,
                                         instance, d);
    int amplitude_id = 0;
    for (const auto* block : *deck->Blocks()) {
      if (status != StatusCode::kOk) return status;
      if (*block->keyword != "AMPLITUDE") continue;
      Amplitude* raw_amplitude = nullptr;
      status = abaqus::BuildAmplitude(deck->SourceFile(), block, ++amplitude_id,
                                      &raw_amplitude, d);
      std::unique_ptr<Amplitude> amplitude(raw_amplitude);
      if (status == StatusCode::kOk)
        status = candidate->RegisterAmplitude(raw_amplitude, d);
      if (status != StatusCode::kOk)
        return Context(status, deck, block, block->line, d);
      amplitude.release();
    }
    int node_dofs = 6;
    if (status == StatusCode::kOk) status = ModelDofs(deck, &node_dofs, d);
    std::map<int, int> node_lines, element_lines;
    if (status == StatusCode::kOk)
      status = RegisterNodesAndMaterials(deck, candidate, node_dofs,
                                         &node_lines, &element_lines, d);
    if (status == StatusCode::kOk)
      status = RegisterSets(deck, candidate, &node_lines, &element_lines, d);
    if (status == StatusCode::kOk)
      status = RegisterSectionsAndElements(deck, candidate, model, d);
    if (status == StatusCode::kOk) status = Conditions(deck, candidate, d);
    if (status == StatusCode::kOk) {
      bool finite_beam = false;
      for (const auto& element : *candidate->Elements())
        finite_beam = finite_beam ||
                      dynamic_cast<const Beam3D*>(element.second) != nullptr;
      finite_beam = finite_beam && !candidate->Steps()->empty() &&
                    candidate->Steps()->begin()->second->HasFiniteGeometry();
      for (const auto* block : *deck->Blocks()) {
        if (*block->keyword != "DLOAD") continue;
        if (finite_beam)
          return Error(deck->SourceFile(), block->line, block->step,
                       "Nonlinear Beam distributed loads unsupported", d,
                       StatusCode::kUnsupported);
        BeamLineLoadProgram* raw_program = nullptr;
        status = BeamLineLoadProgram::Create(deck, candidate, &raw_program, d);
        std::unique_ptr<BeamLineLoadProgram> program(raw_program);
        if (status == StatusCode::kOk)
          status = candidate->AttachLineLoadProgram(raw_program, d);
        if (status == StatusCode::kOk) program.release();
        break;
      }
    }
    if (status != StatusCode::kOk) return status;
    *out = owner.release();
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return Error(deck->SourceFile(), 0, 0, "BuildDomain allocation", d,
                 StatusCode::kAllocationFailure);
  }
}
}  // namespace fesa
