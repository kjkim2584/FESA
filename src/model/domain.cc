#include "fesa/model/domain.h"

#include <cctype>
#include <cstring>
#include <memory>
#include <new>

#include "fesa/loads/beam_line_load.h"
#include "fesa/model/film_property_table.h"
#include "model/model_storage.h"

namespace fesa {
namespace {

template <typename T>
StatusCode Insert(T* object, std::map<int, T*>* objects, Diagnostic* d) {
  if (object == nullptr || !model_storage::ValidId(object->Id())) {
    return model_storage::Invalid("Domain::Register invalid object", d);
  }
  try {
    if (!objects->emplace(object->Id(), object).second) {
      return model_storage::Invalid("Domain::Register duplicate ID", d);
    }
  } catch (const std::bad_alloc&) {
    return model_storage::AllocationFailure("Domain::Register", d);
  }
  return StatusCode::kOk;
}

template <typename T>
StatusCode Find(int id, const std::map<int, T*>* objects, const T** out,
                Diagnostic* d) {
  if (out == nullptr) return model_storage::Invalid("Domain::Find output", d);
  *out = nullptr;
  const auto found = objects->find(id);
  if (found == objects->end()) {
    return model_storage::Invalid("Domain::Find missing ID", d);
  }
  *out = found->second;
  return StatusCode::kOk;
}

// Length-prefix avoids scope/name separator collisions; preserve stored
// spelling.
void NameKey(char kind, const char* scope, const char* name, std::string* out) {
  std::string scope_key(""), name_key("");
  for (const char* p = scope; *p; ++p) {
    if (*p != ' ' && *p != '\t')
      scope_key +=
          *p >= 'A' && *p <= 'Z' ? static_cast<char>(*p - 'A' + 'a') : *p;
  }
  for (const char* p = name; *p; ++p) {
    if (*p != ' ' && *p != '\t')
      name_key +=
          *p >= 'A' && *p <= 'Z' ? static_cast<char>(*p - 'A' + 'a') : *p;
  }
  out->assign(1, kind);
  out->append(std::to_string(scope_key.size()));
  out->push_back(':');
  out->append(scope_key);
  out->append(name_key);
}

template <typename T>
StatusCode InsertNamed(T* object, char kind, const char* scope,
                        std::map<int, T*>* objects,
                        std::map<std::string, int>* names, Diagnostic* d) {
  if (object == nullptr || scope == nullptr) {
    return model_storage::Invalid("Domain::Register named object", d);
  }
  try {
    std::string key("");
    NameKey(kind, scope, object->Name(), &key);
    if (names->find(key) != names->end()) {
      return model_storage::Invalid("Domain::Register duplicate name", d);
    }
    const StatusCode status = Insert(object, objects, d);
    if (status != StatusCode::kOk) return status;
    try {
      names->emplace(key, object->Id());
    } catch (const std::bad_alloc&) {
      objects->erase(object->Id());
      return model_storage::AllocationFailure("Domain::Register name index", d);
    }
  } catch (const std::bad_alloc&) {
    return model_storage::AllocationFailure("Domain::Register name", d);
  }
  return StatusCode::kOk;
}

StatusCode FindName(char kind, const char* scope, const char* name,
                    const std::map<std::string, int>* names, int* id,
                    Diagnostic* d) {
  if (scope == nullptr || name == nullptr) {
    return model_storage::Invalid("Domain::FindName input", d);
  }
  try {
    std::string key("");
    NameKey(kind, scope, name, &key);
    const auto found = names->find(key);
    if (found == names->end()) {
      return model_storage::Invalid("Domain::FindName missing name", d);
    }
    *id = found->second;
  } catch (const std::bad_alloc&) {
    return model_storage::AllocationFailure("Domain::FindName", d);
  }
  return StatusCode::kOk;
}

}  // namespace

Domain::Domain(std::map<std::string, int>* names) : names_(names) {}
Domain::~Domain() {
  delete line_load_program_;
  delete thermal_constants_;
  for (auto entry : thermal_initial_temperatures_) delete entry.second;
  for (auto entry : steps_) delete entry.second;
  for (auto entry : loads_) delete entry.second;
  for (auto entry : boundaries_) delete entry.second;
  for (auto entry : thermal_conditions_) delete entry.second;
  for (auto entry : index_sets_) delete entry.second;
  for (auto entry : elements_) delete entry.second;
  for (auto entry : thermal_elements_) delete entry.second;
  for (auto entry : thermal_solid_sections_) delete entry.second;
  for (auto entry : thermal_shell_sections_) delete entry.second;
  for (auto entry : shell_sections_) delete entry.second;
  for (auto entry : solid_sections_) delete entry.second;
  for (auto entry : beam_sections_) delete entry.second;
  for (auto entry : film_properties_) delete entry.second;
  for (auto entry : amplitudes_) delete entry.second;
  for (auto entry : materials_) delete entry.second;
  for (auto entry : nodes_) delete entry.second;
  delete names_;
  delete[] source_file_;
  delete[] part_name_;
  delete[] assembly_name_;
  delete[] instance_name_;
}
StatusCode Domain::AttachLineLoadProgram(BeamLineLoadProgram* program,
                                         Diagnostic* d) {
  if (!program || line_load_program_ || program->SourceDomain() != this)
    return model_storage::Invalid("Domain line-load program ownership", d);
  line_load_program_ = program;
  return StatusCode::kOk;
}
const BeamLineLoadProgram* Domain::LineLoadProgram() const {
  return line_load_program_;
}
StatusCode Domain::Create(Domain** out, Diagnostic* d) {
  if (out == nullptr) return model_storage::Invalid("Domain::Create", d);
  *out = nullptr;
  try {
    std::unique_ptr<std::map<std::string, int>> names(
        new (std::nothrow) std::map<std::string, int>);
    if (names) *out = new (std::nothrow) Domain(names.get());
    if (*out == nullptr) {
      return model_storage::AllocationFailure("Domain::Create", d);
    }
    names.release();
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return model_storage::AllocationFailure("Domain::Create", d);
  }
}
StatusCode Domain::RegisterNode(Node* object, Diagnostic* d) {
  return Insert(object, &nodes_, d);
}
StatusCode Domain::RegisterThermalElement(ThermalElement* object, Diagnostic* d) {
  if (!object || !object->material() || object->node_count() <= 0 ||
      object->dofs_per_node() <= 0)
    return model_storage::Invalid("Domain thermal element definition", d);
  if (elements_.find(object->Id()) != elements_.end())
    return model_storage::Invalid("Domain duplicate mechanical/thermal ID", d);
  const auto material = materials_.find(object->material()->Id());
  if (material == materials_.end() || material->second != object->material())
    return model_storage::Invalid("Domain thermal element material identity", d);
  for (int i = 0; i < object->node_count(); ++i) {
    const Node* node = object->GetNode(i);
    if (!node || node->temperature_dof_count() != object->dofs_per_node())
      return model_storage::Invalid("Domain thermal element node layout", d);
    const auto registered = nodes_.find(node->Id());
    if (registered == nodes_.end() || registered->second != node)
      return model_storage::Invalid("Domain thermal element node identity", d);
  }
  return Insert(object, &thermal_elements_, d);
}
StatusCode Domain::SetThermalPhysicalConstants(
    const ThermalPhysicalConstants* values, Diagnostic* d) {
  if (thermal_constants_)
    return model_storage::Invalid("Domain duplicate thermal constants", d);
  return CloneThermalPhysicalConstants(values, &thermal_constants_, d);
}
const ThermalPhysicalConstants* Domain::ThermalConstants() const {
  return thermal_constants_;
}
StatusCode Domain::RegisterThermalInitialTemperature(
    ThermalInitialTemperature* object, Diagnostic* d) {
  if (!object) return model_storage::Invalid("Domain thermal IC", d);
  const auto found = nodes_.find(object->Id());
  if (found == nodes_.end() ||
      found->second->temperature_dof_count() == 0 ||
      found->second->temperature_dof_count() != object->values()->Size())
    return model_storage::Invalid("Domain thermal IC node/components", d);
  return Insert(object, &thermal_initial_temperatures_, d);
}
StatusCode Domain::FindThermalInitialTemperature(
    int id, const ThermalInitialTemperature** out, Diagnostic* d) const {
  return Find(id, &thermal_initial_temperatures_, out, d);
}
const std::map<int, ThermalInitialTemperature*>*
Domain::ThermalInitialTemperatures() const {
  return &thermal_initial_temperatures_;
}
StatusCode Domain::RegisterThermalCondition(ThermalConditionDefinition* object,
                                            Diagnostic* d) {
  if (!object) return model_storage::Invalid("Domain thermal condition", d);
  const auto* amplitude = object->specification()->amplitude;
  if (amplitude) {
    const auto found = amplitudes_.find(amplitude->Id());
    if (found == amplitudes_.end() || found->second != amplitude)
      return model_storage::Invalid("Domain thermal amplitude identity", d);
  }
  const auto* spec = object->specification();
  if (spec->film_amplitude) {
    const auto found = amplitudes_.find(spec->film_amplitude->Id());
    if (found == amplitudes_.end() || found->second != spec->film_amplitude)
      return model_storage::Invalid("Domain film amplitude identity", d);
  }
  if (spec->film_property) {
    const auto found = film_properties_.find(spec->film_property->Id());
    if (found == film_properties_.end() || found->second != spec->film_property)
      return model_storage::Invalid("Domain film property identity", d);
  }
  return Insert(object, &thermal_conditions_, d);
}
StatusCode Domain::FindThermalCondition(
    int id, const ThermalConditionDefinition** out, Diagnostic* d) const {
  return Find(id, &thermal_conditions_, out, d);
}
const std::map<int, ThermalConditionDefinition*>*
Domain::ThermalConditions() const {
  return &thermal_conditions_;
}
StatusCode Domain::FindThermalElement(int id, const ThermalElement** out,
                                     Diagnostic* d) const {
  return Find(id, &thermal_elements_, out, d);
}
const std::map<int, ThermalElement*>* Domain::ThermalElements() const {
  return &thermal_elements_;
}
StatusCode Domain::RegisterThermalSolidSection(ThermalSolidSection* object,
                                              Diagnostic* d) {
  if (!object || !object->material())
    return model_storage::Invalid("Domain thermal section definition", d);
  const auto material = materials_.find(object->material()->Id());
  if (material == materials_.end() || material->second != object->material())
    return model_storage::Invalid("Domain thermal section material identity", d);
  return Insert(object, &thermal_solid_sections_, d);
}
StatusCode Domain::FindThermalSolidSection(int id, const ThermalSolidSection** out,
                                          Diagnostic* d) const {
  return Find(id, &thermal_solid_sections_, out, d);
}
StatusCode Domain::RegisterThermalShellSection(ThermalShellSection* object,
                                               Diagnostic* d) {
  if (!object || !object->material())
    return model_storage::Invalid("Domain thermal shell section definition", d);
  const auto material = materials_.find(object->material()->Id());
  if (material == materials_.end() || material->second != object->material())
    return model_storage::Invalid("Domain thermal shell material identity", d);
  return Insert(object, &thermal_shell_sections_, d);
}
StatusCode Domain::FindThermalShellSection(int id,
                                           const ThermalShellSection** out,
                                           Diagnostic* d) const {
  return Find(id, &thermal_shell_sections_, out, d);
}
StatusCode Domain::FindNode(int id, const Node** out,
                                    Diagnostic* d) const {
  return Find(id, &nodes_, out, d);
}
StatusCode Domain::RegisterElement(Element* object, Diagnostic* d) {
  if (object && thermal_elements_.find(object->Id()) != thermal_elements_.end())
    return model_storage::Invalid("Domain duplicate mechanical/thermal ID", d);
  return Insert(object, &elements_, d);
}
StatusCode Domain::FindElement(int id, const Element** out,
                                    Diagnostic* d) const {
  return Find(id, &elements_, out, d);
}
StatusCode Domain::RegisterFilmPropertyTable(FilmPropertyTable* object,
                                             Diagnostic* d) {
  return InsertNamed(object, 'f', "", &film_properties_, names_, d);
}
StatusCode Domain::FindFilmPropertyTable(int id, const FilmPropertyTable** out,
                                         Diagnostic* d) const {
  return Find(id, &film_properties_, out, d);
}
StatusCode Domain::FindFilmPropertyTableByName(const char* name,
                                               const FilmPropertyTable** out,
                                               Diagnostic* d) const {
  if (!out) return model_storage::Invalid("Domain::FindFilmPropertyTable", d);
  *out = nullptr;
  int id = 0;
  const auto status = FindName('f', "", name, names_, &id, d);
  return status == StatusCode::kOk ? FindFilmPropertyTable(id, out, d) : status;
}
StatusCode Domain::RegisterAmplitude(Amplitude* object, Diagnostic* d) {
  return InsertNamed(object, 'a', "", &amplitudes_, names_, d);
}
StatusCode Domain::FindAmplitude(int id, const Amplitude** out,
                                  Diagnostic* d) const {
  return Find(id, &amplitudes_, out, d);
}
StatusCode Domain::FindAmplitudeByName(const char* name, const Amplitude** out,
                                        Diagnostic* d) const {
  if (!out) return model_storage::Invalid("Domain::FindAmplitude", d);
  *out = nullptr;
  int id = 0;
  const auto status = FindName('a', "", name, names_, &id, d);
  return status == StatusCode::kOk ? FindAmplitude(id, out, d) : status;
}
StatusCode Domain::RegisterMaterial(Material* object, Diagnostic* d) {
  return InsertNamed(object, 'm', "", &materials_, names_, d);
}
StatusCode Domain::FindMaterial(int id, const Material** out,
                                    Diagnostic* d) const {
  return Find(id, &materials_, out, d);
}
StatusCode Domain::RegisterBeamSection(BeamSection* object, Diagnostic* d) {
  return Insert(object, &beam_sections_, d);
}
StatusCode Domain::FindBeamSection(int id, const BeamSection** out,
                                    Diagnostic* d) const {
  return Find(id, &beam_sections_, out, d);
}
StatusCode Domain::RegisterShellSection(ShellSection* object, Diagnostic* d) {
  if (!object || !object->material())
    return model_storage::Invalid("Domain::RegisterShellSection definition", d);
  const auto material = materials_.find(object->material()->Id());
  if (material == materials_.end() || material->second != object->material())
    return model_storage::Invalid("Domain::RegisterShellSection material", d);
  return Insert(object, &shell_sections_, d);
}
StatusCode Domain::FindShellSection(int id, const ShellSection** out,
                                    Diagnostic* d) const {
  return Find(id, &shell_sections_, out, d);
}
StatusCode Domain::RegisterSolidSection(SolidSection* object, Diagnostic* d) {
  if (!object || !object->material())
    return model_storage::Invalid("Domain::RegisterSolidSection definition", d);
  const auto material = materials_.find(object->material()->Id());
  if (material == materials_.end() || material->second != object->material())
    return model_storage::Invalid("Domain::RegisterSolidSection material", d);
  return Insert(object, &solid_sections_, d);
}
StatusCode Domain::FindSolidSection(int id, const SolidSection** out,
                                    Diagnostic* d) const {
  return Find(id, &solid_sections_, out, d);
}
StatusCode Domain::RegisterIndexSet(IndexSet* object, Diagnostic* d) {
  if (object == nullptr) {
    return model_storage::Invalid("Domain::RegisterIndexSet", d);
  }
  return InsertNamed(object, object->IsElementSet() ? 'e' : 'n',
                      object->Scope(), &index_sets_, names_, d);
}
StatusCode Domain::FindIndexSet(int id, const IndexSet** out,
                                    Diagnostic* d) const {
  return Find(id, &index_sets_, out, d);
}
StatusCode Domain::RegisterLoadDefinition(LoadDefinition* object,
    Diagnostic* d) {
  return Insert(object, &loads_, d);
}
StatusCode Domain::FindLoadDefinition(int id, const LoadDefinition** out,
                                    Diagnostic* d) const {
  return Find(id, &loads_, out, d);
}
StatusCode Domain::RegisterBoundaryDefinition(BoundaryDefinition* object,
    Diagnostic* d) {
  return Insert(object, &boundaries_, d);
}
StatusCode Domain::FindBoundaryDefinition(int id, const BoundaryDefinition** out,
                                    Diagnostic* d) const {
  return Find(id, &boundaries_, out, d);
}
StatusCode Domain::RegisterStepDefinition(StepDefinition* object,
    Diagnostic* d) {
  if (object != nullptr && !object->Name()[0])
    return Insert(object, &steps_, d);
  return InsertNamed(object, 's', "", &steps_, names_, d);
}
StatusCode Domain::FindStepDefinition(int id, const StepDefinition** out,
                                    Diagnostic* d) const {
  return Find(id, &steps_, out, d);
}
const std::map<int, Node*>* Domain::Nodes() const { return &nodes_; }
const std::map<int, Element*>* Domain::Elements() const { return &elements_; }
StatusCode Domain::FindMaterialByName(const char* name, const Material** out,
                                      Diagnostic* d) const {
  if (out == nullptr) return model_storage::Invalid("Domain::FindMaterial", d);
  *out = nullptr;
  int id = 0;
  const StatusCode status = FindName('m', "", name, names_, &id, d);
  return status == StatusCode::kOk ? FindMaterial(id, out, d) : status;
}
StatusCode Domain::FindStepByName(const char* name, const StepDefinition** out,
                                  Diagnostic* d) const {
  if (out == nullptr) return model_storage::Invalid("Domain::FindStep", d);
  *out = nullptr;
  int id = 0;
  const StatusCode status = FindName('s', "", name, names_, &id, d);
  return status == StatusCode::kOk ? FindStepDefinition(id, out, d) : status;
}
StatusCode Domain::FindIndexSetByName(const char* scope, const char* name,
    bool element_set, const IndexSet** out, Diagnostic* d) const {
  if (out == nullptr) return model_storage::Invalid("Domain::FindIndexSet", d);
  *out = nullptr;
  int id = 0;
  const StatusCode status = FindName(element_set ? 'e' : 'n', scope, name,
                                     names_, &id, d);
  return status == StatusCode::kOk ? FindIndexSet(id, out, d) : status;
}

StatusCode Domain::SetInputIdentity(const char* source, const char* part,
                                    const char* assembly, const char* instance,
                                    Diagnostic* d) {
  if (!source || !part || !assembly || !instance)
    return model_storage::Invalid("Domain::SetInputIdentity", d);
  std::unique_ptr<char[]> a(model_storage::CopyText(source));
  std::unique_ptr<char[]> b(model_storage::CopyText(part));
  std::unique_ptr<char[]> c(model_storage::CopyText(assembly));
  std::unique_ptr<char[]> e(model_storage::CopyText(instance));
  if (!a || !b || !c || !e)
    return model_storage::AllocationFailure("Domain::SetInputIdentity", d);
  delete[] source_file_;
  delete[] part_name_;
  delete[] assembly_name_;
  delete[] instance_name_;
  source_file_ = a.release();
  part_name_ = b.release();
  assembly_name_ = c.release();
  instance_name_ = e.release();
  return StatusCode::kOk;
}
const char* Domain::SourceFile() const {
  return source_file_ ? source_file_ : "";
}
const char* Domain::PartName() const { return part_name_ ? part_name_ : ""; }
const char* Domain::AssemblyName() const {
  return assembly_name_ ? assembly_name_ : "";
}
const char* Domain::InstanceName() const {
  return instance_name_ ? instance_name_ : "";
}
const std::map<int, BoundaryDefinition*>* Domain::Boundaries() const {
  return &boundaries_;
}
const std::map<int, StepDefinition*>* Domain::Steps() const { return &steps_; }

}  // namespace fesa
