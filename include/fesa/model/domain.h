/// @file
/// @brief Linear-static model ownership and immutable input definitions.
#ifndef FESA_MODEL_DOMAIN_H_
#define FESA_MODEL_DOMAIN_H_

#include <map>
#include <string>

#include "fesa/model/amplitude.h"
#include "fesa/model/beam_section.h"
#include "fesa/model/boundary_definition.h"
#include "fesa/model/element.h"
#include "fesa/model/index_set.h"
#include "fesa/model/load_definition.h"
#include "fesa/model/material.h"
#include "fesa/model/node.h"
#include "fesa/model/shell_section.h"
#include "fesa/model/solid_section.h"
#include "fesa/model/step_definition.h"
#include "fesa/model/thermal_condition.h"
#include "fesa/model/thermal_element.h"
#include "fesa/model/thermal_initial_temperature.h"
#include "fesa/model/thermal_section.h"
#include "fesa/model/thermal_physical_constants.h"

namespace fesa {

class BeamLineLoadProgram;
class FilmPropertyTable;

/// Owns objects only after successful registration. Destroy views first.
/// Definitions are immutable during analysis. Before replacing/deleting a
/// Domain,
/// callers must Clear every AnalysisModel borrowing its objects.
class Domain {
 public:
  /// Caller owns *out; failure nulls it, d is optional and borrowed.
  static StatusCode Create(Domain** out, Diagnostic* d);
  ~Domain();
  Domain(const Domain&) = delete;
  Domain& operator=(const Domain&) = delete;
  Domain(Domain&&) = delete;
  Domain& operator=(Domain&&) = delete;
  /// Validates and owns an explicit clone; rejects a second set. Failure
  /// preserves the prior definition and caller ownership. Input borrowed.
  StatusCode SetThermalPhysicalConstants(
      const ThermalPhysicalConstants* values, Diagnostic* d);
  /// Borrowed until destruction, null when no block was supplied.
  const ThermalPhysicalConstants* ThermalConstants() const;
  /// Adopts only on success; node must exist, be thermal and match count.
  /// Duplicate node records (even equal) fail; registry/caller stay unchanged.
  StatusCode RegisterThermalInitialTemperature(
      ThermalInitialTemperature* object, Diagnostic* d);
  /// Borrowed until Domain destruction; missing ID nulls *out, no insertion.
  StatusCode FindThermalInitialTemperature(
      int node_id, const ThermalInitialTemperature** out, Diagnostic* d) const;
  /// Borrowed immutable IC registry; absent nodes have default initial zero.
  const std::map<int, ThermalInitialTemperature*>*
  ThermalInitialTemperatures() const;
  /// Transfers ownership only on success; duplicate/missing/null errors keep
  /// both Domain and caller ownership unchanged. Object ID must be positive.
  StatusCode RegisterNode(Node* object, Diagnostic* d);
  /// Borrows result until Domain destruction. Missing ID nulls *out; no insert.
  StatusCode FindNode(int id, const Node** out, Diagnostic* d) const;
  /// Transfers ownership only on success; duplicate/missing/null errors keep
  /// both Domain and caller ownership unchanged. Object ID must be positive.
  StatusCode RegisterElement(Element* object, Diagnostic* d);
  /// Borrows result until Domain destruction. Missing ID nulls *out; no insert.
  StatusCode FindElement(int id, const Element** out, Diagnostic* d) const;
  /// Adopts on success only. IDs are unique across both element registries.
  /// Nodes/material must be identical registered objects; failure preserves
  /// registries and caller ownership. Register before creating active views.
  StatusCode RegisterThermalElement(ThermalElement* object, Diagnostic* d);
  /// Borrowed until destruction; missing ID nulls *out without insertion.
  StatusCode FindThermalElement(int id, const ThermalElement** out,
                                Diagnostic* d) const;
  /// Borrowed immutable registry; Domain owns all values after registration.
  const std::map<int, ThermalElement*>* ThermalElements() const;
  /// Adopts on success; material must be the identical registered object.
  /// Failure preserves registry/caller ownership. Material outlives section.
  StatusCode RegisterThermalSolidSection(ThermalSolidSection* object,
                                         Diagnostic* d);
  /// Borrowed until destruction; missing ID nulls *out without insertion.
  StatusCode FindThermalSolidSection(int id, const ThermalSolidSection** out,
                                     Diagnostic* d) const;
  /// Adopts on success only; material must be identical to the registered
  /// Domain object. Failure preserves registry/caller ownership. Shell
  /// elements borrow this section and are destroyed before sections.
  StatusCode RegisterThermalShellSection(ThermalShellSection* object,
                                         Diagnostic* d);
  /// Borrows the registered section until Domain destruction; missing ID
  /// nulls *out, without insertion or mutation of the registry.
  StatusCode FindThermalShellSection(int id, const ThermalShellSection** out,
                                     Diagnostic* d) const;
  /// Adopts only on successful registration; optional amplitude must be the
  /// identical registered definition. Failure preserves both owners/registry.
  StatusCode RegisterThermalCondition(ThermalConditionDefinition* object,
                                       Diagnostic* d);
  /// Borrowed until destruction; missing ID nulls *out without insertion.
  StatusCode FindThermalCondition(
      int id, const ThermalConditionDefinition** out, Diagnostic* d) const;
  /// Borrowed registry; model-level and step conditions retain original IDs.
  const std::map<int, ThermalConditionDefinition*>*
  ThermalConditions() const;
  /// Transfers ownership only on success; duplicate/missing/null errors keep
  /// both Domain and caller ownership unchanged. Object ID must be positive.
  StatusCode RegisterMaterial(Material* object, Diagnostic* d);
  /// Adopts only on success; duplicate positive ID/name preserves ownership.
  /// Register before binding active views. Definitions are Domain-owned.
  StatusCode RegisterFilmPropertyTable(FilmPropertyTable* object,
                                       Diagnostic* d);
  /// Borrowed immutable property until Domain destruction; failure nulls out.
  StatusCode FindFilmPropertyTable(int id, const FilmPropertyTable** out,
                                   Diagnostic* d) const;
  /// Case-insensitive name lookup, no insertion; same borrowed-result contract.
  StatusCode FindFilmPropertyTableByName(const char* name,
                                         const FilmPropertyTable** out,
                                         Diagnostic* d) const;
  /// Adopts only on success. Duplicate ID/name leaves caller ownership intact.
  /// Definitions must be registered before any borrowed views are published.
  StatusCode RegisterAmplitude(Amplitude* object, Diagnostic* d);
  /// Borrows until Domain destruction; failure nulls *out without insertion.
  StatusCode FindAmplitude(int id, const Amplitude** out, Diagnostic* d) const;
  /// Case-insensitive name lookup with the same borrowed-result contract.
  StatusCode FindAmplitudeByName(const char* name, const Amplitude** out,
                                 Diagnostic* d) const;
  /// Borrows result until Domain destruction. Missing ID nulls *out; no insert.
  StatusCode FindMaterial(int id, const Material** out, Diagnostic* d) const;
  /// Transfers ownership only on success; duplicate/missing/null errors keep
  /// both Domain and caller ownership unchanged. Object ID must be positive.
  StatusCode RegisterBeamSection(BeamSection* object, Diagnostic* d);
  /// Borrows result until Domain destruction. Missing ID nulls *out; no insert.
  StatusCode FindBeamSection(int id, const BeamSection** out,
      Diagnostic* d) const;
  /// Adopts only a section referencing the identical registered material.
  /// Success transfers ownership; errors preserve Domain and caller ownership.
  /// Register before publishing active views; material outlives section.
  StatusCode RegisterShellSection(ShellSection* object, Diagnostic* d);
  /// Borrowed until Domain destruction; failure nulls *out, d optional.
  StatusCode FindShellSection(int id, const ShellSection** out,
                              Diagnostic* d) const;
  /// Adopts on success only; requires this Domain's identical material pointer.
  /// Errors preserve both ownership and the registry. Material outlives
  /// section.
  StatusCode RegisterSolidSection(SolidSection* object, Diagnostic* d);
  /// Borrows until Domain destruction; missing ID nulls *out without insertion.
  StatusCode FindSolidSection(int id, const SolidSection** out,
                              Diagnostic* d) const;
  /// Transfers ownership only on success; duplicate/missing/null errors keep
  /// both Domain and caller ownership unchanged. Object ID must be positive.
  StatusCode RegisterIndexSet(IndexSet* object, Diagnostic* d);
  /// Borrows result until Domain destruction. Missing ID nulls *out; no insert.
  StatusCode FindIndexSet(int id, const IndexSet** out, Diagnostic* d) const;
  /// Transfers ownership only on success; duplicate/missing/null errors keep
  /// both Domain and caller ownership unchanged. Object ID must be positive.
  StatusCode RegisterLoadDefinition(LoadDefinition* object, Diagnostic* d);
  /// Borrows result until Domain destruction. Missing ID nulls *out; no insert.
  StatusCode FindLoadDefinition(int id, const LoadDefinition** out,
      Diagnostic* d) const;
  /// Transfers ownership only on success; duplicate/missing/null errors keep
  /// both Domain and caller ownership unchanged. Object ID must be positive.
  StatusCode RegisterBoundaryDefinition(BoundaryDefinition* object,
      Diagnostic* d);
  /// Borrows result until Domain destruction. Missing ID nulls *out; no insert.
  StatusCode FindBoundaryDefinition(int id, const BoundaryDefinition** out,
      Diagnostic* d) const;
  /// Transfers ownership only on success; duplicate/missing/null errors keep
  /// both Domain and caller ownership unchanged. Object ID must be positive.
  StatusCode RegisterStepDefinition(StepDefinition* object, Diagnostic* d);
  /// Borrows result until Domain destruction. Missing ID nulls *out; no insert.
  StatusCode FindStepDefinition(int id, const StepDefinition** out,
      Diagnostic* d) const;
  /// Borrowed map for active-view binding; valid while Domain lives.
  const std::map<int, Node*>* Nodes() const;
  /// Borrowed map for active-view binding; valid while Domain lives.
  const std::map<int, Element*>* Elements() const;
  /// Case-insensitive lookup without insertion; borrows returned object.
  StatusCode FindMaterialByName(const char* name, const Material** out,
                                Diagnostic* d) const;
  /// Case-insensitive scope/kind/name lookup; missing result is null.
  StatusCode FindIndexSetByName(const char* scope, const char* name,
                                bool element_set, const IndexSet** out,
                                Diagnostic* d) const;
  /// Case-insensitive lookup, retaining the original stored step name.
  StatusCode FindStepByName(const char* name, const StepDefinition** out,
                            Diagnostic* d) const;

  /// Copies immutable input identity; failure preserves existing strings.
  /// Must be called before publishing the Domain or binding active views.
  StatusCode SetInputIdentity(const char* source, const char* part,
                              const char* assembly, const char* instance,
                              Diagnostic* d);
  /// Borrowed original identity strings, empty before identity is supplied.
  const char* SourceFile() const;
  /// Original Part name.
  const char* PartName() const;
  /// Original Assembly name.
  const char* AssemblyName() const;
  /// Original Instance name.
  const char* InstanceName() const;
  /// Borrowed original definitions, including model boundaries with Step()==0.
  const std::map<int, BoundaryDefinition*>* Boundaries() const;
  /// Borrowed step definitions in input order (positive IDs).
  const std::map<int, StepDefinition*>* Steps() const;
  /// Adopts an immutable line-load program for this Domain only, before active
  /// views are published. Success transfers ownership; null, wrong model or
  /// a second attachment preserves Domain and caller ownership. Views/states
  /// must die first; Domain deletes the program before its borrowed elements.
  StatusCode AttachLineLoadProgram(BeamLineLoadProgram* program, Diagnostic* d);
  /// Borrowed optional program until Domain destruction; null for old models.
  const BeamLineLoadProgram* LineLoadProgram() const;

 private:
  char* source_file_ = nullptr;    ///< Owned UTF-8 source path.
  BeamLineLoadProgram* line_load_program_ = nullptr;  ///< Owned, borrows this.
  char* part_name_ = nullptr;      ///< Owned original identity.
  char* assembly_name_ = nullptr;  ///< Owned original identity.
  char* instance_name_ = nullptr;  ///< Owned original identity.
  explicit Domain(std::map<std::string, int>* names);
  ThermalPhysicalConstants* thermal_constants_ = nullptr;  ///< Owned or null.
  /// Owns immutable initial values, independent of Node result buffers.
  std::map<int, ThermalInitialTemperature*> thermal_initial_temperatures_;
  std::map<int, Node*> nodes_;  ///< Owns registered objects.
  std::map<int, Element*> elements_;  ///< Owns registered objects.
  std::map<int, ThermalElement*> thermal_elements_;  ///< Owns thermal elements.
  /// Owns definitions; borrowed by steps, destroyed before amplitudes.
  std::map<int, ThermalConditionDefinition*> thermal_conditions_;
  std::map<int, ThermalSolidSection*> thermal_solid_sections_;  ///< Owns sections.
  std::map<int, ThermalShellSection*> thermal_shell_sections_;  ///< Owns sections.
  std::map<int, FilmPropertyTable*> film_properties_;  ///< Owns film tables.
  std::map<int, Amplitude*> amplitudes_;  ///< Owns registered definitions.
  std::map<int, Material*> materials_;  ///< Owns registered objects.
  std::map<int, ShellSection*> shell_sections_;  ///< Owns registered sections.
  std::map<int, SolidSection*> solid_sections_;  ///< Owns registered sections.
  std::map<int, BeamSection*> beam_sections_;  ///< Owns registered objects.
  std::map<int, IndexSet*> index_sets_;  ///< Owns registered objects.
  std::map<int, LoadDefinition*> loads_;  ///< Owns registered objects.
  std::map<int, BoundaryDefinition*> boundaries_;  ///< Owns registered objects.
  std::map<int, StepDefinition*> steps_;  ///< Owns registered objects.
  std::map<std::string, int>* names_;  ///< Owned scope/name index.
};

}  // namespace fesa

#endif  // FESA_MODEL_DOMAIN_H_
