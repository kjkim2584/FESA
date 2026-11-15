#include "fesa/analysis/thermal_dof_map.h"

#include <algorithm>
#include <limits>
#include <memory>
#include <new>
#include <vector>

namespace fesa {
namespace {
StatusCode InvalidMap(const char* message, Diagnostic* d,
                      StatusCode code = StatusCode::kInputError) {
  return ReportFailure(code, DiagnosticStage::kValidation, message, d);
}
}  // namespace
ThermalDofMap::~ThermalDofMap() {
  delete[] nodes_;
  delete[] offsets_;
  delete[] equations_;
  delete[] globals_;
  delete[] participates_;
}
StatusCode ThermalDofMap::Create(const Domain* domain,
                                 const ThermalConditions* conditions,
                                 ThermalDofMap** out, Diagnostic* d) {
  if (!out)
    return InvalidMap("Null thermal map output", d,
                      StatusCode::kInvalidArgument);
  *out = nullptr;
  if (!domain || !conditions || conditions->SourceDomain() != domain)
    return InvalidMap("Thermal map source mismatch", d,
                      StatusCode::kInvalidArgument);
  if (!domain->Elements()->empty() || domain->ThermalElements()->empty() ||
      domain->Nodes()->empty())
    return InvalidMap("Thermal map requires thermal-only element connectivity",
                      d);
  const auto maximum = std::numeric_limits<int>::max();
  if (domain->Nodes()->size() >= static_cast<std::size_t>(maximum))
    return InvalidMap("Thermal node count exceeds index range", d);
  try {
    std::unique_ptr<ThermalDofMap> map(new ThermalDofMap);
    map->domain_ = domain;
    map->step_ = conditions->StepNumber();
    map->node_count_ = static_cast<int>(domain->Nodes()->size());
    map->nodes_ = new int[map->node_count_];
    map->offsets_ = new int[map->node_count_ + 1];
    map->participates_ = new bool[map->node_count_];
    std::vector<int> active(0);
    auto status = CollectThermalNodeIds(domain, &active, d);
    if (status != StatusCode::kOk) return status;
    bool solid = false, shell = false;
    for (const auto entry : *domain->ThermalElements()) {
      const auto* element = entry.second;
      shell = shell || element->type() == ThermalElementType::kDs4;
      solid = solid || element->type() != ThermalElementType::kDs4;
      if (shell && solid)
        return InvalidMap("Mixed solid/shell thermal layout", d);
      for (int i = 0; i < element->node_count(); ++i)
        if (element->GetNode(i)->temperature_dof_count() !=
            element->dofs_per_node())
          return InvalidMap("Shared thermal node section-point mismatch", d);
    }
    int row = 0;
    for (const auto entry : *domain->Nodes()) {
      const int count = entry.second->temperature_dof_count();
      if (count <= 0 || count > maximum - map->global_count_)
        return InvalidMap("Invalid thermal node layout/index range", d);
      map->nodes_[row] = entry.first;
      map->offsets_[row] = map->global_count_;
      map->global_count_ += count;
      map->participates_[row] =
          std::binary_search(active.begin(), active.end(), entry.first);
      ++row;
    }
    map->offsets_[row] = map->global_count_;
    map->equations_ = new int[map->global_count_];
    for (int n = 0; n < map->node_count_; ++n)
      for (int g = map->offsets_[n]; g < map->offsets_[n + 1]; ++g)
        map->equations_[g] = map->participates_[n] ? 0 : -2;
    for (const auto* temperature : *conditions->Temperatures()) {
      const int g = map->GlobalDof(temperature->node_id, temperature->label);
      if (g < 0 || !map->Participates(temperature->node_id))
        return InvalidMap("Invalid prescribed thermal component", d);
      map->equations_[g] = -1;
    }
    for (const auto* flux : *conditions->Fluxes())
      if (map->GlobalDof(flux->node_id, flux->label) < 0 ||
          !map->Participates(flux->node_id))
        return InvalidMap("Invalid concentrated thermal component", d);
    for (int g = 0; g < map->global_count_; ++g)
      if (map->equations_[g] == 0) map->equations_[g] = map->free_count_++;
    if (map->free_count_ > 0) {
      map->globals_ = new int[map->free_count_];
      for (int g = 0; g < map->global_count_; ++g)
        if (map->equations_[g] >= 0) map->globals_[map->equations_[g]] = g;
    }
    *out = map.release();
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return InvalidMap("Thermal map allocation", d,
                      StatusCode::kAllocationFailure);
  }
}
int ThermalDofMap::NodeRow(int node) const {
  if (node_count_ == 0) return -1;
  const int* found = std::lower_bound(nodes_, nodes_ + node_count_, node);
  return found == nodes_ + node_count_ || *found != node
             ? -1
             : static_cast<int>(found - nodes_);
}
int ThermalDofMap::NodeCount() const { return node_count_; }
int ThermalDofMap::NodeId(int row) const {
  return row >= 0 && row < node_count_ ? nodes_[row] : -1;
}
int ThermalDofMap::NodeDofCount(int node) const {
  const int row = NodeRow(node);
  return row < 0 ? 0 : offsets_[row + 1] - offsets_[row];
}
int ThermalDofMap::GlobalCount() const { return global_count_; }
int ThermalDofMap::FreeCount() const { return free_count_; }
int ThermalDofMap::GlobalDof(int node, int label) const {
  const int row = NodeRow(node);
  return row < 0 || label < 11 ||
                 label - 11 >= offsets_[row + 1] - offsets_[row]
             ? -1
             : offsets_[row] + label - 11;
}
int ThermalDofMap::Equation(int g) const {
  return g >= 0 && g < global_count_ ? equations_[g] : -2;
}
int ThermalDofMap::GlobalDofForEquation(int e) const {
  return e >= 0 && e < free_count_ ? globals_[e] : -1;
}
bool ThermalDofMap::Participates(int node) const {
  const int row = NodeRow(node);
  return row >= 0 && participates_[row];
}
const Domain* ThermalDofMap::SourceDomain() const { return domain_; }
int ThermalDofMap::StepNumber() const { return step_; }
}  // namespace fesa
