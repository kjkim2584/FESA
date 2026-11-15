#include "fesa/analysis/dof_manager.h"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <memory>
#include <new>
#include <utility>

#include "fesa/elements/beam3d.h"
#include "model/model_storage.h"
#include "model/solid_connectivity.h"

namespace fesa {
DofManager::~DofManager() {
  delete[] internal_element_ids_;
  delete[] node_ids_;
  delete[] node_offsets_;
  delete[] equations_;
  delete[] global_dofs_;
  delete[] participates_;
}
StatusCode DofManager::Create(const Domain* domain,
                              const EffectiveConditions* conditions,
                              DofManager** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("DofManager::Create", d);
  *out = nullptr;
  std::unique_ptr<DofManager> result(new (std::nothrow) DofManager);
  if (!result) return model_storage::AllocationFailure("DofManager::Create", d);
  StatusCode status = result->Build(domain, conditions, d);
  if (status != StatusCode::kOk) return status;
  *out = result.release();
  return StatusCode::kOk;
}
StatusCode DofManager::Build(const Domain* domain,
                             const EffectiveConditions* conditions,
                             Diagnostic* d) {
  if (!domain || !conditions || conditions->SourceDomain() != domain ||
      domain->Nodes()->empty() ||
      domain->Nodes()->size() >= std::numeric_limits<int>::max())
    return model_storage::Invalid("DofManager::Build", d);
  try {
    DofManager next;
    next.node_count_ = static_cast<int>(domain->Nodes()->size());
    next.node_ids_ = new int[next.node_count_];
    next.node_offsets_ = new int[next.node_count_ + 1];
    next.participates_ = new bool[next.node_count_]();
    int node_row = 0;
    for (const auto& entry : *domain->Nodes()) {
      const int components = entry.second->dof_count();
      if (components == 0)
        return model_storage::Invalid("Thermal node in mechanical mapping", d);
      if (next.global_count_ > std::numeric_limits<int>::max() - components)
        return model_storage::Invalid("DofManager physical size overflow", d);
      next.node_ids_[node_row] = entry.first;
      next.node_offsets_[node_row++] = next.global_count_;
      next.global_count_ += components;
    }
    next.node_offsets_[node_row] = next.global_count_;
    next.equations_ = new int[next.global_count_];
    std::fill_n(next.equations_, next.global_count_, -2);
    for (const auto& item : *domain->Elements()) {
      const Element* element = item.second;
      const int nodes = element->node_count();
      const int components = element->dofs_per_node();
      if (nodes <= 0 || components <= 0 ||
          nodes > std::numeric_limits<int>::max() / components ||
          element->dof_count() != nodes * components)
        return model_storage::Invalid("DofManager element layout", d);
      for (int local = 0; local < nodes; ++local) {
        const Node* node = element->GetNode(local);
        const Node* registered = nullptr;
        if (!node || next.NodeDofCount(node->Id()) != components ||
            domain->FindNode(node->Id(), &registered, d) != StatusCode::kOk ||
            node != registered)
          return model_storage::Invalid("DofManager element node/layout", d);
        next.participates_[next.NodeRow(node->Id())] = true;
      }
    }
    const auto topology_status = ValidateSolidConnectivity(domain, d);
    if (topology_status != StatusCode::kOk) return topology_status;
    for (const auto* load : *conditions->Loads()) {
      const int row = next.NodeRow(load->node_id);
      if (row < 0 || next.GlobalDof(load->node_id, load->definition->Dof()) < 0)
        return model_storage::Invalid("DofManager load component", d);
      if (load->contributes_flux || load->change != 0)
        next.participates_[row] = true;
    }
    for (const auto* bc : *conditions->Boundaries()) {
      const int global = next.GlobalDof(bc->node_id, bc->definition->Dof());
      if (global < 0)
        return model_storage::Invalid("DofManager boundary component", d);
      next.participates_[next.NodeRow(bc->node_id)] = true;
      next.equations_[global] = -1;
    }
    for (int row = 0; row < next.node_count_; ++row)
      for (int global = next.node_offsets_[row];
           global < next.node_offsets_[row + 1]; ++global)
        if (next.participates_[row] && next.equations_[global] != -1)
          next.equations_[global] = next.free_count_++;
    if (next.free_count_ > 0) {
      next.global_dofs_ = new int[next.free_count_];
      for (int i = 0; i < next.global_count_; ++i)
        if (next.equations_[i] >= 0) next.global_dofs_[next.equations_[i]] = i;
    }
    if (includes_beam_internal_) {
      if (domain->Elements()->empty() ||
          domain->Elements()->size() >
              static_cast<std::size_t>(
                  (std::numeric_limits<int>::max() - next.free_count_) / 2))
        return model_storage::Invalid("Beam internal equation count", d);
      for (const auto& entry : *domain->Elements()) {
        if (!dynamic_cast<const Beam3D*>(entry.second)) {
          if (d) d->element_id = entry.first;
          return ReportFailure(StatusCode::kUnsupported,
                               DiagnosticStage::kValidation,
                               "Internal equations require exclusive Beam3D",
                               d);
        }
      }
      next.internal_element_count_ =
          static_cast<int>(domain->Elements()->size());
      next.internal_element_ids_ = new int[next.internal_element_count_];
      int element_row = 0;
      for (const auto& entry : *domain->Elements())
        next.internal_element_ids_[element_row++] = entry.first;
    }
    std::swap(global_dofs_, next.global_dofs_);
    std::swap(node_count_, next.node_count_);
    std::swap(global_count_, next.global_count_);
    std::swap(node_offsets_, next.node_offsets_);
    std::swap(free_count_, next.free_count_);
    std::swap(node_ids_, next.node_ids_);
    std::swap(equations_, next.equations_);
    std::swap(participates_, next.participates_);
    std::swap(internal_element_count_, next.internal_element_count_);
    std::swap(internal_element_ids_, next.internal_element_ids_);
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return model_storage::AllocationFailure("DofManager::Build", d);
  }
}
StatusCode DofManager::CreateForNonlinearBeam(
    const Domain* domain, const EffectiveConditions* conditions,
    DofManager** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("Beam DofManager output", d);
  *out = nullptr;
  std::unique_ptr<DofManager> result(new (std::nothrow) DofManager);
  if (!result)
    return model_storage::AllocationFailure("Beam DofManager", d);
  result->includes_beam_internal_ = true;
  const auto status = result->Build(domain, conditions, d);
  if (status != StatusCode::kOk) return status;
  *out = result.release();
  return StatusCode::kOk;
}
int DofManager::EquationCount() const {
  return free_count_ + InternalDofCount();
}
int DofManager::InternalDofCount() const {
  return 2 * internal_element_count_;
}
int DofManager::InternalElementId(int row) const {
  return row >= 0 && row < internal_element_count_
             ? internal_element_ids_[row] : -1;
}
int DofManager::InternalDof(int element_id, int end) const {
  if (end < 0 || end > 1 || internal_element_count_ == 0) return -1;
  const int* begin = internal_element_ids_;
  const int* finish = begin + internal_element_count_;
  const int* found = std::lower_bound(begin, finish, element_id);
  return found == finish || *found != element_id
             ? -1 : 2 * static_cast<int>(found - begin) + end;
}
int DofManager::InternalEquation(int element_id, int end) const {
  const int offset = InternalDof(element_id, end);
  return offset < 0 ? -1 : free_count_ + offset;
}
int DofManager::NodeCount() const { return node_count_; }
int DofManager::GlobalCount() const { return global_count_; }
int DofManager::NodeDofCount(int id) const {
  const int row = NodeRow(id);
  return row < 0 ? 0 : node_offsets_[row + 1] - node_offsets_[row];
}
int DofManager::FreeCount() const { return free_count_; }
int DofManager::NodeId(int row) const {
  return row >= 0 && row < node_count_ ? node_ids_[row] : -1;
}
int DofManager::NodeRow(int id) const {
  if (node_count_ == 0) return -1;
  const int* found = std::lower_bound(node_ids_, node_ids_ + node_count_, id);
  return found == node_ids_ + node_count_ || *found != id
             ? -1
             : static_cast<int>(found - node_ids_);
}
int DofManager::GlobalDof(int node_id, int dof) const {
  const int row = NodeRow(node_id);
  if (row < 0 || dof < 1 || dof > node_offsets_[row + 1] - node_offsets_[row])
    return -1;
  return node_offsets_[row] + dof - 1;
}
int DofManager::Equation(int global_dof) const {
  return global_dof >= 0 && global_dof < global_count_ ? equations_[global_dof]
                                                       : -2;
}
int DofManager::GlobalDofForEquation(int equation) const {
  return equation >= 0 && equation < free_count_ ? global_dofs_[equation] : -1;
}
bool DofManager::Participates(int node_id) const {
  const int row = NodeRow(node_id);
  return row >= 0 && participates_[row];
}
}  // namespace fesa
