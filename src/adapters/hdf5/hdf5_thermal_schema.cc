#include "adapters/hdf5/hdf5_thermal_schema.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <memory>
#include <new>

#include "adapters/hdf5/hdf5_write_support.h"
#include "fesa/analysis/thermal_state.h"
#include "fesa/math/vector.h"
#include "fesa/model/domain.h"
#include "fesa/model/thermal_element.h"

namespace fesa {
namespace hdf5 {
namespace {
StatusCode Invalid(Diagnostic* d) {
  return internal::Failure("Thermal result model identity/layout", d,
                           StatusCode::kInvalidArgument);
}
bool Field(hid_t group, const char* name, hid_t file_type, hid_t memory_type,
           int count, const void* data, const char* unit,
           const char* evaluation, const char* coordinates,
           const char* component, Diagnostic* d) {
  const hsize_t size[] = {static_cast<hsize_t>(count)};
  const internal::DatasetMetadata metadata = {unit, evaluation, coordinates,
                                              &component, 1};
  return internal::Dataset(group, name, file_type, memory_type, 1, size, data,
                           d, &metadata);
}
}  // namespace
StatusCode ValidateThermalResultModel(const Domain* domain, Diagnostic* d) {
  if (!domain || domain->Nodes()->empty() ||
      domain->ThermalElements()->empty() || !domain->Elements()->empty() ||
      domain->Nodes()->size() > INT_MAX / 20 ||
      domain->ThermalElements()->size() > INT_MAX / 8)
    return Invalid(d);
  for (const auto& item : *domain->Nodes()) {
    const Node* node = item.second;
    if (!node || item.first < 1 || item.first > 999999999 ||
        item.first != node->Id() || node->dof_count() != 0 ||
        node->temperature_dof_count() < 1 || node->temperature_dof_count() > 20)
      return Invalid(d);
  }
  const bool shell = domain->ThermalElements()->begin()->second->type() ==
                     ThermalElementType::kDs4;
  for (const auto& item : *domain->ThermalElements()) {
    const ThermalElement* element = item.second;
    if (!element || item.first < 1 || item.first > 999999999 ||
        item.first != element->Id() || !element->material())
      return Invalid(d);
    const Material* material = nullptr;
    if (domain->FindMaterial(element->material()->Id(), &material, d) !=
            StatusCode::kOk ||
        material != element->material())
      return Invalid(d);
    const int count = element->dofs_per_node();
    const bool ds4 = element->type() == ThermalElementType::kDs4 &&
                     element->node_count() == 4 && count >= 1 && count <= 19 &&
                     count % 2 == 1 &&
                     std::strcmp(element->type_name(), "DS4") == 0;
    const bool tetra = element->type() == ThermalElementType::kDc3d4 &&
                       element->node_count() == 4 && count == 1 &&
                       std::strcmp(element->type_name(), "DC3D4") == 0;
    const bool hexa = element->type() == ThermalElementType::kDc3d8 &&
                      element->node_count() == 8 && count == 1 &&
                      std::strcmp(element->type_name(), "DC3D8") == 0;
    if (shell ? !ds4 : !(tetra || hexa)) return Invalid(d);
    if (shell) {
      double thickness = 0;
      const auto status = element->GetShellThickness(&thickness, d);
      if (status != StatusCode::kOk) return status;
      if (!std::isfinite(thickness) || thickness <= 0) return Invalid(d);
    }
    for (int local = 0; local < element->node_count(); ++local) {
      const Node* node = element->GetNode(local);
      const Node* registered = nullptr;
      if (!node || node->temperature_dof_count() != count ||
          domain->FindNode(node->Id(), &registered, d) != StatusCode::kOk ||
          registered != node)
        return Invalid(d);
    }
  }
  return StatusCode::kOk;
}
ThermalSchema::ThermalSchema(const Domain* domain) : domain_(domain) {}
ThermalSchema::~ThermalSchema() {
  delete packed_;
  delete gather_;
  delete nodes_;
}
StatusCode ThermalSchema::Create(const Domain* domain, ThermalSchema** out,
                                 Diagnostic* d) {
  if (!out) return Invalid(d);
  *out = nullptr;
  auto status = ValidateThermalResultModel(domain, d);
  if (status != StatusCode::kOk) return status;
  try {
    std::unique_ptr<ThermalSchema> schema(new ThermalSchema(domain));
    schema->nodes_ = new std::vector<int>(0);
    schema->gather_ = new std::vector<int>(0);
    status = CollectThermalNodeIds(domain, schema->nodes_, d);
    if (status != StatusCode::kOk) return status;
    int nominal = 0;
    for (const auto& item : *domain->Nodes()) {
      const int count = item.second->temperature_dof_count();
      if (std::binary_search(schema->nodes_->begin(), schema->nodes_->end(),
                             item.first))
        for (int i = 0; i < count; ++i) schema->gather_->push_back(nominal + i);
      nominal += count;
    }
    status = Vector::Create(static_cast<int>(schema->gather_->size()),
                            &schema->packed_, d);
    if (status != StatusCode::kOk) return status;
    schema->shell_ = domain->ThermalElements()->begin()->second->type() ==
                     ThermalElementType::kDs4;
    *out = schema.release();
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return internal::Failure("Thermal schema allocation", d,
                             StatusCode::kAllocationFailure);
  }
}
const char* ThermalSchema::Family() const {
  return shell_ ? "heat_transfer_shell" : "heat_transfer_solid";
}
StatusCode ThermalSchema::WriteModel(hid_t model, Diagnostic* d) const {
  using internal::Check;
  using internal::Dataset;
  using internal::Handle;
  const int n = static_cast<int>(nodes_->size());
  const int e = static_cast<int>(domain_->ThermalElements()->size());
  const int count = packed_->Size();
  try {
    std::vector<std::int64_t> offsets(1, 0), connectivity_offsets(1, 0);
    std::vector<int> labels(0), ids(0), connectivity(0), sections(0);
    std::vector<const char*> types(0);
    Vector* raw = nullptr;
    auto status = Vector::Create(3 * n, &raw, d);
    if (status != StatusCode::kOk) return status;
    std::unique_ptr<Vector> coordinates(raw);
    std::unique_ptr<Vector> thickness(nullptr), normalized(nullptr);
    if (shell_) {
      status = Vector::Create(e, &raw, d);
      if (status != StatusCode::kOk) return status;
      thickness.reset(raw);
      status = Vector::Create(count, &raw, d);
      if (status != StatusCode::kOk) return status;
      normalized.reset(raw);
    }
    for (int row = 0; row < n; ++row) {
      const Node* node = nullptr;
      status = domain_->FindNode((*nodes_)[row], &node, d);
      if (status != StatusCode::kOk) return status;
      for (int axis = 0; axis < 3; ++axis) {
        double value = 0;
        status = node->Coordinate(axis, &value, d);
        if (status == StatusCode::kOk)
          status = coordinates->Set(3 * row + axis, value, d);
        if (status != StatusCode::kOk) return status;
      }
      const int ns = node->temperature_dof_count();
      for (int j = 0; j < ns; ++j) {
        if (shell_) {
          status = normalized->Set(static_cast<int>(labels.size()),
                                   ns == 1 ? 0 : -1.0 + 2.0 * j / (ns - 1), d);
          if (status != StatusCode::kOk) return status;
        }
        labels.push_back(11 + j);
      }
      offsets.push_back(static_cast<std::int64_t>(labels.size()));
    }
    int row = 0;
    for (const auto& item : *domain_->ThermalElements()) {
      const ThermalElement* element = item.second;
      ids.push_back(item.first);
      types.push_back(element->type_name());
      for (int local = 0; local < element->node_count(); ++local)
        connectivity.push_back(element->GetNode(local)->Id());
      connectivity_offsets.push_back(
          static_cast<std::int64_t>(connectivity.size()));
      if (shell_) {
        double value = 0;
        status = element->GetShellThickness(&value, d);
        if (status == StatusCode::kOk) status = thickness->Set(row, value, d);
        if (status != StatusCode::kOk) return status;
        sections.push_back(element->dofs_per_node());
      }
      ++row;
    }
    const hsize_t nn[] = {static_cast<hsize_t>(n)};
    const hsize_t xyz[] = {static_cast<hsize_t>(n), 3};
    const hsize_t ee[] = {static_cast<hsize_t>(e)};
    const hsize_t eo[] = {static_cast<hsize_t>(e + 1)};
    const hsize_t cc[] = {static_cast<hsize_t>(connectivity.size())};
    if (!Dataset(model, "node_ids", H5T_STD_I32LE, H5T_NATIVE_INT, 1, nn,
                 nodes_->data(), d) ||
        !Dataset(model, "node_coordinates", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE,
                 2, xyz, coordinates->Data(), d) ||
        !Dataset(model, "element_ids", H5T_STD_I32LE, H5T_NATIVE_INT, 1, ee,
                 ids.data(), d) ||
        !internal::StringDataset(model, "element_types", e, types.data(), d) ||
        !Dataset(model, "connectivity_offsets", H5T_STD_I64LE, H5T_NATIVE_INT64,
                 1, eo, connectivity_offsets.data(), d) ||
        !Dataset(model, "connectivity_node_ids", H5T_STD_I32LE, H5T_NATIVE_INT,
                 1, cc, connectivity.data(), d))
      return StatusCode::kIoFailure;
    Handle thermal(
        H5Gcreate2(model, "thermal", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
        H5Gclose);
    if (!Check(thermal.get(), "Create thermal model group", d) ||
        !internal::ScalarText(thermal.get(), "temperature_unit", "model", d) ||
        !internal::ScalarText(thermal.get(), "reaction_flux_unit",
                              "energy/time (model units)", d) ||
        !Field(thermal.get(), "node_dof_offsets", H5T_STD_I64LE,
               H5T_NATIVE_INT64, n + 1, offsets.data(), "1",
               "initial_node_offset", "index", "node_dof_offset", d) ||
        !Field(thermal.get(), "node_dof_labels", H5T_STD_I32LE, H5T_NATIVE_INT,
               count, labels.data(), "1", "initial_node_dof", "index",
               "temperature_dof_label", d))
      return StatusCode::kIoFailure;
    if (shell_) {
      Handle shell(H5Gcreate2(thermal.get(), "shell", H5P_DEFAULT, H5P_DEFAULT,
                              H5P_DEFAULT),
                   H5Gclose);
      if (!Check(shell.get(), "Create thermal shell group", d) ||
          !Field(shell.get(), "thickness", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, e,
                 thickness->Data(), "L", "initial_element", "scalar",
                 "thickness", d) ||
          !Field(shell.get(), "section_points", H5T_STD_I32LE, H5T_NATIVE_INT,
                 e, sections.data(), "1", "initial_element", "index",
                 "section_points", d) ||
          !Field(shell.get(), "dof_normalized_coordinates", H5T_IEEE_F64LE,
                 H5T_NATIVE_DOUBLE, count, normalized->Data(), "1",
                 "initial_node_dof", "normalized_thickness", "zeta", d) ||
          !shell.Close(d))
        return StatusCode::kIoFailure;
    }
    return thermal.Close(d) ? StatusCode::kOk : StatusCode::kIoFailure;
  } catch (const std::bad_alloc&) {
    return internal::Failure("Thermal model payload allocation", d,
                             StatusCode::kAllocationFailure);
  }
}
StatusCode ThermalSchema::WriteFramePayload(hid_t nodal,
                                            const ThermalState* state,
                                            Diagnostic* d) {
  const Vector* fields[] = {state->CandidateTemperature(),
                            state->CandidateReactionFlux()};
  const char* names[] = {"NT", "RFL"};
  const char* units[] = {"theta", "E/T"};
  for (int i = 0; i < 2; ++i) {
    const auto status =
        packed_->GatherFrom(fields[i], gather_->data(), packed_->Size(), 0, d);
    if (status != StatusCode::kOk) return status;
    if (!Field(nodal, names[i], H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE,
               packed_->Size(), packed_->Data(), units[i],
               "increment_end_node_dof", "scalar", names[i], d))
      return StatusCode::kIoFailure;
  }
  return StatusCode::kOk;
}
}  // namespace hdf5
}  // namespace fesa
