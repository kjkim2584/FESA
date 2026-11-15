#include "fesa/elements/incremental_tetra4.h"

#include <memory>
#include <new>

#include "model/model_storage.h"
namespace fesa {
IncrementalTetra4::IncrementalTetra4(int id, const SolidSection* section)
    : IncrementalSolidElement(id, section, SolidFormulation::kTetra4Full) {}
StatusCode IncrementalTetra4::Create(int id, const Node* const* nodes,
    int node_count, const SolidSection* section, SolidFormulation formulation,
    IncrementalTetra4** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("Incremental tetra4 output", d);
  *out = nullptr;
  if (!model_storage::ValidId(id) || node_count != 4 || !section ||
      !section->material() || formulation != SolidFormulation::kTetra4Full)
    return model_storage::Invalid("Incremental tetra4 definition", d);
  std::unique_ptr<IncrementalTetra4> element(
      new (std::nothrow) IncrementalTetra4(id, section));
  if (!element)
    return model_storage::AllocationFailure("Incremental tetra4 object", d);
  auto status = element->Initialize(nodes, d);
  if (status == StatusCode::kOk) status = element->ValidateGeometry(d);
  if (status != StatusCode::kOk) {
    if (d) d->element_id = id;
    return status;
  }
  *out = element.release();
  return StatusCode::kOk;
}
}  // namespace fesa
