#include "fesa/elements/incremental_tetra10.h"

#include <memory>
#include <new>

#include "model/model_storage.h"
namespace fesa {
IncrementalTetra10::IncrementalTetra10(int id, const SolidSection* section)
    : IncrementalSolidElement(id, section, SolidFormulation::kTetra10Full) {}
StatusCode IncrementalTetra10::Create(int id, const Node* const* nodes,
    int node_count, const SolidSection* section, SolidFormulation formulation,
    IncrementalTetra10** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("Incremental tetra10 output", d);
  *out = nullptr;
  if (!model_storage::ValidId(id) || node_count != 10 || !section ||
      !section->material() || formulation != SolidFormulation::kTetra10Full)
    return model_storage::Invalid("Incremental tetra10 definition", d);
  std::unique_ptr<IncrementalTetra10> element(
      new (std::nothrow) IncrementalTetra10(id, section));
  if (!element)
    return model_storage::AllocationFailure("Incremental tetra10 object", d);
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
