#include "fesa/elements/incremental_wedge6.h"

#include <memory>
#include <new>

#include "model/model_storage.h"
namespace fesa {
IncrementalWedge6::IncrementalWedge6(int id, const SolidSection* section)
    : IncrementalSolidElement(id, section, SolidFormulation::kWedge6) {}
StatusCode IncrementalWedge6::Create(int id, const Node* const* nodes,
    int node_count, const SolidSection* section, SolidFormulation formulation,
    IncrementalWedge6** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("Incremental wedge6 output", d);
  *out = nullptr;
  if (!model_storage::ValidId(id) || node_count != 6 || !section ||
      !section->material() || formulation != SolidFormulation::kWedge6)
    return model_storage::Invalid("Incremental wedge6 definition", d);
  std::unique_ptr<IncrementalWedge6> element(
      new (std::nothrow) IncrementalWedge6(id, section));
  if (!element)
    return model_storage::AllocationFailure("Incremental wedge6 object", d);
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
