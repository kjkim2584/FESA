#include "fesa/elements/wedge6.h"

#include <memory>
#include <new>

#include "model/model_storage.h"
namespace fesa {
Wedge6::Wedge6(int id, const SolidSection* section)
    : SolidElement(id, 6, 2, section, SolidFormulation::kWedge6) {}
StatusCode Wedge6::Create(int id, const Node* const* nodes, int node_count,
                          const SolidSection* section,
                          SolidFormulation formulation, Wedge6** out,
                          Diagnostic* d) {
  if (!out) return model_storage::Invalid("Wedge6 output", d);
  *out = nullptr;
  if (!model_storage::ValidId(id) || node_count != 6 || !section ||
      !section->material() || formulation != SolidFormulation::kWedge6)
    return model_storage::Invalid("Wedge6 definition", d);
  std::unique_ptr<Wedge6> element(new (std::nothrow) Wedge6(id, section));
  if (!element) return model_storage::AllocationFailure("Wedge6 object", d);
  auto status = element->Initialize(nodes, d);
  if (status != StatusCode::kOk) return status;
  status = element->ValidateGeometry(d);
  if (status != StatusCode::kOk) {
    if (d) d->element_id = id;
    return status;
  }
  *out = element.release();
  return StatusCode::kOk;
}
}  // namespace fesa
