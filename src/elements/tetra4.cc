#include "fesa/elements/tetra4.h"

#include <memory>
#include <new>

#include "model/model_storage.h"
namespace fesa {
Tetra4::Tetra4(int id, const SolidSection* section)
    : SolidElement(id, 4, 1, section, SolidFormulation::kTetra4Full) {}
StatusCode Tetra4::Create(int id, const Node* const* nodes, int node_count,
                          const SolidSection* section,
                          SolidFormulation formulation, Tetra4** out,
                          Diagnostic* d) {
  if (!out) return model_storage::Invalid("Tetra4 output", d);
  *out = nullptr;
  if (!model_storage::ValidId(id) || node_count != 4 || !section ||
      !section->material() || formulation != SolidFormulation::kTetra4Full)
    return model_storage::Invalid("Tetra4 definition", d);
  std::unique_ptr<Tetra4> element(new (std::nothrow) Tetra4(id, section));
  if (!element) return model_storage::AllocationFailure("Tetra4 object", d);
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
