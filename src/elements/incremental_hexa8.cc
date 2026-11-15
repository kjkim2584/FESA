#include "fesa/elements/incremental_hexa8.h"

#include <memory>
#include <new>

#include "model/model_storage.h"
namespace fesa {
IncrementalHexa8::IncrementalHexa8(int id, const SolidSection* section)
    : IncrementalSolidElement(id, section, SolidFormulation::kHexa8Full) {}
StatusCode IncrementalHexa8::Create(int id, const Node* const* nodes,
                                    int node_count, const SolidSection* section,
                                    SolidFormulation formulation,
                                    IncrementalHexa8** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("Incremental hexa8 output", d);
  *out = nullptr;
  if (!model_storage::ValidId(id) || node_count != 8 || !section ||
      !section->material() || formulation != SolidFormulation::kHexa8Full)
    return model_storage::Invalid("Incremental hexa8 definition", d);
  std::unique_ptr<IncrementalHexa8> element(new (std::nothrow)
                                                IncrementalHexa8(id, section));
  if (!element)
    return model_storage::AllocationFailure("Incremental hexa8 object", d);
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
