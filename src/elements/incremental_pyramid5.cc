#include "fesa/elements/incremental_pyramid5.h"

#include <memory>
#include <new>

#include "model/model_storage.h"
namespace fesa {
IncrementalPyramid5::IncrementalPyramid5(int id, const SolidSection* section)
    : IncrementalSolidElement(id, section, SolidFormulation::kPyramid5Full) {}
StatusCode IncrementalPyramid5::Create(int id, const Node* const* nodes,
    int node_count, const SolidSection* section, SolidFormulation formulation,
    IncrementalPyramid5** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("Incremental pyramid5 output", d);
  *out = nullptr;
  if (!model_storage::ValidId(id) || node_count != 5 || !section ||
      !section->material() || formulation != SolidFormulation::kPyramid5Full)
    return model_storage::Invalid("Incremental pyramid5 definition", d);
  std::unique_ptr<IncrementalPyramid5> element(
      new (std::nothrow) IncrementalPyramid5(id, section));
  if (!element)
    return model_storage::AllocationFailure("Incremental pyramid5 object", d);
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
