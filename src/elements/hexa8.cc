#include "fesa/elements/hexa8.h"

#include <memory>
#include <new>

#include "model/model_storage.h"
namespace fesa {
Hexa8::Hexa8(int id, const SolidSection* section, SolidFormulation formulation)
    : SolidElement(id, 8,
                   formulation == SolidFormulation::kHexa8Reduced ? 1 : 8,
                   section, formulation) {}
StatusCode Hexa8::Create(int id, const Node* const* nodes, int node_count,
                         const SolidSection* section,
                         SolidFormulation formulation, Hexa8** out,
                         Diagnostic* d) {
  if (!out) return model_storage::Invalid("Hexa8 output", d);
  *out = nullptr;
  if (!model_storage::ValidId(id) || node_count != 8 || !section ||
      !section->material() ||
      (formulation != SolidFormulation::kHexa8Full &&
       formulation != SolidFormulation::kHexa8Reduced))
    return model_storage::Invalid("Hexa8 definition", d);
  std::unique_ptr<Hexa8> element(new (std::nothrow)
                                     Hexa8(id, section, formulation));
  if (!element) return model_storage::AllocationFailure("Hexa8 object", d);
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
