#include "fesa/elements/pyramid5.h"

#include <memory>
#include <new>

#include "model/model_storage.h"
namespace fesa {
Pyramid5::Pyramid5(int id, const SolidSection* section,
                   SolidFormulation formulation)
    : SolidElement(id, 5,
                   formulation == SolidFormulation::kPyramid5Full ? 5 : 1,
                   section, formulation) {}
StatusCode Pyramid5::Create(int id, const Node* const* nodes, int node_count,
                            const SolidSection* section,
                            SolidFormulation formulation, Pyramid5** out,
                            Diagnostic* d) {
  if (!out) return model_storage::Invalid("Pyramid5 output", d);
  *out = nullptr;
  if (!model_storage::ValidId(id) || node_count != 5 || !section ||
      !section->material() ||
      (formulation != SolidFormulation::kPyramid5Full &&
       formulation != SolidFormulation::kPyramid5Reduced))
    return model_storage::Invalid("Pyramid5 definition", d);
  std::unique_ptr<Pyramid5> element(new (std::nothrow)
                                        Pyramid5(id, section, formulation));
  if (!element) return model_storage::AllocationFailure("Pyramid5 object", d);
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
