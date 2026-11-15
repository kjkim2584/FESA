#include "fesa/elements/tetra10.h"

#include <memory>
#include <new>

#include "model/model_storage.h"
namespace fesa {
Tetra10::Tetra10(int id, const SolidSection* section,
                 SolidFormulation formulation)
    : SolidElement(id, 10, 4, section, formulation) {}
StatusCode Tetra10::Create(int id, const Node* const* nodes, int node_count,
                           const SolidSection* section,
                           SolidFormulation formulation, Tetra10** out,
                           Diagnostic* d) {
  if (!out) return model_storage::Invalid("Tetra10 output", d);
  *out = nullptr;
  if (formulation == SolidFormulation::kTetra10Reduced)
    return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                         "FESA_C3D10R retired; use C3D10", d);
  if (!model_storage::ValidId(id) || node_count != 10 || !section ||
      !section->material() ||
      formulation != SolidFormulation::kTetra10Full)
    return model_storage::Invalid("Tetra10 definition", d);
  std::unique_ptr<Tetra10> element(new (std::nothrow)
                                       Tetra10(id, section, formulation));
  if (!element) return model_storage::AllocationFailure("Tetra10 object", d);
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
