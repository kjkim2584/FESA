#include "elements/solid_element_read_view.h"
namespace fesa {
StatusCode SolidElementReadView::Bind(const Element* element,
    SolidElementReadView* out, Diagnostic* d) {
  if (!out || !element)
    return ReportFailure(StatusCode::kInvalidArgument,
        DiagnosticStage::kValidation, "Solid read view input", d);
  const auto* legacy = dynamic_cast<const SolidElement*>(element);
  const auto* incremental =
      dynamic_cast<const IncrementalSolidElement*>(element);
  if (!legacy && !incremental)
    return ReportFailure(StatusCode::kUnsupported,
        DiagnosticStage::kValidation, "Solid read view family", d);
  out->legacy_ = legacy;
  out->incremental_ = incremental;
  return StatusCode::kOk;
}
const Element* SolidElementReadView::element() const {
  return legacy_ ? static_cast<const Element*>(legacy_)
                 : static_cast<const Element*>(incremental_);
}
SolidFormulation SolidElementReadView::formulation() const {
  return legacy_ ? legacy_->formulation() : incremental_->formulation();
}
int SolidElementReadView::node_count() const {
  return legacy_ ? legacy_->node_count() : incremental_->node_count();
}
int SolidElementReadView::integration_point_count() const {
  return legacy_ ? legacy_->integration_point_count() : incremental_->integration_point_count();
}
const SolidSection* SolidElementReadView::section() const {
  return legacy_ ? legacy_->section() : incremental_->section();
}
const ElasticMaterial* SolidElementReadView::material() const {
  return legacy_ ? legacy_->material() : incremental_->material();
}
const Vector* SolidElementReadView::internal_force() const {
  return legacy_ ? legacy_->internal_force() : incremental_->internal_force();
}
const Matrix* SolidElementReadView::ip_stress() const {
  return legacy_ ? legacy_->ip_stress() : incremental_->ip_stress();
}
const Matrix* SolidElementReadView::ip_strain() const {
  return legacy_ ? legacy_->ip_strain() : incremental_->ip_strain();
}
const Matrix* SolidElementReadView::nodal_stress() const {
  return legacy_ ? legacy_->nodal_stress() : incremental_->nodal_stress();
}
const Matrix* SolidElementReadView::nodal_strain() const {
  return legacy_ ? legacy_->nodal_strain() : incremental_->nodal_strain();
}
const Node* SolidElementReadView::GetNode(int local) const {
  return element()->GetNode(local);
}
}  // namespace fesa
