#include "fesa/model/thermal_condition.h"

#include <cmath>
#include <memory>
#include <new>

#include "model/model_storage.h"

namespace fesa {
const char* NormalizedThermalFaceLabel(int code) {
  if (code == -1) return "SNEG";
  if (code == -2) return "SPOS";
  static constexpr const char* labels[] = {"BF", "S1", "S2", "S3",
                                           "S4", "S5", "S6"};
  return code >= 0 && code <= 6 ? labels[code] : "";
}
ThermalConditionDefinition::ThermalConditionDefinition(
    ThermalConditionSpec* spec)
    : spec_(spec) {}
ThermalConditionDefinition::~ThermalConditionDefinition() {
  delete[] spec_->scope;
  delete[] spec_->name;
  delete[] spec_->face_label;
  delete[] spec_->source_file;
  delete[] spec_->expanded_element_ids;
  delete[] spec_->expanded_face_codes;
  delete spec_;
}
StatusCode ThermalConditionDefinition::Create(const ThermalConditionSpec* s,
                                              ThermalConditionDefinition** out,
                                              Diagnostic* d) {
  if (!out) return model_storage::Invalid("Thermal condition output", d);
  *out = nullptr;
  if (!s) return model_storage::Invalid("Thermal condition specification", d);
  const bool film = s->kind == ThermalConditionKind::kFilm;
  const bool radiation = s->kind == ThermalConditionKind::kRadiation;
  const bool distributed = radiation || film ||
                           s->kind == ThermalConditionKind::kBodyFlux ||
                           s->kind == ThermalConditionKind::kSurfaceFlux;
  if (s->kind != ThermalConditionKind::kTemperature &&
      s->kind != ThermalConditionKind::kConcentratedFlux && !distributed)
    return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                         "M2 thermal condition gate pending", d);
  const bool node = s->target == ThermalConditionTarget::kNode ||
                    s->target == ThermalConditionTarget::kElement;
  const bool set = s->target == ThermalConditionTarget::kNodeSet ||
                   s->target == ThermalConditionTarget::kElementSet ||
                   s->target == ThermalConditionTarget::kSurface;
  const bool element_target =
      s->target == ThermalConditionTarget::kElement ||
      s->target == ThermalConditionTarget::kElementSet ||
      s->target == ThermalConditionTarget::kSurface;
  if (!model_storage::ValidId(s->id) || !s->scope || !s->name ||
      !s->source_file || !*s->source_file || !s->face_label ||
      s->source_line <= 0 || s->step < 0 || (!node && !set) ||
      distributed != element_target ||
      (node && (!model_storage::ValidId(s->target_id) || *s->name)) ||
      (set && (s->target_id != 0 || !*s->name)) || s->first < 11 ||
      s->last < s->first || s->last > 30 || !std::isfinite(s->magnitude) ||
      (s->operation != ConditionOperation::kMod &&
       s->operation != ConditionOperation::kNew))
    return ReportFailure(StatusCode::kInputError, DiagnosticStage::kValidation,
                         "Invalid thermal condition row/label", d);
  if ((!distributed && *s->face_label) ||
      (!radiation && (s->emissivity != 0 || s->ambient_temperature != 0)) ||
      (!film && (s->film_property || s->film_amplitude || s->coefficient != 0 ||
                 s->sink_temperature != 0)))
    return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                         "M2 thermal condition properties gate pending", d);
  if (radiation && (!std::isfinite(s->ambient_temperature) ||
                    !std::isfinite(s->emissivity) || s->emissivity < 0 ||
                    s->emissivity > 1 || s->magnitude != 0))
    return ReportFailure(StatusCode::kInputError, DiagnosticStage::kValidation,
                         "Invalid radiation ambient/emissivity", d);
  if (film &&
      (!std::isfinite(s->sink_temperature) || !std::isfinite(s->coefficient) ||
       s->coefficient < 0 || s->magnitude != 0))
    return ReportFailure(StatusCode::kInputError, DiagnosticStage::kValidation,
                         "Invalid film sink/coefficient", d);
  if (distributed &&
      (s->expanded_count <= 0 || !s->expanded_element_ids ||
       !s->expanded_face_codes || !*s->face_label ||
       s->family != (s->target == ThermalConditionTarget::kSurface
                         ? (radiation ? ThermalLoadFamily::kSurfaceRadiation
                            : film    ? ThermalLoadFamily::kSurfaceFilm
                                      : ThermalLoadFamily::kSurfaceFlux)
                         : (radiation ? ThermalLoadFamily::kElementRadiation
                            : film    ? ThermalLoadFamily::kElementFilm
                                      : ThermalLoadFamily::kElementFlux))))
    return model_storage::Invalid("Thermal distributed expansion", d);
  std::unique_ptr<int[]> elements, faces;
  if (distributed) {
    elements.reset(new (std::nothrow) int[s->expanded_count]);
    faces.reset(new (std::nothrow) int[s->expanded_count]);
    if (!elements || !faces)
      return model_storage::AllocationFailure("Thermal face expansion", d);
    for (int i = 0; i < s->expanded_count; ++i) {
      const int face = s->expanded_face_codes[i];
      if (!model_storage::ValidId(s->expanded_element_ids[i]) ||
          (s->kind == ThermalConditionKind::kBodyFlux
               ? face != 0
               : face == 0 || face < -2 || face > 6))
        return model_storage::Invalid("Thermal distributed expansion entry", d);
      elements[i] = s->expanded_element_ids[i];
      faces[i] = face;
    }
  }
  std::unique_ptr<char[]> scope(model_storage::CopyText(s->scope));
  std::unique_ptr<char[]> name(model_storage::CopyText(s->name));
  std::unique_ptr<char[]> face(model_storage::CopyText(s->face_label));
  std::unique_ptr<char[]> file(model_storage::CopyText(s->source_file));
  std::unique_ptr<ThermalConditionSpec> copy(new (std::nothrow)
                                                 ThermalConditionSpec);
  if (!scope || !name || !face || !file || !copy)
    return model_storage::AllocationFailure("Thermal condition storage", d);
  copy->id = s->id;
  copy->kind = s->kind;
  copy->family = s->family;
  copy->expanded_count = distributed ? s->expanded_count : 0;
  copy->expanded_element_ids = elements.get();
  copy->expanded_face_codes = faces.get();
  copy->target = s->target;
  copy->target_id = s->target_id;
  copy->first = s->first;
  copy->last = s->last;
  copy->magnitude = s->magnitude;
  copy->sink_temperature = s->sink_temperature;
  copy->coefficient = s->coefficient;
  copy->emissivity = s->emissivity;
  copy->ambient_temperature = s->ambient_temperature;
  copy->film_property = s->film_property;
  copy->operation = s->operation;
  copy->amplitude = s->amplitude;
  copy->film_amplitude = s->film_amplitude;
  copy->source_line = s->source_line;
  copy->step = s->step;
  copy->scope = scope.get();
  copy->name = name.get();
  copy->face_label = face.get();
  copy->source_file = file.get();
  auto* result = new (std::nothrow) ThermalConditionDefinition(copy.get());
  if (!result) return model_storage::AllocationFailure("Thermal condition", d);
  copy.release();
  scope.release();
  name.release();
  face.release();
  file.release();
  elements.release();
  faces.release();
  *out = result;
  return StatusCode::kOk;
}
int ThermalConditionDefinition::Id() const { return spec_->id; }
const ThermalConditionSpec* ThermalConditionDefinition::specification() const {
  return spec_;
}
}  // namespace fesa
