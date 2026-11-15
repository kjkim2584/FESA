/// @file
/// @brief Original thermal condition identity and owned input spelling.
#ifndef FESA_MODEL_THERMAL_CONDITION_H_
#define FESA_MODEL_THERMAL_CONDITION_H_

#include "fesa/model/load_definition.h"

namespace fesa {
class FilmPropertyTable;
/// Thermal condition kinds supported by the thermal path.
enum class ThermalConditionKind {
  kTemperature,
  kConcentratedFlux,
  kBodyFlux,
  kSurfaceFlux,
  kFilm,
  kRadiation
};
/// Preserves original target identity before set/surface expansion.
enum class ThermalConditionTarget {
  kNode,
  kNodeSet,
  kElement,
  kElementSet,
  kSurface
};
/// Independent keyword histories; spatial aliases do not merge OP=NEW groups.
enum class ThermalLoadFamily {
  kElementFlux,
  kSurfaceFlux,
  kElementFilm,
  kSurfaceFilm,
  kElementRadiation,
  kSurfaceRadiation
};
/// Borrowed static spelling for validated expansion code BF0/S1..6/NEG-1/POS-2.
/// Unknown code returns an empty spelling; no state or storage ownership.
const char* NormalizedThermalFaceLabel(int code);
/// Borrowed input definition for one original row. No resource ownership.
/// Class pointers/strings/expansion arrays must remain valid through Create.
/// Film properties borrow Domain; radiation uses scalar ambient/emissivity.
struct ThermalConditionSpec {
  ThermalConditionSpec() = default;
  ThermalConditionSpec(const ThermalConditionSpec&) = delete;
  ThermalConditionSpec& operator=(const ThermalConditionSpec&) = delete;
  ThermalConditionSpec(ThermalConditionSpec&&) = delete;
  ThermalConditionSpec& operator=(ThermalConditionSpec&&) = delete;
  int id = 0;
  ThermalConditionKind kind = ThermalConditionKind::kTemperature;
  ThermalConditionTarget target = ThermalConditionTarget::kNode;
  int target_id = 0;  ///< Direct node/element ID, zero for a named target.
  const char* scope = "";
  const char* name = "";
  int first = 11;
  int last = 11;
  /// Temperature[theta], concentrated power[E/T], or BF/face power density.
  double magnitude = 0;
  const char* face_label = "";
  ThermalLoadFamily family = ThermalLoadFamily::kElementFlux;
  /// Paired topology expansion borrowed during Create, then copied/owned by
  /// the immutable definition. Original target/label above remain the key.
  /// Codes: BF0, S1..S6=1..6, SNEG=-1, SPOS=-2. No Node/element pointers.
  int expanded_count = 0;
  const int* expanded_element_ids = nullptr;
  const int* expanded_face_codes = nullptr;
  double sink_temperature = 0;
  double coefficient = 0;
  double emissivity = 0;
  double ambient_temperature = 0;  ///< Radiation reference ambient [theta].
  const FilmPropertyTable* film_property = nullptr;
  ConditionOperation operation = ConditionOperation::kMod;
  /// Temperature/CFLUX/flux amplitude, film SINK AMPLITUDE or radiation AMBIENT
  /// AMPLITUDE.
  const Amplitude* amplitude = nullptr;
  const Amplitude* film_amplitude = nullptr;  ///< Scalar h only; table ignores.
  const char* source_file = "";
  int source_line = 0;  ///< One-based input line.
  int step = 0;         ///< Zero denotes a model-level definition.
};
/// Immutable original row; owns spec, four text buffers and two topology
/// arrays. Amplitude/film property are borrowed from Domain and must outlive
/// this.
class ThermalConditionDefinition final {
 public:
  /// Borrows spec for this call, copies text/scalars and expansion arrays.
  /// Caller owns success; failure nulls *out and cleans partial allocations.
  /// Nodal labels11..30, distributed BF/face flux and film/radiation are
  /// supported. Physical constants are checked by step validation.
  static StatusCode Create(const ThermalConditionSpec* spec,
                           ThermalConditionDefinition** out, Diagnostic* d);
  ~ThermalConditionDefinition();
  ThermalConditionDefinition(const ThermalConditionDefinition&) = delete;
  ThermalConditionDefinition& operator=(const ThermalConditionDefinition&) =
      delete;
  ThermalConditionDefinition(ThermalConditionDefinition&&) = delete;
  ThermalConditionDefinition& operator=(ThermalConditionDefinition&&) = delete;
  /// Original positive definition ID, distinct from target identity.
  int Id() const;
  /// Borrowed immutable spec/text/topology until destruction. No evaluation.
  const ThermalConditionSpec* specification() const;

 private:
  explicit ThermalConditionDefinition(ThermalConditionSpec* spec);
  ThermalConditionSpec*
      spec_;  ///< Owned; owns copied scope/name/face/file and topology arrays.
};
}  // namespace fesa

#endif  // FESA_MODEL_THERMAL_CONDITION_H_
