/// @file
/// @brief Abaqus-first convergence settings; no convergence calculation.
#ifndef FESA_IO_CONTROL_SETTINGS_H_
#define FESA_IO_CONTROL_SETTINGS_H_

#include <string>

#include "fesa/analysis/increment_controls.h"
#include "fesa/analysis/newton_controls.h"
#include "fesa/io/input_deck.h"
#include "fesa/io/run_options.h"

namespace fesa {
/// Scalars from docs/Linear static/abaqus_controls.md. Ratios are dimensionless.
/// q0/qu have field force or moment units. Defaults are automatic scales.
struct FieldControls {
  FieldControls() = default;
  FieldControls(const FieldControls&) = delete;
  FieldControls& operator=(const FieldControls&) = delete;
  FieldControls(FieldControls&&) = delete;
  FieldControls& operator=(FieldControls&&) = delete;
  double rn = .005;         ///< General residual ratio.
  double cn = .01;          ///< General correction ratio; 1 disables the check.
  double q0 = .01;          ///< First-step default only when not explicit.
  double qu = 0;            ///< Ignored unless qu_explicit is true.
  double rp = .02;          ///< Parsed/logged; unused by one direct solve.
  double epsilon = 1e-5;    ///< Zero-flux threshold.
  double cepsilon = .001;   ///< Zero-flux correction ratio.
  double rl = 1e-8;         ///< Linear residual ratio.
  double cf = 1;            ///< Field conversion coefficient; 0 disables.
  double epsilon_l = 1e-5;  ///< Active contribution threshold; 0 includes all.
  double epsilon_d = 1e-8;  ///< Displacement/Lc threshold; unused for rotation.
  bool q0_explicit = false;  ///< Otherwise use previous committed raw mean.
  bool qu_explicit = false;  ///< Otherwise use automatic average.
};
/// Source selection applies to the whole validated deck.
enum class ControlSource { kAbaqus, kConfig };
/// Global mass request. Default defers to each supported element's policy.
enum class MassSelection { kDefault, kConsistent, kLumped };
/// Whole-analysis selection. Automatic uses incremental nonlinear Solid and
/// preserves existing linear/Beam/Shell construction. kStvk remains for legacy
/// linear construction and historical code; nonlinear execution is deprecated.
enum class SolidConstitutiveModel { kStvk, kIncrementalElastic, kAutomatic };
/// Owns field/time settings and source path; borrowed access expires at
/// destruction. ApplyStepControls mutates only settings, never the deck or
/// committed flux.
struct ControlSettings {
  ControlSettings() = default;
  ~ControlSettings();
  ControlSettings(const ControlSettings&) = delete;
  ControlSettings& operator=(const ControlSettings&) = delete;
  ControlSettings(ControlSettings&&) = delete;
  ControlSettings& operator=(ControlSettings&&) = delete;
  FieldControls* displacement = nullptr;  ///< Owned, global X/Y/Z force field.
  FieldControls* rotation = nullptr;      ///< Owned, global X/Y/Z moment field.
  FieldControls* temperature = nullptr;   ///< Owned thermal field; null mechanically.
  NewtonControls* newton = nullptr;  ///< Owned; optional in manual linear use.
  IncrementControls* increments = nullptr;  ///< Owned; same lifetime.
  std::string* source_path = nullptr;     ///< Owned INP or selected cfg path.
  ControlSource source = ControlSource::kAbaqus;  ///< Whole-deck source.
  std::string* config_path = nullptr;  ///< Owned validated cfg path, if selected.
  MassSelection mass = MassSelection::kDefault;  ///< Independent of CONTROLS.
  int mass_line = 0;  ///< Physical cfg line; 0 means omitted/default generation.
  SolidConstitutiveModel solid_constitutive_model =
      SolidConstitutiveModel::kAutomatic;
  ///< Independent of step CONTROLS/reset and convergence source selection.
  int solid_constitutive_model_line = 0;  ///< Cfg line; zero means omitted.
  bool config_created = false;  ///< Missing cfg created this call.
};
/// Selects only after ReadInputDeck succeeds. Always validates explicit
/// cwd-relative cfg or exe/solver.cfg, creating missing files without parents.
/// With any CONTROLS, cfg convergence values are validated but not applied;
/// Abaqus defaults start the step sequence. Mass is always selected from cfg.
/// Caller owns *out; failure nulls it. Success preserves d. Failure d borrows
/// options->config_path or executable_directory (kDefaultConfigDirectory);
/// these bases must outlive diagnostic reporting.
StatusCode SelectInitialControls(const InputDeck* deck,
                                 const RunOptions* options,
                                 const char* executable_directory,
                                 ControlSettings** out, Diagnostic* d);
/// Applies ordered blocks for the requested positive step. Absent blocks
/// propagate; FIELD blanks reset that field; TIME blanks reset time defaults.
/// RESET restores field/time defaults. TIME requires both owned time pointers.
/// Caller supplies steps in execution order. Failure preserves all settings;
/// d borrows deck.SourceFile(). Automatic q0 remains a policy, never a frozen
/// analysis mean. Success preserves d.
StatusCode ApplyStepControls(const InputDeck* deck, int step_number,
                             ControlSettings* settings, Diagnostic* d);
/// Deep clone of field settings, optional time settings and source metadata.
/// Caller owns *out; failure nulls it and preserves source. Source path must be
/// non-null.
StatusCode CloneControlSettings(const ControlSettings* source,
                                ControlSettings** out, Diagnostic* d);
}  // namespace fesa
#endif  // FESA_IO_CONTROL_SETTINGS_H_
