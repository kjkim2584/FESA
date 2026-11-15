/// @file
/// @brief Owned HHT time settings independent of input deck lifetime.
#ifndef FESA_ANALYSIS_DYNAMIC_STEP_SETTINGS_H_
#define FESA_ANALYSIS_DYNAMIC_STEP_SETTINGS_H_

#include <string>

#include "fesa/analysis/hht_integrator.h"
#include "fesa/model/dynamic_step_parameters.h"

namespace fesa {
/// Solver-owned immutable settings; outlives borrowing analysis/integrator.
/// Owns scalar input copies, source path and HHT settings, no numerical
/// buffers.
class DynamicStepSettings {
 public:
  /// Resolves DIRECT h/T/min/max using the documented DIRECT contract.
  /// Borrows inputs only for this call. Copies source_file (UTF-8) and scalars.
  /// source_line/step are positive original DYNAMIC line/one-based step.
  /// Caller owns *out; creation-only output is null on failure. Inputs
  /// unchanged. Does not validate element capabilities or initialize physical
  /// state.
  static StatusCode CreateDirect(const DynamicStepParameters* parameters,
                                 const char* source_file, int source_line,
                                 int step, DynamicStepSettings** out,
                                 Diagnostic* d);
  /// Resolves AGGRESSIVE min/max defaults and first bounded nominal.
  /// Same ownership/failure/source contract as CreateDirect. Original inputs
  /// remain available in Parameters; resolved limits have separate getters.
  /// @see docs/HHT/automatic_increment_control.md.
  static StatusCode CreateAutomatic(const DynamicStepParameters* parameters,
                                    const char* source_file, int source_line,
                                    int step, DynamicStepSettings** out,
                                    Diagnostic* d);
  ~DynamicStepSettings();
  DynamicStepSettings(const DynamicStepSettings&) = delete;
  DynamicStepSettings& operator=(const DynamicStepSettings&) = delete;
  DynamicStepSettings(DynamicStepSettings&&) = delete;
  DynamicStepSettings& operator=(DynamicStepSettings&&) = delete;
  /// Borrowed original scalars/flags, valid until destruction.
  const DynamicStepParameters* Parameters() const;
  /// Borrowed immutable coefficients; must outlive the HHT integrator.
  const HhtParameters* IntegratorParameters() const;
  /// DIRECT min(h,T,nonzero max); automatic min(h,resolved max), in T units.
  /// Automatic endpoint duration is clipped separately to preserve valid
  /// nominal metadata even when T<minimum.
  double NominalIncrement() const;
  /// True for convergence-history AGGRESSIVE; false for DIRECT.
  bool IsAutomatic() const;
  /// Resolved automatic lower/upper limits in T units; DIRECT original limits.
  double MinimumIncrement() const;
  double MaximumIncrement() const;
  /// Borrowed owned copy; does not borrow InputDeck.
  const char* SourceFile() const;
  /// Original DYNAMIC keyword line.
  int SourceLine() const;
  /// One-based step number.
  int Step() const;

 private:
  DynamicStepSettings() = default;
  /// Shared factory implementation; all inputs borrowed for the call only.
  static StatusCode Create(const DynamicStepParameters* parameters,
                           const char* source_file, int source_line, int step,
                           bool automatic, DynamicStepSettings** out,
                           Diagnostic* d);
  DynamicStepParameters* parameters_ = nullptr;  ///< Owned scalar input copy.
  HhtParameters* hht_ = nullptr;                 ///< Owned coefficients.
  std::string* source_file_ = nullptr;           ///< Owned path.
  double nominal_ = 0;
  double minimum_ = 0, maximum_ = 0;
  bool automatic_ = false;
  int source_line_ = 0;
  int step_ = 0;
};
}  // namespace fesa
#endif  // FESA_ANALYSIS_DYNAMIC_STEP_SETTINGS_H_
