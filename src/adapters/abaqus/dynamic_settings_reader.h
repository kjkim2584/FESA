/// @file
/// @brief Resolves validated Abaqus HHT settings without starting analysis.
#ifndef FESA_ADAPTERS_ABAQUS_DYNAMIC_SETTINGS_READER_H_
#define FESA_ADAPTERS_ABAQUS_DYNAMIC_SETTINGS_READER_H_
#include "fesa/analysis/dynamic_step_settings.h"
#include "fesa/io/input_deck.h"
namespace fesa {
namespace abaqus {
/// Selects a one-based step from a validated deck, preserving original scalars.
/// DIRECT only; automatic/static procedures are unsupported here. Empty min/max
/// are absent, explicit 0 remains explicitly specified. No state/file changes.
/// Caller owns *out on success; creation-only output is null on failure.
/// Settings own their source path; failure diagnostic borrows deck's path.
StatusCode ReadDirectDynamicSettings(const InputDeck* deck, int step,
                                     DynamicStepSettings** out, Diagnostic* d);
/// Selects DIRECT when explicit, otherwise FESA AGGRESSIVE. Same borrowed
/// input/owned output/diagnostic lifetime contract as the DIRECT-only reader.
/// Preserves original min/max flags and resolves procedure-specific defaults.
StatusCode ReadDynamicSettings(const InputDeck* deck, int step,
                               DynamicStepSettings** out, Diagnostic* d);
}  // namespace abaqus
}  // namespace fesa
#endif  // FESA_ADAPTERS_ABAQUS_DYNAMIC_SETTINGS_READER_H_
