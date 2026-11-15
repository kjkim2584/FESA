/// @file
/// @brief Separate approved thermal lexical whitelist and block validation.
#ifndef FESA_ADAPTERS_ABAQUS_THERMAL_READER_H_
#define FESA_ADAPTERS_ABAQUS_THERMAL_READER_H_
#include <vector>

#include "fesa/io/input_deck.h"
#include "fesa/model/heat_transfer_step.h"
namespace fesa {
namespace abaqus {
/// Borrowed option list, null for an unsupported thermal keyword.
const char* ThermalAllowedOptions(const char* keyword);
/// Annotates scopes/steps and validates the thermal subset before publication.
/// Blocks remain caller-owned; failed annotations must be discarded with deck.
StatusCode ValidateThermalDeck(const char* source,
                               std::vector<InputBlock*>* blocks, Diagnostic* d);
/// Reads approved time settings and STEP overrides; input borrowed, out scalar
/// settings changed only on success. Does not resolve conditions or nodes.
StatusCode ReadThermalStep(const char* source, const InputBlock* step,
                           const InputBlock* heat,
                           HeatTransferStepParameters* out, Diagnostic* d);
}  // namespace abaqus
}  // namespace fesa
#endif  // FESA_ADAPTERS_ABAQUS_THERMAL_READER_H_
