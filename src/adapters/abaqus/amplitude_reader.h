/// @file
/// @brief Shared TABULAR row validation and immutable model construction.
#ifndef FESA_ADAPTERS_ABAQUS_AMPLITUDE_READER_H_
#define FESA_ADAPTERS_ABAQUS_AMPLITUDE_READER_H_
#include "fesa/io/input_deck.h"
#include "fesa/model/amplitude.h"
namespace fesa {
namespace abaqus {
/// Borrows a lexical AMPLITUDE block with validated options/scope. Validates
/// finite ordered pairs and physical row widths. Writes count only on success.
/// Failure d borrows source, which must outlive reporting; block unchanged.
StatusCode ValidateAmplitudeRows(const char* source, const InputBlock* block,
                                  int* count, Diagnostic* d);
/// Revalidates rows and copies the table into a caller-owned Amplitude.
/// Options/scope must have passed ReadInputDeck. Failure nulls *out and frees
/// partial storage; inputs borrowed only for this call. d borrows source.
StatusCode BuildAmplitude(const char* source, const InputBlock* block, int id,
                           Amplitude** out, Diagnostic* d);
}  // namespace abaqus
}  // namespace fesa
#endif  // FESA_ADAPTERS_ABAQUS_AMPLITUDE_READER_H_
