/// @file
/// @brief Thermal model preparation independent of mechanical model creation.
#ifndef FESA_ADAPTERS_ABAQUS_THERMAL_MODEL_BUILDER_H_
#define FESA_ADAPTERS_ABAQUS_THERMAL_MODEL_BUILDER_H_
#include <cstddef>
#include <vector>

#include "fesa/io/input_deck.h"
#include "fesa/model/domain.h"
namespace fesa {
/// Reads one logical temperature IC using the original target resolver.
/// Borrows immutable inputs; participating is sorted CollectThermalNodeIds.
/// First row carries target and at most seven values; continuation rows hold
/// eight values each (last holds the remainder). UNIFORM has one value.
/// Success swaps caller-owned ids/values and writes next row/first source line.
/// Any failure preserves all outputs. No deck/model/file mutation or retention.
/// Empty participating targets with seven values are explicitly unsupported;
/// shorter rows are consumed alone for the existing empty-IC warning.
StatusCode ReadThermalInitialRow(const InputDeck* deck, const Domain* domain,
                                 const InputBlock* block, std::size_t first_row,
                                 const std::vector<int>* participating,
                                 std::vector<int>* ids,
                                 std::vector<double>* values,
                                 std::size_t* next_row, int* source_line,
                                 Diagnostic* d);
/// Counts one validated IC/BC/CFLUX row through the original Target resolver.
/// participating must be the sorted CollectThermalNodeIds result for domain.
/// Borrows inputs only for this call, retains none; no model/file mutation.
/// Success writes total/included. Any failure preserves both scalar outputs,
/// including allocation failure. Diagnostic borrows the deck source lifetime.
StatusCode CountThermalTargetParticipation(
    const InputDeck* deck, const Domain* domain, const InputBlock* block,
    const InputLine* row, const std::vector<int>* participating, int* total,
    int* included, Diagnostic* d);
/// Builds caller-owned thermal Domain from validated deck; nulls *out on
/// failure and releases partial objects. No outputs or references are written.
/// Failure diagnostic borrows deck source, which must outlive reporting.
StatusCode BuildThermalDomain(const InputDeck* deck, Domain** out,
                              Diagnostic* d);
}  // namespace fesa
#endif  // FESA_ADAPTERS_ABAQUS_THERMAL_MODEL_BUILDER_H_
