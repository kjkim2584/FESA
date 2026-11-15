/// @file
/// @brief Resolves source-preserving input into an independently owned Domain.
#ifndef FESA_ADAPTERS_ABAQUS_MODEL_BUILDER_H_
#define FESA_ADAPTERS_ABAQUS_MODEL_BUILDER_H_
#include <map>

#include "fesa/io/control_settings.h"
#include "fesa/io/input_deck.h"
#include "fesa/model/domain.h"
namespace fesa {
namespace abaqus {
/// Shared single-instance index-set expansion after node/element IDs validate.
/// Adopts created sets into candidate Domain only; failure discards candidate.
StatusCode BuildInputIndexSets(const InputDeck* deck, Domain* domain,
    const std::map<int, int>* node_lines,
    const std::map<int, int>* element_lines, Diagnostic* d);
}  // namespace abaqus
/// Resolves an omitted selection without creating a Domain. Option lookup may
/// allocate; allocation failure returns kAllocationFailure without exceptions.
/// Nonlinear Solid uses incremental elastic; other models retain their linear
/// construction. Explicit nonlinear StVK is rejected. Failure preserves *out;
/// d borrows deck's source path. out is required and may alias the input scalar.
StatusCode ResolveSolidConstitutiveModel(const InputDeck* deck,
    SolidConstitutiveModel requested, SolidConstitutiveModel* out, Diagnostic* d);
/// Resolves references and validates static/linear DIRECT model data without
/// solver/output IO. Initial velocity remains in InputDeck for separate
/// resolution. Model success does not certify dynamic execution; FESASolver
/// gates unsupported execution. Caller owns *out; failure nulls it and releases
/// only the candidate Domain. Existing models are untouched. out is
/// creation-only. d borrows deck's path; deck must outlive diagnostic
/// reporting. Success leaves d unchanged.
/// Constitutive selection is fixed before allocation. Incremental mode rejects
/// unsupported TYPE/step before returning any Domain; failure nulls out.
/// Automatic preserves linear solids; nonlinear StVK is no longer supported.
StatusCode BuildDomain(const InputDeck* deck, Domain** out, Diagnostic* d,
    SolidConstitutiveModel model = SolidConstitutiveModel::kAutomatic);
}  // namespace fesa
#endif  // FESA_ADAPTERS_ABAQUS_MODEL_BUILDER_H_
