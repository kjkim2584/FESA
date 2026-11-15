/// @file
/// @brief Applies validated ordered CONTROLS rows without reparsing the deck.
#ifndef FESA_ADAPTERS_ABAQUS_CONTROL_READER_H_
#define FESA_ADAPTERS_ABAQUS_CONTROL_READER_H_
#include "fesa/io/control_settings.h"
namespace fesa {
namespace abaqus {
/// Number/order follows FIELD rows: 8 first-row and 3 second-row scalars.
constexpr int kControlCount = 11;
/// Static case-sensitive cfg suffix, index must be in [0,kControlCount).
const char* ControlSymbol(int index);
/// Borrowed address of one scalar, same index range; field must be nonnull.
double* ControlValue(FieldControls* field, int index);
/// Explicit scalar configuration copy; no analysis vectors/flux are copied.
void CopyFieldControls(const FieldControls* from, FieldControls* to);
/// Copies scalar iteration counts explicitly; inputs borrowed, no allocation.
void CopyNewtonControls(const NewtonControls* from, NewtonControls* to);
/// Copies scalar increment settings explicitly, same ownership/lifetime.
void CopyIncrementControls(const IncrementControls* from,
                           IncrementControls* to);
/// Reads one TIME block from defaults, validates supported positions/domains.
/// Caller owns temporary outputs and publishes only on success. d optional.
/// dynamic rejects explicit static extrapolation DE and stability DT positions;
/// no HHT predictor/half-step accuracy behavior is implied by scalar parsing.
StatusCode ReadTimeControlsBlock(const char* source, const InputBlock* block,
                                 NewtonControls* counts,
                                 IncrementControls* increments, Diagnostic* d,
                                 bool dynamic = false);
/// Applies ordered TIME/RESET blocks to borrowed temporary settings only.
StatusCode ReadStepTimeControls(const InputDeck* deck, int step,
                                NewtonControls* counts,
                                IncrementControls* increments, Diagnostic* d);
/// Applies only this step's blocks into borrowed temporary fields. Callers
/// publish the pair only on success, so failure cannot alter live settings.
/// Diagnostic borrows deck's original path. May throw bad_alloc internally.
StatusCode ReadStepControls(const InputDeck* deck, int step,
                            FieldControls* displacement,
                            FieldControls* rotation, Diagnostic* d);
/// Applies validated thermal FIELD/RESET blocks to a temperature field only.
/// Borrows deck/out; failure preserves out. No mechanical field is modified.
/// q0/qu use power units. epsilon_d is not a thermal input or criterion.
StatusCode ReadThermalFieldControls(const InputDeck* deck, int step,
                                    FieldControls* out, Diagnostic* d);
}  // namespace abaqus
}  // namespace fesa
#endif  // FESA_ADAPTERS_ABAQUS_CONTROL_READER_H_
