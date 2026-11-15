/// @file
/// @brief Shared linear/nonlinear F11 interpolation of converted tensors.
#ifndef FESA_ELEMENTS_SOLID_RECOVERY_H_
#define FESA_ELEMENTS_SOLID_RECOVERY_H_

#include "fesa/elements/solid_reference_data.h"

namespace fesa {
/// Borrows all objects, validates stiffness-purpose reference, ng x6 IP and
/// n x6 nodal shapes and pairwise disjoint storage. Does not read numerical
/// payload or mutate any matrix. d is optional failure metadata.
StatusCode ValidateSolidTensorRecovery(const SolidReferenceData* data,
    const Matrix* ip_stress, const Matrix* ip_strain,
    const Matrix* nodal_stress, const Matrix* nodal_strain, Diagnostic* d);
/// Borrows reference and finite IP tensors for this call. Caller owns disjoint
/// n x6 nodal outputs. Preserves component order, units, basis and time;
/// performs no stress/strain transformation. F11 spaces follow the reference.
/// One call-local LU serves twelve RHS, except single-IP constant replication.
/// Preflight failure preserves outputs; computation failure invalidates both
/// nodal outputs. No shared-node averaging, state commit or retained cache.
StatusCode ExtrapolateSolidTensors(const SolidReferenceData* data,
    const Matrix* ip_stress, const Matrix* ip_strain,
    Matrix* nodal_stress, Matrix* nodal_strain, Diagnostic* d);
}  // namespace fesa
#endif  // FESA_ELEMENTS_SOLID_RECOVERY_H_
