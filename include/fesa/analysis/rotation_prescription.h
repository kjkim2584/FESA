/// @file
/// @brief Fixed-step Abaqus-style accumulated rotation-command path.
#ifndef FESA_ANALYSIS_ROTATION_PRESCRIPTION_H_
#define FESA_ANALYSIS_ROTATION_PRESCRIPTION_H_
#include "fesa/math/vector.h"
namespace fesa {
/// q_target=Exp(lambda*(command_end-command_start))*q_step_start.
/// Borrowed unit wxyz q_step_start[4], global rad commands[3]. Continued
/// constraints retain previous INPUT command; first addition uses principal
/// Log(q_step_start) as command_start, the explicit FESA supplement.
/// Snapshots remain fixed throughout this step and every cutback.
/// Requires0<=committed_lambda<=lambda<=1. Guards the known unwrapped increment
/// norm((lambda-committed_lambda)*Delta)<pi with relative-scalar256epsilon.
/// Thus one-step2pi/4pi fails even if endpoint q is identity; many smaller
/// accepted increments may make a full turn. This helper retains quaternion
/// cover. RotationState separately guards actual committed-relative endpoint.
/// Inputs/outputs distinct, readonly inputs may alias each other. Caller owns
/// out[4]; required branch flag cleared on entry. Preflight/branch failure
/// preserves out; later backend failure invalidates out. No state is retained
/// or modified, d optional. Branch returns kEquilibriumFailure, others retain
/// invalid-input/resource/library status. Never add commands to outputUR.
StatusCode EvaluateRotationPrescription(const Vector* q_step_start,
                                        const Vector* command_start,
                                        const Vector* command_end,
                                        double lambda, double committed_lambda,
                                        Vector* out, bool* branch_exceeded,
                                        Diagnostic* d);
}  // namespace fesa
#endif  // FESA_ANALYSIS_ROTATION_PRESCRIPTION_H_
