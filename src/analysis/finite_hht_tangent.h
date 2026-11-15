/// @file
/// @brief H6 chain rule for local finite-rotation HHT correction tangents.
#ifndef FESA_ANALYSIS_FINITE_HHT_TANGENT_H_
#define FESA_ANALYSIS_FINITE_HHT_TANGENT_H_

#include "fesa/analysis/time_integrator.h"
#include "fesa/model/element_dynamic.h"

namespace fesa {
/// Writes K=(1+alpha)*(Kt+GD+CD*H)+GI+CI*H+MI*B in full general storage.
/// Borrows local connectivity-ordered6n-square structural/physical derivatives,
/// n-by-9 R2 angular blocks, and positive Cartesian c1/c0 [1/T,1/T^2]. H/B
/// combine these node diagonal blocks in spatial [du,dw] order. alpha is finite
/// in [-1/3,0]; all coefficients/blocks must refer to the same trial and dt.
/// Forces in physical are ignored; all five derivative pointers are required.
/// The caller decides supported elements/damping and sparse MatrixType; no
/// symmetry is imposed. No state, result recovery, solve or retained scratch.
/// Input-only aliasing is allowed. out is distinct from every input matrix.
/// Failure preserves out until final CopyFrom (whose backend failure may
/// invalidate it). All numerical scratch is call-local. @see
/// docs/HHT/rotational_integrator.md and implementation_contract.md.
StatusCode ComposeFiniteHhtTangent(
    double alpha, const TranslationalLinearization* translation,
    const RotationalLinearization* rotation, const Matrix* structural,
    const ElementDynamicResult* physical, Matrix* out, Diagnostic* d);
/// Beam retained-eta counterpart of ComposeFiniteHhtTangent. Borrows two
/// physical nodes plus eta_a/eta_b, exactly14-square blocks and2-by-9 R2
/// blocks. Coordinate columns are [L,rad,L,rad,L,L]; force rows
/// [F,FL,F,FL,F,F]. H/B use Cartesian c1/c0 on translation and eta; physical
/// rotations use R2. Same call-only ownership, alias, failure and scratch
/// contract as above.
StatusCode ComposeBeamFiniteHhtTangent(
    double alpha, const TranslationalLinearization* translation,
    const RotationalLinearization* rotation, const Matrix* structural,
    const ElementDynamicResult* physical, Matrix* out, Diagnostic* d);
}  // namespace fesa
#endif  // FESA_ANALYSIS_FINITE_HHT_TANGENT_H_
