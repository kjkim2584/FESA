/// @file
/// @brief H5 element contribution for linear Solid/MITC4 dynamics.
#ifndef FESA_ANALYSIS_LINEAR_DYNAMIC_ELEMENT_H_
#define FESA_ANALYSIS_LINEAR_DYNAMIC_ELEMENT_H_

#include "fesa/elements/element_workspace.h"
#include "fesa/model/element.h"

namespace fesa {
class AnalysisModel;
class DofManager;
class EffectiveConditions;
class StaticWorkspace;
/// Resolves analysis kDefault per element (row-sum when supported, otherwise
/// consistent), or validates an explicit request. Actual Element calls always
/// receive a physical mode. Borrowed inputs; output scalar preserved on failure.
StatusCode ResolveDynamicMassMode(const Element* element, MassMode requested,
                                  MassMode* out, Diagnostic* d);
/// Reads constant Rayleigh scalars of Solid/MITC4 materials, independent of
/// kinematics. Null/another element returns false. Borrows without computation.
bool HasMaterialRayleighDamping(const Element* element);
/// Requests CD from the physical E0 role at zero motion. All scratch is local;
/// caller owns matching square out. Input/geometry failures preserve or
/// invalidate output according to the underlying E0 contract. No retained data.
StatusCode ComputeLinearDynamicDamping(const Element* element,
                                       const ElementWorkspace* workspace,
                                       MassMode mode, Matrix* out,
                                       Diagnostic* d);
/// Validates one supported linear Solid(3N) or MITC4(6N) family model;
/// different supported Solid formulations/local sizes may coexist.
/// including density, mass capability, node layout and linear step kinematics.
/// Borrows model; no mutation/retained pointers. Writes shell only on success.
StatusCode ValidateLinearDynamicModel(const AnalysisModel* model, MassMode mode,
                                      bool* shell, Diagnostic* d);
/// Validates exclusive worker contexts. Solid requires distinct nonnull entries;
/// MITC4 permits null contexts/array because its linear operators ignore them.
/// Nonnull entries remain distinct. Positive count required. No mutation.
StatusCode ValidateLinearDynamicWorkers(bool shell, ElementWorkspace** workers,
                                        int count, Diagnostic* d);
/// Borrowed physical global arrays for one joined H5 reduced assembly call.
/// External is already (1+alpha)*p_next-alpha*p_previous [F, FL]; no reactions.
/// Motion arrays [L, rad] follow DofManager order; no pointer retained after
/// join.
struct HhtSolidAssemblyInput {
  const Vector* committed = nullptr;
  const Vector* predictor = nullptr;
  const Vector* external = nullptr;
  const Vector* previous_velocity = nullptr;  ///< Borrowed v_old [L/T,rad/T].
  const Vector* predicted_velocity = nullptr; ///< Borrowed a_new=0 predictor.
  double alpha = -.05;  ///< Dimensionless [-1/3,0].
  double c0 = 0;        ///< Positive acceleration derivative [1/T^2].
  double c1 = 0;        ///< Positive velocity derivative [1/T] for damping.
  MassMode mode = MassMode::kConsistent;
};
/// Checks supported linear Solid/MITC4 family, scalars, workers, and global
/// shapes/finite data/no workspace alias, including fully fixed models.
/// Borrows all pointers; no mutation or retained pointer on any outcome.
StatusCode ValidateHhtSolidAssembly(const AnalysisModel* model,
                                    const DofManager* dofs,
                                    const EffectiveConditions* conditions,
                                    StaticWorkspace* workspace,
                                    const HhtSolidAssemblyInput* dynamic,
                                    ElementWorkspace** workers,
                                    int worker_count, Diagnostic* d);
/// Borrows element/context and actual element-sized Cartesian vectors for this
/// call. Computes Ke=c0*M+(1+alpha)*(K+c1*C) and
/// be=alpha*K*committed+c0*M*predictor-Ke*prescribed [F, FL].
/// For nonzero Rayleigh, adds C*((1+alpha)*(c1*predictor-predicted_velocity)
/// +alpha*previous_velocity) to be. Both velocities are required then; otherwise
/// they may be null. C is requested from E0 and excludes element stabilization.
/// Motion components are [L, rad] in element node order. Matrix rows use
/// [F, FL] and columns
/// [L, rad]; Solid has translation only. MITC4 uses small additive rotations.
/// prescribed contains fixed values and zeros on free components; caller
/// establishes that mapping. External force is NOT included: caller assembles
/// (1+alpha)*p_next-alpha*p_previous separately. c0=1/(beta*dt^2)>0 [1/T^2],
/// alpha in [-1/3,0]; no finite rotation. Inputs remain unchanged.
/// Solid context scratch is overwritten; MITC4 permits null context. Caller
/// owns Ke/be; rhs must not alias an input. Invalid arguments preserve outputs;
/// subsequent numerical/allocation failures invalidate both outputs. No
/// state/result mutation, global assembly, caching, convergence decision or
/// commit.
/// @see docs/HHT/implementation_contract.md and formulation.md H4/H5.
StatusCode EvaluateHhtSolidContribution(
    const Element* element, ElementWorkspace* workspace, MassMode mode,
    double alpha, double c0, const Vector* committed, const Vector* predictor,
    const Vector* prescribed, Matrix* effective, Vector* rhs, Diagnostic* d,
    const Vector* previous_velocity = nullptr,
    const Vector* predicted_velocity = nullptr, double c1 = 0);
}  // namespace fesa
#endif  // FESA_ANALYSIS_LINEAR_DYNAMIC_ELEMENT_H_
