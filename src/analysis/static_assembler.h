/// @file
/// @brief Direct reduced assembly, E1; no full global stiffness allocation.
#ifndef FESA_ANALYSIS_STATIC_ASSEMBLER_H_
#define FESA_ANALYSIS_STATIC_ASSEMBLER_H_
#include <cstddef>

#include "analysis/linear_dynamic_element.h"
#include "fesa/analysis/analysis_model.h"
#include "fesa/analysis/dof_manager.h"
#include "fesa/analysis/equation_system.h"
#include "fesa/analysis/static_workspace.h"
namespace fesa {
/// Allocation-free conservative LP64/CSR-byte count update. Starts at n_free
/// diagonal entries and adds m*(m+1)/2 for each element's m free components.
/// Rejects before overflow; leaves *count unchanged on failure. Huge
/// overlapping models can fail even if vendor-coalesced CSR would fit.
StatusCode AccumulateAssemblyCount(int free_count, int element_free_count,
                                   std::size_t* count, Diagnostic* d);
/// Applies the same conservative count rule before file truncation and
/// assembly. No stiffness evaluation/singularity check; zero free DOFs bypass
/// safely. Optional layout outputs change only on success; default no-element
/// layout is 12 and beam-only. No extra numerical calculation is performed.
StatusCode ValidateAssemblyCounts(const AnalysisModel* model,
                                  const DofManager* dofs, Diagnostic* d,
                                  int* maximum_element_dofs = nullptr,
                                  bool* beam_only = nullptr);
/// Borrows inputs and overwrites workspace scratch; failed scratch is invalid.
/// Caller deletes successful *out; failure nulls it and preserves
/// model/results. FreeCount()==0 succeeds with null *out, creating no reduced
/// math objects.
/// Optional workers are borrowed through the join; null creates temporary
/// workers. Each worker writes disjoint chunk ranges; only the coordinator
/// scatters. Optional dynamic input selects H5 undamped linear C3D8
/// contributions and already-weighted external force. Borrows all arrays
/// through join, distinct from workspace buffers; caller owns and discards
/// failed scratch. Fully fixed dynamics still validate support/inputs before
/// null-system success.
StatusCode AssembleReducedSystem(
    const AnalysisModel* model, const DofManager* dofs,
    const EffectiveConditions* conditions, StaticWorkspace* workspace,
    EquationSystem** out, Diagnostic* d, ElementWorkspace** workers = nullptr,
    int worker_count = 1, const HhtSolidAssemblyInput* dynamic = nullptr);
}  // namespace fesa
#endif  // FESA_ANALYSIS_STATIC_ASSEMBLER_H_
