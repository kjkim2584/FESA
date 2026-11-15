/// @file
/// @brief IE2-IE4 Cauchy update and exact candidate PK1 Jacobian.
#ifndef FESA_ELEMENTS_INCREMENTAL_ELASTIC_POINT_H_
#define FESA_ELEMENTS_INCREMENTAL_ELASTIC_POINT_H_
#include "elements/incremental_elastic_workspace.h"
#include "fesa/model/elastic_material.h"
#include "fesa/model/element_trial.h"
namespace fesa {
/// Borrows fixed committed Fn/σn and candidate F3×3 in global axes. Stress/PK1
/// are disjoint caller-owned3×3 in F/L²; optional dp_df is9×9 row-major pairs
/// 3*i+A,3*k+B. All input/output/scratch storage is disjoint. Workspace and
/// failure are required. Preflight preserves outputs; computational failure
/// invalidates requested outputs. No official result/history/commit change.
/// Committed invalidity is fatal with kNone; candidate geometry singularity
/// reports kGeometry. Backend/resource errors never become geometry errors.
/// Null dp_df skips IE4. @see incremental_elastic_formulation.md IE2-IE4.
StatusCode EvaluateIncrementalElasticPoint(const Matrix* committed_f,
    const Matrix* candidate_f, const Matrix* committed_stress,
    const ElasticMaterial* material, IncrementalElasticWorkspace* workspace,
    Matrix* stress, Matrix* pk1, Matrix* dp_df,
    TrialEvaluationFailure* failure, Diagnostic* d);
}  // namespace fesa
#endif  // FESA_ELEMENTS_INCREMENTAL_ELASTIC_POINT_H_
