#include "elements/incremental_solid_recovery.h"

#include "elements/solid_recovery.h"
#include "elements/incremental_hexa8_workspace.h"
namespace fesa {
namespace {
thread_local void (*recovery_observer)() = nullptr;
}
void SetIncrementalRecoveryObserver(void (*observer)()) {
  recovery_observer = observer;
}
StatusCode RecoverIncrementalSolidPoints(
    const SolidReferenceData* data, const Matrix* const* candidate_f,
    const Matrix* candidate_stress, IncrementalElasticWorkspace* workspace,
    Matrix* ip_stress, Matrix* ip_strain, Matrix* nodal_stress,
    Matrix* nodal_strain, TrialEvaluationFailure* failure, Diagnostic* d) {
  if (failure) *failure = TrialEvaluationFailure::kNone;
  if (data && data->formulation() == SolidFormulation::kHexa8Full &&
      (!workspace || !workspace->hexa8_workspace() || candidate_f !=
       workspace->hexa8_workspace()->modified_f_data()))
    return ReportFailure(StatusCode::kInvalidArgument,DiagnosticStage::kValidation,
                         "Incremental Hexa8 recovery requires modified F",d);
  if (!failure || !workspace || !data || !candidate_f ||
      data->node_count() != workspace->node_count() ||
      data->point_count() != workspace->point_count() || !candidate_stress ||
      candidate_stress->Rows() != data->point_count() ||
      candidate_stress->Columns() != 6 || candidate_stress == ip_stress ||
      candidate_stress == ip_strain || candidate_stress == nodal_stress ||
      candidate_stress == nodal_strain)
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation,
                         "Incremental recovery input", d);
  auto status = ValidateSolidTensorRecovery(data, ip_stress, ip_strain,
                                            nodal_stress, nodal_strain, d);
  if (status != StatusCode::kOk) return status;
  for (int g = 0; g < data->point_count(); ++g) {
    const auto* f = candidate_f[g];
    if (!f || f->Rows() != 3 || f->Columns() != 3 || f == ip_stress ||
        f == ip_strain || f == nodal_stress || f == nodal_strain)
      return ReportFailure(StatusCode::kInvalidArgument,
                           DiagnosticStage::kValidation,
                           "Incremental recovery point F", d);
  }
  if (recovery_observer) recovery_observer();
  auto* left_cauchy_green = workspace->matrix(1);
  auto* logarithm = workspace->matrix(2);
  for (int g = 0; g < data->point_count(); ++g) {
    status = Matrix::Multiply(candidate_f[g], false, candidate_f[g], true,
                              left_cauchy_green, d);
    bool not_positive = false;
    if (status == StatusCode::kOk)
      status = Matrix::SymmetricPositiveDefiniteLog(left_cauchy_green, logarithm,
                                                    &not_positive, d);
    if (not_positive) *failure = TrialEvaluationFailure::kGeometry;
    if (status == StatusCode::kOk) status = logarithm->Scale(.5, d);
    const int packed[] = {0, 4, 8, 1, 2, 5};
    if (status == StatusCode::kOk)
      status = logarithm->GatherTo(packed, 6, workspace->tensor(), 0, d);
    if (status == StatusCode::kOk)
      status = ip_strain->SetRow(g, workspace->tensor(), d);
    if (status != StatusCode::kOk) {
      if (d) {
        const double* point = data->natural_coordinates()->Data() + 3 * g;
        d->geometry_point_available = true;
        d->natural_xi = point[0];
        d->natural_eta = point[1];
        d->natural_zeta = point[2];
      }
      return status;
    }
  }
  status = ip_stress->CopyFrom(candidate_stress, d);
  if (status == StatusCode::kOk)
    status = ExtrapolateSolidTensors(data, ip_stress, ip_strain, nodal_stress,
                                     nodal_strain, d);
  return status;
}
}  // namespace fesa
