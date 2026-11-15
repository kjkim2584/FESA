#include "elements/beam3d_dynamic_metrics.h"

#include <cmath>
#include <memory>
#include <new>

namespace fesa {
namespace {

#define FESA_METRIC_TRY(expression)                \
  do {                                            \
    const StatusCode status = (expression);        \
    if (status != StatusCode::kOk) return status;   \
  } while (false)

constexpr BeamKinematicField kFields[] = {
    BeamKinematicField::kPosition, BeamKinematicField::kTangent,
    BeamKinematicField::kNormal1, BeamKinematicField::kNormal2,
    BeamKinematicField::kStrain, BeamKinematicField::kStretch};

StatusCode Invalid(const char* operation, Diagnostic* d) {
  return ReportFailure(StatusCode::kInvalidArgument,
                       DiagnosticStage::kValidation, operation, d);
}

bool Shape(const Matrix* matrix, int rows, int columns) {
  return matrix && matrix->Rows() == rows && matrix->Columns() == columns;
}

bool PointAlias(const BeamKinematicWorkspace* point, const Matrix* output) {
  if (!point || !output) return false;
  for (const auto field : kFields)
    if (output == point->Value(field) || output == point->Variation(field) ||
        output == point->Directional(field) || output == point->Mixed(field))
      return true;
  return false;
}

StatusCode CreateMatrix(int rows, int columns, std::unique_ptr<Matrix>* owner,
                        Diagnostic* d) {
  Matrix* raw = nullptr;
  const auto status = Matrix::Create(rows, columns, &raw, d);
  owner->reset(raw);
  return status;
}

StatusCode ValidateOutputs(const BeamKinematicInput* pose,
                           const BeamKinematicWorkspace* point,
                           const Matrix* mass, const Matrix* dm,
                           const Matrix* damping, const Matrix* dc,
                           Diagnostic* d) {
  if (!Shape(mass, 14, 14) || !Shape(damping, 14, 14) ||
      ((dm == nullptr) != (dc == nullptr)) ||
      (dm && (!Shape(dm, 14, 14) || !Shape(dc, 14, 14))))
    return Invalid("BeamDynamicMetrics output shape/pair", d);
  const Matrix* outputs[] = {mass, dm, damping, dc};
  for (int i = 0; i < 4; ++i) {
    if (!outputs[i]) continue;
    for (int j = 0; j < i; ++j)
      if (outputs[i] == outputs[j])
        return Invalid("BeamDynamicMetrics output alias", d);
    if (outputs[i] == pose->positions || outputs[i] == pose->tangents ||
        outputs[i] == pose->normals || PointAlias(point, outputs[i]))
      return Invalid("BeamDynamicMetrics borrowed input alias", d);
  }
  return StatusCode::kOk;
}

// Bz is the Cauchy axial/fibre rate packet, not the Green strain B.
// Ordered K1 mixed rows directly differentiate the fixed spatial basis.
StatusCode SpatialRows(const BeamKinematicWorkspace* point, Matrix* bz,
                       Matrix* dbz, Matrix* bt, Matrix* dbt, Matrix* row,
                       Matrix* drow, Diagnostic* d) {
  const auto stretch = BeamKinematicField::kStretch;
  const auto strain = BeamKinematicField::kStrain;
  const double lambda = point->Value(stretch)->Data()[0];
  const Matrix* bl = point->Variation(stretch);
  const Matrix* dbl = point->Mixed(stretch);
  const double dl = dbz ? point->Directional(stretch)->Data()[0] : 0;
  const double inverse = 1 / lambda;
  FESA_METRIC_TRY(row->CopyFrom(bl, d));
  FESA_METRIC_TRY(row->Scale(inverse, d));
  FESA_METRIC_TRY(bz->CopyBlockFrom(row, 0, 0, 1, 14, 0, 0, d));
  if (dbz) {
    FESA_METRIC_TRY(drow->CopyFrom(dbl, d));
    FESA_METRIC_TRY(drow->Scale(inverse, d));
    FESA_METRIC_TRY(drow->AddScaled(-dl * inverse * inverse, bl, d));
    FESA_METRIC_TRY(dbz->CopyBlockFrom(drow, 0, 0, 1, 14, 0, 0, d));
  }
  for (int i = 1; i < 3; ++i) {
    const double curvature = point->Value(strain)->Data()[i];
    const double inverse2 = inverse * inverse;
    const double inverse3 = inverse2 * inverse;
    FESA_METRIC_TRY(row->CopyBlockFrom(point->Variation(strain), i, 0, 1,
                                      14, 0, 0, d));
    if (dbz) {
      const double dk = point->Directional(strain)->Data()[i];
      FESA_METRIC_TRY(drow->CopyBlockFrom(point->Mixed(strain), i, 0, 1,
                                         14, 0, 0, d));
      FESA_METRIC_TRY(drow->Scale(inverse2, d));
      FESA_METRIC_TRY(drow->AddScaled(-2 * dl * inverse3, row, d));
      FESA_METRIC_TRY(drow->AddScaled(
          -2 * dk * inverse3 + 6 * curvature * dl * inverse3 * inverse, bl, d));
      FESA_METRIC_TRY(drow->AddScaled(-2 * curvature * inverse3, dbl, d));
      FESA_METRIC_TRY(dbz->CopyBlockFrom(drow, 0, 0, 1, 14, i, 0, d));
    }
    FESA_METRIC_TRY(row->Scale(inverse2, d));
    FESA_METRIC_TRY(row->AddScaled(-2 * curvature * inverse3, bl, d));
    FESA_METRIC_TRY(bz->CopyBlockFrom(row, 0, 0, 1, 14, i, 0, d));
  }
  FESA_METRIC_TRY(bt->CopyBlockFrom(point->Variation(strain), 3, 0, 1,
                                   14, 0, 0, d));
  if (dbt)
    FESA_METRIC_TRY(dbt->CopyBlockFrom(point->Mixed(strain), 3, 0, 1,
                                      14, 0, 0, d));
  return StatusCode::kOk;
}

// Accumulate every lambda-prefactor and B-product derivative separately.
StatusCode AddSpatialMetric(const BeamKinematicWorkspace* point,
                            const Matrix* coefficients, double torsion,
                            double weight, const Matrix* bz, const Matrix* dbz,
                            const Matrix* bt, const Matrix* dbt, Matrix* cbz,
                            Matrix* cdbz, Matrix* damping, Matrix* dc,
                            Diagnostic* d) {
  const double lambda = point->Value(BeamKinematicField::kStretch)->Data()[0];
  FESA_METRIC_TRY(Matrix::Multiply(coefficients, false, bz, false, cbz, d));
  if (damping) {
    FESA_METRIC_TRY(Matrix::MultiplyAdd(bz, true, cbz, false, weight * lambda,
                                        1, damping, d));
    FESA_METRIC_TRY(Matrix::MultiplyAdd(
        bt, true, bt, false, weight * torsion / lambda, 1, damping, d));
  }
  if (!dc) return StatusCode::kOk;
  const double dl = point->Directional(BeamKinematicField::kStretch)->Data()[0];
  FESA_METRIC_TRY(Matrix::Multiply(coefficients, false, dbz, false, cdbz, d));
  FESA_METRIC_TRY(
      Matrix::MultiplyAdd(bz, true, cbz, false, weight * dl, 1, dc, d));
  FESA_METRIC_TRY(
      Matrix::MultiplyAdd(dbz, true, cbz, false, weight * lambda, 1, dc, d));
  FESA_METRIC_TRY(
      Matrix::MultiplyAdd(bz, true, cdbz, false, weight * lambda, 1, dc, d));
  FESA_METRIC_TRY(Matrix::MultiplyAdd(dbt, true, bt, false,
                                      weight * torsion / lambda, 1, dc, d));
  FESA_METRIC_TRY(Matrix::MultiplyAdd(bt, true, dbt, false,
                                      weight * torsion / lambda, 1, dc, d));
  FESA_METRIC_TRY(Matrix::MultiplyAdd(
      bt, true, bt, false, -weight * torsion * dl / (lambda * lambda), 1, dc,
      d));
  return StatusCode::kOk;
}

}  // namespace

StatusCode BeamKineticMetricWorkspace::Create(
    BeamKineticMetricWorkspace** out, Diagnostic* d) {
  if (!out)
    return ReportFailure(StatusCode::kInvalidArgument,
                         DiagnosticStage::kValidation,
                         "BeamKineticMetricWorkspace output", d);
  *out = nullptr;
  auto* created = new (std::nothrow) BeamKineticMetricWorkspace;
  if (!created)
    return ReportFailure(StatusCode::kAllocationFailure,
                         DiagnosticStage::kAllocation,
                         "BeamKineticMetricWorkspace allocation", d);
  StatusCode status = Matrix::Create(1, 14, &created->spin_, d);
  if (status == StatusCode::kOk)
    status = Matrix::Create(1, 14, &created->spin_directional_, d);
  if (status != StatusCode::kOk) {
    delete created;
    return status;
  }
  *out = created;
  return StatusCode::kOk;
}

BeamKineticMetricWorkspace::~BeamKineticMetricWorkspace() {
  delete spin_directional_;
  delete spin_;
}

StatusCode ComputeBeamPointKineticMetric(const BeamKinematicWorkspace* point,
                                         double rho_area, double polar,
                                         BeamKineticMetricWorkspace* scratch,
                                         Matrix* mass, Matrix* dm,
                                         Diagnostic* d,
                                         BeamKineticMetricRequest request) {
  const bool need_mass =
      request == BeamKineticMetricRequest::kMassAndOptionalDirection;
  if ((!need_mass && request != BeamKineticMetricRequest::kDirectionOnly) ||
      !point || !scratch ||
      (need_mass ? !Shape(mass, 14, 14) : mass != nullptr) ||
      (!need_mass && !dm) || (dm && !Shape(dm, 14, 14)) || mass == dm ||
      PointAlias(point, mass) || PointAlias(point, dm) ||
      !std::isfinite(rho_area) || rho_area <= 0 || !std::isfinite(polar) ||
      polar <= 0)
    return Invalid("BeamPointKineticMetric shape/alias/coefficients", d);
  const auto position = BeamKinematicField::kPosition;
  const auto normal1 = BeamKinematicField::kNormal1;
  const auto normal2 = BeamKinematicField::kNormal2;
  const Matrix* jr = point->Variation(position);
  const Matrix* jn = point->Variation(normal1);
  const Matrix* n2 = point->Value(normal2);
  if (!Shape(jr, 3, 14) || !Shape(jn, 3, 14) || !Shape(n2, 3, 1) ||
      (dm && (!Shape(point->Mixed(position), 3, 14) ||
              !Shape(point->Mixed(normal1), 3, 14) ||
              !Shape(point->Directional(normal2), 3, 1))))
    return Invalid("BeamPointKineticMetric invalid current point", d);
  Matrix* spin = scratch->spin_;
  Matrix* dspin = scratch->spin_directional_;
  FESA_METRIC_TRY(Matrix::Multiply(n2, true, jn, false, spin, d));
  if (need_mass) {
    FESA_METRIC_TRY(
        Matrix::MultiplyAdd(jr, true, jr, false, rho_area, 0, mass, d));
    FESA_METRIC_TRY(
        Matrix::MultiplyAdd(spin, true, spin, false, polar, 1, mass, d));
  }
  if (!dm) return StatusCode::kOk;
  FESA_METRIC_TRY(
      Matrix::Multiply(point->Directional(normal2), true, jn, false, dspin, d));
  FESA_METRIC_TRY(Matrix::MultiplyAdd(n2, true, point->Mixed(normal1), false, 1,
                                      1, dspin, d));
  FESA_METRIC_TRY(Matrix::MultiplyAdd(point->Mixed(position), true, jr, false,
                                      rho_area, 0, dm, d));
  FESA_METRIC_TRY(Matrix::MultiplyAdd(jr, true, point->Mixed(position), false,
                                      rho_area, 1, dm, d));
  FESA_METRIC_TRY(
      Matrix::MultiplyAdd(dspin, true, spin, false, polar, 1, dm, d));
  FESA_METRIC_TRY(
      Matrix::MultiplyAdd(spin, true, dspin, false, polar, 1, dm, d));
  return StatusCode::kOk;
}

StatusCode EvaluateBeamDynamicMetrics(const BeamKinematicInput* pose,
                                      const BeamSection* section,
                                      BeamKinematicWorkspace* workspace,
                                      Matrix* mass, Matrix* dm, Matrix* damping,
                                      Matrix* dc,
                                      TrialEvaluationFailure* failure,
                                      Diagnostic* d) {
  if (failure) *failure = TrialEvaluationFailure::kNone;
  if (!pose || !section || !workspace || !failure)
    return Invalid("BeamDynamicMetrics required inputs", d);
  FESA_METRIC_TRY(ValidateOutputs(pose, workspace, mass, dm, damping, dc, d));
  const ElasticMaterial* material = section->GetMaterial();
  if (!material) return Invalid("BeamDynamicMetrics elastic material", d);
  double density = 0;
  FESA_METRIC_TRY(material->Density(&density, d));
  const double rho_area = density * section->Area();
  const double polar = density * (section->I11() + section->I22());
  const double young = material->young();
  const double shear = young / (2 * (1 + material->poisson()));
  const double coefficients[] = {young * section->Area(),
                                 young * section->I11(),
                                 young * section->I22(),
                                 shear * section->TorsionJ(), rho_area, polar};
  for (double value : coefficients)
    if (!std::isfinite(value) || value <= 0)
      return ReportFailure(StatusCode::kInputError,
                           DiagnosticStage::kValidation,
                           "BeamDynamicMetrics section coefficient range", d);
  const double alpha = material->rayleigh_alpha();
  const double beta = material->rayleigh_beta();
  if (!std::isfinite(alpha) || !std::isfinite(beta) || alpha < 0 || beta < 0)
    return ReportFailure(StatusCode::kInputError, DiagnosticStage::kValidation,
                         "BeamDynamicMetrics Rayleigh coefficient range", d);
  const double root = std::sqrt(6. / 5);
  const double nodes[] = {-std::sqrt((3 + 2 * root) / 7),
                          -std::sqrt((3 - 2 * root) / 7),
                          std::sqrt((3 - 2 * root) / 7),
                          std::sqrt((3 + 2 * root) / 7)};
  const double weights[] = {(18 - std::sqrt(30.)) / 36,
                            (18 + std::sqrt(30.)) / 36,
                            (18 + std::sqrt(30.)) / 36,
                            (18 - std::sqrt(30.)) / 36};
  // K1's argument/triad preflight must finish before any output is written.
  FESA_METRIC_TRY(workspace->Evaluate(pose, (nodes[0] + 1) / 2, failure, d,
      dm ? BeamKinematicDerivative::kFull : BeamKinematicDerivative::kFirst));
  BeamKineticMetricWorkspace* raw = nullptr;
  FESA_METRIC_TRY(BeamKineticMetricWorkspace::Create(&raw, d));
  std::unique_ptr<BeamKineticMetricWorkspace> kinetic(raw);
  std::unique_ptr<Matrix> point_mass, point_dm, cb, bz, dbz, bt, dbt;
  std::unique_ptr<Matrix> row, drow, cbz, cdbz;
  FESA_METRIC_TRY(CreateMatrix(14, 14, &point_mass, d));
  if (dm) FESA_METRIC_TRY(CreateMatrix(14, 14, &point_dm, d));
  if (beta != 0) {
    FESA_METRIC_TRY(CreateMatrix(3, 3, &cb, d));
    FESA_METRIC_TRY(CreateMatrix(3, 14, &bz, d));
    FESA_METRIC_TRY(CreateMatrix(3, 14, &cbz, d));
    FESA_METRIC_TRY(CreateMatrix(1, 14, &bt, d));
    FESA_METRIC_TRY(CreateMatrix(1, 14, &row, d));
    if (dm) {
      FESA_METRIC_TRY(CreateMatrix(3, 14, &dbz, d));
      FESA_METRIC_TRY(CreateMatrix(3, 14, &cdbz, d));
      FESA_METRIC_TRY(CreateMatrix(1, 14, &dbt, d));
      FESA_METRIC_TRY(CreateMatrix(1, 14, &drow, d));
    }
    for (int i = 0; i < 3; ++i)
      FESA_METRIC_TRY(cb->Set(i, i, coefficients[i], d));
  }
  // Outputs may contain nonfinite debris from an earlier failed attempt.
  // point_mass is freshly zeroed: overwrite without reading old outputs.
  FESA_METRIC_TRY(mass->CopyFrom(point_mass.get(), d));
  FESA_METRIC_TRY(damping->CopyFrom(point_mass.get(), d));
  if (dm) {
    FESA_METRIC_TRY(dm->CopyFrom(point_mass.get(), d));
    FESA_METRIC_TRY(dc->CopyFrom(point_mass.get(), d));
  }
  for (int i = 0; i < 4; ++i) {
    if (i != 0)
      FESA_METRIC_TRY(workspace->Evaluate(
          pose, (nodes[i] + 1) / 2, failure, d,
          dm ? BeamKinematicDerivative::kFull : BeamKinematicDerivative::kFirst));
    const double weight = pose->length * weights[i] / 2;
    FESA_METRIC_TRY(ComputeBeamPointKineticMetric(
        workspace, rho_area, polar, kinetic.get(), point_mass.get(),
        point_dm.get(), d));
    FESA_METRIC_TRY(mass->AddScaled(weight, point_mass.get(), d));
    if (dm) FESA_METRIC_TRY(dm->AddScaled(weight, point_dm.get(), d));
    if (beta == 0) continue;
    FESA_METRIC_TRY(SpatialRows(workspace, bz.get(), dbz.get(), bt.get(),
                                dbt.get(), row.get(), drow.get(), d));
    FESA_METRIC_TRY(AddSpatialMetric(
        workspace, cb.get(), coefficients[3], beta * weight, bz.get(),
        dbz.get(), bt.get(), dbt.get(), cbz.get(), cdbz.get(), damping, dc, d));
  }
  FESA_METRIC_TRY(damping->AddScaled(alpha, mass, d));
  if (dc) FESA_METRIC_TRY(dc->AddScaled(alpha, dm, d));
  return StatusCode::kOk;
}

struct BeamRayleighWorkspace::Storage {
  enum Buffer {
    kZero,
    kMass,
    kDm,
    kPointMass,
    kPointDm,
    kCoefficients,
    kBz,
    kDbz,
    kBt,
    kDbt,
    kRow,
    kDrow,
    kCbz,
    kCdbz,
    kBufferCount
  };
  BeamKinematicWorkspace* points[4] = {};
  BeamKineticMetricWorkspace* kinetic = nullptr;
  Matrix* buffer[kBufferCount] = {};
  double xi[4] = {}, weights[4] = {};
  double coefficients[6] = {};  // EA,EI11,EI22,GJ,rho*A,rho*(I11+I22).
  double alpha = 0, beta = 0;
  bool active = false;
  Storage() = default;
  Storage(const Storage&) = delete;
  Storage& operator=(const Storage&) = delete;
  Storage(Storage&&) = delete;
  Storage& operator=(Storage&&) = delete;
  ~Storage() {
    for (auto* point : points) delete point;
    delete kinetic;
    for (auto* matrix : buffer) delete matrix;
  }
  StatusCode Initialize(Diagnostic* d) {
    for (int i = 0; i < 4; ++i)
      FESA_METRIC_TRY(BeamKinematicWorkspace::Create(&points[i], d));
    FESA_METRIC_TRY(BeamKineticMetricWorkspace::Create(&kinetic, d));
    const int rows[] = {14, 14, 14, 14, 14, 3, 3, 3, 1, 1, 1, 1, 3, 3};
    const int columns[] = {14, 14, 14, 14, 14, 3,  14,
                           14, 14, 14, 14, 14, 14, 14};
    for (int i = 0; i < kBufferCount; ++i)
      FESA_METRIC_TRY(Matrix::Create(rows[i], columns[i], &buffer[i], d));
    const double root = std::sqrt(6. / 5);
    const double nodes[] = {
        -std::sqrt((3 + 2 * root) / 7), -std::sqrt((3 - 2 * root) / 7),
        std::sqrt((3 - 2 * root) / 7), std::sqrt((3 + 2 * root) / 7)};
    for (int i = 0; i < 4; ++i) xi[i] = (nodes[i] + 1) / 2;
    weights[0] = weights[3] = (18 - std::sqrt(30.)) / 36;
    weights[1] = weights[2] = (18 + std::sqrt(30.)) / 36;
    return StatusCode::kOk;
  }
  void Invalidate() {
    active = false;
    for (auto* point : points)
      if (point) point->Invalidate();
  }
  StatusCode Material(const BeamSection* section, Diagnostic* d) {
    const auto* material = section->GetMaterial();
    if (!material) return Invalid("BeamDynamicMetrics elastic material", d);
    double density = 0;
    FESA_METRIC_TRY(material->Density(&density, d));
    const double young = material->young();
    const double shear = young / (2 * (1 + material->poisson()));
    coefficients[0] = young * section->Area();
    coefficients[1] = young * section->I11();
    coefficients[2] = young * section->I22();
    coefficients[3] = shear * section->TorsionJ();
    coefficients[4] = density * section->Area();
    coefficients[5] = density * (section->I11() + section->I22());
    for (double value : coefficients)
      if (!std::isfinite(value) || value <= 0)
        return ReportFailure(StatusCode::kInputError,
                             DiagnosticStage::kValidation,
                             "BeamDynamicMetrics section coefficient range", d);
    alpha = material->rayleigh_alpha();
    beta = material->rayleigh_beta();
    if (!std::isfinite(alpha) || !std::isfinite(beta) || alpha < 0 || beta < 0)
      return ReportFailure(StatusCode::kInputError,
                           DiagnosticStage::kValidation,
                           "BeamDynamicMetrics Rayleigh coefficient range", d);
    if (beta != 0)
      for (int i = 0; i < 3; ++i)
        FESA_METRIC_TRY(buffer[kCoefficients]->Set(i, i, coefficients[i], d));
    return StatusCode::kOk;
  }
  StatusCode Base(const BeamKinematicInput* pose, const BeamSection* section,
                  const Matrix* current_mass, Matrix* damping,
                  TrialEvaluationFailure* failure, Diagnostic* d) {
    if (failure) *failure = TrialEvaluationFailure::kNone;
    if (!pose || !section || !failure || !Shape(damping, 14, 14) ||
        (current_mass &&
         (!Shape(current_mass, 14, 14) || current_mass == damping)))
      return Invalid("BeamRayleighWorkspace base inputs/shape/alias", d);
    FESA_METRIC_TRY(Material(section, d));
    // K1 preflight precedes caller output writes, as in the complete metric.
    FESA_METRIC_TRY(points[0]->Evaluate(pose, xi[0], failure, d,
                                        BeamKinematicDerivative::kFirst));
    FESA_METRIC_TRY(damping->CopyFrom(buffer[kZero], d));
    const bool integrate_mass = alpha != 0 && !current_mass;
    if (integrate_mass)
      FESA_METRIC_TRY(buffer[kMass]->CopyFrom(buffer[kZero], d));
    for (int i = 0; i < 4; ++i) {
      if (i != 0)
        FESA_METRIC_TRY(points[i]->Evaluate(pose, xi[i], failure, d,
                                            BeamKinematicDerivative::kFirst));
      const double weight = pose->length * weights[i] / 2;
      if (integrate_mass) {
        FESA_METRIC_TRY(ComputeBeamPointKineticMetric(
            points[i], coefficients[4], coefficients[5], kinetic,
            buffer[kPointMass], nullptr, d));
        FESA_METRIC_TRY(
            buffer[kMass]->AddScaled(weight, buffer[kPointMass], d));
      }
      if (beta == 0) continue;
      FESA_METRIC_TRY(SpatialRows(points[i], buffer[kBz], nullptr, buffer[kBt],
                                  nullptr, buffer[kRow], nullptr, d));
      FESA_METRIC_TRY(AddSpatialMetric(
          points[i], buffer[kCoefficients], coefficients[3], beta * weight,
          buffer[kBz], nullptr, buffer[kBt], nullptr, buffer[kCbz], nullptr,
          damping, nullptr, d));
    }
    if (alpha != 0)
      FESA_METRIC_TRY(damping->AddScaled(
          alpha, current_mass ? current_mass : buffer[kMass], d));
    active = true;
    return StatusCode::kOk;
  }
  StatusCode Direction(const BeamKinematicInput* pose, Matrix* dc,
                       TrialEvaluationFailure* failure, Diagnostic* d) {
    if (failure) *failure = TrialEvaluationFailure::kNone;
    if (!active || !pose || !failure || !Shape(dc, 14, 14))
      return Invalid("BeamRayleighWorkspace active candidate/direction output",
                     d);
    const auto derivative = alpha != 0 || beta != 0
                                ? BeamKinematicDerivative::kFull
                                : BeamKinematicDerivative::kFirst;
    FESA_METRIC_TRY(points[0]->Evaluate(pose, xi[0], failure, d, derivative));
    FESA_METRIC_TRY(dc->CopyFrom(buffer[kZero], d));
    if (alpha != 0) FESA_METRIC_TRY(buffer[kDm]->CopyFrom(buffer[kZero], d));
    for (int i = 0; i < 4; ++i) {
      if (i != 0)
        FESA_METRIC_TRY(
            points[i]->Evaluate(pose, xi[i], failure, d, derivative));
      const double weight = pose->length * weights[i] / 2;
      if (alpha != 0) {
        FESA_METRIC_TRY(ComputeBeamPointKineticMetric(
            points[i], coefficients[4], coefficients[5], kinetic, nullptr,
            buffer[kPointDm], d, BeamKineticMetricRequest::kDirectionOnly));
        FESA_METRIC_TRY(buffer[kDm]->AddScaled(weight, buffer[kPointDm], d));
      }
      if (beta == 0) continue;
      FESA_METRIC_TRY(SpatialRows(points[i], buffer[kBz], buffer[kDbz],
                                  buffer[kBt], buffer[kDbt], buffer[kRow],
                                  buffer[kDrow], d));
      FESA_METRIC_TRY(AddSpatialMetric(
          points[i], buffer[kCoefficients], coefficients[3], beta * weight,
          buffer[kBz], buffer[kDbz], buffer[kBt], buffer[kDbt], buffer[kCbz],
          buffer[kCdbz], nullptr, dc, d));
    }
    if (alpha != 0) FESA_METRIC_TRY(dc->AddScaled(alpha, buffer[kDm], d));
    return StatusCode::kOk;
  }
};
StatusCode BeamRayleighWorkspace::Create(BeamRayleighWorkspace** out,
                                         Diagnostic* d) {
  if (!out) return Invalid("BeamRayleighWorkspace output", d);
  *out = nullptr;
  std::unique_ptr<BeamRayleighWorkspace> work(new (std::nothrow)
                                                  BeamRayleighWorkspace);
  if (!work)
    return ReportFailure(StatusCode::kAllocationFailure,
                         DiagnosticStage::kAllocation,
                         "BeamRayleighWorkspace allocation", d);
  work->storage_ = new (std::nothrow) Storage;
  if (!work->storage_)
    return ReportFailure(StatusCode::kAllocationFailure,
                         DiagnosticStage::kAllocation,
                         "BeamRayleighWorkspace storage allocation", d);
  const auto status = work->storage_->Initialize(d);
  if (status != StatusCode::kOk) return status;
  *out = work.release();
  return StatusCode::kOk;
}
BeamRayleighWorkspace::~BeamRayleighWorkspace() { delete storage_; }
StatusCode BeamRayleighWorkspace::EvaluateBase(const BeamKinematicInput* pose,
                                               const BeamSection* section,
                                               const Matrix* current_mass,
                                               Matrix* damping,
                                               TrialEvaluationFailure* failure,
                                               Diagnostic* d) {
  storage_->Invalidate();
  const auto status =
      storage_->Base(pose, section, current_mass, damping, failure, d);
  if (status != StatusCode::kOk) storage_->Invalidate();
  return status;
}
StatusCode BeamRayleighWorkspace::EvaluateDirection(
    const BeamKinematicInput* pose, Matrix* dc, TrialEvaluationFailure* failure,
    Diagnostic* d) {
  const auto status = storage_->Direction(pose, dc, failure, d);
  if (status != StatusCode::kOk) storage_->Invalidate();
  return status;
}

#undef FESA_METRIC_TRY

}  // namespace fesa
