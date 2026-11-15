#include "elements/beam3d_inertia.h"

#include <cmath>

#include "elements/beam3d_dynamic_metrics.h"
#include "elements/beam3d_kinematic_graph.h"

namespace fesa {
namespace {

#define FESA_INERTIA_TRY(expression)              \
  do {                                            \
    const StatusCode status = (expression);       \
    if (status != StatusCode::kOk) return status; \
  } while (false)

constexpr auto kPosition = BeamKinematicField::kPosition;
constexpr auto kTangent = BeamKinematicField::kTangent;
constexpr auto kNormal1 = BeamKinematicField::kNormal1;
constexpr auto kNormal2 = BeamKinematicField::kNormal2;
constexpr BeamKinematicField kFields[] = {kPosition,
                                          kTangent,
                                          kNormal1,
                                          kNormal2,
                                          BeamKinematicField::kStrain,
                                          BeamKinematicField::kStretch};

// Fixed buffer slots. Every allocation is completed before the Gauss loop.
enum Slot {
  kJr,
  kJt,
  kJn1,
  kT,
  kN2,
  kSpin,
  kSkew,
  kVirtualSpin,
  kRAcceleration,
  kMomentumRate,
  kScalar,
  kSpinDerivative,
  kSpinRateDerivative,
  kRAccelerationDerivative,
  kMomentumRateDerivative,
  kDSpin,
  kDSkew,
  kDVirtualSpin,
  kColumn,
  kSelector,
  kPointG,
  kPointC,
  kPointM,
  kForce,
  kG,
  kC,
  kM,
  kSlotCount
};

struct Scratch {
  Matrix* matrix[kSlotCount] = {};
  Vector* direction = nullptr;
  BeamKineticMetricWorkspace* kinetic = nullptr;
  Scratch() = default;
  Scratch(const Scratch&) = delete;
  Scratch& operator=(const Scratch&) = delete;
  Scratch(Scratch&&) = delete;
  Scratch& operator=(Scratch&&) = delete;
  ~Scratch() {
    delete kinetic;
    delete direction;
    for (Matrix* value : matrix) delete value;
  }
  StatusCode Initialize(bool full, Diagnostic* d) {
    constexpr int shapes[kSlotCount][2] = {
        {3, 14},  {3, 14},  {3, 14}, {3, 1},   {3, 1},   {1, 14}, {3, 3},
        {3, 14},  {3, 1},   {3, 1},  {1, 1},   {1, 14},  {1, 14}, {3, 14},
        {3, 14},  {1, 14},  {3, 3},  {3, 14},  {14, 1},  {14, 1}, {14, 14},
        {14, 14}, {14, 14}, {14, 1}, {14, 14}, {14, 14}, {14, 14}};
    for (int i = 0; i < kSlotCount; ++i) {
      if (!full && i >= kSpinDerivative && i != kForce) continue;
      FESA_INERTIA_TRY(
          Matrix::Create(shapes[i][0], shapes[i][1], &matrix[i], d));
    }
    FESA_INERTIA_TRY(Vector::Create(14, &direction, d));
    if (full) FESA_INERTIA_TRY(BeamKineticMetricWorkspace::Create(&kinetic, d));
    return StatusCode::kOk;
  }
};

StatusCode Invalid(const char* operation, Diagnostic* d) {
  return ReportFailure(StatusCode::kInvalidArgument,
                       DiagnosticStage::kValidation, operation, d);
}

bool Shape(const Matrix* matrix) {
  return matrix && matrix->Rows() == 14 && matrix->Columns() == 14;
}

bool MotionValid(const Vector* vector) {
  return vector && vector->Size() == 14 &&
         beam3d_graph::Finite(vector->Data(), 14);
}

bool GetterAlias(const BeamKinematicWorkspace* current,
                 const BeamTimeKinematicWorkspace* time, const Matrix* out) {
  if (!out) return false;
  for (auto field : kFields) {
    if (out == current->Value(field) || out == current->Variation(field) ||
        out == current->Directional(field) || out == current->Mixed(field))
      return true;
    for (int order = 0; order < 3; ++order)
      if (out == time->Coefficient(field, order) ||
          out == time->Derivative(field, order))
        return true;
  }
  return false;
}

StatusCode Validate(const BeamKinematicInput* pose, const Vector* velocity,
                    const Vector* acceleration,
                    const BeamKinematicWorkspace* current,
                    const BeamTimeKinematicWorkspace* time,
                    DynamicEvaluationRequest request, const Vector* force,
                    const Matrix* g, const Matrix* c, const Matrix* m,
                    Diagnostic* d) {
  if (!beam3d_graph::Valid(pose, .5) || !MotionValid(velocity) ||
      !MotionValid(acceleration) || !current || !time || !force ||
      force->Size() != 14 ||
      (request != DynamicEvaluationRequest::kForces &&
       request != DynamicEvaluationRequest::kForcesAndLinearization))
    return Invalid("BeamInertia required inputs/shape/finite", d);
  const bool full =
      request == DynamicEvaluationRequest::kForcesAndLinearization;
  if ((full && (!Shape(g) || !Shape(c) || !Shape(m))) ||
      (!full && (g || c || m)) || force == velocity || force == acceleration ||
      force == pose->internal_position || force == pose->direction)
    return Invalid("BeamInertia output request/shape/alias", d);
  const Matrix* outputs[] = {g, c, m};
  for (int i = 0; i < 3; ++i) {
    if (!outputs[i]) continue;
    for (int j = 0; j < i; ++j)
      if (outputs[i] == outputs[j])
        return Invalid("BeamInertia output alias", d);
    if (outputs[i] == pose->positions || outputs[i] == pose->tangents ||
        outputs[i] == pose->normals || GetterAlias(current, time, outputs[i]))
      return Invalid("BeamInertia borrowed output alias", d);
  }
  return StatusCode::kOk;
}

// Scalar skew coefficients define cross(v,.); products remain common math.
StatusCode Skew(const Matrix* vector, Matrix* out, Diagnostic* d) {
  const double* v = vector->Data();
  FESA_INERTIA_TRY(out->Set(0, 1, -v[2], d));
  FESA_INERTIA_TRY(out->Set(0, 2, v[1], d));
  FESA_INERTIA_TRY(out->Set(1, 0, v[2], d));
  FESA_INERTIA_TRY(out->Set(1, 2, -v[0], d));
  FESA_INERTIA_TRY(out->Set(2, 0, -v[1], d));
  FESA_INERTIA_TRY(out->Set(2, 1, v[0], d));
  return StatusCode::kOk;
}

// Snapshot all current quantities before any direction evaluation expires K1.
StatusCode CurrentSpin(const BeamKinematicWorkspace* point, Scratch* work,
                       Diagnostic* d) {
  Matrix** b = work->matrix;
  FESA_INERTIA_TRY(b[kJr]->CopyFrom(point->Variation(kPosition), d));
  FESA_INERTIA_TRY(b[kJt]->CopyFrom(point->Variation(kTangent), d));
  FESA_INERTIA_TRY(b[kJn1]->CopyFrom(point->Variation(kNormal1), d));
  FESA_INERTIA_TRY(b[kT]->CopyFrom(point->Value(kTangent), d));
  FESA_INERTIA_TRY(b[kN2]->CopyFrom(point->Value(kNormal2), d));
  FESA_INERTIA_TRY(Matrix::Multiply(b[kN2], true, b[kJn1], false, b[kSpin], d));
  FESA_INERTIA_TRY(Skew(b[kT], b[kSkew], d));
  FESA_INERTIA_TRY(
      Matrix::Multiply(b[kSkew], false, b[kJt], false, b[kVirtualSpin], d));
  FESA_INERTIA_TRY(Matrix::MultiplyAdd(b[kT], false, b[kSpin], false, 1, 1,
                                       b[kVirtualSpin], d));
  return StatusCode::kOk;
}

// n_i coefficients are Taylor coefficients; coefficient2=.5*second derivative.
// The full product rule keeps both transverse momentum terms.
StatusCode TimeMotion(const BeamTimeKinematicWorkspace* time, double polar,
                      bool differentiate, Scratch* work, Diagnostic* d) {
  Matrix** b = work->matrix;
  const Matrix* t0 = time->Coefficient(kTangent, 0);
  const Matrix* t1 = time->Coefficient(kTangent, 1);
  const Matrix* n20 = time->Coefficient(kNormal2, 0);
  const Matrix* n21 = time->Coefficient(kNormal2, 1);
  const Matrix* n11 = time->Coefficient(kNormal1, 1);
  const Matrix* n12 = time->Coefficient(kNormal1, 2);
  FESA_INERTIA_TRY(
      b[kRAcceleration]->CopyFrom(time->Coefficient(kPosition, 2), d));
  FESA_INERTIA_TRY(b[kRAcceleration]->Scale(2, d));
  FESA_INERTIA_TRY(Matrix::Multiply(n20, true, n11, false, b[kScalar], d));
  const double spin = b[kScalar]->Data()[0];
  FESA_INERTIA_TRY(Matrix::Multiply(n21, true, n11, false, b[kScalar], d));
  FESA_INERTIA_TRY(
      Matrix::MultiplyAdd(n20, true, n12, false, 2, 1, b[kScalar], d));
  const double spin_rate = b[kScalar]->Data()[0];
  FESA_INERTIA_TRY(b[kMomentumRate]->CopyFrom(t0, d));
  FESA_INERTIA_TRY(b[kMomentumRate]->Scale(polar * spin_rate, d));
  FESA_INERTIA_TRY(b[kMomentumRate]->AddScaled(polar * spin, t1, d));
  if (!differentiate) return StatusCode::kOk;
  const Matrix* dn20 = time->Derivative(kNormal2, 0);
  const Matrix* dn21 = time->Derivative(kNormal2, 1);
  const Matrix* dn11 = time->Derivative(kNormal1, 1);
  const Matrix* dn12 = time->Derivative(kNormal1, 2);
  FESA_INERTIA_TRY(
      b[kRAccelerationDerivative]->CopyFrom(time->Derivative(kPosition, 2), d));
  FESA_INERTIA_TRY(b[kRAccelerationDerivative]->Scale(2, d));
  FESA_INERTIA_TRY(
      Matrix::Multiply(n11, true, dn20, false, b[kSpinDerivative], d));
  FESA_INERTIA_TRY(
      Matrix::MultiplyAdd(n20, true, dn11, false, 1, 1, b[kSpinDerivative], d));
  FESA_INERTIA_TRY(
      Matrix::Multiply(n11, true, dn21, false, b[kSpinRateDerivative], d));
  FESA_INERTIA_TRY(Matrix::MultiplyAdd(n21, true, dn11, false, 1, 1,
                                       b[kSpinRateDerivative], d));
  FESA_INERTIA_TRY(Matrix::MultiplyAdd(n12, true, dn20, false, 2, 1,
                                       b[kSpinRateDerivative], d));
  FESA_INERTIA_TRY(Matrix::MultiplyAdd(n20, true, dn12, false, 2, 1,
                                       b[kSpinRateDerivative], d));
  FESA_INERTIA_TRY(Matrix::MultiplyAdd(t0, false, b[kSpinRateDerivative], false,
                                       polar, 0, b[kMomentumRateDerivative],
                                       d));
  FESA_INERTIA_TRY(b[kMomentumRateDerivative]->AddScaled(
      polar * spin_rate, time->Derivative(kTangent, 0), d));
  FESA_INERTIA_TRY(Matrix::MultiplyAdd(t1, false, b[kSpinDerivative], false,
                                       polar, 1, b[kMomentumRateDerivative],
                                       d));
  FESA_INERTIA_TRY(b[kMomentumRateDerivative]->AddScaled(
      polar * spin, time->Derivative(kTangent, 1), d));
  return StatusCode::kOk;
}

StatusCode Partial(double rho_area, Scratch* work, Matrix* out, Diagnostic* d) {
  Matrix** b = work->matrix;
  FESA_INERTIA_TRY(Matrix::MultiplyAdd(
      b[kJr], true, b[kRAccelerationDerivative], false, rho_area, 0, out, d));
  FESA_INERTIA_TRY(Matrix::MultiplyAdd(
      b[kVirtualSpin], true, b[kMomentumRateDerivative], false, 1, 1, out, d));
  return StatusCode::kOk;
}

// D Btheta = [Dt]x Jt + [t]x DJt + Dt*w + t*Dw. Mixed is the
// ordered derivative of the fixed spatial virtual basis, with no half skew.
StatusCode GeometryColumn(const BeamKinematicWorkspace* point, double rho_area,
                          Scratch* work, Diagnostic* d) {
  Matrix** b = work->matrix;
  const Matrix* dt = point->Directional(kTangent);
  FESA_INERTIA_TRY(Matrix::Multiply(point->Directional(kNormal2), true, b[kJn1],
                                    false, b[kDSpin], d));
  FESA_INERTIA_TRY(Matrix::MultiplyAdd(b[kN2], true, point->Mixed(kNormal1),
                                       false, 1, 1, b[kDSpin], d));
  FESA_INERTIA_TRY(Skew(dt, b[kDSkew], d));
  FESA_INERTIA_TRY(
      Matrix::Multiply(b[kDSkew], false, b[kJt], false, b[kDVirtualSpin], d));
  FESA_INERTIA_TRY(Matrix::MultiplyAdd(b[kSkew], false, point->Mixed(kTangent),
                                       false, 1, 1, b[kDVirtualSpin], d));
  FESA_INERTIA_TRY(Matrix::MultiplyAdd(dt, false, b[kSpin], false, 1, 1,
                                       b[kDVirtualSpin], d));
  FESA_INERTIA_TRY(Matrix::MultiplyAdd(b[kT], false, b[kDSpin], false, 1, 1,
                                       b[kDVirtualSpin], d));
  FESA_INERTIA_TRY(Matrix::MultiplyAdd(point->Mixed(kPosition), true,
                                       b[kRAcceleration], false, rho_area, 0,
                                       b[kColumn], d));
  FESA_INERTIA_TRY(Matrix::MultiplyAdd(b[kDVirtualSpin], true, b[kMomentumRate],
                                       false, 1, 1, b[kColumn], d));
  FESA_INERTIA_TRY(Matrix::MultiplyAdd(b[kColumn], false, b[kSelector], true, 1,
                                       1, b[kPointG], d));
  return StatusCode::kOk;
}

}  // namespace

StatusCode EvaluateBeamInertia(
    const BeamKinematicInput* pose, const BeamSection* section,
    const Vector* velocity14, const Vector* acceleration14,
    BeamKinematicWorkspace* current, BeamTimeKinematicWorkspace* time,
    DynamicEvaluationRequest request, Vector* force14, Matrix* gi14,
    Matrix* ci14, Matrix* mi14, TrialEvaluationFailure* failure,
    Diagnostic* d) {
  if (failure) *failure = TrialEvaluationFailure::kNone;
  if (!section || !failure || !section->GetMaterial())
    return Invalid("BeamInertia required section/failure/material", d);
  FESA_INERTIA_TRY(Validate(pose, velocity14, acceleration14, current, time,
                            request, force14, gi14, ci14, mi14, d));
  double density = 0;
  FESA_INERTIA_TRY(section->GetMaterial()->Density(&density, d));
  const double rho_area = density * section->Area();
  const double polar = density * (section->I11() + section->I22());
  if (!std::isfinite(rho_area) || rho_area <= 0 || !std::isfinite(polar) ||
      polar <= 0)
    return ReportFailure(StatusCode::kInputError, DiagnosticStage::kValidation,
                         "BeamInertia density section coefficient range", d);
  const bool full =
      request == DynamicEvaluationRequest::kForcesAndLinearization;
  Scratch scratch;
  FESA_INERTIA_TRY(scratch.Initialize(full, d));
  Matrix** b = scratch.matrix;
  BeamKinematicInput local;
  local.length = pose->length;
  local.positions = pose->positions;
  local.relative_positions = pose->relative_positions;
  local.tangents = pose->tangents;
  local.normals = pose->normals;
  local.internal_position = pose->internal_position;
  local.direction = scratch.direction;
  const double root = std::sqrt(6. / 5);
  const double nodes[] = {
      -std::sqrt((3 + 2 * root) / 7), -std::sqrt((3 - 2 * root) / 7),
      std::sqrt((3 - 2 * root) / 7), std::sqrt((3 + 2 * root) / 7)};
  const double weights[] = {
      (18 - std::sqrt(30.)) / 36, (18 + std::sqrt(30.)) / 36,
      (18 + std::sqrt(30.)) / 36, (18 - std::sqrt(30.)) / 36};
  for (int point = 0; point < 4; ++point) {
    const double xi = (1 + nodes[point]) / 2;
    const double weight = local.length * weights[point] / 2;
    FESA_INERTIA_TRY(current->Evaluate(&local, xi, failure, d,
                                       BeamKinematicDerivative::kFirst));
    FESA_INERTIA_TRY(CurrentSpin(current, &scratch, d));
    FESA_INERTIA_TRY(time->Evaluate(
        &local, velocity14, acceleration14,
        full ? BeamTimeDerivative::kDisplacement : BeamTimeDerivative::kNone,
        xi, failure, d));
    FESA_INERTIA_TRY(TimeMotion(time, polar, full, &scratch, d));
    FESA_INERTIA_TRY(Matrix::MultiplyAdd(b[kJr], true, b[kRAcceleration], false,
                                         weight * rho_area, point == 0 ? 0 : 1,
                                         b[kForce], d));
    FESA_INERTIA_TRY(Matrix::MultiplyAdd(b[kVirtualSpin], true,
                                         b[kMomentumRate], false, weight, 1,
                                         b[kForce], d));
    if (!full) continue;
    FESA_INERTIA_TRY(Partial(rho_area, &scratch, b[kPointG], d));
    for (int column = 0; column < 14; ++column) {
      FESA_INERTIA_TRY(scratch.direction->Set(column, 1, d));
      FESA_INERTIA_TRY(b[kSelector]->Set(column, 0, 1, d));
      FESA_INERTIA_TRY(current->Evaluate(&local, xi, failure, d));
      FESA_INERTIA_TRY(GeometryColumn(current, rho_area, &scratch, d));
      FESA_INERTIA_TRY(scratch.direction->Set(column, 0, d));
      FESA_INERTIA_TRY(b[kSelector]->Set(column, 0, 0, d));
    }
    FESA_INERTIA_TRY(b[kG]->AddScaled(weight, b[kPointG], d));
    FESA_INERTIA_TRY(time->Evaluate(&local, velocity14, acceleration14,
                                    BeamTimeDerivative::kVelocity, xi, failure,
                                    d));
    FESA_INERTIA_TRY(TimeMotion(time, polar, true, &scratch, d));
    FESA_INERTIA_TRY(Partial(rho_area, &scratch, b[kPointC], d));
    FESA_INERTIA_TRY(b[kC]->AddScaled(weight, b[kPointC], d));
    FESA_INERTIA_TRY(ComputeBeamPointKineticMetric(
        current, rho_area, polar, scratch.kinetic, b[kPointM], nullptr, d));
    FESA_INERTIA_TRY(b[kM]->AddScaled(weight, b[kPointM], d));
  }
  const int indices[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
  FESA_INERTIA_TRY(b[kForce]->GatherTo(indices, 14, force14, 0, d));
  if (full) {
    FESA_INERTIA_TRY(gi14->CopyFrom(b[kG], d));
    FESA_INERTIA_TRY(ci14->CopyFrom(b[kC], d));
    FESA_INERTIA_TRY(mi14->CopyFrom(b[kM], d));
  }
  return StatusCode::kOk;
}

#undef FESA_INERTIA_TRY

}  // namespace fesa
