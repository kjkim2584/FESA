#include "fesa/elements/mitc4_results.h"

#include <cmath>
#include <memory>
#include <new>

#include "elements/mitc4_geometry.h"
#include "elements/mitc4_initial.h"
#include "elements/mitc4_linear_ops.h"
#include "elements/mitc4_nonlinear_recovery.h"
#include "elements/mitc4_trial.h"
#include "model/model_storage.h"

namespace fesa {
namespace {

StatusCode CreateMatrix(int rows, int columns, std::unique_ptr<Matrix>* owner,
                        Diagnostic* d) {
  Matrix* raw = nullptr;
  const auto status = Matrix::Create(rows, columns, &raw, d);
  owner->reset(raw);
  return status;
}

StatusCode ValidateInput(const Node* const* nodes, int count,
                         const ShellSection* section, const Vector* u,
                         const Vector* force, double drilling, Diagnostic* d) {
  if (!nodes || count != 4 || !section || !u || u == force || u->Size() != 24 ||
      !std::isfinite(drilling) || !(drilling > 0.0)) {
    return model_storage::Invalid("MITC4 linear results input/shape/alias", d);
  }
  for (int i = 0; i < 24; ++i) {
    if (!std::isfinite(u->Data()[i]))
      return model_storage::Invalid("MITC4 linear results nonfinite input", d);
  }
  for (int i = 0; i < 4; ++i) {
    if (!nodes[i] || nodes[i]->dof_count() != 6)
      return model_storage::Invalid("MITC4 linear results node components", d);
    for (int j = 0; j < i; ++j) {
      if (nodes[i] == nodes[j] || nodes[i]->Id() == nodes[j]->Id())
        return model_storage::Invalid("MITC4 linear results duplicate node", d);
    }
  }
  return StatusCode::kOk;
}

StatusCode FillTransform(const Matrix* frame, Matrix* transform,
                         Diagnostic* d) {
  for (int block = 0; block < 8; ++block) {
    const auto status =
        transform->CopyBlockFrom(frame, 0, 0, 3, 3, 3 * block, 3 * block, d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}

StatusCode FillStressSelector(Matrix* selector, Diagnostic* d) {
  // The constitutive order is xx,yy,xy,yz,xz. Official S is
  // 11,22,33,12,13,23 with tensor shear; the 33 row remains zero.
  constexpr int kOutput[5] = {0, 1, 3, 5, 4};
  for (int component = 0; component < 5; ++component) {
    const auto status = selector->Set(kOutput[component], component, 1.0, d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}

StatusCode FillDirectStrain(const Matrix* b, const Matrix* gradients, double r,
                            double s, Matrix* direct, Diagnostic* d) {
  // The three in-plane components equal the linear MITC B rows. LE's
  // transverse engineering shears use the DI gradient, without tying.
  auto status = direct->CopyBlockFrom(b, 0, 0, 2, 24, 0, 0, d);
  if (status != StatusCode::kOk) return status;
  status = direct->CopyBlockFrom(b, 2, 0, 1, 24, 3, 0, d);
  if (status != StatusCode::kOk) return status;
  constexpr double kR[4] = {-1.0, 1.0, 1.0, -1.0};
  constexpr double kS[4] = {-1.0, -1.0, 1.0, 1.0};
  for (int node = 0; node < 4; ++node) {
    const double shape = (1.0 + kR[node] * r) * (1.0 + kS[node] * s) / 4.0;
    status = direct->Set(4, 6 * node + 2, gradients->Data()[node], d);
    if (status != StatusCode::kOk) return status;
    status = direct->Set(4, 6 * node + 4, shape, d);
    if (status != StatusCode::kOk) return status;
    status = direct->Set(5, 6 * node + 2, gradients->Data()[4 + node], d);
    if (status != StatusCode::kOk) return status;
    status = direct->Set(5, 6 * node + 3, -shape, d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}

StatusCode FillExtrapolation(Matrix* out, Diagnostic* d) {
  constexpr double kR[4] = {-1.0, 1.0, 1.0, -1.0};
  constexpr double kS[4] = {-1.0, -1.0, 1.0, 1.0};
  const double root_three = std::sqrt(3.0);
  for (int corner = 0; corner < 4; ++corner) {
    for (int plane = 0; plane < 4; ++plane) {
      // Bilinear interpolation at the four Gauss points, evaluated at a
      // connectivity corner. This is not a shared-node averaging operator.
      const double weight = (1.0 + kR[plane] * root_three * kR[corner]) *
                            (1.0 + kS[plane] * root_three * kS[corner]) / 4.0;
      const auto status = out->Set(corner, plane, weight, d);
      if (status != StatusCode::kOk) return status;
    }
  }
  return StatusCode::kOk;
}

StatusCode PublishSamples(const Matrix* samples, const Matrix* extrapolation,
                          Matrix* layer, Matrix* extrapolated, Matrix* ip,
                          Matrix* surface, Matrix* nodal, Diagnostic* d) {
  for (int plane = 0; plane < 4; ++plane) {
    for (int thickness = 0; thickness < 2; ++thickness) {
      const auto status = ip->CopyBlockFrom(samples, 4 * thickness + plane, 0,
                                            1, 6, 2 * plane + thickness, 0, d);
      if (status != StatusCode::kOk) return status;
    }
  }
  auto status = surface->CopyBlockFrom(samples, 8, 0, 8, 6, 0, 0, d);
  if (status != StatusCode::kOk) return status;
  for (int thickness = 0; thickness < 4; ++thickness) {
    status = layer->CopyBlockFrom(samples, 4 * thickness, 0, 4, 6, 0, 0, d);
    if (status != StatusCode::kOk) return status;
    status =
        Matrix::Multiply(extrapolation, false, layer, false, extrapolated, d);
    if (status != StatusCode::kOk) return status;
    status =
        nodal->CopyBlockFrom(extrapolated, 0, 0, 4, 6, 4 * thickness, 0, d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}

}  // namespace

StatusCode Mitc4Results::Create(Mitc4Results** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("Mitc4Results::Create output", d);
  *out = nullptr;
  std::unique_ptr<Mitc4Results> owner(new (std::nothrow) Mitc4Results());
  if (!owner) return model_storage::AllocationFailure("Mitc4Results", d);
  auto status = Vector::Create(24, &owner->internal_force_, d);
  if (status != StatusCode::kOk) return status;
  Matrix** outputs[7] = {&owner->axes_,           &owner->ip_stress_,
                         &owner->ip_strain_,      &owner->surface_stress_,
                         &owner->surface_strain_, &owner->nodal_stress_,
                         &owner->nodal_strain_};
  const int rows[7] = {3, 8, 8, 8, 8, 16, 16};
  for (int i = 0; i < 7; ++i) {
    status = Matrix::Create(rows[i], i == 0 ? 3 : 6, outputs[i], d);
    if (status != StatusCode::kOk) return status;
  }
  *out = owner.release();
  return StatusCode::kOk;
}

Mitc4Results::~Mitc4Results() {
  delete nodal_strain_;
  delete nodal_stress_;
  delete surface_strain_;
  delete surface_stress_;
  delete ip_strain_;
  delete ip_stress_;
  delete axes_;
  delete internal_force_;
}

StatusCode Mitc4Results::RecoverLinear(const Node* const* nodes, int count,
                                       const ShellSection* section,
                                       const Vector* u, double drilling,
                                       Diagnostic* d) {
  auto status =
      ValidateInput(nodes, count, section, u, internal_force_, drilling, d);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<Matrix> matrices[15];
  const int rows[15] = {4, 3, 24, 4, 5, 2, 2, 5, 6, 6, 16, 16, 4, 4, 4};
  const int columns[15] = {3, 3, 24, 24, 5, 4, 24, 24, 24, 5, 6, 6, 4, 6, 6};
  // These owners are call-local, without a workspace or persistent cache.
  for (int i = 0; i < 15; ++i) {
    status = CreateMatrix(rows[i], columns[i], &matrices[i], d);
    if (status != StatusCode::kOk) return status;
  }
  Matrix* coordinates = matrices[0].get();
  Matrix* frame = matrices[1].get();
  Matrix* transform = matrices[2].get();
  Matrix* ties = matrices[3].get();
  Matrix* constitutive = matrices[4].get();
  Matrix* gradients = matrices[5].get();
  Matrix* shears = matrices[6].get();
  Matrix* b = matrices[7].get();
  Matrix* direct = matrices[8].get();
  Matrix* stress_selector = matrices[9].get();
  Matrix* stress_samples = matrices[10].get();
  Matrix* strain_samples = matrices[11].get();
  Matrix* extrapolation = matrices[12].get();
  Matrix* layer = matrices[13].get();
  Matrix* extrapolated = matrices[14].get();
  double diameter = 0.0;
  double length = 0.0;
  status = ValidateMitc4Geometry(nodes, count, coordinates, frame, &diameter,
                                 &length, d);
  if (status != StatusCode::kOk) return status;
  status = FillTransform(frame, transform, d);
  if (status != StatusCode::kOk) return status;
  status = BuildMitc4LinearTying(coordinates, ties, d);
  if (status != StatusCode::kOk) return status;
  status = BuildMitc4LinearConstitutive(section->material(), constitutive, d);
  if (status != StatusCode::kOk) return status;
  status = FillStressSelector(stress_selector, d);
  if (status != StatusCode::kOk) return status;
  status = FillExtrapolation(extrapolation, d);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<Vector> vectors[5];
  const int sizes[5] = {24, 5, 5, 6, 6};
  for (int i = 0; i < 5; ++i) {
    Vector* raw = nullptr;
    status = Vector::Create(sizes[i], &raw, d);
    vectors[i].reset(raw);
    if (status != StatusCode::kOk) return status;
  }
  Vector* local_u = vectors[0].get();
  Vector* assumed_strain = vectors[1].get();
  Vector* physical_stress = vectors[2].get();
  Vector* output_stress = vectors[3].get();
  Vector* output_strain = vectors[4].get();
  status = Matrix::MultiplyVector(transform, false, u, local_u, d);
  if (status != StatusCode::kOk) return status;
  // Direct integration includes frozen drilling exactly once and does not
  // create K or M. The result validity/commit belongs to the caller.
  status = ComputeMitc4LinearInternalForce(nodes, count, section, u, drilling,
                                           internal_force_, d);
  if (status != StatusCode::kOk) return status;
  status = axes_->CopyFrom(frame, d);
  if (status != StatusCode::kOk) return status;
  const double a = 1.0 / std::sqrt(3.0);
  const double plane_r[4] = {-a, a, a, -a};
  const double plane_s[4] = {-a, -a, a, a};
  const double zeta[4] = {-a, a, -1.0, 1.0};
  for (int plane = 0; plane < 4; ++plane) {
    double determinant = 0.0;
    status = ComputeMitc4LinearPointOperators(coordinates, ties, plane_r[plane],
                                              plane_s[plane], gradients, shears,
                                              &determinant, d);
    if (status != StatusCode::kOk) return status;
    for (int thickness = 0; thickness < 4; ++thickness) {
      const double z = (section->thickness() / 2.0) * zeta[thickness];
      status = FillMitc4LinearB(gradients, shears, z, b, d);
      if (status != StatusCode::kOk) return status;
      status = Matrix::MultiplyVector(b, false, local_u, assumed_strain, d);
      if (status != StatusCode::kOk) return status;
      status = Matrix::MultiplyVector(constitutive, false, assumed_strain,
                                      physical_stress, d);
      if (status != StatusCode::kOk) return status;
      status = Matrix::MultiplyVector(stress_selector, false, physical_stress,
                                      output_stress, d);
      if (status != StatusCode::kOk) return status;
      status = FillDirectStrain(b, gradients, plane_r[plane], plane_s[plane],
                                direct, d);
      if (status != StatusCode::kOk) return status;
      status = Matrix::MultiplyVector(direct, false, local_u, output_strain, d);
      if (status != StatusCode::kOk) return status;
      const int row = 4 * thickness + plane;
      status = stress_samples->SetRow(row, output_stress, d);
      if (status != StatusCode::kOk) return status;
      status = strain_samples->SetRow(row, output_strain, d);
      if (status != StatusCode::kOk) return status;
    }
  }
  status = PublishSamples(stress_samples, extrapolation, layer, extrapolated,
                          ip_stress_, surface_stress_, nodal_stress_, d);
  if (status != StatusCode::kOk) return status;
  return PublishSamples(strain_samples, extrapolation, layer, extrapolated,
                        ip_strain_, surface_strain_, nodal_strain_, d);
}

StatusCode Mitc4Results::RecoverNonlinear(const Node* const* nodes, int count,
                                          const ShellSection* section,
                                          const ElementTrialState* state,
                                          double drilling,
                                          TrialEvaluationFailure* failure,
                                          Diagnostic* d) {
  if (!failure) return model_storage::Invalid("MITC4 recovery failure flag", d);
  *failure = TrialEvaluationFailure::kNone;
  if (!std::isfinite(drilling) || !(drilling > 0))
    return model_storage::Invalid("MITC4 recovery fixed drilling coefficient",
                                  d);
  auto status = EvaluateMitc4Trial(
      nodes, count, section, state, ElementEvaluationRequest::kInternalForce,
      drilling, internal_force_, nullptr, failure, d);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<Matrix> stress, strain, extrapolation, layer, extrapolated;
  status = CreateMatrix(16, 6, &stress, d);
  if (status == StatusCode::kOk) status = CreateMatrix(16, 6, &strain, d);
  if (status == StatusCode::kOk) status = CreateMatrix(4, 4, &extrapolation, d);
  if (status == StatusCode::kOk) status = CreateMatrix(4, 6, &layer, d);
  if (status == StatusCode::kOk) status = CreateMatrix(4, 6, &extrapolated, d);
  if (status != StatusCode::kOk) return status;
  status = RecoverMitc4FiniteSamples(nodes, count, section, state->translations,
                                     state->orientations, stress.get(),
                                     strain.get(), axes_, failure, d);
  if (status == StatusCode::kOk)
    status = FillExtrapolation(extrapolation.get(), d);
  if (status == StatusCode::kOk)
    status = PublishSamples(stress.get(), extrapolation.get(), layer.get(),
                            extrapolated.get(), ip_stress_, surface_stress_,
                            nodal_stress_, d);
  if (status != StatusCode::kOk) return status;
  return PublishSamples(strain.get(), extrapolation.get(), layer.get(),
                        extrapolated.get(), ip_strain_, surface_strain_,
                        nodal_strain_, d);
}
const Vector* Mitc4Results::internal_force() const { return internal_force_; }
const Matrix* Mitc4Results::output_axes() const { return axes_; }
const Matrix* Mitc4Results::ip_stress() const { return ip_stress_; }
const Matrix* Mitc4Results::ip_strain() const { return ip_strain_; }
const Matrix* Mitc4Results::surface_stress() const { return surface_stress_; }
const Matrix* Mitc4Results::surface_strain() const { return surface_strain_; }
const Matrix* Mitc4Results::nodal_stress() const { return nodal_stress_; }
const Matrix* Mitc4Results::nodal_strain() const { return nodal_strain_; }

}  // namespace fesa
