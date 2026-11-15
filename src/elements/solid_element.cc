#include "fesa/elements/solid_element.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <initializer_list>
#include <limits>
#include <memory>
#include <new>

#include "elements/hexa_hourglass.h"
#include "elements/hexa_interpolation.h"
#include "elements/pyramid_interpolation.h"
#include "elements/solid_access.h"
#include "elements/solid_nonlinear_ops.h"
#include "elements/solid_recovery.h"
#include "elements/solid_trial_geometry.h"
#include "elements/wedge_interpolation.h"
#include "fesa/elements/element_workspace.h"
#include "fesa/math/dense_factorization.h"
#include "model/model_storage.h"

namespace fesa {
namespace {
thread_local void (*evaluation_observer)(SolidEvaluationEvent) = nullptr;
void Observe(SolidEvaluationEvent event) {
  if (evaluation_observer) evaluation_observer(event);
}
StatusCode InvalidCall(const char* operation, Diagnostic* d) {
  return ReportFailure(StatusCode::kInvalidArgument,
                       DiagnosticStage::kValidation, operation, d);
}
StatusCode InvalidGeometry(Diagnostic* d) {
  return ReportFailure(
      StatusCode::kInputError, DiagnosticStage::kValidation,
      "Tetra nonpositive or nearly singular normalized Jacobian", d);
}
StatusCode AtNaturalPoint(StatusCode status, const double* location,
                          Diagnostic* d) {
  if (status != StatusCode::kOk && d) {
    d->geometry_point_available = true;
    d->natural_xi = location[0];
    d->natural_eta = location[1];
    d->natural_zeta = location[2];
  }
  return status;
}
// User-adopted empirical Standard convergence length, not a public Abaqus
// general formula. Only the initial center mapping is needed; no LU or D.
StatusCode CenterJacobianLength(const SolidElement* element, double* out,
                                Diagnostic* d) {
  const bool wedge = element->formulation() == SolidFormulation::kWedge6;
  const int nodes = wedge ? 6 : 8;
  Matrix *coordinates = nullptr, *gradient = nullptr, *jacobian = nullptr;
  Vector *anchor = nullptr, *row = nullptr, *lengths = nullptr,
         *weights = nullptr;
  auto status = Matrix::Create(nodes, 3, &coordinates, d);
  std::unique_ptr<Matrix> co(coordinates);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Create(3, nodes, &gradient, d);
  std::unique_ptr<Matrix> go(gradient);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Create(3, 3, &jacobian, d);
  std::unique_ptr<Matrix> jo(jacobian);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(3, &anchor, d);
  std::unique_ptr<Vector> ao(anchor);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(3, &row, d);
  std::unique_ptr<Vector> ro(row);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(3, &lengths, d);
  std::unique_ptr<Vector> lo(lengths);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(3, &weights, d);
  std::unique_ptr<Vector> wo(weights);
  if (status != StatusCode::kOk) return status;
  const int signs[8][3] = {{-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1},
                           {-1, -1, 1},  {1, -1, 1},  {1, 1, 1},  {-1, 1, 1}};
  const int triangle_derivative[2][3] = {{-1, 1, 0}, {-1, 0, 1}};
  for (int axis = 0; axis < 3; ++axis) {
    double value = 0;
    status = element->GetNode(0)->Coordinate(axis, &value, d);
    if (status != StatusCode::kOk) return status;
    status = anchor->Set(axis, value, d);
    if (status != StatusCode::kOk) return status;
    status = weights->Set(axis, wedge && axis < 2 ? 1. / 3 : 2. / 3, d);
    if (status != StatusCode::kOk) return status;
  }
  for (int node = 0; node < nodes; ++node) {
    for (int axis = 0; axis < 3; ++axis) {
      double value = 0;
      status = element->GetNode(node)->Coordinate(axis, &value, d);
      if (status != StatusCode::kOk) return status;
      status = row->Set(axis, value, d);
      if (status != StatusCode::kOk) return status;
      // F4 at (1/3,1/3,0) or F3 at the Hexa center; only dN is needed.
      const double natural =
          wedge ? (axis < 2 ? triangle_derivative[axis][node % 3] / 2.
                            : (node < 3 ? -1. : 1.) / 6)
                : signs[node][axis] / 8.;
      status = gradient->Set(axis, node, natural, d);
      if (status != StatusCode::kOk) return status;
    }
    status = row->AddScaled(-1, anchor, d);
    if (status != StatusCode::kOk) return status;
    status = coordinates->SetRow(node, row, d);
    if (status != StatusCode::kOk) return status;
  }
  status = Matrix::Multiply(coordinates, true, gradient, true, jacobian, d);
  if (status != StatusCode::kOk) return status;
  for (int axis = 0; axis < 3; ++axis) {
    const int indices[] = {axis, 3 + axis, 6 + axis};
    status = jacobian->GatherTo(indices, 3, row, 0, d);
    if (status != StatusCode::kOk) return status;
    double norm = 0;
    status = row->Norm2(&norm, d);
    if (status != StatusCode::kOk) return status;
    status = lengths->Set(axis, norm, d);
    if (status != StatusCode::kOk) return status;
  }
  double length = 0;
  status = lengths->Dot(weights, &length, d);
  if (status != StatusCode::kOk) return status;
  if (!(length > 0) || !std::isfinite(length))
    return ReportFailure(StatusCode::kLibraryFailure,
                         DiagnosticStage::kComputation,
                         "Solid center-J characteristic length range", d);
  *out = length;
  return StatusCode::kOk;
}
// User-adopted empirical Pyramid length: four base samples at a,b=+/-0.51.
// Only bilinear base tangents and apex distances are needed; no apex gradient.
StatusCode PyramidLength(const SolidElement* element, double* out,
                         Diagnostic* d) {
  Matrix *coordinates = nullptr, *coefficients = nullptr, *values = nullptr;
  Vector *anchor = nullptr, *apex = nullptr, *row = nullptr, *lengths = nullptr,
         *weights = nullptr;
  auto status = Matrix::Create(4, 3, &coordinates, d);
  std::unique_ptr<Matrix> co(coordinates);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Create(12, 4, &coefficients, d);
  std::unique_ptr<Matrix> cfo(coefficients);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Create(12, 3, &values, d);
  std::unique_ptr<Matrix> vo(values);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(3, &anchor, d);
  std::unique_ptr<Vector> ao(anchor);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(3, &apex, d);
  std::unique_ptr<Vector> apo(apex);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(3, &row, d);
  std::unique_ptr<Vector> ro(row);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(12, &lengths, d);
  std::unique_ptr<Vector> lo(lengths);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(12, &weights, d);
  std::unique_ptr<Vector> wo(weights);
  if (status != StatusCode::kOk) return status;
  for (int axis = 0; axis < 3; ++axis) {
    double value = 0;
    status = element->GetNode(0)->Coordinate(axis, &value, d);
    if (status != StatusCode::kOk) return status;
    status = anchor->Set(axis, value, d);
    if (status != StatusCode::kOk) return status;
    status = element->GetNode(4)->Coordinate(axis, &value, d);
    if (status != StatusCode::kOk) return status;
    status = apex->Set(axis, value, d);
    if (status != StatusCode::kOk) return status;
  }
  status = apex->AddScaled(-1, anchor, d);
  if (status != StatusCode::kOk) return status;
  // Base node zero is the translation anchor, hence its relative row is zero.
  for (int node = 1; node < 4; ++node) {
    for (int axis = 0; axis < 3; ++axis) {
      double value = 0;
      status = element->GetNode(node)->Coordinate(axis, &value, d);
      if (status != StatusCode::kOk) return status;
      status = row->Set(axis, value, d);
      if (status != StatusCode::kOk) return status;
    }
    status = row->AddScaled(-1, anchor, d);
    if (status != StatusCode::kOk) return status;
    status = coordinates->SetRow(node, row, d);
    if (status != StatusCode::kOk) return status;
  }
  const int signs[4][2] = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}};
  for (int point = 0; point < 4; ++point) {
    const double a = .51 * signs[point][0];
    const double b = .51 * signs[point][1];
    for (int node = 0; node < 4; ++node) {
      const double coefficients_at_point[] = {
          signs[node][0] * (1 + signs[node][1] * b) / 4,
          signs[node][1] * (1 + signs[node][0] * a) / 4,
          (1 + signs[node][0] * a) * (1 + signs[node][1] * b) / 4};
      for (int component = 0; component < 3; ++component) {
        status = coefficients->Set(3 * point + component, node,
                                   coefficients_at_point[component], d);
        if (status != StatusCode::kOk) return status;
      }
    }
  }
  status = Matrix::Multiply(coefficients, false, coordinates, false, values, d);
  if (status != StatusCode::kOk) return status;
  for (int sample = 0; sample < 12; ++sample) {
    const int indices[] = {3 * sample, 3 * sample + 1, 3 * sample + 2};
    status = values->GatherTo(indices, 3, row, 0, d);
    if (status != StatusCode::kOk) return status;
    if (sample % 3 == 2) {
      status = row->AddScaled(-1, apex, d);
      if (status != StatusCode::kOk) return status;
    }
    double norm = 0;
    status = row->Norm2(&norm, d);
    if (status != StatusCode::kOk) return status;
    status = lengths->Set(sample, norm, d);
    if (status != StatusCode::kOk) return status;
    status = weights->Set(sample, sample % 3 == 2 ? 1. / 12 : .273 / 4, d);
    if (status != StatusCode::kOk) return status;
  }
  double length = 0;
  status = lengths->Dot(weights, &length, d);
  if (status != StatusCode::kOk) return status;
  if (!(length > 0) || !std::isfinite(length))
    return ReportFailure(StatusCode::kLibraryFailure,
                         DiagnosticStage::kComputation,
                         "Pyramid characteristic length range", d);
  *out = length;
  return StatusCode::kOk;
}
StatusCode CornerGradient(Matrix** out, Diagnostic* d) {
  auto status = Matrix::Create(3, 4, out, d);
  if (status != StatusCode::kOk) return status;
  for (int axis = 0; axis < 3; ++axis)
    for (int node = 0; node < 4; ++node) {
      status = (*out)->Set(axis, node,
                           node == 0          ? -1
                           : node == axis + 1 ? 1
                                              : 0,
                           d);
      if (status != StatusCode::kOk) return status;
    }
  return StatusCode::kOk;
}
// All buffers and the LU have one evaluation lifetime. Geometry scale h is
// solely the guard scale from geometry_contract.md, never convergence Lc.
class SolidGeometry {
 public:
  SolidGeometry() = default;
  ~SolidGeometry() {
    delete factor_;
    delete coordinates_;
  }
  SolidGeometry(const SolidGeometry&) = delete;
  SolidGeometry& operator=(const SolidGeometry&) = delete;
  StatusCode Initialize(const SolidElement* element, const Matrix* natural,
                        Diagnostic* d) {
    Observe(SolidEvaluationEvent::kGeometry);
    const int count = element->node_count() == 10 ? 4 : element->node_count();
    Matrix *coordinates = nullptr, *ones = nullptr, *anchor = nullptr,
           *translation = nullptr;
    Vector *first = nullptr, *second = nullptr;
    auto status = Matrix::Create(count, 3, &coordinates, d);
    std::unique_ptr<Matrix> coordinates_owner(coordinates);
    if (status != StatusCode::kOk) return status;
    status = Vector::Create(3, &first, d);
    std::unique_ptr<Vector> first_owner(first);
    if (status != StatusCode::kOk) return status;
    status = Vector::Create(3, &second, d);
    std::unique_ptr<Vector> second_owner(second);
    if (status != StatusCode::kOk) return status;
    for (int node = 0; node < count; ++node) {
      for (int axis = 0; axis < 3; ++axis) {
        double value = 0;
        status = element->GetNode(node)->Coordinate(axis, &value, d);
        if (status != StatusCode::kOk) return status;
        status = coordinates->Set(node, axis, value, d);
        if (status != StatusCode::kOk) return status;
      }
    }
    for (int i = 0; i < count; ++i) {
      const int indices_i[] = {3 * i, 3 * i + 1, 3 * i + 2};
      for (int j = i + 1; j < count; ++j) {
        const int indices_j[] = {3 * j, 3 * j + 1, 3 * j + 2};
        status = coordinates->GatherTo(indices_i, 3, first, 0, d);
        if (status != StatusCode::kOk) return status;
        status = coordinates->GatherTo(indices_j, 3, second, 0, d);
        if (status != StatusCode::kOk) return status;
        status = first->AddScaled(-1, second, d);
        if (status != StatusCode::kOk) return status;
        double length = 0;
        status = first->Norm2(&length, d);
        if (status != StatusCode::kOk) return status;
        scale_ = std::max(scale_, length);
      }
    }
    if (!(scale_ > 0)) return InvalidGeometry(d);
    if (!std::isfinite(1 / scale_))
      return ReportFailure(StatusCode::kLibraryFailure,
                           DiagnosticStage::kComputation,
                           "Tetra geometry reciprocal scale range", d);
    status = Matrix::Create(count, 1, &ones, d);
    std::unique_ptr<Matrix> ones_owner(ones);
    if (status != StatusCode::kOk) return status;
    for (int i = 0; i < count; ++i) {
      status = ones->Set(i, 0, 1, d);
      if (status != StatusCode::kOk) return status;
    }
    status = Matrix::Create(1, 3, &anchor, d);
    std::unique_ptr<Matrix> anchor_owner(anchor);
    if (status != StatusCode::kOk) return status;
    status = anchor->CopyBlockFrom(coordinates, 0, 0, 1, 3, 0, 0, d);
    if (status != StatusCode::kOk) return status;
    status = Matrix::Create(count, 3, &translation, d);
    std::unique_ptr<Matrix> translation_owner(translation);
    if (status != StatusCode::kOk) return status;
    status = Matrix::Multiply(ones, false, anchor, false, translation, d);
    if (status != StatusCode::kOk) return status;
    status = coordinates->AddScaled(-1, translation, d);
    if (status != StatusCode::kOk) return status;
    status = coordinates->Scale(1 / scale_, d);
    if (status != StatusCode::kOk) return status;
    coordinates_ = coordinates_owner.release();
    return Factor(natural, d);
  }
  StatusCode SetPoint(const Matrix* natural, Diagnostic* d) {
    return coordinates_->Rows() != 4 ? Factor(natural, d) : StatusCode::kOk;
  }
  StatusCode Factor(const Matrix* natural, Diagnostic* d) {
    Matrix* jacobian = nullptr;
    auto status = Matrix::Create(3, 3, &jacobian, d);
    std::unique_ptr<Matrix> jacobian_owner(jacobian);
    if (status != StatusCode::kOk) return status;
    // Straight quadratic geometry uses corners, not quadratic displacement
    // gradients.
    Matrix* corner = nullptr;
    std::unique_ptr<Matrix> corner_owner;
    if (natural->Columns() != coordinates_->Rows()) {
      status = CornerGradient(&corner, d);
      corner_owner.reset(corner);
      if (status != StatusCode::kOk) return status;
    }
    const Matrix* geometry_gradient = corner ? corner : natural;
    // J rows are physical axes, columns natural axes (F6).
    status = Matrix::Multiply(coordinates_, true, geometry_gradient, true,
                              jacobian, d);
    if (status != StatusCode::kOk) return status;
    delete factor_;
    factor_ = nullptr;
    Observe(SolidEvaluationEvent::kJacobianFactor);
    status = DenseLuFactorization::Create(jacobian, &factor_, d);
    if (status == StatusCode::kSingular) return InvalidGeometry(d);
    if (status != StatusCode::kOk) return status;
    double normalized = 0;
    status = factor_->Determinant(&normalized, d);
    if (status != StatusCode::kOk) return status;
    if (!(normalized > 256 * std::numeric_limits<double>::epsilon()))
      return InvalidGeometry(d);
    determinant_ = normalized * scale_ * scale_ * scale_;
    if (!(determinant_ > 0) || !std::isfinite(determinant_))
      return ReportFailure(StatusCode::kLibraryFailure,
                           DiagnosticStage::kComputation,
                           "Tetra physical determinant range", d);
    return StatusCode::kOk;
  }
  StatusCode Gradient(const Matrix* natural, Matrix* physical,
                      Diagnostic* d) const {
    auto status = factor_->Solve(natural, true, physical, d);
    if (status != StatusCode::kOk) return status;
    return physical->Scale(1 / scale_, d);
  }
  double determinant() const { return determinant_; }
  double scale() const { return scale_; }
  const Matrix* normalized_coordinates() const { return coordinates_; }

 private:
  Matrix* coordinates_ =
      nullptr;  // Owned normalized coordinates, call lifetime.
  DenseLuFactorization* factor_ = nullptr;
  double scale_ = 0;
  double determinant_ = 0;
};
StatusCode CheckRule(const SolidElement* element,
                     const ElementWorkspace* context, bool mass,
                     const SolidReferenceData** out, Diagnostic* d) {
  if (!context) return InvalidCall("Solid evaluation context", d);
  const auto* data = context->reference(
      element->formulation(),
      mass ? IntegrationPurpose::kMass : IntegrationPurpose::kStiffness);
  const int nodes = element->node_count();
  const bool reduced =
      element->formulation() == SolidFormulation::kHexa8Reduced;
  const bool hexa =
      element->formulation() == SolidFormulation::kHexa8Full || reduced;
  const bool wedge = element->formulation() == SolidFormulation::kWedge6;
  const bool pyramid =
      element->formulation() == SolidFormulation::kPyramid5Full ||
      element->formulation() == SolidFormulation::kPyramid5Reduced;
  const int points = mass ? (pyramid       ? 8
                             : wedge       ? 18
                             : hexa        ? 27
                             : nodes == 10 ? 15
                                           : 4)
                          : element->integration_point_count();
  const char* rule = pyramid ? (points == 1   ? "PYR1"
                                : points == 5 ? "PYR5_FELIPPA"
                                              : "PYR8")
                     : wedge ? (mass ? "WEDGE18" : "WEDGE2")
                     : hexa  ? (mass      ? "HEX27"
                                : reduced ? "UNIFORM_STRAIN_STANDARD_STIFFNESS"
                                          : "HEX8_BBAR")
                     : points == 15 ? "TET15"
                     : points == 4  ? "TET4"
                                    : "TET1";
  const auto purpose =
      mass ? IntegrationPurpose::kMass : IntegrationPurpose::kStiffness;
  if (!data || data->purpose() != purpose || data->node_count() != nodes ||
      data->point_count() != points || std::strcmp(data->rule_id(), rule) != 0)
    return InvalidCall("Solid mismatched purpose reference data", d);
  if (reduced && !mass) {
    const auto* geometry = context->reference(element->formulation(),
                                              IntegrationPurpose::kGeometry);
    if (!geometry || geometry->purpose() != IntegrationPurpose::kGeometry ||
        geometry->node_count() != 8 || geometry->point_count() != 27 ||
        std::strcmp(geometry->rule_id(), "HEX27") != 0)
      return InvalidCall("Uniform Hexa missing geometry reference data", d);
  }
  *out = data;
  return StatusCode::kOk;
}
// Scalar placement of the isotropic F8 constitutive coefficients. The matrix
// is local to this call; it is never cached in Material or the element.
StatusCode Constitutive(const ElasticMaterial* material, Matrix* out,
                        Diagnostic* d) {
  Observe(SolidEvaluationEvent::kConstitutive);
  const double e = material->young(), nu = material->poisson();
  const double mu = (e / 2) / (1 + nu);
  const double lambda = ((e / (1 + nu)) / (1 - 2 * nu)) * nu;
  for (int row = 0; row < 6; ++row) {
    for (int col = 0; col < 6; ++col) {
      const double value = row < 3 && col < 3
                               ? lambda + (row == col ? 2 * mu : 0)
                           : row == col ? mu
                                        : 0;
      if (!std::isfinite(value))
        return InvalidCall("Solid elastic coefficient range", d);
      auto status = out->Set(row, col, value, d);
      if (status != StatusCode::kOk) return status;
    }
  }
  return StatusCode::kOk;
}
StatusCode StrainDisplacement(const Matrix* gradient, Matrix* b,
                              Diagnostic* d) {
  // F7 engineering order [xx,yy,zz,xy,yz,xz]. Each entry is a scalar
  // derivative coefficient; products and accumulation are external math.
  const int nodes = gradient->Columns();
  for (int node = 0; node < nodes; ++node) {
    const double x = gradient->Data()[node], y = gradient->Data()[nodes + node],
                 z = gradient->Data()[2 * nodes + node];
    const double coefficients[6][3] = {{x, 0, 0}, {0, y, 0}, {0, 0, z},
                                       {y, x, 0}, {0, z, y}, {z, 0, x}};
    for (int row = 0; row < 6; ++row)
      for (int axis = 0; axis < 3; ++axis) {
        auto status = b->Set(row, 3 * node + axis, coefficients[row][axis], d);
        if (status != StatusCode::kOk) return status;
      }
  }
  return StatusCode::kOk;
}
// One call's numerical buffers; no new workspace policy or retained caches.
class SolidEvaluation {
 public:
  explicit SolidEvaluation(SolidGeometry* geometry) : geometry_(geometry) {}
  ~SolidEvaluation() {
    delete gradient_;
    delete b_;
    delete constitutive_;
    delete db_;
    delete strain_;
    delete stress_;
    delete force_;
    delete hourglass_;
    delete point_gradients_;
    delete measures_;
    delete delta_;
    delete correction_;
    delete correction_row_;
  }
  SolidEvaluation(const SolidEvaluation&) = delete;
  SolidEvaluation& operator=(const SolidEvaluation&) = delete;
  StatusCode Initialize(const SolidElement* element,
                        const SolidReferenceData* data,
                        const ElementWorkspace* context, bool stiffness,
                        Diagnostic* d, bool include_hourglass = true) {
    reduced_ = element->formulation() == SolidFormulation::kHexa8Reduced;
    selective_ = element->formulation() == SolidFormulation::kHexa8Full;
    const auto* geometry_data =
        reduced_ ? context->reference(element->formulation(),
                                      IntegrationPurpose::kGeometry)
                 : data;
    auto status =
        geometry_->Initialize(element, geometry_data->natural_gradients(0), d);
    if (status != StatusCode::kOk)
      return AtNaturalPoint(status,
                            geometry_data->natural_coordinates()->Data(), d);
    status = Matrix::Create(3, element->node_count(), &gradient_, d);
    if (status != StatusCode::kOk) return status;
    if (selective_) {
      status = Matrix::Create(24, 8, &point_gradients_, d);
      if (status != StatusCode::kOk) return status;
      status = Vector::Create(8, &measures_, d);
      if (status != StatusCode::kOk) return status;
      status = Matrix::Create(3, 8, &delta_, d);
      if (status != StatusCode::kOk) return status;
      status = Matrix::Create(6, 24, &correction_, d);
      if (status != StatusCode::kOk) return status;
      status = Vector::Create(24, &correction_row_, d);
      if (status != StatusCode::kOk) return status;
    }
    status =
        (reduced_ || selective_)
            ? AverageGradient(geometry_data, d)
            : geometry_->Gradient(data->natural_gradients(0), gradient_, d);
    if (status != StatusCode::kOk) return status;
    status = Matrix::Create(6, element->dof_count(), &b_, d);
    if (status != StatusCode::kOk) return status;
    status = selective_ ? SelectivePoint(0, d)
                        : StrainDisplacement(gradient_, b_, d);
    if (status != StatusCode::kOk) return status;
    status = Matrix::Create(6, 6, &constitutive_, d);
    if (status != StatusCode::kOk) return status;
    status = Constitutive(element->material(), constitutive_, d);
    if (status != StatusCode::kOk) return status;
    if (stiffness) {
      status = Matrix::Create(6, element->dof_count(), &db_, d);
      if (status != StatusCode::kOk) return status;
    } else {
      status = Vector::Create(6, &strain_, d);
      if (status != StatusCode::kOk) return status;
      status = Vector::Create(6, &stress_, d);
      if (status != StatusCode::kOk) return status;
      status = Vector::Create(element->dof_count(), &force_, d);
      if (status != StatusCode::kOk) return status;
    }
    if (reduced_ && include_hourglass) {
      Observe(SolidEvaluationEvent::kHourglass);
      return HexaHourglass::Create(
          geometry_->normalized_coordinates(), geometry_->scale(), gradient_,
          volume_, element->material(), stiffness, &hourglass_, d);
    }
    return StatusCode::kOk;
  }
  StatusCode Point(const SolidReferenceData* data, int point, Diagnostic* d) {
    if (selective_)
      return AtNaturalPoint(SelectivePoint(point, d),
                            data->natural_coordinates()->Data() + 3 * point, d);
    auto status = geometry_->SetPoint(data->natural_gradients(point), d);
    if (status != StatusCode::kOk)
      return AtNaturalPoint(status,
                            data->natural_coordinates()->Data() + 3 * point, d);
    status = geometry_->Gradient(data->natural_gradients(point), gradient_, d);
    return status == StatusCode::kOk ? StrainDisplacement(gradient_, b_, d)
                                     : status;
  }
  StatusCode Force(const Vector* u, double measure, Diagnostic* d) {
    auto status = Matrix::MultiplyVector(b_, false, u, strain_, d);
    if (status != StatusCode::kOk) return status;
    status = Matrix::MultiplyVector(constitutive_, false, strain_, stress_, d);
    if (status != StatusCode::kOk) return status;
    status = Matrix::MultiplyVector(b_, true, stress_, force_, d);
    if (status != StatusCode::kOk) return status;
    status = force_->Scale(measure, d);
    if (status != StatusCode::kOk) return status;
    return hourglass_ ? hourglass_->AddForce(u, force_, d) : StatusCode::kOk;
  }
  double Measure(const SolidReferenceData* data, int point) const {
    return reduced_ ? volume_
           : selective_
               ? measures_->Data()[point]
               : data->weights()->Data()[point] * geometry_->determinant();
  }
  StatusCode AddHourglassStiffness(Matrix* out, Diagnostic* d) {
    return hourglass_ ? hourglass_->AddStiffness(out, d) : StatusCode::kOk;
  }
  Matrix* b() const { return b_; }
  Matrix* constitutive() const { return constitutive_; }
  Matrix* db() const { return db_; }
  Vector* strain() const { return strain_; }
  Vector* stress() const { return stress_; }
  Vector* force() const { return force_; }

 private:
  // The eight physical gradients/measures were evaluated once in Initialize.
  // Only volumetric rows change; shear rows retain the local gradient.
  StatusCode SelectivePoint(int point, Diagnostic* d) {
    auto status =
        delta_->CopyBlockFrom(point_gradients_, 3 * point, 0, 3, 8, 0, 0, d);
    if (status != StatusCode::kOk) return status;
    status = StrainDisplacement(delta_, b_, d);
    if (status != StatusCode::kOk) return status;
    status = delta_->Scale(-1, d);
    if (status != StatusCode::kOk) return status;
    status = delta_->AddScaled(1, gradient_, d);
    if (status != StatusCode::kOk) return status;
    status = delta_->Scale(1. / 3, d);
    if (status != StatusCode::kOk) return status;
    int indices[24] = {};
    for (int node = 0; node < 8; ++node)
      for (int axis = 0; axis < 3; ++axis)
        indices[3 * node + axis] = 8 * axis + node;
    status = delta_->GatherTo(indices, 24, correction_row_, 0, d);
    if (status != StatusCode::kOk) return status;
    for (int row = 0; row < 3; ++row) {
      status = correction_->SetRow(row, correction_row_, d);
      if (status != StatusCode::kOk) return status;
    }
    return b_->AddScaled(1, correction_, d);
  }
  StatusCode AverageGradient(const SolidReferenceData* data, Diagnostic* d) {
    Matrix* current = nullptr;
    Vector* weights = measures_;
    auto status = Matrix::Create(3, 8, &current, d);
    std::unique_ptr<Matrix> co(current);
    if (status != StatusCode::kOk) return status;
    if (!weights) status = Vector::Create(data->point_count(), &weights, d);
    std::unique_ptr<Vector> wo(selective_ ? nullptr : weights);
    if (status != StatusCode::kOk) return status;
    for (int p = 0; p < data->point_count(); ++p) {
      if (p > 0) {
        status = geometry_->SetPoint(data->natural_gradients(p), d);
        if (status != StatusCode::kOk)
          return AtNaturalPoint(status,
                                data->natural_coordinates()->Data() + 3 * p, d);
      }
      status = geometry_->Gradient(data->natural_gradients(p), current, d);
      if (status != StatusCode::kOk) return status;
      if (selective_) {
        status =
            point_gradients_->CopyBlockFrom(current, 0, 0, 3, 8, 3 * p, 0, d);
        if (status != StatusCode::kOk) return status;
      }
      const double measure =
          data->weights()->Data()[p] * geometry_->determinant();
      status = weights->Set(p, measure, d);
      if (status != StatusCode::kOk) return status;
      status = gradient_->AddScaled(measure, current, d);
      if (status != StatusCode::kOk) return status;
    }
    status = weights->AbsoluteSum(&volume_, d);
    if (status != StatusCode::kOk) return status;
    if (!(volume_ > 0) || !std::isfinite(1 / volume_))
      return ReportFailure(StatusCode::kLibraryFailure,
                           DiagnosticStage::kComputation,
                           "Uniform Hexa volume range", d);
    return gradient_->Scale(1 / volume_, d);
  }
  HexaHourglass* hourglass_ = nullptr;  // Owned only in reduced evaluations.
  bool reduced_ = false;
  bool selective_ = false;
  Matrix *point_gradients_ = nullptr, *delta_ = nullptr, *correction_ = nullptr;
  Vector *measures_ = nullptr, *correction_row_ = nullptr;
  double volume_ = 0;
  Matrix *gradient_ = nullptr, *b_ = nullptr, *constitutive_ = nullptr,
         *db_ = nullptr;
  Vector *strain_ = nullptr, *stress_ = nullptr, *force_ = nullptr;
  SolidGeometry* geometry_;  // Borrowed local geometry, outlives this helper.
};
// Call-local owner shared by trial and recovery, never retained by Element.
class SolidTrialData {
 public:
  SolidTrialData() = default;
  SolidTrialData(const SolidTrialData&) = delete;
  SolidTrialData& operator=(const SolidTrialData&) = delete;
  ~SolidTrialData() {
    for (const auto* gradient : gradients_) delete gradient;
    delete weights_;
  }
  StatusCode Initialize(const SolidElement* element,
                        const SolidReferenceData* data, Diagnostic* d) {
    // M1 uses the same anchored/scaled initial Jacobian and LU as linear
    // evaluation. There is no B-bar average here: F-bar needs each actual G/q.
    SolidGeometry geometry;
    auto status = geometry.Initialize(element, data->natural_gradients(0), d);
    if (status != StatusCode::kOk) return status;
    status = Vector::Create(data->point_count(), &weights_, d);
    if (status != StatusCode::kOk) return status;
    for (int point = 0; point < data->point_count(); ++point) {
      const auto* natural = data->natural_gradients(point);
      if (point > 0) {
        status = geometry.SetPoint(natural, d);
        if (status != StatusCode::kOk) return status;
      }
      Matrix* gradient = nullptr;
      status = Matrix::Create(3, element->node_count(), &gradient, d);
      gradients_[point] = gradient;
      if (status != StatusCode::kOk) return status;
      status = geometry.Gradient(natural, gradient, d);
      if (status != StatusCode::kOk) return status;
      const double measure =
          data->weights()->Data()[point] * geometry.determinant();
      if (!(measure > 0) || !std::isfinite(measure))
        return ReportFailure(StatusCode::kLibraryFailure,
                             DiagnosticStage::kComputation,
                             "Solid nonlinear integration measure range", d);
      status = weights_->Set(point, measure, d);
      if (status != StatusCode::kOk) return status;
    }
    return StatusCode::kOk;
  }
  const Matrix* const* gradients() const { return gradients_; }
  const Vector* weights() const { return weights_; }

 private:
  const Matrix* gradients_[8] = {};  // Owned physical gradients, call lifetime.
  Vector* weights_ = nullptr;       // Owned initial measures, call lifetime.
};
StatusCode ValidateSolidTrialInput(const SolidElement* element,
    const ElementTrialState* state, const ElementWorkspace* workspace,
    const SolidReferenceData** data, Diagnostic* d) {
  if (!state || !state->translations || state->orientations ||
      state->committed_orientations || state->body_spin ||
      state->translations->Size() != element->dof_count())
    return InvalidCall("Solid nonlinear trial state", d);
  for (int i = 0; i < element->dof_count(); ++i)
    if (!std::isfinite(state->translations->Data()[i]))
      return InvalidCall("Solid nonlinear nonfinite translation", d);
  return CheckRule(element, workspace, false, data, d);
}
}  // namespace
void SetSolidEvaluationObserver(void (*observer)(SolidEvaluationEvent)) {
  evaluation_observer = observer;
}
void ObserveSolidEvaluation(SolidEvaluationEvent event) {
  Observe(event);
}
SolidElement::SolidElement(int id, int nodes, int points,
                           const SolidSection* section,
                           SolidFormulation formulation)
    : Element(id),
      node_count_(nodes),
      point_count_(points),
      formulation_(formulation),
      section_(section),
      material_(section->material()) {}
StatusCode SolidElement::Initialize(const Node* const* nodes, Diagnostic* d) {
  if (!nodes) return model_storage::Invalid("Solid nodes", d);
  for (int i = 0; i < node_count_; ++i) {
    if (!nodes[i] || nodes[i]->dof_count() != 3)
      return model_storage::Invalid("Solid node layout", d);
    for (int j = 0; j < i; ++j)
      if (nodes[i]->Id() == nodes[j]->Id())
        return model_storage::Invalid("Solid repeated node ID", d);
  }
  nodes_ = new (std::nothrow) const Node*[node_count_];
  if (!nodes_) return model_storage::AllocationFailure("Solid connectivity", d);
  for (int i = 0; i < node_count_; ++i) nodes_[i] = nodes[i];
  auto status = Vector::Create(3 * node_count_, &internal_force_, d);
  if (status != StatusCode::kOk) return status;
  for (auto** store : {&ip_stress_, &ip_strain_}) {
    status = Matrix::Create(point_count_, 6, store, d);
    if (status != StatusCode::kOk) return status;
  }
  for (auto** store : {&nodal_stress_, &nodal_strain_}) {
    status = Matrix::Create(node_count_, 6, store, d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}
StatusCode SolidElement::ValidateGeometry(Diagnostic* d) const {
  const bool wedge = formulation_ == SolidFormulation::kWedge6;
  const bool pyramid = formulation_ == SolidFormulation::kPyramid5Full ||
                       formulation_ == SolidFormulation::kPyramid5Reduced;
  if (formulation_ == SolidFormulation::kHexa8Full ||
      formulation_ == SolidFormulation::kHexa8Reduced || wedge || pyramid) {
    Matrix *shape = nullptr, *natural = nullptr;
    auto status = Matrix::Create(1, node_count_, &shape, d);
    std::unique_ptr<Matrix> so(shape);
    if (status != StatusCode::kOk) return status;
    status = Matrix::Create(3, node_count_, &natural, d);
    std::unique_ptr<Matrix> no(natural);
    if (status != StatusCode::kOk) return status;
    double points[62][3] = {};
    int count = 0;
    auto add = [&](double x, double y, double z) {
      for (int i = 0; i < count; ++i)
        if (points[i][0] == x && points[i][1] == y && points[i][2] == z) return;
      points[count][0] = x;
      points[count][1] = y;
      points[count++][2] = z;
    };
    if (wedge) {
      const double triangle[][2] = {{0, 0},   {1, 0},  {0, 1},          {.5, 0},
                                    {.5, .5}, {0, .5}, {1. / 3, 1. / 3}};
      for (int z = -1; z <= 1; ++z)
        for (const auto& point : triangle) add(point[0], point[1], z);
    } else if (pyramid) {
      for (int y = -1; y <= 1; ++y)
        for (int x = -1; x <= 1; ++x) add(x, y, 0);
      for (double y : {-.5, .5})
        for (double x : {-.5, .5}) add(x, y, .5);
      add(0, 0, .25);  // PYR1; apex gradient is deliberately never requested.
    } else {
      for (int z = -1; z <= 1; ++z)
        for (int y = -1; y <= 1; ++y)
          for (int x = -1; x <= 1; ++x) add(x, y, z);
    }
    for (auto purpose :
         {IntegrationPurpose::kStiffness, IntegrationPurpose::kMass}) {
      SolidReferenceData* raw = nullptr;
      status =
          SolidReferenceData::Create(pyramid ? SolidFormulation::kPyramid5Full
                                     : wedge ? SolidFormulation::kWedge6
                                             : SolidFormulation::kHexa8Full,
                                     purpose, &raw, d);
      std::unique_ptr<SolidReferenceData> data(raw);
      if (status != StatusCode::kOk) return status;
      for (int i = 0; i < data->point_count(); ++i) {
        const double* p = data->natural_coordinates()->Data() + 3 * i;
        add(p[0], p[1], p[2]);
      }
    }
    SolidGeometry geometry;
    for (int i = 0; i < count; ++i) {
      status = pyramid
                   ? FillPyramidInterpolation(points[i][0], points[i][1],
                                              points[i][2], shape, natural, d)
               : wedge ? FillWedgeInterpolation(points[i][0], points[i][1],
                                                points[i][2], shape, natural, d)
                       : FillHexaInterpolation(points[i][0], points[i][1],
                                               points[i][2], shape, natural, d);
      if (status != StatusCode::kOk) return status;
      status = i == 0 ? geometry.Initialize(this, natural, d)
                      : geometry.SetPoint(natural, d);
      if (status != StatusCode::kOk)
        return AtNaturalPoint(status, points[i], d);
    }
    return StatusCode::kOk;
  }
  Matrix* natural = nullptr;
  auto status = Matrix::Create(3, 4, &natural, d);
  std::unique_ptr<Matrix> owner(natural);
  if (status != StatusCode::kOk) return status;
  for (int axis = 0; axis < 3; ++axis)
    for (int node = 0; node < 4; ++node) {
      status = natural->Set(axis, node,
                            node == 0          ? -1
                            : node == axis + 1 ? 1
                                               : 0,
                            d);
      if (status != StatusCode::kOk) return status;
    }
  SolidGeometry geometry;
  status = geometry.Initialize(this, natural, d);
  if (status != StatusCode::kOk) {
    const double location[] = {0, 0, 0};
    return AtNaturalPoint(status, location, d);
  }
  if (node_count_ != 10) return status;
  Vector *first = nullptr, *second = nullptr, *error = nullptr;
  status = Vector::Create(3, &first, d);
  std::unique_ptr<Vector> first_owner(first);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(3, &second, d);
  std::unique_ptr<Vector> second_owner(second);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(3, &error, d);
  std::unique_ptr<Vector> error_owner(error);
  if (status != StatusCode::kOk) return status;
  const int edges[6][2] = {{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}};
  for (int edge = 0; edge < 6; ++edge) {
    for (int axis = 0; axis < 3; ++axis) {
      double a = 0, b = 0, m = 0;
      status = GetNode(edges[edge][0])->Coordinate(axis, &a, d);
      if (status != StatusCode::kOk) return status;
      status = GetNode(edges[edge][1])->Coordinate(axis, &b, d);
      if (status != StatusCode::kOk) return status;
      status = GetNode(4 + edge)->Coordinate(axis, &m, d);
      if (status != StatusCode::kOk) return status;
      status = first->Set(axis, a, d);
      if (status != StatusCode::kOk) return status;
      status = second->Set(axis, b, d);
      if (status != StatusCode::kOk) return status;
      status = error->Set(axis, m, d);
      if (status != StatusCode::kOk) return status;
    }
    status = first->Scale(.5, d);
    if (status != StatusCode::kOk) return status;
    status = first->AddScaled(.5, second, d);
    if (status != StatusCode::kOk) return status;
    status = error->AddScaled(-1, first, d);
    if (status != StatusCode::kOk) return status;
    double deviation = 0;
    status = error->Norm2(&deviation, d);
    if (status != StatusCode::kOk) return status;
    if (deviation / geometry.scale() > 1e-10) {
      if (d) d->node_id = GetNode(4 + edge)->Id();
      return ReportFailure(
          StatusCode::kInputError, DiagnosticStage::kValidation,
          "Tetra10 curved midside outside straight tolerance", d);
    }
  }
  return StatusCode::kOk;
}
SolidElement::~SolidElement() {
  delete[] nodes_;
  delete internal_force_;
  delete ip_stress_;
  delete ip_strain_;
  delete nodal_stress_;
  delete nodal_strain_;
}
ElementFamily SolidElement::family() const { return ElementFamily::kSolid; }
const char* SolidElement::type_name() const {
  if (formulation_ == SolidFormulation::kPyramid5Full) return "C3D5";
  if (formulation_ == SolidFormulation::kPyramid5Reduced) return "FESA_C3D5R";
  if (formulation_ == SolidFormulation::kWedge6) return "C3D6";
  if (formulation_ == SolidFormulation::kHexa8Reduced) return "C3D8R";
  return formulation_ == SolidFormulation::kHexa8Full     ? "C3D8"
         : formulation_ == SolidFormulation::kTetra10Full ? "C3D10"
                                                         : "C3D4";
}
SolidFormulation SolidElement::formulation() const { return formulation_; }
int SolidElement::node_count() const { return node_count_; }
int SolidElement::dofs_per_node() const { return 3; }
int SolidElement::dof_count() const { return 3 * node_count_; }
const Node* SolidElement::GetNode(int i) const {
  return i >= 0 && i < node_count_ ? nodes_[i] : nullptr;
}
StatusCode SolidElement::CharacteristicLength(double* out,
                                              Diagnostic* d) const {
  if (!out) return InvalidCall("Solid characteristic length output", d);
  if (formulation_ == SolidFormulation::kHexa8Full ||
      formulation_ == SolidFormulation::kHexa8Reduced ||
      formulation_ == SolidFormulation::kWedge6)
    return CenterJacobianLength(this, out, d);
  if (formulation_ == SolidFormulation::kPyramid5Full ||
      formulation_ == SolidFormulation::kPyramid5Reduced)
    return PyramidLength(this, out, d);
  if (formulation_ != SolidFormulation::kTetra4Full &&
      formulation_ != SolidFormulation::kTetra10Full)
    return Element::CharacteristicLength(out, d);
  Vector *anchor = nullptr, *edge = nullptr, *lengths = nullptr,
         *weights = nullptr;
  auto status = Vector::Create(3, &anchor, d);
  std::unique_ptr<Vector> anchor_owner(anchor);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(3, &edge, d);
  std::unique_ptr<Vector> edge_owner(edge);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(3, &lengths, d);
  std::unique_ptr<Vector> lengths_owner(lengths);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(3, &weights, d);
  std::unique_ptr<Vector> weights_owner(weights);
  if (status != StatusCode::kOk) return status;
  for (int axis = 0; axis < 3; ++axis) {
    double value = 0;
    status = GetNode(0)->Coordinate(axis, &value, d);
    if (status != StatusCode::kOk) return status;
    status = anchor->Set(axis, value, d);
    if (status != StatusCode::kOk) return status;
    status = weights->Set(axis, axis == 2 ? 2. / 3 : 1. / 3, d);
    if (status != StatusCode::kOk) return status;
  }
  for (int node = 1; node < 4; ++node) {
    for (int axis = 0; axis < 3; ++axis) {
      double value = 0;
      status = GetNode(node)->Coordinate(axis, &value, d);
      if (status != StatusCode::kOk) return status;
      status = edge->Set(axis, value, d);
      if (status != StatusCode::kOk) return status;
    }
    status = edge->AddScaled(-1, anchor, d);
    if (status != StatusCode::kOk) return status;
    double length = 0;
    status = edge->Norm2(&length, d);
    if (status != StatusCode::kOk) return status;
    status = lengths->Set(node - 1, length, d);
    if (status != StatusCode::kOk) return status;
  }
  double length = 0;
  status = lengths->Dot(weights, &length, d);
  if (status != StatusCode::kOk) return status;
  if (!(length > 0) || !std::isfinite(length))
    return ReportFailure(StatusCode::kLibraryFailure,
                         DiagnosticStage::kComputation,
                         "Solid characteristic length range", d);
  *out = length;
  return StatusCode::kOk;
}
StatusCode SolidElement::ComputeStiffness(ElementWorkspace* context,
                                          Matrix* out, Diagnostic* d) const {
  return ComputeStiffnessReadOnly(context, out, d);
}
StatusCode SolidElement::ComputeStiffnessReadOnly(
    const ElementWorkspace* context, Matrix* out, Diagnostic* d,
    bool include_hourglass) const {
  const int size = dof_count();
  if (!out || out->Rows() != size || out->Columns() != size)
    return InvalidCall("Tetra stiffness shape", d);
  const SolidReferenceData* data = nullptr;
  auto status = CheckRule(this, context, false, &data, d);
  if (status != StatusCode::kOk) return status;
  Observe(SolidEvaluationEvent::kStiffness);
  SolidGeometry geometry;
  SolidEvaluation work(&geometry);
  status = work.Initialize(this, data, context, true, d, include_hourglass);
  if (status != StatusCode::kOk) return status;
  for (int p = 0; p < data->point_count(); ++p) {
    if (p > 0) {
      status = work.Point(data, p, d);
      if (status != StatusCode::kOk) return status;
    }
    status = Matrix::Multiply(work.constitutive(), false, work.b(), false,
                              work.db(), d);
    if (status != StatusCode::kOk) return status;
    status = Matrix::MultiplyAdd(work.b(), true, work.db(), false,
                                 work.Measure(data, p), p == 0 ? 0 : 1, out, d);
    if (status != StatusCode::kOk) return status;
  }
  return work.AddHourglassStiffness(out, d);
}
StatusCode SolidElement::ComputeInternalForce(const Vector* u,
                                              ElementWorkspace* context,
                                              Vector* out,
                                              Diagnostic* d) const {
  const int size = dof_count();
  if (!u || !out || u == out || u->Size() != size || out->Size() != size)
    return InvalidCall("Tetra force shape/alias", d);
  const SolidReferenceData* data = nullptr;
  auto status = CheckRule(this, context, false, &data, d);
  if (status != StatusCode::kOk) return status;
  Observe(SolidEvaluationEvent::kForce);
  SolidGeometry geometry;
  SolidEvaluation work(&geometry);
  status = work.Initialize(this, data, context, false, d);
  if (status != StatusCode::kOk) return status;
  for (int p = 0; p < data->point_count(); ++p) {
    if (p > 0) {
      status = work.Point(data, p, d);
      if (status != StatusCode::kOk) return status;
    }
    status = work.Force(u, work.Measure(data, p), d);
    if (status != StatusCode::kOk) return status;
    status = p == 0 ? out->CopyFrom(work.force(), d)
                    : out->AddScaled(1, work.force(), d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}
bool SolidElement::SupportsMass(MassMode mode) const {
  return mode == MassMode::kConsistent ||
         (node_count_ != 10 && mode == MassMode::kRowSum);
}
StatusCode SolidElement::ComputeMass(MassMode mode, ElementWorkspace* context,
                                     Matrix* out, Diagnostic* d) const {
  return ComputeMassReadOnly(mode, context, out, d);
}
StatusCode SolidElement::ComputeMassReadOnly(MassMode mode,
                                             const ElementWorkspace* context,
                                             Matrix* out, Diagnostic* d) const {
  if (!SupportsMass(mode))
    return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                         "Tetra mass mode", d);
  const int size = dof_count();
  if (!out || out->Rows() != size || out->Columns() != size)
    return InvalidCall("Tetra mass shape", d);
  const SolidReferenceData* data = nullptr;
  auto status = CheckRule(this, context, true, &data, d);
  if (status != StatusCode::kOk) return status;
  double density = 0;
  status = material_->Density(&density, d);
  if (status != StatusCode::kOk) return status;
  Observe(SolidEvaluationEvent::kMass);
  SolidGeometry geometry;
  status = geometry.Initialize(this, data->natural_gradients(0), d);
  if (status != StatusCode::kOk)
    return AtNaturalPoint(status, data->natural_coordinates()->Data(), d);
  Matrix *h = nullptr, *consistent = nullptr;
  status = Matrix::Create(3, size, &h, d);
  std::unique_ptr<Matrix> h_owner(h);
  if (status != StatusCode::kOk) return status;
  if (mode == MassMode::kRowSum) {
    status = Matrix::Create(size, size, &consistent, d);
    if (status != StatusCode::kOk) return status;
  }
  std::unique_ptr<Matrix> consistent_owner(consistent);
  Matrix* target = consistent ? consistent : out;
  for (int p = 0; p < data->point_count(); ++p) {
    if (p > 0) {
      status = geometry.SetPoint(data->natural_gradients(p), d);
      if (status != StatusCode::kOk)
        return AtNaturalPoint(status,
                              data->natural_coordinates()->Data() + 3 * p, d);
    }
    for (int node = 0; node < node_count_; ++node)
      for (int axis = 0; axis < 3; ++axis) {
        status =
            h->Set(axis, 3 * node + axis,
                   data->shape_values()->Data()[node_count_ * p + node], d);
        if (status != StatusCode::kOk) return status;
      }
    const double weight =
        density * data->weights()->Data()[p] * geometry.determinant();
    status = Matrix::MultiplyAdd(h, true, h, false, weight, p == 0 ? 0 : 1,
                                 target, d);
    if (status != StatusCode::kOk) return status;
  }
  return consistent ? LumpFromConsistentMass(consistent, out, d)
                    : StatusCode::kOk;
}
bool SolidElement::SupportsDynamicEvaluation(
    DynamicKinematics kinematics) const {
  if (kinematics != DynamicKinematics::kLinear) return false;
  switch (formulation_) {
    case SolidFormulation::kTetra4Full:
    case SolidFormulation::kTetra10Full:
    case SolidFormulation::kHexa8Full:
    case SolidFormulation::kHexa8Reduced:
    case SolidFormulation::kWedge6:
    case SolidFormulation::kPyramid5Full:
    case SolidFormulation::kPyramid5Reduced:
      return true;
    default:
      return false;
  }
}
StatusCode SolidElement::EvaluateDynamic(const ElementMotionInput* input,
                                         MassMode mass_mode,
                                         const ElementWorkspace* workspace,
                                         DynamicEvaluationRequest request,
                                         ElementDynamicResult* out,
                                         TrialEvaluationFailure* failure,
                                         Diagnostic* d) const {
  if (failure) *failure = TrialEvaluationFailure::kNone;
  if (!SupportsDynamicEvaluation(DynamicKinematics::kLinear) ||
      (input && input->kinematics != DynamicKinematics::kLinear))
    return Element::EvaluateDynamic(input, mass_mode, workspace, request, out,
                                    failure, d);
  if (!failure) return InvalidCall("Solid dynamic failure output", d);
  auto status =
      ValidateElementDynamicEvaluation(node_count_, 3, input, request, out, d);
  if (status != StatusCode::kOk) return status;
  if (!SupportsMass(mass_mode))
    return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                         "Solid dynamic mass mode", d);
  Matrix* mass = nullptr;
  status = Matrix::Create(dof_count(), dof_count(), &mass, d);
  std::unique_ptr<Matrix> mass_owner(mass);
  if (status != StatusCode::kOk) return status;
  status = ComputeMassReadOnly(mass_mode, workspace, mass, d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::MultiplyVector(mass, false, input->acceleration,
                                  out->inertia_force, d);
  if (status != StatusCode::kOk) return status;
  Matrix* raw_damping = nullptr;
  status = Matrix::Create(dof_count(), dof_count(), &raw_damping, d);
  std::unique_ptr<Matrix> damping(raw_damping);
  if (status != StatusCode::kOk) return status;
  if (material_->rayleigh_beta() != 0) {
    status = ComputeStiffnessReadOnly(workspace, damping.get(), d, false);
    if (status == StatusCode::kOk)
      status = damping->Scale(material_->rayleigh_beta(), d);
  }
  if (status == StatusCode::kOk && material_->rayleigh_alpha() != 0)
    status = damping->AddScaled(material_->rayleigh_alpha(), mass, d);
  if (status == StatusCode::kOk)
    status = Matrix::MultiplyVector(damping.get(), false, input->velocity,
                                    out->damping_force, d);
  if (status != StatusCode::kOk || request == DynamicEvaluationRequest::kForces)
    return status;
  status = out->inertia_acceleration->CopyFrom(mass, d);
  if (status != StatusCode::kOk) return status;
  Matrix* zero = nullptr;
  status = Matrix::Create(dof_count(), dof_count(), &zero, d);
  std::unique_ptr<Matrix> zero_owner(zero);
  if (status != StatusCode::kOk) return status;
  for (auto* derivative : {out->inertia_displacement, out->inertia_velocity,
                           out->damping_displacement}) {
    status = derivative->CopyFrom(zero, d);
    if (status != StatusCode::kOk) return status;
  }
  return out->damping_velocity->CopyFrom(damping.get(), d);
}
StatusCode SolidElement::LumpFromConsistentMass(const Matrix* consistent,
                                                Matrix* out,
                                                Diagnostic* d) const {
  if (!SupportsMass(MassMode::kRowSum))
    return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                         "Tetra10 row-sum mass unsupported", d);
  const int size = dof_count();
  if (!consistent || !out || consistent == out || consistent->Rows() != size ||
      consistent->Columns() != size || out->Rows() != size ||
      out->Columns() != size)
    return InvalidCall("Tetra lumped mass shape/alias", d);
  Vector* sums = nullptr;
  auto status = Vector::Create(size, &sums, d);
  std::unique_ptr<Vector> owner(sums);
  if (status != StatusCode::kOk) return status;
  status = consistent->RowSums(sums, d);
  if (status != StatusCode::kOk) return status;
  return Matrix::Diagonal(sums, out, d);
}
StatusCode SolidElement::UpdateResults(const Vector* u,
                                       ElementWorkspace* context,
                                       Diagnostic* d) {
  if (!u || u->Size() != dof_count())
    return InvalidCall("Tetra official displacement shape", d);
  const SolidReferenceData* data = nullptr;
  auto status = CheckRule(this, context, false, &data, d);
  if (status != StatusCode::kOk) return status;
  SolidGeometry geometry;
  SolidEvaluation work(&geometry);
  status = work.Initialize(this, data, context, false, d);
  if (status != StatusCode::kOk) return status;
  Vector* ordered = nullptr;
  status = Vector::Create(6, &ordered, d);
  std::unique_ptr<Vector> ordered_owner(ordered);
  if (status != StatusCode::kOk) return status;
  Matrix* conversion = nullptr;
  status = Matrix::Create(6, 6, &conversion, d);
  std::unique_ptr<Matrix> conversion_owner(conversion);
  if (status != StatusCode::kOk) return status;
  const int order[] = {0, 1, 2, 3, 5, 4};
  for (int i = 0; i < 6; ++i) {
    status = conversion->Set(i, order[i], i < 3 ? 1 : .5, d);
    if (status != StatusCode::kOk) return status;
  }
  for (int p = 0; p < data->point_count(); ++p) {
    if (p > 0) {
      status = work.Point(data, p, d);
      if (status != StatusCode::kOk) return status;
    }
    status = work.Force(u, work.Measure(data, p), d);
    if (status != StatusCode::kOk) return status;
    status = p == 0 ? internal_force_->CopyFrom(work.force(), d)
                    : internal_force_->AddScaled(1, work.force(), d);
    if (status != StatusCode::kOk) return status;
    status = ordered->GatherFrom(work.stress(), order, 6, 0, d);
    if (status != StatusCode::kOk) return status;
    status = ip_stress_->SetRow(p, ordered, d);
    if (status != StatusCode::kOk) return status;
    status =
        Matrix::MultiplyVector(conversion, false, work.strain(), ordered, d);
    if (status != StatusCode::kOk) return status;
    status = ip_strain_->SetRow(p, ordered, d);
    if (status != StatusCode::kOk) return status;
  }
  return ExtrapolateSolidTensors(data, ip_stress_, ip_strain_, nodal_stress_,
                                  nodal_strain_, d);
}

bool SolidElement::SupportsTrialEvaluation() const {
  return formulation_ == SolidFormulation::kTetra4Full ||
         formulation_ == SolidFormulation::kTetra10Full ||
         formulation_ == SolidFormulation::kHexa8Full ||
         formulation_ == SolidFormulation::kWedge6 ||
         formulation_ == SolidFormulation::kPyramid5Full;
}
StatusCode SolidElement::EvaluateTrial(
    const ElementTrialState* state, const ElementWorkspace* workspace,
    ElementEvaluationRequest request, Vector* force, Matrix* tangent,
    TrialEvaluationFailure* failure, Diagnostic* d) const {
  if (failure) *failure = TrialEvaluationFailure::kNone;
  if (!SupportsTrialEvaluation())
    return Element::EvaluateTrial(state, workspace, request, force, tangent,
                                   failure, d);
  const bool need_force = request == ElementEvaluationRequest::kInternalForce ||
                          request == ElementEvaluationRequest::kBoth;
  const bool need_tangent = request == ElementEvaluationRequest::kTangent ||
                            request == ElementEvaluationRequest::kBoth;
  const SolidReferenceData* data = nullptr;
  auto status = ValidateSolidTrialInput(this, state, workspace, &data, d);
  if (status != StatusCode::kOk) return status;
  if (!failure || (!need_force && !need_tangent) ||
      (need_force ? !force || force->Size() != dof_count()
                  : force != nullptr) ||
      (need_tangent ? !tangent || tangent->Rows() != dof_count() ||
                         tangent->Columns() != dof_count()
                    : tangent != nullptr) || force == state->translations)
    return InvalidCall("Solid nonlinear trial input/request/output", d);
  SolidTrialData work;
  status = work.Initialize(this, data, d);
  if (status != StatusCode::kOk) return status;
  return EvaluateSolidNonlinearPoints(formulation_, work.gradients(),
      work.weights(), state->translations, material_, request, nullptr, force,
      tangent,
      failure, d);
}

StatusCode SolidElement::UpdateTrialResults(
    const ElementTrialState* state, const ElementWorkspace* workspace,
    TrialEvaluationFailure* failure, Diagnostic* d) {
  if (failure) *failure = TrialEvaluationFailure::kNone;
  if (!SupportsTrialEvaluation())
    return Element::UpdateTrialResults(state, workspace, failure, d);
  const SolidReferenceData* data = nullptr;
  auto status = ValidateSolidTrialInput(this, state, workspace, &data, d);
  if (status != StatusCode::kOk) return status;
  if (!failure || state->translations == internal_force_)
    return InvalidCall("Solid nonlinear recovery input/failure", d);
  SolidTrialData work;
  status = work.Initialize(this, data, d);
  if (status != StatusCode::kOk) return status;
  return RecoverSolidNonlinearPoints(data, work.gradients(), work.weights(),
      state->translations, material_, internal_force_, ip_stress_, ip_strain_,
      nodal_stress_, nodal_strain_, failure, d);
}

bool SolidElement::SupportsTrialGeometryCheck() const {
  return formulation_ == SolidFormulation::kTetra4Full ||
         formulation_ == SolidFormulation::kTetra10Full ||
         formulation_ == SolidFormulation::kHexa8Full ||
         formulation_ == SolidFormulation::kWedge6 ||
         formulation_ == SolidFormulation::kPyramid5Full;
}
StatusCode SolidElement::CheckTrialGeometry(const ElementTrialState* state,
    const ElementWorkspace* workspace, TrialEvaluationFailure* failure,
    Diagnostic* d) const {
  if (!SupportsTrialGeometryCheck())
    return Element::CheckTrialGeometry(state, workspace, failure, d);
  const auto status = CheckSolidTrialGeometry(this, state, workspace, failure, d);
  if (status != StatusCode::kOk && d) d->element_id = Id();
  return status;
}

int SolidElement::integration_point_count() const { return point_count_; }
const Vector* SolidElement::internal_force() const { return internal_force_; }
const Matrix* SolidElement::ip_stress() const { return ip_stress_; }
const Matrix* SolidElement::ip_strain() const { return ip_strain_; }
const Matrix* SolidElement::nodal_stress() const { return nodal_stress_; }
const Matrix* SolidElement::nodal_strain() const { return nodal_strain_; }
const SolidSection* SolidElement::section() const { return section_; }
const ElasticMaterial* SolidElement::material() const { return material_; }
}  // namespace fesa
