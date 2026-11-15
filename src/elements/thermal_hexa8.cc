#include "fesa/elements/thermal_hexa8.h"

#include <cmath>
#include <memory>
#include <new>

#include "fesa/math/dense_factorization.h"
#include "elements/thermal_face_geometry.h"
#include "model/model_storage.h"

namespace fesa {
namespace {
constexpr int kSigns[8][3] = {{-1, -1, -1}, {1, -1, -1}, {1, 1, -1},
                              {-1, 1, -1},  {-1, -1, 1}, {1, -1, 1},
                              {1, 1, 1},    {-1, 1, 1}};

// Anchors geometry at node 1 through BLAS to remove the global origin.
// Caller owns success; failure cleans partial storage and nulls *out.
StatusCode Coordinates(const ThermalHexa8* element, Matrix** out,
                       Diagnostic* d) {
  *out = nullptr;
  Matrix* raw = nullptr;
  auto status = Matrix::Create(8, 3, &raw, d);
  std::unique_ptr<Matrix> coordinates(raw);
  if (status != StatusCode::kOk) return status;
  for (int i = 0; i < 8; ++i) {
    for (int axis = 0; axis < 3; ++axis) {
      double value = 0;
      status = element->GetNode(i)->Coordinate(axis, &value, d);
      if (status == StatusCode::kOk)
        status = coordinates->Set(i, axis, value, d);
      if (status != StatusCode::kOk) return status;
    }
  }
  status = Matrix::Create(8, 1, &raw, d);
  std::unique_ptr<Matrix> ones(raw);
  if (status != StatusCode::kOk) return status;
  for (int i = 0; i < 8; ++i) {
    status = ones->Set(i, 0, 1, d);
    if (status != StatusCode::kOk) return status;
  }
  status = Matrix::Create(1, 3, &raw, d);
  std::unique_ptr<Matrix> anchor(raw);
  if (status != StatusCode::kOk) return status;
  status = anchor->CopyBlockFrom(coordinates.get(), 0, 0, 1, 3, 0, 0, d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Create(8, 3, &raw, d);
  std::unique_ptr<Matrix> translation(raw);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Multiply(ones.get(), false, anchor.get(), false,
                            translation.get(), d);
  if (status == StatusCode::kOk)
    status = coordinates->AddScaled(-1, translation.get(), d);
  if (status != StatusCode::kOk) return status;
  *out = coordinates.release();
  return StatusCode::kOk;
}

// Formulation 3.1: J[alpha,i]=dx_i/dxi_alpha; J*grad_x(N)=grad_xi(N).
// Station p is nearest node p. Null gradient skips its solve. LU determinant
// and physical gradients use the common external-library math adapters.
StatusCode MapPoint(const Matrix* coordinates, int point, Matrix* gradient,
                    double* determinant, Diagnostic* d) {
  const double gauss = 1 / std::sqrt(3.0);
  const double xi = gauss * kSigns[point][0];
  const double eta = gauss * kSigns[point][1];
  const double zeta = gauss * kSigns[point][2];
  Matrix* raw = nullptr;
  auto status = Matrix::Create(3, 8, &raw, d);
  std::unique_ptr<Matrix> natural(raw);
  if (status != StatusCode::kOk) return status;
  for (int a = 0; a < 8; ++a) {
    const int sx = kSigns[a][0], sy = kSigns[a][1], sz = kSigns[a][2];
    const double values[3] = {.125 * sx * (1 + sy * eta) * (1 + sz * zeta),
                              .125 * sy * (1 + sx * xi) * (1 + sz * zeta),
                              .125 * sz * (1 + sx * xi) * (1 + sy * eta)};
    for (int axis = 0; axis < 3; ++axis) {
      status = natural->Set(axis, a, values[axis], d);
      if (status != StatusCode::kOk) return status;
    }
  }
  status = Matrix::Create(3, 3, &raw, d);
  std::unique_ptr<Matrix> jacobian(raw);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Multiply(natural.get(), false, coordinates, false,
                            jacobian.get(), d);
  if (status != StatusCode::kOk) return status;
  DenseLuFactorization* factor = nullptr;
  status = DenseLuFactorization::Create(jacobian.get(), &factor, d);
  std::unique_ptr<DenseLuFactorization> owner(factor);
  if (status == StatusCode::kSingular)
    return ReportFailure(StatusCode::kInputError, DiagnosticStage::kValidation,
                         "DC3D8 singular Jacobian", d);
  if (status != StatusCode::kOk) return status;
  double det = 0;
  status = factor->Determinant(&det, d);
  if (status != StatusCode::kOk) return status;
  if (det <= 0)
    return ReportFailure(StatusCode::kInputError, DiagnosticStage::kValidation,
                         "DC3D8 nonpositive Jacobian", d);
  if (gradient) status = factor->Solve(natural.get(), false, gradient, d);
  if (status == StatusCode::kOk) *determinant = det;
  return status;
}
}  // namespace

ThermalHexa8::ThermalHexa8(int id, const ThermalMaterial* material)
    : ThermalElement(id), material_(material) {}
StatusCode ThermalHexa8::Create(int id, const Node* const* nodes, int count,
                                const ThermalSolidSection* section,
                                ThermalHexa8** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("ThermalHexa8 output", d);
  *out = nullptr;
  if (!model_storage::ValidId(id) || !nodes || count != 8 || !section ||
      !section->material())
    return model_storage::Invalid("ThermalHexa8 definition", d);
  for (int i = 0; i < 8; ++i) {
    if (!nodes[i] || nodes[i]->temperature_dof_count() != 1)
      return model_storage::Invalid("ThermalHexa8 node layout", d);
    for (int j = 0; j < i; ++j)
      if (nodes[i]->Id() == nodes[j]->Id())
        return model_storage::Invalid("ThermalHexa8 repeated node", d);
  }
  std::unique_ptr<ThermalHexa8> element(
      new (std::nothrow) ThermalHexa8(id, section->material()));
  if (!element) return model_storage::AllocationFailure("ThermalHexa8", d);
  for (int i = 0; i < 8; ++i) element->nodes_[i] = nodes[i];
  const auto status = element->ValidateGeometry(d);
  if (status != StatusCode::kOk) return status;
  *out = element.release();
  return StatusCode::kOk;
}
ThermalElementType ThermalHexa8::type() const {
  return ThermalElementType::kDc3d8;
}
const char* ThermalHexa8::type_name() const { return "DC3D8"; }
int ThermalHexa8::node_count() const { return 8; }
int ThermalHexa8::dofs_per_node() const { return 1; }
const Node* ThermalHexa8::GetNode(int local) const {
  return local >= 0 && local < 8 ? nodes_[local] : nullptr;
}
const ThermalMaterial* ThermalHexa8::material() const { return material_; }
StatusCode ThermalHexa8::ValidateGeometry(Diagnostic* d) const {
  Matrix* raw = nullptr;
  auto status = Coordinates(this, &raw, d);
  std::unique_ptr<Matrix> coordinates(raw);
  for (int point = 0; point < 8 && status == StatusCode::kOk; ++point) {
    double determinant = 0;
    status = MapPoint(coordinates.get(), point, nullptr, &determinant, d);
  }
  if (status != StatusCode::kOk && d) d->element_id = Id();
  return status;
}
StatusCode ThermalHexa8::ComputeConduction(const Vector* theta, Matrix* out,
                                           Diagnostic* d) const {
  if (!out || out->Rows() != 8 || out->Columns() != 8 ||
      (theta && theta->Size() != 8))
    return model_storage::Invalid("ThermalHexa8 conduction sizes", d);
  if (theta)
    for (int i = 0; i < 8; ++i)
      if (!std::isfinite(theta->Data()[i]))
        return model_storage::Invalid("ThermalHexa8 finite temperatures", d);
  if (material_->IsConductivityTemperatureDependent() && !theta)
    return model_storage::Invalid("DC3D8 material requires temperatures", d);
  double k = 0;
  auto status = material_->Conductivity(0, &k, d);
  if (status != StatusCode::kOk) return status;
  Matrix* raw = nullptr;
  status = Coordinates(this, &raw, d);
  std::unique_ptr<Matrix> coordinates(raw);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Create(3, 8, &raw, d);
  std::unique_ptr<Matrix> gradient(raw);
  if (status != StatusCode::kOk) return status;
  Vector* vector = nullptr;
  std::unique_ptr<Vector> shape;
  if (material_->IsConductivityTemperatureDependent()) {
    status = Vector::Create(8, &vector, d);
    if (status != StatusCode::kOk) return status;
    shape.reset(vector);
  }
  for (int point = 0; point < 8; ++point) {
    double determinant = 0;
    status =
        MapPoint(coordinates.get(), point, gradient.get(), &determinant, d);
    if (status != StatusCode::kOk) return status;
    if (shape) {
      // Formulation 8.3: interpolate temperature, then evaluate k at Gauss g.
      const double gauss = 1 / std::sqrt(3.0);
      for (int a = 0; a < 8; ++a) {
        const double value = .125 *
                             (1 + kSigns[a][0] * kSigns[point][0] * gauss) *
                             (1 + kSigns[a][1] * kSigns[point][1] * gauss) *
                             (1 + kSigns[a][2] * kSigns[point][2] * gauss);
        status = shape->Set(a, value, d);
        if (status != StatusCode::kOk) return status;
      }
      double temperature = 0;
      status = shape->Dot(theta, &temperature, d);
      if (status == StatusCode::kOk)
        status = material_->Conductivity(temperature, &k, d);
      if (status != StatusCode::kOk) return status;
    }
    const double coefficient = k * determinant;
    if (!std::isfinite(coefficient) || coefficient <= 0)
      return ReportFailure(StatusCode::kLibraryFailure,
                           DiagnosticStage::kComputation,
                           "DC3D8 conduction coefficient range", d);
    status = Matrix::MultiplyAdd(gradient.get(), true, gradient.get(), false,
                                 coefficient, point == 0 ? 0 : 1, out, d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}
StatusCode ThermalHexa8::ComputeCapacityWeights(Vector* out,
                                                Diagnostic* d) const {
  if (!out || out->Size() != 8)
    return model_storage::Invalid("ThermalHexa8 capacity sizes", d);
  Matrix* raw = nullptr;
  auto status = Coordinates(this, &raw, d);
  std::unique_ptr<Matrix> coordinates(raw);
  if (status != StatusCode::kOk) return status;
  for (int point = 0; point < 8; ++point) {
    double determinant = 0;
    status = MapPoint(coordinates.get(), point, nullptr, &determinant, d);
    if (status == StatusCode::kOk) status = out->Set(point, determinant, d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}
StatusCode ThermalHexa8::ComputeSurfaceWeights(const char* face_label,
                                               Vector* out,
                                               Diagnostic* d) const {
  if (!face_label || face_label[0] != 'S' || face_label[1] < '1' ||
      face_label[1] > '6' || face_label[2] != '\0')
    return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                         "Unsupported ThermalHexa8 flux face", d);
  // Abaqus ELM-3D Table 4; zero-based local node indices.
  constexpr int faces[6][4] = {{0, 1, 2, 3}, {4, 7, 6, 5}, {0, 4, 5, 1},
                               {1, 5, 6, 2}, {2, 6, 7, 3}, {3, 7, 4, 0}};
  return thermal_face_geometry::AssignWeights(this, faces[face_label[1] - '1'],
                                              4, 0, out, d);
}
}  // namespace fesa
