#include "fesa/elements/thermal_tetra4.h"

#include <cmath>
#include <memory>
#include <new>

#include "fesa/math/dense_factorization.h"
#include "elements/thermal_face_geometry.h"
#include "model/model_storage.h"

namespace fesa {
ThermalTetra4::ThermalTetra4(int id, const ThermalMaterial* material)
    : ThermalElement(id), material_(material) {}
StatusCode ThermalTetra4::Create(int id, const Node* const* nodes, int count,
                                 const ThermalSolidSection* section,
                                 ThermalTetra4** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("ThermalTetra4 output", d);
  *out = nullptr;
  if (!model_storage::ValidId(id) || !nodes || count != 4 || !section ||
      !section->material())
    return model_storage::Invalid("ThermalTetra4 definition", d);
  for (int i = 0; i < 4; ++i) {
    if (!nodes[i] || nodes[i]->temperature_dof_count() != 1)
      return model_storage::Invalid("ThermalTetra4 node layout", d);
    for (int j = 0; j < i; ++j)
      if (nodes[i]->Id() == nodes[j]->Id())
        return model_storage::Invalid("ThermalTetra4 repeated node", d);
  }
  std::unique_ptr<ThermalTetra4> element(
      new (std::nothrow) ThermalTetra4(id, section->material()));
  if (!element) return model_storage::AllocationFailure("ThermalTetra4", d);
  for (int i = 0; i < 4; ++i) element->nodes_[i] = nodes[i];
  const auto status = element->ValidateGeometry(d);
  if (status != StatusCode::kOk) return status;
  *out = element.release();
  return StatusCode::kOk;
}
ThermalElementType ThermalTetra4::type() const {
  return ThermalElementType::kDc3d4;
}
const char* ThermalTetra4::type_name() const { return "DC3D4"; }
int ThermalTetra4::node_count() const { return 4; }
int ThermalTetra4::dofs_per_node() const { return 1; }
const Node* ThermalTetra4::GetNode(int local) const {
  return local >= 0 && local < 4 ? nodes_[local] : nullptr;
}
const ThermalMaterial* ThermalTetra4::material() const { return material_; }
StatusCode ThermalTetra4::Geometry(Matrix* gradient, double* volume,
                                   Diagnostic* d) const {
  Matrix* raw = nullptr;
  auto status = Matrix::Create(4, 3, &raw, d);
  std::unique_ptr<Matrix> coordinates(raw);
  if (status != StatusCode::kOk) return status;
  for (int i = 0; i < 4; ++i) {
    for (int axis = 0; axis < 3; ++axis) {
      double value = 0;
      status = nodes_[i]->Coordinate(axis, &value, d);
      if (status == StatusCode::kOk)
        status = coordinates->Set(i, axis, value, d);
      if (status != StatusCode::kOk) return status;
    }
  }
  status = Matrix::Create(3, 4, &raw, d);
  std::unique_ptr<Matrix> natural(raw);
  if (status != StatusCode::kOk) return status;
  // Formulation 3.2, DN rows: derivatives of (1-r-s-t,r,s,t).
  for (int axis = 0; axis < 3; ++axis) {
    status = natural->Set(axis, 0, -1, d);
    if (status == StatusCode::kOk) status = natural->Set(axis, axis + 1, 1, d);
    if (status != StatusCode::kOk) return status;
  }
  status = Matrix::Create(3, 3, &raw, d);
  std::unique_ptr<Matrix> jacobian(raw);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Multiply(natural.get(), false, coordinates.get(), false,
                            jacobian.get(), d);
  if (status != StatusCode::kOk) return status;
  DenseLuFactorization* factor = nullptr;
  status = DenseLuFactorization::Create(jacobian.get(), &factor, d);
  std::unique_ptr<DenseLuFactorization> owner(factor);
  if (status == StatusCode::kSingular)
    return ReportFailure(StatusCode::kInputError, DiagnosticStage::kValidation,
                         "DC3D4 singular Jacobian", d);
  if (status != StatusCode::kOk) return status;
  double determinant = 0;
  status = factor->Determinant(&determinant, d);
  if (status != StatusCode::kOk) return status;
  if (determinant <= 0)
    return ReportFailure(StatusCode::kInputError, DiagnosticStage::kValidation,
                         "DC3D4 nonpositive volume", d);
  const double v = determinant / 6;
  if (!(v > 0) || !std::isfinite(v))
    return ReportFailure(StatusCode::kLibraryFailure,
                         DiagnosticStage::kComputation, "DC3D4 volume range",
                         d);
  if (gradient) status = factor->Solve(natural.get(), false, gradient, d);
  if (status == StatusCode::kOk) *volume = v;
  return status;
}
StatusCode ThermalTetra4::ValidateGeometry(Diagnostic* d) const {
  double volume = 0;
  const auto status = Geometry(nullptr, &volume, d);
  if (status != StatusCode::kOk && d) d->element_id = Id();
  return status;
}
StatusCode ThermalTetra4::ComputeConduction(const Vector* theta, Matrix* out,
                                            Diagnostic* d) const {
  if (!out || out->Rows() != 4 || out->Columns() != 4 ||
      (theta && theta->Size() != 4))
    return model_storage::Invalid("ThermalTetra4 conduction sizes", d);
  if (theta)
    for (int i = 0; i < 4; ++i)
      if (!std::isfinite(theta->Data()[i]))
        return model_storage::Invalid("ThermalTetra4 finite temperatures", d);
  if (material_->IsConductivityTemperatureDependent() && !theta)
    return model_storage::Invalid("DC3D4 material requires temperatures", d);
  double k = 0;
  auto status = material_->Conductivity(0, &k, d);
  if (status != StatusCode::kOk) return status;
  if (material_->IsConductivityTemperatureDependent()) {
    Vector* raw_values = nullptr;
    status = Vector::Create(4, &raw_values, d);
    std::unique_ptr<Vector> values(raw_values);
    if (status != StatusCode::kOk) return status;
    // Formulation 8.3: average nodal conductivities, not nodal temperatures.
    for (int a = 0; a < 4; ++a) {
      double value = 0;
      status = material_->Conductivity(theta->Data()[a], &value, d);
      if (status == StatusCode::kOk) status = values->Set(a, value, d);
      if (status != StatusCode::kOk) return status;
    }
    status = values->AbsoluteSum(&k, d);
    if (status != StatusCode::kOk) return status;
    k /= 4;
  }
  Matrix* raw = nullptr;
  status = Matrix::Create(3, 4, &raw, d);
  std::unique_ptr<Matrix> gradient(raw);
  if (status != StatusCode::kOk) return status;
  double volume = 0;
  status = Geometry(gradient.get(), &volume, d);
  if (status != StatusCode::kOk) return status;
  const double coefficient = k * volume;
  if (!std::isfinite(coefficient) || coefficient <= 0)
    return ReportFailure(StatusCode::kLibraryFailure,
                         DiagnosticStage::kComputation,
                         "DC3D4 conduction coefficient range", d);
  return Matrix::MultiplyAdd(gradient.get(), true, gradient.get(), false,
                             coefficient, 0, out, d);
}
StatusCode ThermalTetra4::ComputeCapacityWeights(Vector* out,
                                                 Diagnostic* d) const {
  if (!out || out->Size() != 4)
    return model_storage::Invalid("ThermalTetra4 capacity sizes", d);
  double volume = 0;
  auto status = Geometry(nullptr, &volume, d);
  if (status != StatusCode::kOk) return status;
  const double weight = volume / 4;
  if (!(weight > 0))
    return ReportFailure(StatusCode::kLibraryFailure,
                         DiagnosticStage::kComputation,
                         "DC3D4 capacity coefficient range", d);
  for (int i = 0; i < 4; ++i) {
    status = out->Set(i, weight, d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}
StatusCode ThermalTetra4::ComputeSurfaceWeights(const char* face_label,
                                                Vector* out,
                                                Diagnostic* d) const {
  if (!face_label || face_label[0] != 'S' || face_label[1] < '1' ||
      face_label[1] > '4' || face_label[2] != '\0')
    return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                         "Unsupported ThermalTetra4 flux face", d);
  // Abaqus ELM-3D Table 1; zero-based local node indices.
  constexpr int faces[4][3] = {{0, 1, 2}, {0, 3, 1}, {1, 3, 2}, {2, 3, 0}};
  return thermal_face_geometry::AssignWeights(this, faces[face_label[1] - '1'],
                                              3, 0, out, d);
}
}  // namespace fesa
