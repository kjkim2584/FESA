#include "fesa/elements/thermal_shell4.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <limits>
#include <memory>
#include <new>

#include "elements/thermal_face_geometry.h"
#include "model/model_storage.h"

namespace fesa {
namespace {
StatusCode Coordinates(const ThermalShell4* element, Matrix** out,
                       double* scale, Diagnostic* d) {
  const Node* nodes[4] = {};
  for (int i = 0; i < 4; ++i) nodes[i] = element->GetNode(i);
  return thermal_face_geometry::CreateCoordinates(nodes, out, scale, d);
}
StatusCode MapPoint(const Matrix* xy, int point, Matrix* gradient,
                    double* shape, double* determinant, Diagnostic* d) {
  return thermal_face_geometry::MapPoint(xy, point, gradient, shape,
                                         determinant, d);
}
// Scalar coefficient evaluation only; all dense algebra stays in Matrix.
// Each call has at most six numerator factors and one nonzero denominator.
// Normalizing bounded mantissas avoids losing compensating physical scales.
// Coefficient underflow is allowed; the final K/volume range gates decide
// whether a physically required normal result was lost (formulation 3.3).
bool ScaledProduct(const double* factors, int count, double denominator,
                   double* out) {
  if (!std::isfinite(denominator) || denominator == 0) return false;
  double mantissa = 1;
  std::int64_t exponent = 0;
  for (int i = 0; i <= count; ++i) {
    const double factor = i < count ? factors[i] : denominator;
    if (!std::isfinite(factor)) return false;
    int factor_exponent = 0, shift = 0;
    const double part = std::frexp(factor, &factor_exponent);
    if (i < count) {
      mantissa *= part;
      exponent += factor_exponent;
    } else {
      mantissa /= part;
      exponent -= factor_exponent;
    }
    mantissa = std::frexp(mantissa, &shift);
    exponent += shift;
  }
  // Seven finite double factors bound the exponent well inside int range.
  const double value = std::ldexp(mantissa, static_cast<int>(exponent));
  if (!std::isfinite(value)) return false;
  *out = value;
  return true;
}
double SimpsonRatio(int i, int ns) {
  return ns == 1
             ? 1
             : (i == 0 || i == ns - 1 ? 1 : (i % 2 ? 4 : 2)) / (3.0 * (ns - 1));
}
}  // namespace
ThermalShell4::ThermalShell4(int id, const ThermalShellSection* section)
    : ThermalElement(id), section_(section) {}
StatusCode ThermalShell4::Create(int id, const Node* const* nodes, int count,
                                 const ThermalShellSection* section,
                                 ThermalShell4** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("DS4 output", d);
  *out = nullptr;
  if (!model_storage::ValidId(id) || !nodes || count != 4 || !section)
    return model_storage::Invalid("DS4 definition", d);
  for (int a = 0; a < 4; ++a) {
    if (!nodes[a] ||
        nodes[a]->temperature_dof_count() != section->section_points())
      return model_storage::Invalid("DS4 node layout", d);
    for (int b = 0; b < a; ++b)
      if (nodes[a]->Id() == nodes[b]->Id())
        return model_storage::Invalid("DS4 repeated node", d);
  }
  std::unique_ptr<ThermalShell4> element(new (std::nothrow)
                                             ThermalShell4(id, section));
  if (!element) return model_storage::AllocationFailure("DS4", d);
  for (int a = 0; a < 4; ++a) element->nodes_[a] = nodes[a];
  const auto status = element->ValidateGeometry(d);
  if (status == StatusCode::kOk) *out = element.release();
  return status;
}
ThermalElementType ThermalShell4::type() const {
  return ThermalElementType::kDs4;
}
const char* ThermalShell4::type_name() const { return "DS4"; }
int ThermalShell4::node_count() const { return 4; }
int ThermalShell4::dofs_per_node() const { return section_->section_points(); }
const Node* ThermalShell4::GetNode(int i) const {
  return i >= 0 && i < 4 ? nodes_[i] : nullptr;
}
const ThermalMaterial* ThermalShell4::material() const {
  return section_->material();
}
StatusCode ThermalShell4::GetShellThickness(double* out, Diagnostic* d) const {
  if (!out) return model_storage::Invalid("DS4 thickness output", d);
  *out = section_->thickness();
  return StatusCode::kOk;
}
StatusCode ThermalShell4::ValidateGeometry(Diagnostic* d) const {
  Matrix* raw = nullptr;
  double length = 0;
  auto status = Coordinates(this, &raw, &length, d);
  std::unique_ptr<Matrix> xy(raw);
  for (int g = 0; g < 4 && status == StatusCode::kOk; ++g) {
    double det = 0;
    status = MapPoint(raw, g, nullptr, nullptr, &det, d);
  }
  if (status != StatusCode::kOk && d) d->element_id = Id();
  return status;
}
StatusCode ThermalShell4::ComputeConduction(const Vector* theta, Matrix* out,
                                            Diagnostic* d) const {
  const int ns = dofs_per_node(), n = 4 * ns;
  if (!out || out->Rows() != n || out->Columns() != n ||
      (theta && theta->Size() != n))
    return model_storage::Invalid("DS4 conduction sizes", d);
  if (theta)
    for (int i = 0; i < n; ++i)
      if (!std::isfinite(theta->Data()[i]))
        return model_storage::Invalid("DS4 finite temperatures", d);
  const bool variable = material()->IsConductivityTemperatureDependent();
  if (variable && !theta)
    return model_storage::Invalid("DS4 material requires temperatures", d);
  double k = 0;
  auto status = material()->Conductivity(0, &k, d);
  if (status != StatusCode::kOk) return status;
  Matrix* raw = nullptr;
  double length = 0;
  status = Coordinates(this, &raw, &length, d);
  std::unique_ptr<Matrix> xy(raw);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Create(2, 4, &raw, d);
  std::unique_ptr<Matrix> gradient(raw);
  if (status != StatusCode::kOk) return status;
  const double h = section_->thickness();
  for (int g = 0; g < 4; ++g) {
    double shape[4] = {}, det = 0;
    status = MapPoint(xy.get(), g, gradient.get(), shape, &det, d);
    if (status != StatusCode::kOk) return status;
    const int thickness_stations = variable ? 3 : 2;
    const int rows = 2 * ns + thickness_stations * ((ns - 1) / 2);
    status = Matrix::Create(rows, n, &raw, d);
    std::unique_ptr<Matrix> b(raw);
    if (status != StatusCode::kOk) return status;
    double station_k[19] = {};
    for (int i = 0; i < ns; ++i) {
      station_k[i] = k;
      if (variable) {
        // Task16: natural-sign pairing g -> node g, not physical nearest.
        status = material()->Conductivity(theta->Data()[g * ns + i],
                                          &station_k[i], d);
        if (status != StatusCode::kOk) return status;
      }
    }
    for (int i = 0; i < ns; ++i)
      for (int a = 0; a < 4; ++a)
        for (int axis = 0; axis < 2; ++axis) {
          const double factors[5] = {std::sqrt(det),
                                     gradient->Data()[axis * 4 + a],
                                     std::sqrt(station_k[i]), std::sqrt(h),
                                     std::sqrt(SimpsonRatio(i, ns))};
          double coefficient = 0;
          if (!ScaledProduct(factors, 5, 1, &coefficient))
            return ReportFailure(
                StatusCode::kLibraryFailure, DiagnosticStage::kComputation,
                "DS4 coefficient outside representable range", d);
          status = b->Set(2 * i + axis, a * ns + i, coefficient, d);
          if (status != StatusCode::kOk) return status;
        }
    if (ns > 1) {
      for (int m = 0; m < (ns - 1) / 2; ++m)
        for (int q = 0; q < thickness_stations; ++q) {
          const double xi =
              variable ? q - 1.0 : (q == 0 ? -1 : 1) / std::sqrt(3.0);
          const double weight = variable ? (q == 1 ? 4. / 3 : 1. / 3) : 1.;
          const double root_k = std::sqrt(variable ? station_k[2 * m + q] : k);
          const double derivative[3] = {xi - .5, -2 * xi, xi + .5};
          for (int a = 0; a < 4; ++a)
            for (int r = 0; r < 3; ++r) {
              const double factors[6] = {
                  std::sqrt(det * weight),
                  length,
                  shape[a],
                  derivative[r],
                  root_k,
                  std::sqrt(static_cast<double>(ns - 1))};
              double coefficient = 0;
              if (!ScaledProduct(factors, 6, std::sqrt(h), &coefficient))
                return ReportFailure(
                    StatusCode::kLibraryFailure, DiagnosticStage::kComputation,
                    "DS4 coefficient outside representable range", d);
              status = b->Set(2 * ns + thickness_stations * m + q,
                              a * ns + 2 * m + r, coefficient, d);
              if (status != StatusCode::kOk) return status;
            }
        }
    }
    status = Matrix::MultiplyAdd(b.get(), true, b.get(), false, 1,
                                 g == 0 ? 0 : 1, out, d);
    if (status != StatusCode::kOk) return status;
  }
  for (int i = 0; i < n * n; ++i)
    if (!std::isfinite(out->Data()[i]))
      return ReportFailure(StatusCode::kLibraryFailure,
                           DiagnosticStage::kComputation,
                           "DS4 conduction outside representable range", d);
  // Diagonals are strictly positive. Subnormal diagonals cannot retain the
  // normal-double relative precision required by the DS4 range contract.
  for (int i = 0; i < n; ++i)
    if (!(out->Data()[i * n + i] > 0) || !std::isnormal(out->Data()[i * n + i]))
      return ReportFailure(StatusCode::kLibraryFailure,
                           DiagnosticStage::kComputation,
                           "DS4 conduction outside representable range", d);
  return StatusCode::kOk;
}
StatusCode ThermalShell4::ComputeCapacityWeights(Vector* out,
                                                 Diagnostic* d) const {
  const int ns = dofs_per_node(), n = 4 * ns;
  if (!out || out->Size() != n)
    return model_storage::Invalid("DS4 capacity sizes", d);
  Matrix* raw = nullptr;
  double length = 0;
  auto status = Coordinates(this, &raw, &length, d);
  std::unique_ptr<Matrix> xy(raw);
  if (status != StatusCode::kOk) return status;
  Vector* vector = nullptr;
  status = Vector::Create(n, &vector, d);
  std::unique_ptr<Vector> station(vector);
  if (status != StatusCode::kOk) return status;
  for (int g = 0; g < 4; ++g) {
    double shape[4] = {}, det = 0;
    status = MapPoint(xy.get(), g, nullptr, shape, &det, d);
    if (status != StatusCode::kOk) return status;
    for (int a = 0; a < 4; ++a)
      for (int i = 0; i < ns; ++i) {
        const double factors[6] = {det,
                                   length,
                                   length,
                                   shape[a],
                                   section_->thickness(),
                                   SimpsonRatio(i, ns)};
        double weight = 0;
        if (!ScaledProduct(factors, 6, 1, &weight) || !std::isnormal(weight) ||
            !(weight > 0))
          return ReportFailure(StatusCode::kLibraryFailure,
                               DiagnosticStage::kComputation,
                               "DS4 capacity outside representable range", d);
        status = station->Set(a * ns + i, weight, d);
        if (status != StatusCode::kOk) return status;
      }
    status = g == 0 ? out->CopyFrom(station.get(), d)
                    : out->AddScaled(1, station.get(), d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}
StatusCode ThermalShell4::ComputeSurfaceWeights(const char* face_label,
                                                Vector* out,
                                                Diagnostic* d) const {
  if (!face_label || (std::strcmp(face_label, "SNEG") != 0 &&
                      std::strcmp(face_label, "SPOS") != 0))
    return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                         "Unsupported DS4 flux face", d);
  const int nodes[4] = {0, 1, 2, 3};
  const int point =
      std::strcmp(face_label, "SNEG") == 0 ? 0 : dofs_per_node() - 1;
  return thermal_face_geometry::AssignWeights(this, nodes, 4, point, out, d);
}
}  // namespace fesa
