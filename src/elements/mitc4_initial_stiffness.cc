#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>

#include "elements/mitc4_geometry.h"
#include "elements/mitc4_initial.h"
#include "elements/mitc4_linear_ops.h"
#include "fesa/math/dense_factorization.h"
#include "model/model_storage.h"

namespace fesa {
namespace {

StatusCode NumericalFailure(const char* operation, Diagnostic* d) {
  return ReportFailure(StatusCode::kLibraryFailure,
                       DiagnosticStage::kComputation, operation, d);
}

bool Shape(const Matrix* matrix, int rows, int columns) {
  return matrix && matrix->Rows() == rows && matrix->Columns() == columns;
}

bool Finite(const Matrix* matrix) {
  if (!matrix) return false;
  for (int i = 0; i < matrix->Rows() * matrix->Columns(); ++i)
    if (!std::isfinite(matrix->Data()[i])) return false;
  return true;
}

// Scalar shape-function coefficients follow the node order in the initial
// kernel contract. Derivative row order is r,s.
StatusCode ShapeDerivatives(double r, double s, Matrix* derivatives,
                            Diagnostic* d) {
  constexpr double kR[4] = {-1.0, 1.0, 1.0, -1.0};
  constexpr double kS[4] = {-1.0, -1.0, 1.0, 1.0};
  for (int i = 0; i < 4; ++i) {
    auto status = derivatives->Set(0, i, kR[i] * (1.0 + kS[i] * s) / 4.0, d);
    if (status != StatusCode::kOk) return status;
    status = derivatives->Set(1, i, kS[i] * (1.0 + kR[i] * r) / 4.0, d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}

// Each row is a natural engineering shear evaluated at its own tying point.
// No Cartesian shear averaging or reduced quadrature is used.
StatusCode TyingRows(const Matrix* coordinates, Matrix* rows, Diagnostic* d) {
  Matrix* raw = nullptr;
  auto status = Matrix::Create(2, 4, &raw, d);
  std::unique_ptr<Matrix> derivatives(raw);
  if (status != StatusCode::kOk) return status;
  raw = nullptr;
  status = Matrix::Create(2, 3, &raw, d);
  std::unique_ptr<Matrix> tangents(raw);
  if (status != StatusCode::kOk) return status;
  constexpr double kR[4] = {0.0, 0.0, -1.0, 1.0};
  constexpr double kS[4] = {-1.0, 1.0, 0.0, 0.0};
  constexpr double kNodeR[4] = {-1.0, 1.0, 1.0, -1.0};
  constexpr double kNodeS[4] = {-1.0, -1.0, 1.0, 1.0};
  for (int point = 0; point < 4; ++point) {
    status = ShapeDerivatives(kR[point], kS[point], derivatives.get(), d);
    if (status != StatusCode::kOk) return status;
    status = Matrix::Multiply(derivatives.get(), false, coordinates, false,
                              tangents.get(), d);
    if (status != StatusCode::kOk) return status;
    const int direction = point / 2;
    const double x_derivative = tangents->Data()[3 * direction];
    const double y_derivative = tangents->Data()[3 * direction + 1];
    for (int i = 0; i < 4; ++i) {
      const double shape =
          (1.0 + kNodeR[i] * kR[point]) * (1.0 + kNodeS[i] * kS[point]) / 4.0;
      status = rows->Set(point, 6 * i + 2,
                         derivatives->Data()[4 * direction + i], d);
      if (status != StatusCode::kOk) return status;
      status = rows->Set(point, 6 * i + 3, -shape * y_derivative, d);
      if (status != StatusCode::kOk) return status;
      status = rows->Set(point, 6 * i + 4, shape * x_derivative, d);
      if (status != StatusCode::kOk) return status;
    }
  }
  return StatusCode::kOk;
}

StatusCode Constitutive(const ElasticMaterial* material, Matrix* out,
                        Diagnostic* d) {
  const double young = material->young();
  const double poisson = material->poisson();
  const double normal = young / (1.0 - poisson * poisson);
  const double shear = young / (2.0 * (1.0 + poisson));
  const double diagonal[5] = {normal, normal, shear, (5.0 / 6.0) * shear,
                              (5.0 / 6.0) * shear};
  for (int i = 0; i < 5; ++i) {
    if (!std::isfinite(diagonal[i]) || diagonal[i] <= 0.0)
      return NumericalFailure("MITC4 initial constitutive range", d);
    const auto status = out->Set(i, i, diagonal[i], d);
    if (status != StatusCode::kOk) return status;
  }
  auto status = out->Set(0, 1, poisson * normal, d);
  if (status != StatusCode::kOk) return status;
  return out->Set(1, 0, poisson * normal, d);
}

// At initial director n=(0,0,1), delta d=(theta_y,-theta_x,0).
// Entries below are scalar B coefficients, not a hand-written matrix product.
StatusCode MembraneB(const Matrix* gradients, double z, Matrix* out,
                     Diagnostic* d) {
  for (int i = 0; i < 4; ++i) {
    const double dx = gradients->Data()[i];
    const double dy = gradients->Data()[4 + i];
    const int row[8] = {0, 0, 1, 1, 2, 2, 2, 2};
    const int column[8] = {0, 4, 1, 3, 0, 1, 3, 4};
    const double value[8] = {dx, z * dx, dy, -z * dy, dy, dx, -z * dx, z * dy};
    for (int j = 0; j < 8; ++j) {
      const auto status = out->Set(row[j], 6 * i + column[j], value[j], d);
      if (status != StatusCode::kOk) return status;
    }
  }
  return StatusCode::kOk;
}

// The two transverse shears are independent of z in the initial flat limit.
// Integrate both z points at this in-plane Gauss location.
StatusCode IntegratePoint(const Matrix* coordinates, const Matrix* ties,
                          const Matrix* constitutive, double r, double s,
                          double thickness, Matrix* stiffness,
                          const Vector* displacement, Vector* force,
                          Diagnostic* d) {
  std::unique_ptr<Matrix> owners[3];
  constexpr int kRows[3] = {2, 2, 5};
  constexpr int kColumns[3] = {4, 24, 24};
  for (int i = 0; i < 3; ++i) {
    Matrix* raw = nullptr;
    const auto status = Matrix::Create(kRows[i], kColumns[i], &raw, d);
    owners[i].reset(raw);
    if (status != StatusCode::kOk) return status;
  }
  Matrix* gradients = owners[0].get();
  Matrix* cartesian_shears = owners[1].get();
  Matrix* strain = owners[2].get();
  // D*B is needed only for K. A force request owns only the three vectors
  // strain(5), stress(5), and the point's nodal contribution(24).
  std::unique_ptr<Matrix> stress;
  std::unique_ptr<Vector> vectors[3];
  if (stiffness) {
    Matrix* raw = nullptr;
    const auto code = Matrix::Create(5, 24, &raw, d);
    stress.reset(raw);
    if (code != StatusCode::kOk) return code;
  } else {
    constexpr int kSizes[3] = {5, 5, 24};
    for (int i = 0; i < 3; ++i) {
      Vector* raw = nullptr;
      const auto code = Vector::Create(kSizes[i], &raw, d);
      vectors[i].reset(raw);
      if (code != StatusCode::kOk) return code;
    }
  }
  double determinant = 0.0;
  auto status = ComputeMitc4LinearPointOperators(
      coordinates, ties, r, s, gradients, cartesian_shears, &determinant, d);
  if (status != StatusCode::kOk) return status;
  const double weight = determinant * (thickness / 2.0);
  if (!std::isfinite(weight) || weight <= 0.0)
    return NumericalFailure("MITC4 initial integration weight range", d);
  const double gauss = 1.0 / std::sqrt(3.0);
  for (double zeta : {-gauss, gauss}) {
    const double z = (thickness / 2.0) * zeta;
    if (!std::isfinite(z) || z == 0.0)
      return NumericalFailure("MITC4 initial thickness coordinate range", d);
    status = FillMitc4LinearB(gradients, cartesian_shears, z, strain, d);
    if (status != StatusCode::kOk) return status;
    if (stiffness) {
      status =
          Matrix::Multiply(constitutive, false, strain, false, stress.get(), d);
      if (status != StatusCode::kOk) return status;
      status = Matrix::MultiplyAdd(strain, true, stress.get(), false, weight,
                                   1.0, stiffness, d);
    } else {
      status = Matrix::MultiplyVector(strain, false, displacement,
                                      vectors[0].get(), d);
      if (status != StatusCode::kOk) return status;
      status = Matrix::MultiplyVector(constitutive, false, vectors[0].get(),
                                      vectors[1].get(), d);
      if (status != StatusCode::kOk) return status;
      status = Matrix::MultiplyVector(strain, true, vectors[1].get(),
                                      vectors[2].get(), d);
      if (status != StatusCode::kOk) return status;
      status = force->AddScaled(weight, vectors[2].get(), d);
    }
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}

StatusCode AddDrilling(Matrix* local, double* drilling, Diagnostic* d) {
  double minimum = std::numeric_limits<double>::max();
  bool found = false;
  for (int i = 0; i < 4; ++i)
    for (int component : {3, 4}) {
      const int dof = 6 * i + component;
      const double diagonal = local->Data()[24 * dof + dof];
      if (!std::isfinite(diagonal))
        return NumericalFailure("MITC4 initial rotation diagonal range", d);
      if (diagonal > 0.0) {
        found = true;
        minimum = std::min(minimum, diagonal);
      }
    }
  const double coefficient = minimum * 1.0e-3;
  if (!found || !std::isfinite(coefficient) || coefficient <= 0.0)
    return NumericalFailure("MITC4 initial drilling coefficient range", d);
  for (int i = 0; i < 4; ++i) {
    const auto status = local->Set(6 * i + 5, 6 * i + 5, coefficient, d);
    if (status != StatusCode::kOk) return status;
  }
  *drilling = coefficient;
  return StatusCode::kOk;
}

StatusCode BuildTransform(const Matrix* frame, Matrix* transform,
                          Diagnostic* d) {
  for (int block = 0; block < 8; ++block) {
    const auto status =
        transform->CopyBlockFrom(frame, 0, 0, 3, 3, 3 * block, 3 * block, d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}

}  // namespace

StatusCode BuildMitc4LinearTying(const Matrix* coordinates, Matrix* out,
                                 Diagnostic* d) {
  if (!Shape(coordinates, 4, 3) || !Shape(out, 4, 24) || coordinates == out ||
      !Finite(coordinates))
    return model_storage::Invalid("MITC4 tying arguments", d);
  Matrix* raw = nullptr;
  auto status = Matrix::Create(4, 24, &raw, d);
  std::unique_ptr<Matrix> result(raw);
  if (status != StatusCode::kOk) return status;
  status = TyingRows(coordinates, result.get(), d);
  if (status != StatusCode::kOk) return status;
  return out->CopyFrom(result.get(), d);
}

StatusCode BuildMitc4LinearConstitutive(const ElasticMaterial* material,
                                        Matrix* out, Diagnostic* d) {
  if (!material || !Shape(out, 5, 5))
    return model_storage::Invalid("MITC4 constitutive arguments", d);
  Matrix* raw = nullptr;
  auto status = Matrix::Create(5, 5, &raw, d);
  std::unique_ptr<Matrix> result(raw);
  if (status != StatusCode::kOk) return status;
  status = Constitutive(material, result.get(), d);
  if (status != StatusCode::kOk) return status;
  return out->CopyFrom(result.get(), d);
}

StatusCode ComputeMitc4LinearPointOperators(const Matrix* coordinates,
                                            const Matrix* ties, double r,
                                            double s, Matrix* gradients,
                                            Matrix* shears, double* determinant,
                                            Diagnostic* d) {
  if (!Shape(coordinates, 4, 3) || !Shape(ties, 4, 24) ||
      !Shape(gradients, 2, 4) || !Shape(shears, 2, 24) || !determinant ||
      !std::isfinite(r) || !std::isfinite(s) || r < -1.0 || r > 1.0 ||
      s < -1.0 || s > 1.0 || !Finite(coordinates) || !Finite(ties))
    return model_storage::Invalid("MITC4 point operator arguments", d);
  const Matrix* matrices[4] = {coordinates, ties, gradients, shears};
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < i; ++j)
      if (matrices[i] == matrices[j])
        return model_storage::Invalid("MITC4 point operator aliases", d);
  std::unique_ptr<Matrix> owners[5];
  constexpr int kRows[5] = {2, 2, 2, 2, 2};
  constexpr int kColumns[5] = {4, 3, 2, 4, 24};
  for (int i = 0; i < 5; ++i) {
    Matrix* raw = nullptr;
    const auto status = Matrix::Create(kRows[i], kColumns[i], &raw, d);
    owners[i].reset(raw);
    if (status != StatusCode::kOk) return status;
  }
  Matrix* derivatives = owners[0].get();
  Matrix* tangents = owners[1].get();
  Matrix* jacobian = owners[2].get();
  Matrix* interpolation = owners[3].get();
  Matrix* natural_shears = owners[4].get();
  auto status = ShapeDerivatives(r, s, derivatives, d);
  if (status != StatusCode::kOk) return status;
  status =
      Matrix::Multiply(derivatives, false, coordinates, false, tangents, d);
  if (status != StatusCode::kOk) return status;
  status = jacobian->CopyBlockFrom(tangents, 0, 0, 2, 2, 0, 0, d);
  if (status != StatusCode::kOk) return status;
  DenseLuFactorization* raw_factor = nullptr;
  status = DenseLuFactorization::Create(jacobian, &raw_factor, d);
  std::unique_ptr<DenseLuFactorization> factor(raw_factor);
  if (status != StatusCode::kOk) return status;
  double value = 0.0;
  status = factor->Determinant(&value, d);
  if (status != StatusCode::kOk) return status;
  if (!std::isfinite(value) || value <= 0.0)
    return NumericalFailure("MITC4 initial positive Jacobian", d);
  status = factor->Solve(derivatives, false, gradients, d);
  if (status != StatusCode::kOk) return status;
  const double weights[4] = {(1.0 - s) / 2.0, (1.0 + s) / 2.0, (1.0 - r) / 2.0,
                             (1.0 + r) / 2.0};
  for (int i = 0; i < 4; ++i) {
    status = interpolation->Set(i / 2, i, weights[i], d);
    if (status != StatusCode::kOk) return status;
  }
  status =
      Matrix::Multiply(interpolation, false, ties, false, natural_shears, d);
  if (status != StatusCode::kOk) return status;
  status = factor->Solve(natural_shears, false, shears, d);
  if (status != StatusCode::kOk) return status;
  *determinant = value;
  return StatusCode::kOk;
}

StatusCode FillMitc4LinearB(const Matrix* gradients, const Matrix* shears,
                            double z, Matrix* out, Diagnostic* d) {
  if (!Shape(gradients, 2, 4) || !Shape(shears, 2, 24) || !Shape(out, 5, 24) ||
      gradients == shears || gradients == out || shears == out ||
      !std::isfinite(z) || !Finite(gradients) || !Finite(shears))
    return model_storage::Invalid("MITC4 strain operator arguments", d);
  Matrix* raw = nullptr;
  auto status = Matrix::Create(5, 24, &raw, d);
  std::unique_ptr<Matrix> result(raw);
  if (status != StatusCode::kOk) return status;
  status = MembraneB(gradients, z, result.get(), d);
  if (status != StatusCode::kOk) return status;
  // Engineering order is xx,yy,xy,yz,xz; point operators return xz,yz.
  status = result->CopyBlockFrom(shears, 1, 0, 1, 24, 3, 0, d);
  if (status != StatusCode::kOk) return status;
  status = result->CopyBlockFrom(shears, 0, 0, 1, 24, 4, 0, d);
  if (status != StatusCode::kOk) return status;
  return out->CopyFrom(result.get(), d);
}

namespace {

StatusCode InitialStiffness(const Node* const* nodes, int count,
                            const ShellSection* section, Matrix* stiffness,
                            Matrix* axes, double* drilling, Diagnostic* d,
                            bool elastic_only = false) {
  if (!nodes || count != 4 || !section || !section->material() || !stiffness ||
      !drilling || stiffness->Rows() != 24 || stiffness->Columns() != 24)
    return model_storage::Invalid("MITC4 initial stiffness arguments", d);
  std::unique_ptr<Matrix> owners[7];
  constexpr int kRows[7] = {4, 3, 4, 5, 24, 24, 24};
  constexpr int kColumns[7] = {3, 3, 24, 5, 24, 24, 24};
  for (int i = 0; i < 7; ++i) {
    Matrix* raw = nullptr;
    const auto status = Matrix::Create(kRows[i], kColumns[i], &raw, d);
    owners[i].reset(raw);
    if (status != StatusCode::kOk) return status;
  }
  Matrix* coordinates = owners[0].get();
  Matrix* frame = owners[1].get();
  Matrix* ties = owners[2].get();
  Matrix* constitutive = owners[3].get();
  Matrix* local = owners[4].get();
  Matrix* transform = owners[5].get();
  Matrix* product = owners[6].get();
  double diameter = 0.0;
  double length = 0.0;
  auto status = ValidateMitc4Geometry(nodes, count, coordinates, frame,
                                      &diameter, &length, d);
  if (status != StatusCode::kOk) return status;
  status = BuildMitc4LinearConstitutive(section->material(), constitutive, d);
  if (status != StatusCode::kOk) return status;
  if (elastic_only) {
    // RD12: beta damping acts on membrane/bending, not transverse shear.
    status = constitutive->Set(3, 3, 0, d);
    if (status == StatusCode::kOk) status = constitutive->Set(4, 4, 0, d);
    if (status != StatusCode::kOk) return status;
  }
  status = BuildMitc4LinearTying(coordinates, ties, d);
  if (status != StatusCode::kOk) return status;
  const double gauss = 1.0 / std::sqrt(3.0);
  for (double r : {-gauss, gauss})
    for (double s : {-gauss, gauss}) {
      status = IntegratePoint(coordinates, ties, constitutive, r, s,
                              section->thickness(), local, nullptr, nullptr, d);
      if (status != StatusCode::kOk) return status;
    }
  double coefficient = 0.0;
  if (!elastic_only) status = AddDrilling(local, &coefficient, d);
  if (status != StatusCode::kOk) return status;
  status = BuildTransform(frame, transform, d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Multiply(local, false, transform, false, product, d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Multiply(transform, true, product, false, stiffness, d);
  if (status != StatusCode::kOk) return status;
  if (axes) {
    status = axes->CopyFrom(frame, d);
    if (status != StatusCode::kOk) return status;
  }
  *drilling = coefficient;
  return StatusCode::kOk;
}

}  // namespace

StatusCode ComputeMitc4InitialStiffness(const Node* const* nodes, int count,
                                        const ShellSection* section,
                                        Matrix* stiffness, double* drilling,
                                        Diagnostic* d) {
  return InitialStiffness(nodes, count, section, stiffness, nullptr, drilling,
                          d);
}

StatusCode ComputeMitc4InitialElasticStiffness(const Node* const* nodes,
                                               int count,
                                               const ShellSection* section,
                                               Matrix* out, Diagnostic* d) {
  double unused_drilling = 0;
  return InitialStiffness(nodes, count, section, out, nullptr,
                          &unused_drilling, d, true);
}

StatusCode ComputeMitc4InitialProperties(const Node* const* nodes, int count,
                                         const ShellSection* section,
                                         Matrix* axes, double* drilling,
                                         Diagnostic* d) {
  if (!nodes || count != 4 || !section || !section->material() ||
      !Shape(axes, 3, 3) || !drilling)
    return model_storage::Invalid("MITC4 initial properties arguments", d);
  Matrix* raw = nullptr;
  const auto status = Matrix::Create(24, 24, &raw, d);
  std::unique_ptr<Matrix> stiffness(raw);
  if (status != StatusCode::kOk) return status;
  return InitialStiffness(nodes, count, section, stiffness.get(), axes,
                          drilling, d);
}

StatusCode ComputeMitc4LinearInternalForce(const Node* const* nodes, int count,
                                           const ShellSection* section,
                                           const Vector* u, double drilling,
                                           Vector* force, Diagnostic* d) {
  if (!nodes || count != 4 || !section || !section->material() || !u ||
      !force || u == force || u->Size() != 24 || force->Size() != 24 ||
      !std::isfinite(drilling) || drilling <= 0.0)
    return model_storage::Invalid("MITC4 linear force arguments", d);
  for (int i = 0; i < 24; ++i)
    if (!std::isfinite(u->Data()[i]))
      return model_storage::Invalid("MITC4 linear force finite displacement",
                                    d);
  std::unique_ptr<Matrix> owners[5];
  constexpr int kRows[5] = {4, 3, 4, 5, 24};
  constexpr int kColumns[5] = {3, 3, 24, 5, 24};
  for (int i = 0; i < 5; ++i) {
    Matrix* raw = nullptr;
    const auto status = Matrix::Create(kRows[i], kColumns[i], &raw, d);
    owners[i].reset(raw);
    if (status != StatusCode::kOk) return status;
  }
  std::unique_ptr<Vector> vectors[3];
  constexpr int kSizes[3] = {24, 24, 4};
  for (int i = 0; i < 3; ++i) {
    Vector* raw = nullptr;
    const auto status = Vector::Create(kSizes[i], &raw, d);
    vectors[i].reset(raw);
    if (status != StatusCode::kOk) return status;
  }
  Matrix* coordinates = owners[0].get();
  Matrix* frame = owners[1].get();
  Matrix* ties = owners[2].get();
  Matrix* constitutive = owners[3].get();
  Matrix* transform = owners[4].get();
  Vector* displacement = vectors[0].get();
  Vector* local_force = vectors[1].get();
  Vector* drilling_moments = vectors[2].get();
  double diameter = 0.0;
  double length = 0.0;
  auto status = ValidateMitc4Geometry(nodes, count, coordinates, frame,
                                      &diameter, &length, d);
  if (status != StatusCode::kOk) return status;
  status = BuildMitc4LinearConstitutive(section->material(), constitutive, d);
  if (status != StatusCode::kOk) return status;
  status = BuildMitc4LinearTying(coordinates, ties, d);
  if (status != StatusCode::kOk) return status;
  status = BuildTransform(frame, transform, d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::MultiplyVector(transform, false, u, displacement, d);
  if (status != StatusCode::kOk) return status;
  const double gauss = 1.0 / std::sqrt(3.0);
  for (double r : {-gauss, gauss})
    for (double s : {-gauss, gauss}) {
      status = IntegratePoint(coordinates, ties, constitutive, r, s,
                              section->thickness(), nullptr, displacement,
                              local_force, d);
      if (status != StatusCode::kOk) return status;
    }
  // Frozen local drilling springs are independent of constitutive quadrature.
  const int indices[4] = {5, 11, 17, 23};
  status = drilling_moments->GatherFrom(displacement, indices, 4, 0, d);
  if (status != StatusCode::kOk) return status;
  status = drilling_moments->Scale(drilling, d);
  if (status != StatusCode::kOk) return status;
  status = local_force->ScatterAdd(drilling_moments, indices, 4, d);
  if (status != StatusCode::kOk) return status;
  return Matrix::MultiplyVector(transform, true, local_force, force, d);
}

}  // namespace fesa
