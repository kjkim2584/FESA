#include "fesa/elements/solid_reference_data.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <new>

#include "elements/hexa_interpolation.h"
#include "elements/pyramid_interpolation.h"
#include "elements/tetra_interpolation.h"
#include "elements/wedge_interpolation.h"
#include "model/model_storage.h"
namespace fesa {
StatusCode FillPyramidInterpolation(double xi, double eta, double zeta,
                                    Matrix* shape, Matrix* gradient,
                                    Diagnostic* d) {
  if (!shape || shape->Rows() != 1 || shape->Columns() != 5 ||
      (gradient && (gradient->Rows() != 3 || gradient->Columns() != 5)) ||
      !std::isfinite(xi) || !std::isfinite(eta) || !std::isfinite(zeta) ||
      zeta < 0 || zeta > 1 || std::abs(xi) > 1 - zeta ||
      std::abs(eta) > 1 - zeta || (zeta == 1 && gradient))
    return model_storage::Invalid("Pyramid interpolation arguments", d);
  const int signs[4][2] = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}};
  const double q = 1 - zeta;
  // F5 apex limit is shape-only; never evaluate its path-dependent gradient.
  const double a = q == 0 ? 0 : xi / q, b = q == 0 ? 0 : eta / q;
  for (int node = 0; node < 5; ++node) {
    const double s = node == 4 ? 0 : signs[node][0];
    const double t = node == 4 ? 0 : signs[node][1];
    auto status = shape->Set(
        0, node, node == 4 ? zeta : q * (1 + s * a) * (1 + t * b) / 4, d);
    if (status != StatusCode::kOk) return status;
    if (gradient) {
      const double derivatives[] = {node == 4 ? 0 : s * (1 + t * b) / 4,
                                    node == 4 ? 0 : t * (1 + s * a) / 4,
                                    node == 4 ? 1 : (-1 + s * t * a * b) / 4};
      for (int axis = 0; axis < 3; ++axis) {
        status = gradient->Set(axis, node, derivatives[axis], d);
        if (status != StatusCode::kOk) return status;
      }
    }
  }
  return StatusCode::kOk;
}
StatusCode FillWedgeInterpolation(double xi, double eta, double zeta,
                                  Matrix* shape, Matrix* gradient,
                                  Diagnostic* d) {
  if (!shape || !gradient || shape->Rows() != 1 || shape->Columns() != 6 ||
      gradient->Rows() != 3 || gradient->Columns() != 6 || !std::isfinite(xi) ||
      !std::isfinite(eta) || !std::isfinite(zeta) || xi < 0 || eta < 0 ||
      xi + eta > 1 || std::abs(zeta) > 1)
    return model_storage::Invalid("Wedge interpolation arguments", d);
  const double l[] = {1 - xi - eta, xi, eta};
  const int derivatives[2][3] = {{-1, 1, 0}, {-1, 0, 1}};
  for (int node = 0; node < 6; ++node) {
    const int local = node % 3;
    const double sign = node < 3 ? -1 : 1;
    auto status = shape->Set(0, node, l[local] * (1 + sign * zeta) / 2, d);
    if (status != StatusCode::kOk) return status;
    for (int axis = 0; axis < 3; ++axis) {
      status = gradient->Set(
          axis, node,
          axis < 2 ? derivatives[axis][local] * (1 + sign * zeta) / 2
                   : sign * l[local] / 2,
          d);
      if (status != StatusCode::kOk) return status;
    }
  }
  return StatusCode::kOk;
}
StatusCode FillHexaInterpolation(double xi, double eta, double zeta,
                                 Matrix* shape, Matrix* gradient,
                                 Diagnostic* d) {
  if (!shape || !gradient || shape->Rows() != 1 || shape->Columns() != 8 ||
      gradient->Rows() != 3 || gradient->Columns() != 8 || !std::isfinite(xi) ||
      !std::isfinite(eta) || !std::isfinite(zeta) || std::abs(xi) > 1 ||
      std::abs(eta) > 1 || std::abs(zeta) > 1)
    return model_storage::Invalid("Hexa interpolation arguments", d);
  const int signs[8][3] = {{-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1},
                           {-1, -1, 1},  {1, -1, 1},  {1, 1, 1},  {-1, 1, 1}};
  for (int node = 0; node < 8; ++node) {
    const double factors[] = {1 + signs[node][0] * xi, 1 + signs[node][1] * eta,
                              1 + signs[node][2] * zeta};
    auto status =
        shape->Set(0, node, factors[0] * factors[1] * factors[2] / 8, d);
    if (status != StatusCode::kOk) return status;
    for (int axis = 0; axis < 3; ++axis) {
      status = gradient->Set(axis, node,
                             signs[node][axis] * factors[(axis + 1) % 3] *
                                 factors[(axis + 2) % 3] / 8,
                             d);
      if (status != StatusCode::kOk) return status;
    }
  }
  return StatusCode::kOk;
}
StatusCode FillTetraInterpolation(int nodes, double xi, double eta, double zeta,
                                  Matrix* shape, Matrix* gradient,
                                  Diagnostic* d) {
  if ((nodes != 4 && nodes != 10) || !shape || !gradient ||
      shape->Rows() != 1 || shape->Columns() != nodes ||
      gradient->Rows() != 3 || gradient->Columns() != nodes ||
      !std::isfinite(xi) || !std::isfinite(eta) || !std::isfinite(zeta) ||
      !std::isfinite(1 - xi - eta - zeta))
    return model_storage::Invalid("Tetra interpolation arguments", d);
  const double l[] = {1 - xi - eta - zeta, xi, eta, zeta};
  const double dl[3][4] = {{-1, 1, 0, 0}, {-1, 0, 1, 0}, {-1, 0, 0, 1}};
  const int edges[6][2] = {{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}};
  for (int node = 0; node < nodes; ++node) {
    const double value =
        nodes == 4 ? l[node]
        : node < 4 ? l[node] * (2 * l[node] - 1)
                   : 4 * l[edges[node - 4][0]] * l[edges[node - 4][1]];
    auto status = shape->Set(0, node, value, d);
    if (status != StatusCode::kOk) return status;
    for (int axis = 0; axis < 3; ++axis) {
      const double derivative =
          nodes == 4 ? dl[axis][node]
          : node < 4
              ? (4 * l[node] - 1) * dl[axis][node]
              : 4 * (l[edges[node - 4][1]] * dl[axis][edges[node - 4][0]] +
                     l[edges[node - 4][0]] * dl[axis][edges[node - 4][1]]);
      status = gradient->Set(axis, node, derivative, d);
      if (status != StatusCode::kOk) return status;
    }
  }
  return StatusCode::kOk;
}
SolidReferenceData::~SolidReferenceData() {
  if (gradients_)
    for (int p = 0; p < point_count_; ++p) delete gradients_[p];
  delete[] gradients_;
  delete shapes_;
  delete weights_;
  delete coordinates_;
}
StatusCode SolidReferenceData::InitializeTrialGeometry(Diagnostic* d) {
  double points[62][3] = {};
  int count = 0;
  auto add = [&](double x, double y, double z) {
    for (int i = 0; i < count; ++i)
      if (points[i][0] == x && points[i][1] == y && points[i][2] == z) return;
    points[count][0] = x;
    points[count][1] = y;
    points[count++][2] = z;
  };
  if (node_count_ == 4) {
    add(.25, .25, .25);
  } else if (node_count_ == 10) {
    const double nodes[10][3] = {{0,0,0},{1,0,0},{0,1,0},{0,0,1},
        {.5,0,0},{.5,.5,0},{0,.5,0},{0,0,.5},{.5,0,.5},{0,.5,.5}};
    for (const auto& node : nodes) add(node[0], node[1], node[2]);
    add(.25, .25, .25);
  } else if (node_count_ == 8) {
    for (int z = -1; z <= 1; ++z)
      for (int y = -1; y <= 1; ++y)
        for (int x = -1; x <= 1; ++x) add(x, y, z);
  } else if (node_count_ == 6) {
    const double triangle[7][2] = {{0,0},{1,0},{0,1},{.5,0},
                                  {.5,.5},{0,.5},{1./3,1./3}};
    for (int z = -1; z <= 1; ++z)
      for (const auto& p : triangle) add(p[0], p[1], z);
  } else {
    for (int y = -1; y <= 1; ++y)
      for (int x = -1; x <= 1; ++x) add(x, y, 0);
    for (double y : {-.5, .5})
      for (double x : {-.5, .5}) add(x, y, .5);
    add(0, 0, .25);
  }
  if (node_count_ != 4) {
    for (auto purpose : {IntegrationPurpose::kStiffness,
                         IntegrationPurpose::kMass}) {
      SolidReferenceData* raw = nullptr;
      auto status = Create(formulation_, purpose, &raw, d);
      std::unique_ptr<SolidReferenceData> quadrature(raw);
      if (status != StatusCode::kOk) return status;
      for (int i = 0; i < raw->point_count(); ++i) {
        const auto* p = raw->natural_coordinates()->Data() + 3 * i;
        add(p[0], p[1], p[2]);
      }
    }
  }
  point_count_ = count;
  rule_id_ = "Q23_Q53_TRIAL_GEOMETRY";
  auto status = Matrix::Create(count, 3, &coordinates_, d);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(count, &weights_, d);  // Nonintegration samples: zero.
  if (status != StatusCode::kOk) return status;
  status = Matrix::Create(count, node_count_, &shapes_, d);
  if (status != StatusCode::kOk) return status;
  gradients_ = new (std::nothrow) Matrix*[count] {};
  if (!gradients_)
    return model_storage::AllocationFailure("Trial geometry gradients", d);
  Matrix* raw = nullptr;
  status = Matrix::Create(1, node_count_, &raw, d);
  std::unique_ptr<Matrix> row(raw);
  if (status != StatusCode::kOk) return status;
  for (int i = 0; i < count; ++i) {
    for (int axis = 0; axis < 3; ++axis) {
      status = coordinates_->Set(i, axis, points[i][axis], d);
      if (status != StatusCode::kOk) return status;
    }
    status = Matrix::Create(3, node_count_, &gradients_[i], d);
    if (status != StatusCode::kOk) return status;
    const double* p = points[i];
    status = node_count_ == 8
        ? FillHexaInterpolation(p[0], p[1], p[2], raw, gradients_[i], d)
        : node_count_ == 6
        ? FillWedgeInterpolation(p[0], p[1], p[2], raw, gradients_[i], d)
        : node_count_ == 5
        ? FillPyramidInterpolation(p[0], p[1], p[2], raw, gradients_[i], d)
        : FillTetraInterpolation(node_count_, p[0], p[1], p[2], raw,
                                  gradients_[i], d);
    if (status != StatusCode::kOk) return status;
    status = shapes_->CopyBlockFrom(raw, 0, 0, 1, node_count_, i, 0, d);
    if (status != StatusCode::kOk) return status;
  }
  return StatusCode::kOk;
}
StatusCode SolidReferenceData::Create(SolidFormulation formulation,
                                      IntegrationPurpose purpose,
                                      SolidReferenceData** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("SolidReferenceData output", d);
  *out = nullptr;
  if (formulation == SolidFormulation::kTetra10Reduced)
    return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                         "FESA_C3D10R retired; use C3D10", d);
  const bool reduced_hexa = formulation == SolidFormulation::kHexa8Reduced;
  const bool hexa = formulation == SolidFormulation::kHexa8Full || reduced_hexa;
  const bool wedge = formulation == SolidFormulation::kWedge6;
  const bool reduced_pyramid =
      formulation == SolidFormulation::kPyramid5Reduced;
  const bool pyramid =
      formulation == SolidFormulation::kPyramid5Full || reduced_pyramid;
  const bool quadratic = formulation == SolidFormulation::kTetra10Full;
  if ((!hexa && !wedge && !pyramid && !quadratic &&
       formulation != SolidFormulation::kTetra4Full) ||
      (purpose != IntegrationPurpose::kStiffness &&
       purpose != IntegrationPurpose::kMass &&
       !(purpose == IntegrationPurpose::kTrialGeometry &&
         !reduced_hexa && !reduced_pyramid) &&
       !((hexa || wedge || pyramid) &&
         purpose == IntegrationPurpose::kGeometry)))
    return ReportFailure(StatusCode::kUnsupported, DiagnosticStage::kValidation,
                         "Solid reference rule not implemented", d);
  std::unique_ptr<SolidReferenceData> data(new (std::nothrow)
                                               SolidReferenceData);
  if (!data)
    return model_storage::AllocationFailure("SolidReferenceData object", d);
  data->formulation_ = formulation;
  data->purpose_ = purpose;
  data->node_count_ = pyramid ? 5 : wedge ? 6 : hexa ? 8 : quadratic ? 10 : 4;
  if (purpose == IntegrationPurpose::kTrialGeometry) {
    const auto status = data->InitializeTrialGeometry(d);
    if (status != StatusCode::kOk) return status;
    *out = data.release();
    return StatusCode::kOk;
  }
  data->point_count_ =
      pyramid ? (purpose == IntegrationPurpose::kStiffness
                     ? (reduced_pyramid ? 1 : 5)
                     : 8)
      : wedge ? (purpose == IntegrationPurpose::kStiffness ? 2 : 18)
      : hexa
          ? (purpose == IntegrationPurpose::kStiffness ? (reduced_hexa ? 1 : 8)
                                                       : 27)
      : purpose == IntegrationPurpose::kMass          ? (quadratic ? 15 : 4)
      : formulation == SolidFormulation::kTetra10Full ? 4
                                                      : 1;
  data->rule_id_ =
      pyramid ? (data->point_count_ == 1   ? "PYR1"
                 : data->point_count_ == 5 ? "PYR5_FELIPPA"
                                           : "PYR8")
      : wedge ? (data->point_count_ == 2 ? "WEDGE2" : "WEDGE18")
      : hexa  ? (data->point_count_ == 1   ? "UNIFORM_STRAIN_STANDARD_STIFFNESS"
                 : data->point_count_ == 8 ? "HEX8_BBAR"
                                           : "HEX27")
      : data->point_count_ == 15 ? "TET15"
      : data->point_count_ == 4  ? "TET4"
                                 : "TET1";
  double points[27][4] = {}, weights[27] = {};
  int count = 0;
  if (pyramid) {
    if (data->point_count_ == 1) {
      points[0][3] = .25;
      weights[count++] = 4. / 3;
    } else if (data->point_count_ == 5) {
      // Felippa O05: source weights sum to one; these include volume 4/3.
      const double a = std::sqrt(32. / 135);
      const int signs[4][2] = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}};
      for (int p = 0; p < 4; ++p) {
        points[count][1] = signs[p][0] * a;
        points[count][2] = signs[p][1] * a;
        points[count][3] = 1. / 6;
        weights[count++] = 9. / 32;
      }
      points[count][3] = 7. / 10;
      weights[count++] = 5. / 24;
    } else {
      const double r = std::sqrt(10.);
      const double z[] = {(5 - r) / 15, (5 + r) / 15};
      const double w[] = {1. / 6 + r / 48, 1. / 6 - r / 48};
      for (int layer = 0; layer < 2; ++layer)
        for (int eta = -1; eta <= 1; eta += 2)
          for (int xi = -1; xi <= 1; xi += 2) {
            points[count][1] = xi * (1 - z[layer]) / std::sqrt(3.);
            points[count][2] = eta * (1 - z[layer]) / std::sqrt(3.);
            points[count][3] = z[layer];
            weights[count++] = w[layer];
          }
    }
  } else if (wedge) {
    const bool stiffness = purpose == IntegrationPurpose::kStiffness;
    double triangle[6][3] = {{1. / 3, 1. / 3, 1. / 3}},
           triangle_weights[6] = {.5};
    int triangle_count = 1;
    if (!stiffness) {
      triangle_count = 0;
      const double representatives[2][3] = {
          {.108103018168070, .445948490915965, .445948490915965},
          {.816847572980459, .091576213509771, .091576213509771}};
      const double w[] = {.223381589678011 / 2, .109951743655322 / 2};
      for (int orbit = 0; orbit < 2; ++orbit) {
        double l[] = {representatives[orbit][0], representatives[orbit][1],
                      representatives[orbit][2]};
        std::sort(l, l + 3);
        do {
          for (int axis = 0; axis < 3; ++axis)
            triangle[triangle_count][axis] = l[axis];
          triangle_weights[triangle_count++] = w[orbit];
        } while (std::next_permutation(l, l + 3));
      }
    }
    int order[6] = {};
    for (int i = 0; i < triangle_count; ++i) order[i] = i;
    std::sort(order, order + triangle_count, [&triangle](int a, int b) {
      return std::lexicographical_compare(triangle[a], triangle[a] + 3,
                                          triangle[b], triangle[b] + 3);
    });
    const double a = stiffness ? 1 / std::sqrt(3.) : std::sqrt(3. / 5);
    const double z[] = {-a, stiffness ? a : 0, a},
                 wz[] = {stiffness ? 1. : 5. / 9, stiffness ? 1. : 8. / 9,
                         5. / 9};
    for (int axis = 0; axis < (stiffness ? 2 : 3); ++axis)
      for (int t = 0; t < triangle_count; ++t) {
        const int index = order[t];
        for (int j = 0; j < 3; ++j) points[count][j] = triangle[index][j];
        points[count][3] = z[axis];
        weights[count++] = triangle_weights[index] * wz[axis];
      }
  } else if (hexa) {
    const int degree = data->point_count_ == 1   ? 1
                       : data->point_count_ == 8 ? 2
                                                 : 3;
    const double a = degree == 1   ? 0
                     : degree == 2 ? 1 / std::sqrt(3.)
                                   : std::sqrt(3. / 5);
    const double x[] = {-a, degree == 2 ? a : 0, a};
    const double w[] = {degree == 1   ? 2.
                        : degree == 2 ? 1.
                                      : 5. / 9,
                        degree == 2 ? 1. : 8. / 9, 5. / 9};
    for (int z = 0; z < degree; ++z)
      for (int y = 0; y < degree; ++y)
        for (int xi = 0; xi < degree; ++xi) {
          points[count][1] = x[xi];
          points[count][2] = x[y];
          points[count][3] = x[z];
          weights[count++] = w[xi] * w[y] * w[z];
        }
  } else if (data->point_count_ == 15) {
    const double representatives[4][4] = {
        {.25, .25, .25, .25},
        {0, 1. / 3, 1. / 3, 1. / 3},
        {8. / 11, 1. / 11, 1. / 11, 1. / 11},
        {.0665501535736642813, .0665501535736642813, .433449846426335728,
         .433449846426335728}};
    const double w[] = {.0302836780970891856, .00602678571428571597,
                        .0116452490860289742, .0109491415613864534};
    for (int orbit = 0; orbit < 4; ++orbit) {
      double l[4] = {};
      for (int i = 0; i < 4; ++i) l[i] = representatives[orbit][i];
      std::sort(l, l + 4);
      do {
        for (int i = 0; i < 4; ++i) points[count][i] = l[i];
        weights[count++] = w[orbit];
      } while (std::next_permutation(l, l + 4));
    }
  } else if (data->point_count_ == 4) {
    const double a = (5 + 3 * std::sqrt(5.)) / 20, b = (5 - std::sqrt(5.)) / 20;
    for (int p = 0; p < 4; ++p) {
      for (int i = 0; i < 4; ++i) points[p][i] = i == 3 - p ? a : b;
      weights[count++] = 1. / 24;
    }
  } else {
    for (int i = 0; i < 4; ++i) points[0][i] = .25;
    weights[count++] = 1. / 6;
  }
  int order[27] = {};
  for (int i = 0; i < count; ++i) order[i] = i;
  if (!hexa && !wedge && !pyramid)
    std::sort(order, order + count, [&points](int a, int b) {
      return std::lexicographical_compare(points[a], points[a] + 4, points[b],
                                          points[b] + 4);
    });
  auto status = Matrix::Create(count, 3, &data->coordinates_, d);
  if (status != StatusCode::kOk) return status;
  status = Vector::Create(count, &data->weights_, d);
  if (status != StatusCode::kOk) return status;
  status = Matrix::Create(count, data->node_count_, &data->shapes_, d);
  if (status != StatusCode::kOk) return status;
  data->gradients_ = new (std::nothrow) Matrix* [count] {};
  if (!data->gradients_)
    return model_storage::AllocationFailure("Solid reference gradients", d);
  Matrix* row = nullptr;
  status = Matrix::Create(1, data->node_count_, &row, d);
  std::unique_ptr<Matrix> row_owner(row);
  if (status != StatusCode::kOk) return status;
  for (int p = 0; p < count; ++p) {
    const int index = order[p];
    status = data->weights_->Set(p, weights[index], d);
    if (status != StatusCode::kOk) return status;
    for (int axis = 0; axis < 3; ++axis) {
      status = data->coordinates_->Set(p, axis, points[index][axis + 1], d);
      if (status != StatusCode::kOk) return status;
    }
    status = Matrix::Create(3, data->node_count_, &data->gradients_[p], d);
    if (status != StatusCode::kOk) return status;
    status =
        pyramid ? FillPyramidInterpolation(points[index][1], points[index][2],
                                           points[index][3], row,
                                           data->gradients_[p], d)
        : wedge ? FillWedgeInterpolation(points[index][1], points[index][2],
                                         points[index][3], row,
                                         data->gradients_[p], d)
        : hexa ? FillHexaInterpolation(points[index][1], points[index][2],
                                       points[index][3], row,
                                       data->gradients_[p], d)
               : FillTetraInterpolation(data->node_count_, points[index][1],
                                        points[index][2], points[index][3], row,
                                        data->gradients_[p], d);
    if (status != StatusCode::kOk) return status;
    status =
        data->shapes_->CopyBlockFrom(row, 0, 0, 1, data->node_count_, p, 0, d);
    if (status != StatusCode::kOk) return status;
  }
  *out = data.release();
  return StatusCode::kOk;
}
SolidFormulation SolidReferenceData::formulation() const {
  return formulation_;
}
IntegrationPurpose SolidReferenceData::purpose() const { return purpose_; }
int SolidReferenceData::node_count() const { return node_count_; }
int SolidReferenceData::point_count() const { return point_count_; }
const char* SolidReferenceData::rule_id() const { return rule_id_; }
const Matrix* SolidReferenceData::natural_coordinates() const {
  return coordinates_;
}
const Vector* SolidReferenceData::weights() const { return weights_; }
const Matrix* SolidReferenceData::shape_values() const { return shapes_; }
const Matrix* SolidReferenceData::natural_gradients(int p) const {
  return p >= 0 && p < point_count_ ? gradients_[p] : nullptr;
}
}  // namespace fesa
