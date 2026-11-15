/// @file
/// @brief Immutable elasticity, optional density and Rayleigh coefficients.
#ifndef FESA_MODEL_ELASTIC_MATERIAL_H_
#define FESA_MODEL_ELASTIC_MATERIAL_H_

#include "fesa/model/material.h"

namespace fesa {

/// Owns E, nu, optional rho and Rayleigh scalars; elements form matrices.
class ElasticMaterial final : public Material {
 public:
  /// Copies name. Caller owns success; failure nulls *out and frees resources.
  /// Requires positive ID/E, -1<nu<.5, and finite positive rho when present.
  /// Absent density ignores density, including nonfinite values.
  static StatusCode Create(int id, const char* name, double young,
                           double poisson, bool has_density, double density,
                           ElasticMaterial** out, Diagnostic* d);
  /// Same ownership/failure contract as Create; stores finite nonnegative
  /// Rayleigh coefficients alpha [1/T], beta [T], independent of HHT alpha.
  /// No damping matrix is built; elements consume the immutable definition.
  static StatusCode Create(int id, const char* name, double young,
                           double poisson, bool has_density, double density,
                           double rayleigh_alpha, double rayleigh_beta,
                           ElasticMaterial** out, Diagnostic* d);
  /// Mass proportional Rayleigh coefficient, 1/time; default zero.
  double rayleigh_alpha() const;
  /// Elastic stiffness proportional coefficient, time; default zero.
  /// Nonlinear element exclusions follow docs/HHT/mitc4_rayleigh.md.
  double rayleigh_beta() const;
  ~ElasticMaterial() override = default;
  ElasticMaterial(const ElasticMaterial&) = delete;
  ElasticMaterial& operator=(const ElasticMaterial&) = delete;
  ElasticMaterial(ElasticMaterial&&) = delete;
  ElasticMaterial& operator=(ElasticMaterial&&) = delete;
  /// Borrowed this pointer, valid through destruction.
  const ElasticMaterial* AsElastic() const override;
  /// E in force/length squared.
  double young() const;
  /// Dimensionless Poisson ratio.
  double poisson() const;
  /// True only when density was explicitly provided.
  bool has_density() const;
  /// Writes mass/length cubed; missing density returns kInputError.
  /// Invalid output or missing density preserves the caller's scalar.
  StatusCode Density(double* out, Diagnostic* d) const;

 private:
  ElasticMaterial(int id, char* name, double young, double poisson,
                  bool has_density, double density, double rayleigh_alpha,
                  double rayleigh_beta);
  double young_;
  double poisson_;
  bool has_density_;
  double density_;
  double rayleigh_alpha_;
  double rayleigh_beta_;
};

}  // namespace fesa

#endif  // FESA_MODEL_ELASTIC_MATERIAL_H_
