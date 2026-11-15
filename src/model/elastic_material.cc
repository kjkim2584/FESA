#include "fesa/model/elastic_material.h"

#include <cmath>
#include <memory>
#include <new>

#include "model/model_storage.h"

namespace fesa {

ElasticMaterial::ElasticMaterial(int id, char* name, double young,
                                 double poisson, bool has_density,
                                 double density, double rayleigh_alpha,
                                 double rayleigh_beta)
    : Material(id, name),
      young_(young),
      poisson_(poisson),
      has_density_(has_density),
      density_(has_density ? density : 0),
      rayleigh_alpha_(rayleigh_alpha),
      rayleigh_beta_(rayleigh_beta) {}
StatusCode ElasticMaterial::Create(int id, const char* name, double young,
                                   double poisson, bool has_density,
                                   double density, ElasticMaterial** out,
                                   Diagnostic* d) {
  return Create(id, name, young, poisson, has_density, density, 0, 0, out, d);
}
StatusCode ElasticMaterial::Create(int id, const char* name, double young,
                                   double poisson, bool has_density,
                                   double density, double rayleigh_alpha,
                                   double rayleigh_beta, ElasticMaterial** out,
                                   Diagnostic* d) {
  if (!out) return model_storage::Invalid("ElasticMaterial::Create output", d);
  *out = nullptr;
  if (!model_storage::ValidId(id) || !name || !name[0] ||
      !std::isfinite(young) || young <= 0 || !std::isfinite(poisson) ||
      poisson <= -1 || poisson >= .5 ||
      (has_density && (!std::isfinite(density) || density <= 0)) ||
      !std::isfinite(rayleigh_alpha) || rayleigh_alpha < 0 ||
      !std::isfinite(rayleigh_beta) || rayleigh_beta < 0)
    return model_storage::Invalid("ElasticMaterial::Create properties", d);
  std::unique_ptr<char[]> copy(model_storage::CopyText(name));
  if (copy)
    *out = new (std::nothrow)
        ElasticMaterial(id, copy.get(), young, poisson, has_density, density,
                        rayleigh_alpha, rayleigh_beta);
  if (!*out)
    return model_storage::AllocationFailure("ElasticMaterial::Create", d);
  copy.release();
  return StatusCode::kOk;
}
const ElasticMaterial* ElasticMaterial::AsElastic() const { return this; }
double ElasticMaterial::rayleigh_alpha() const { return rayleigh_alpha_; }
double ElasticMaterial::rayleigh_beta() const { return rayleigh_beta_; }
double ElasticMaterial::young() const { return young_; }
double ElasticMaterial::poisson() const { return poisson_; }
bool ElasticMaterial::has_density() const { return has_density_; }
StatusCode ElasticMaterial::Density(double* out, Diagnostic* d) const {
  if (!out) return model_storage::Invalid("ElasticMaterial::Density output", d);
  if (!has_density_)
    return ReportFailure(StatusCode::kInputError, DiagnosticStage::kValidation,
                         "ElasticMaterial density is missing", d);
  *out = density_;
  return StatusCode::kOk;
}

}  // namespace fesa
