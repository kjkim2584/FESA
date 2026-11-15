#include "fesa/model/thermal_section.h"

#include <cmath>
#include <new>

#include "model/model_storage.h"

namespace fesa {
ThermalSolidSection::ThermalSolidSection(int id,
                                         const ThermalMaterial* material)
    : id_(id), material_(material) {}
StatusCode ThermalSolidSection::Create(int id, const ThermalMaterial* material,
                                       ThermalSolidSection** out,
                                       Diagnostic* d) {
  if (!out) return model_storage::Invalid("ThermalSolidSection output", d);
  *out = nullptr;
  if (!model_storage::ValidId(id) || !material)
    return model_storage::Invalid("ThermalSolidSection definition", d);
  *out = new (std::nothrow) ThermalSolidSection(id, material);
  return *out ? StatusCode::kOk
              : model_storage::AllocationFailure("ThermalSolidSection", d);
}
int ThermalSolidSection::Id() const { return id_; }
const ThermalMaterial* ThermalSolidSection::material() const {
  return material_;
}
}  // namespace fesa

namespace fesa {
ThermalShellSection::ThermalShellSection(int id, double h, int ns,
                                         const ThermalMaterial* material)
    : id_(id), thickness_(h), section_points_(ns), material_(material) {}
StatusCode ThermalShellSection::Create(int id, double h, int ns,
                                       const ThermalMaterial* material,
                                       ThermalShellSection** out,
                                       Diagnostic* d) {
  if (!out) return model_storage::Invalid("ThermalShellSection output", d);
  *out = nullptr;
  if (!model_storage::ValidId(id) || !(h > 0) || !std::isfinite(h) || ns < 1 ||
      ns > 19 || ns % 2 == 0 || !material)
    return model_storage::Invalid("ThermalShellSection definition", d);
  *out = new (std::nothrow) ThermalShellSection(id, h, ns, material);
  return *out ? StatusCode::kOk
              : model_storage::AllocationFailure("ThermalShellSection", d);
}
int ThermalShellSection::Id() const { return id_; }
double ThermalShellSection::thickness() const { return thickness_; }
int ThermalShellSection::section_points() const { return section_points_; }
const ThermalMaterial* ThermalShellSection::material() const {
  return material_;
}
}  // namespace fesa
