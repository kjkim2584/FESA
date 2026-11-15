#include "fesa/model/shell_section.h"

#include <cmath>
#include <new>

#include "model/model_storage.h"

namespace fesa {
ShellSection::ShellSection(int id, double thickness,
                           const ElasticMaterial* material)
    : id_(id), thickness_(thickness), material_(material) {}
StatusCode ShellSection::Create(int id, double thickness,
                                const ElasticMaterial* material,
                                ShellSection** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("ShellSection::Create output", d);
  *out = nullptr;
  if (!model_storage::ValidId(id) || !material || !std::isfinite(thickness) ||
      thickness <= 0)
    return model_storage::Invalid("ShellSection::Create definition", d);
  *out = new (std::nothrow) ShellSection(id, thickness, material);
  if (!*out) return model_storage::AllocationFailure("ShellSection::Create", d);
  return StatusCode::kOk;
}
int ShellSection::Id() const { return id_; }
double ShellSection::thickness() const { return thickness_; }
const ElasticMaterial* ShellSection::material() const { return material_; }
}  // namespace fesa
