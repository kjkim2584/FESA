#include "fesa/model/solid_section.h"

#include <new>

#include "model/model_storage.h"

namespace fesa {

SolidSection::SolidSection(int id, const ElasticMaterial* material)
    : id_(id), material_(material) {}
StatusCode SolidSection::Create(int id, const ElasticMaterial* material,
                                SolidSection** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("SolidSection::Create output", d);
  *out = nullptr;
  if (!model_storage::ValidId(id) || !material)
    return model_storage::Invalid("SolidSection::Create definition", d);
  *out = new (std::nothrow) SolidSection(id, material);
  if (!*out) return model_storage::AllocationFailure("SolidSection::Create", d);
  return StatusCode::kOk;
}
int SolidSection::Id() const { return id_; }
const ElasticMaterial* SolidSection::material() const { return material_; }

}  // namespace fesa
