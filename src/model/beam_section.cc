#include "fesa/model/beam_section.h"

#include <cmath>
#include <memory>
#include <new>

#include "model/model_storage.h"

namespace fesa {

BeamSection::BeamSection(int id, const ElasticMaterial* material, double area,
                         double i11, double i12, double i22, double torsion_j,
                         Vector* direction)
    : id_(id),
      material_(material),
      area_(area),
      i11_(i11),
      i12_(i12),
      i22_(i22),
      torsion_j_(torsion_j),
      approximate_n1_(direction) {}
BeamSection::~BeamSection() { delete approximate_n1_; }
StatusCode BeamSection::Create(int id, const ElasticMaterial* material,
                               double area, double i11, double i12, double i22,
                               double torsion_j, const Vector* approximate_n1,
                               BeamSection** out, Diagnostic* d) {
  if (out == nullptr) return model_storage::Invalid("BeamSection::Create", d);
  *out = nullptr;
  if (!model_storage::ValidId(id) || material == nullptr
      || approximate_n1 == nullptr || approximate_n1->Size() != 3
      || !std::isfinite(area) || area <= 0 || !std::isfinite(i11) || i11 <= 0
      || i12 != 0 || !std::isfinite(i22) || i22 <= 0
      || !std::isfinite(torsion_j) || torsion_j <= 0) {
    return model_storage::Invalid("BeamSection::Create", d);
  }
  Vector* raw = nullptr;
  StatusCode status = Vector::Create(3, &raw, d);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<Vector> direction(raw);
  status = direction->CopyFrom(approximate_n1, d);
  if (status != StatusCode::kOk) return status;
  *out = new (std::nothrow) BeamSection(id, material, area, i11, i12, i22,
                                      torsion_j, direction.get());
  if (*out == nullptr) {
    return model_storage::AllocationFailure("BeamSection::Create", d);
  }
  direction.release();
  return StatusCode::kOk;
}
int BeamSection::Id() const { return id_; }
const ElasticMaterial* BeamSection::GetMaterial() const { return material_; }
const Vector* BeamSection::ApproximateN1() const { return approximate_n1_; }
double BeamSection::Area() const { return area_; }
double BeamSection::I11() const { return i11_; }
double BeamSection::I12() const { return i12_; }
double BeamSection::I22() const { return i22_; }
double BeamSection::TorsionJ() const { return torsion_j_; }

}  // namespace fesa
