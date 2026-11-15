#include "fesa/analysis/residual_report.h"

#include <memory>
#include <new>

#include "model/model_storage.h"
namespace fesa {
StatusCode ResidualReport::Create(int global_count, ResidualReport** out,
                                  Diagnostic* d) {
  if (!out) return model_storage::Invalid("ResidualReport::Create", d);
  *out = nullptr;
  if (global_count <= 0 || global_count % 3 != 0)
    return model_storage::Invalid("ResidualReport::Create", d);
  std::unique_ptr<ResidualReport> result(new (std::nothrow) ResidualReport);
  if (!result) return model_storage::AllocationFailure("ResidualReport", d);
  result->displacement_ = new (std::nothrow) FieldResidual;
  result->rotation_ = new (std::nothrow) FieldResidual;
  if (!result->displacement_ || !result->rotation_)
    return model_storage::AllocationFailure("ResidualReport", d);
  const StatusCode status = Vector::Create(global_count, &result->reaction_, d);
  if (status != StatusCode::kOk) return status;
  *out = result.release();
  return StatusCode::kOk;
}
ResidualReport::~ResidualReport() {
  delete displacement_;
  delete rotation_;
  delete reaction_;
}
const FieldResidual* ResidualReport::Displacement() const {
  return displacement_;
}
const FieldResidual* ResidualReport::Rotation() const { return rotation_; }
const Vector* ResidualReport::Reaction() const { return reaction_; }
double ResidualReport::CharacteristicLength() const { return length_; }
}  // namespace fesa
