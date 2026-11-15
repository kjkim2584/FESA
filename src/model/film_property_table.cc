#include "fesa/model/film_property_table.h"

#include <memory>
#include <new>

#include "model/model_storage.h"

namespace fesa {
FilmPropertyTable::FilmPropertyTable(int id, char* name,
                                     TemperatureTable* h_of_theta)
    : id_(id), name_(name), h_of_theta_(h_of_theta) {}
FilmPropertyTable::~FilmPropertyTable() {
  delete[] name_;
  delete h_of_theta_;
}
StatusCode FilmPropertyTable::Create(int id, const char* name,
                                     const TemperatureTable* h_of_theta,
                                     FilmPropertyTable** out, Diagnostic* d) {
  if (!out) return model_storage::Invalid("FilmPropertyTable output", d);
  *out = nullptr;
  if (!model_storage::ValidId(id) || !name || !*name || !h_of_theta ||
      !h_of_theta->HasNonnegativeValues())
    return model_storage::Invalid("FilmPropertyTable properties", d);
  std::unique_ptr<char[]> owned_name(model_storage::CopyText(name));
  if (!owned_name)
    return model_storage::AllocationFailure("FilmPropertyTable name", d);
  std::unique_ptr<FilmPropertyTable> result(
      new (std::nothrow) FilmPropertyTable(id, owned_name.get(), nullptr));
  if (!result) return model_storage::AllocationFailure("FilmPropertyTable", d);
  owned_name.release();
  auto status = h_of_theta->Clone(&result->h_of_theta_, d);
  if (status != StatusCode::kOk) return status;
  *out = result.release();
  return StatusCode::kOk;
}
int FilmPropertyTable::Id() const { return id_; }
const char* FilmPropertyTable::Name() const { return name_; }
StatusCode FilmPropertyTable::Evaluate(double theta, double* value,
                                       double* slope, Diagnostic* d) const {
  return h_of_theta_->Evaluate(theta, value, slope, d);
}
bool FilmPropertyTable::HasHeatExchange() const {
  return !h_of_theta_->IsZero();
}
}  // namespace fesa
