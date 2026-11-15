#include "fesa/model/nonlinear_static_parameters.h"

#include <cmath>

#include "model/model_storage.h"

namespace fesa {
StatusCode ValidateNonlinearStaticParameters(const NonlinearStaticParameters* p,
                                             Diagnostic* d) {
  if (!p || !std::isfinite(p->period) || p->period <= 0 ||
      !std::isfinite(p->initial_increment) || p->initial_increment <= 0 ||
      !std::isfinite(p->minimum_increment) || p->minimum_increment <= 0 ||
      p->minimum_increment > p->initial_increment || p->increment_limit <= 0 ||
      (p->has_maximum && (!std::isfinite(p->maximum_increment) ||
                          p->maximum_increment < p->initial_increment))) {
    return model_storage::Invalid("Nonlinear STATIC parameters", d);
  }
  return StatusCode::kOk;
}
}  // namespace fesa
