#include "fesa/base/diagnostic.h"

namespace fesa {

void Diagnostic::Record(StatusCode value, DiagnosticStage failed_stage,
                        const char* failed_operation, int native_code) {
  active_flux_failure = false;
  geometry_point_available = false;
  geometry_value_available = false;
  geometry_quantity = GeometryQuantity::kNone;
  code = value;
  stage = failed_stage;
  operation = failed_operation;
  vendor_code = native_code;
}

StatusCode ReportFailure(StatusCode code, DiagnosticStage stage,
                         const char* operation, Diagnostic* d) {
  if (d != nullptr) d->Record(code, stage, operation);
  return code;
}

}  // namespace fesa
