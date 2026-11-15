/// @file
/// @brief Explicit scalar/borrowed metadata transfer across joined error paths.
#ifndef FESA_BASE_DIAGNOSTIC_COPY_H_
#define FESA_BASE_DIAGNOSTIC_COPY_H_

#include "fesa/base/diagnostic.h"

namespace fesa {
/// Copies all diagnostic metadata without allocation or physical mutation.
/// Null is a no-op, self-transfer is allowed. Source pointer fields must refer
/// to stable model/controller owners that outlive reporting through target;
/// this never transfers ownership or borrows worker-local string storage.
inline void CopyDiagnosticMetadata(const Diagnostic* source,
                                   Diagnostic* target) noexcept {
  if (!target || !source) return;
  target->code = source->code;
  target->stage = source->stage;
  target->operation = source->operation;
  target->source_file = source->source_file;
  target->path_kind = source->path_kind;
  target->source_line = source->source_line;
  target->first_source_line = source->first_source_line;
  target->step = source->step;
  target->node_id = source->node_id;
  target->element_id = source->element_id;
  target->condition_dof = source->condition_dof;
  target->duplicate_boundary = source->duplicate_boundary;
  target->boundary_dof = source->boundary_dof;
  target->first_boundary_value = source->first_boundary_value;
  target->boundary_value = source->boundary_value;
  target->active_flux_failure = source->active_flux_failure;
  target->flux_field = source->flux_field;
  target->flux_count = source->flux_count;
  target->active_flux_count = source->active_flux_count;
  target->flux_epsilon_l = source->flux_epsilon_l;
  target->flux_history_max = source->flux_history_max;
  target->flux_threshold = source->flux_threshold;
  target->flux_current_max = source->flux_current_max;
  target->geometry_point_available = source->geometry_point_available;
  target->natural_xi = source->natural_xi;
  target->natural_eta = source->natural_eta;
  target->natural_zeta = source->natural_zeta;
  target->geometry_value_available = source->geometry_value_available;
  target->geometry_quantity = source->geometry_quantity;
  target->geometry_value = source->geometry_value;
  target->geometry_limit = source->geometry_limit;
  target->vendor_code = source->vendor_code;
}
}  // namespace fesa
#endif  // FESA_BASE_DIAGNOSTIC_COPY_H_
