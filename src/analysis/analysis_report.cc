#include "fesa/analysis/analysis_report.h"

namespace fesa {

void AnalysisReport::Reset() noexcept {
  status = calculation_status = StatusCode::kOk;
  calculation_finished = residual_valid = step_completed = false;
  residual = nullptr;
  displacement_field = rotation_field = nullptr;
  heat_flux_field = nullptr;
  logging_diagnostic = nullptr;
  characteristic_length = 0;
  equation_count = 0;
  step_number = increment_number = attempt_count = iteration_count = 0;
  linear_solver_calls = 0;
  last_committed_step = last_committed_increment = 0;
  last_committed_frame = 0;
  last_committed_total_time = 0;
  assembly_seconds = solve_seconds = equilibrium_recovery_seconds = 0;
  element_dispatch_seconds = calculation_seconds = output_seconds = 0;
  observed_element_team_max = 0;
}

}  // namespace fesa
