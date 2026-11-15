/// @file
/// @brief Analysis-owned scalar time proposals shared by step controllers.
#ifndef FESA_ANALYSIS_INCREMENT_CLOCK_H_
#define FESA_ANALYSIS_INCREMENT_CLOCK_H_

namespace fesa {
/// Analysis-owned metadata, no physical state. Next nominal is an intended
/// increment, not subtraction of rounded consecutive time values.
struct IncrementClock {
  double time = 0, compensation = 0;
  double nominal = 1;
  double last_increment = 0;  ///< Last accepted logical controller increment.
  int successful = 0, quick_successes = 0;  ///< quick capped at2.
};
/// Candidate tied to one exact committed clock. Caller never edits fields.
struct IncrementTarget {
  double base_time = 0, base_compensation = 0;
  int base_successful = 0;
  double nominal = 0;  ///< Requested increment after maximum cap, >=minimum.
  double size = 0;  ///< Logical h, or remaining time for endpoint clip/stretch.
  double time = 0, compensation = 0;
  bool completes_step = false, use_extrapolation = false;
  double extrapolation_ratio = 0;  ///< NA when use_extrapolation=false.
};
}  // namespace fesa
#endif  // FESA_ANALYSIS_INCREMENT_CLOCK_H_
