/// @file
/// @brief Immutable nodal initial temperature, independent of result buffers.
#ifndef FESA_MODEL_THERMAL_INITIAL_TEMPERATURE_H_
#define FESA_MODEL_THERMAL_INITIAL_TEMPERATURE_H_

#include "fesa/math/vector.h"

namespace fesa {
/// Domain-owned input record. Values follow labels 11..10+count in model
/// temperature units. Source spelling and values are owned explicit copies.
class ThermalInitialTemperature final {
 public:
  /// Borrows inputs for this call. Positive node/line, nonempty source and
  /// 1..20 finite values required. Caller owns success; failure nulls *out.
  static StatusCode Create(int node_id, const Vector* values,
                           const char* source_file, int source_line,
                           ThermalInitialTemperature** out, Diagnostic* d);
  ~ThermalInitialTemperature();
  ThermalInitialTemperature(const ThermalInitialTemperature&) = delete;
  ThermalInitialTemperature& operator=(const ThermalInitialTemperature&) =
      delete;
  ThermalInitialTemperature(ThermalInitialTemperature&&) = delete;
  ThermalInitialTemperature& operator=(ThermalInitialTemperature&&) = delete;
  /// Node ID used as the Domain registry key.
  int Id() const;
  /// Borrowed immutable values, valid until destruction.
  const Vector* values() const;
  /// Borrowed owned spelling, valid until destruction.
  const char* source_file() const;
  /// One-based input source line.
  int source_line() const;

 private:
  ThermalInitialTemperature(int node_id, Vector* values, char* source_file,
                            int source_line);
  int node_id_;
  Vector* values_;     ///< Owned; never aliases a Node result.
  char* source_file_;  ///< Owned copied input spelling.
  int source_line_;
};
}  // namespace fesa

#endif  // FESA_MODEL_THERMAL_INITIAL_TEMPERATURE_H_
