/// @file
/// @brief Immutable Domain-owned film coefficient table h(surface temperature).
#ifndef FESA_MODEL_FILM_PROPERTY_TABLE_H_
#define FESA_MODEL_FILM_PROPERTY_TABLE_H_

#include "fesa/model/temperature_table.h"

namespace fesa {
/// Owns a deep table copy and name. Conditions borrow this until Domain dies.
/// Scalar/table identity is retained even for constant h; no trial state.
class FilmPropertyTable final {
 public:
  /// Borrows name/table only for call. Requires positive ID, nonempty name,
  /// nonnegative h. Caller owns success; failure nulls *out/cleans resources.
  static StatusCode Create(int id, const char* name,
                           const TemperatureTable* h_of_theta,
                           FilmPropertyTable** out, Diagnostic* d);
  ~FilmPropertyTable();
  FilmPropertyTable(const FilmPropertyTable&) = delete;
  FilmPropertyTable& operator=(const FilmPropertyTable&) = delete;
  FilmPropertyTable(FilmPropertyTable&&) = delete;
  FilmPropertyTable& operator=(FilmPropertyTable&&) = delete;
  int Id() const;            ///< Original positive property ID.
  const char* Name() const;  ///< Borrowed owned spelling until destruction.
  /// Pure finite-temperature evaluation of h[E/T/L^2/theta] and dh/dtheta.
  /// Distinct required outputs are preserved on failure; no input retention.
  StatusCode Evaluate(double theta, double* value, double* slope,
                      Diagnostic* d) const;
  /// False only for an all-zero table; potential steady exchange anchor.
  bool HasHeatExchange() const;

 private:
  FilmPropertyTable(int id, char* name, TemperatureTable* h_of_theta);
  int id_;
  char* name_;                    ///< Owned null-terminated name.
  TemperatureTable* h_of_theta_;  ///< Owned immutable cloned h table.
};
}  // namespace fesa
#endif  // FESA_MODEL_FILM_PROPERTY_TABLE_H_
