/// @file
/// @brief Adapter-owned thermal output selection and numeric metadata.
#ifndef FESA_ADAPTERS_HDF5_HDF5_THERMAL_SCHEMA_H_
#define FESA_ADAPTERS_HDF5_HDF5_THERMAL_SCHEMA_H_
#include <vector>

#include "fesa/base/diagnostic.h"
#include "hdf5.h"
namespace fesa {
class Domain;
class ThermalState;
class Vector;
namespace hdf5 {
/// Validates output shape/identities without opening a file or computing FEM.
/// Borrows Domain/d for the call, retains nothing and changes no model state.
StatusCode ValidateThermalResultModel(const Domain* domain, Diagnostic* d);
/// Owns participating node selection, nominal gather indices and scratch.
/// Domain remains immutable/alive until this schema is destroyed.
class ThermalSchema {
 public:
  /// Caller owns *out; failure nulls it and releases partial storage.
  static StatusCode Create(const Domain* domain, ThermalSchema** out,
                           Diagnostic* d);
  ~ThermalSchema();
  ThermalSchema(const ThermalSchema&) = delete;
  ThermalSchema& operator=(const ThermalSchema&) = delete;
  ThermalSchema(ThermalSchema&&) = delete;
  ThermalSchema& operator=(ThermalSchema&&) = delete;
  /// Immutable root family name, borrowed until schema destruction.
  const char* Family() const;
  /// Writes model metadata once. Failure may leave partial native output;
  /// caller aborts writer. Does not mutate/retain additional model objects.
  StatusCode WriteModel(hid_t model, Diagnostic* d) const;
  /// Gathers candidate NT/RFL via common Vector backend. Modifies only scratch
  /// and native output; never confirms/commits/mutates borrowed state.
  /// Caller validates nominal layout and result tags before creating a frame.
  StatusCode WriteFramePayload(hid_t nodal, const ThermalState* state,
                               Diagnostic* d);

 private:
  explicit ThermalSchema(const Domain* domain);
  const Domain* domain_;                ///< Borrowed immutable definitions.
  std::vector<int>* nodes_ = nullptr;   ///< Owned participating sorted IDs.
  std::vector<int>* gather_ = nullptr;  ///< Owned nominal indices.
  Vector* packed_ = nullptr;            ///< Owned reusable candidate scratch.
  bool shell_ = false;
};
}  // namespace hdf5
}  // namespace fesa
#endif  // FESA_ADAPTERS_HDF5_HDF5_THERMAL_SCHEMA_H_
