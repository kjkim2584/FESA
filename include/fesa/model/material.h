/// @file
/// @brief Linear-static model ownership and immutable input definitions.
#ifndef FESA_MODEL_MATERIAL_H_
#define FESA_MODEL_MATERIAL_H_

#include "fesa/base/diagnostic.h"

namespace fesa {

class ElasticMaterial;
class ThermalMaterial;

/// Immutable material identity; Domain owns registered definitions.
/// Owns its original name. Constitutive properties belong to derived types.
class Material {
 public:
  /// Releases owned name; permits deletion through Domain base pointers.
  virtual ~Material();
  /// Borrowed elastic definition, or null for a different material kind.
  virtual const ElasticMaterial* AsElastic() const;
  /// Borrowed thermal definition, or null for a different material kind.
  virtual const ThermalMaterial* AsThermal() const;
  Material(const Material&) = delete;
  Material& operator=(const Material&) = delete;
  Material(Material&&) = delete;
  Material& operator=(Material&&) = delete;
  /// Original/internal positive ID.
  int Id() const;
  /// Borrowed original name, valid until destruction.
  const char* Name() const;

 protected:
  /// Adopts caller-provided owned name.
  Material(int id, char* name);

 private:
  int id_;
  char* name_;  ///< Owned original spelling.
};

}  // namespace fesa

#endif  // FESA_MODEL_MATERIAL_H_
