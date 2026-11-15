/// @file
/// @brief Immutable solid section connecting an element to its material.
#ifndef FESA_MODEL_SOLID_SECTION_H_
#define FESA_MODEL_SOLID_SECTION_H_

#include "fesa/model/elastic_material.h"

namespace fesa {

/// Borrows material, which must outlive this section. No geometry is stored.
class SolidSection {
 public:
  /// Caller owns success; failure nulls *out. Positive ID/material required.
  /// Borrowed inputs are unchanged. Domain adopts only on registration success.
  static StatusCode Create(int id, const ElasticMaterial* material,
                           SolidSection** out, Diagnostic* d);
  ~SolidSection() = default;
  SolidSection(const SolidSection&) = delete;
  SolidSection& operator=(const SolidSection&) = delete;
  SolidSection(SolidSection&&) = delete;
  SolidSection& operator=(SolidSection&&) = delete;
  /// Original positive section ID.
  int Id() const;
  /// Borrowed material, valid while owning Domain lives.
  const ElasticMaterial* material() const;

 private:
  SolidSection(int id, const ElasticMaterial* material);
  int id_;
  const ElasticMaterial* material_;  ///< Borrowed; Domain owns material.
};

}  // namespace fesa

#endif  // FESA_MODEL_SOLID_SECTION_H_
