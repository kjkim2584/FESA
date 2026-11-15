/// @file
/// @brief Immutable constant-thickness homogeneous shell section.
#ifndef FESA_MODEL_SHELL_SECTION_H_
#define FESA_MODEL_SHELL_SECTION_H_

#include "fesa/model/elastic_material.h"

namespace fesa {
/// Borrows material through destruction. Contains no geometry, D or K cache.
/// Thickness uses the model length unit; constitutive integration belongs to
/// the shell element. Missing density is legal for static analyses.
class ShellSection {
 public:
  /// Positive ID, finite positive thickness and nonnull material required.
  /// Borrows material without changing it. Caller owns successful *out;
  /// failure nulls *out and releases partial resources. d is optional.
  /// Domain membership is checked by registration, not by this factory.
  static StatusCode Create(int id, double thickness,
                           const ElasticMaterial* material, ShellSection** out,
                           Diagnostic* d);
  ~ShellSection() = default;
  ShellSection(const ShellSection&) = delete;
  ShellSection& operator=(const ShellSection&) = delete;
  ShellSection(ShellSection&&) = delete;
  ShellSection& operator=(ShellSection&&) = delete;
  /// Original positive section ID.
  int Id() const;
  /// Constant thickness in model length units.
  double thickness() const;
  /// Identical borrowed material; must outlive this section.
  const ElasticMaterial* material() const;

 private:
  ShellSection(int id, double thickness, const ElasticMaterial* material);
  int id_;
  double thickness_;
  const ElasticMaterial* material_;  ///< Borrowed; never deleted here.
};
}  // namespace fesa
#endif  // FESA_MODEL_SHELL_SECTION_H_
