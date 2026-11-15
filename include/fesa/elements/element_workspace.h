/// @file
/// @brief Caller-owned Beam buffers or borrowed solid reference context.
#ifndef FESA_ELEMENTS_ELEMENT_WORKSPACE_H_
#define FESA_ELEMENTS_ELEMENT_WORKSPACE_H_

#include "fesa/math/matrix.h"

namespace fesa {
class SolidReferenceData;
class IncrementalElasticWorkspace;
class BeamKinematicWorkspace;
class BeamTimeKinematicWorkspace;
class BeamRayleighWorkspace;
enum class SolidFormulation;
enum class IntegrationPurpose;

/// One owner per serial worker; never stored by an Element. Calls overwrite
/// scratch data, including on failure. No scratch data is a committed result.
class ElementWorkspace {
 public:
  /// Creates buffers for element_dofs=12 only. Caller deletes *out; failure
  /// nulls it and releases partial storage. d is optional and borrowed.
  static StatusCode Create(int element_dofs, ElementWorkspace** out,
                           Diagnostic* d);
  /// Borrows immutable solid reference data through this context's lifetime;
  /// at least one pointer is required. Allocates no Beam numerical storage.
  /// Caller owns success; failure nulls output and cleans partial resources.
  static StatusCode CreateForSolid(const SolidReferenceData* stiffness,
                                   const SolidReferenceData* mass,
                                   const SolidReferenceData* geometry,
                                   ElementWorkspace** out, Diagnostic* d);
  /// Borrows the immutable array AND its nonnull unique-purpose tables until
  /// context destruction. Allocates no Beam numerical storage or pointer copy.
  /// Caller owns successful context; failure nulls output and cleans resources.
  static StatusCode CreateForSolidCatalog(
      const SolidReferenceData* const* references, int count,
      ElementWorkspace** out, Diagnostic* d);
  /// Borrows the catalog through destruction and owns exact-sized scratch for
  /// each stiffness formulation. Caller deletes success; failure nulls out.
  static StatusCode CreateForIncrementalSolidCatalog(
      const SolidReferenceData* const* references, int count,
      ElementWorkspace** out, Diagnostic* d);
  /// Borrowed mutable single-worker scratch; null for other constitutive modes.
  IncrementalElasticWorkspace* incremental_solid_workspace() const;
  /// Exact stiffness-formulation scratch, or null when absent. Borrowed until
  /// worker destruction; one worker must not evaluate concurrent element calls.
  IncrementalElasticWorkspace* incremental_solid_workspace(
      SolidFormulation formulation) const;
  /// Exact catalog lookup; missing pair returns null. Legacy direct contexts
  /// return their corresponding purpose pointer; element validates rule/shape.
  const SolidReferenceData* reference(SolidFormulation formulation,
                                      IntegrationPurpose purpose) const;
  /// Borrowed optional purpose data, valid while its Analysis/caller lives.
  const SolidReferenceData* stiffness_reference() const;
  /// Borrowed optional consistent-mass reference data.
  const SolidReferenceData* mass_reference() const;
  /// Borrowed optional uniform-strain geometry data.
  const SolidReferenceData* geometry_reference() const;
  ~ElementWorkspace();
  ElementWorkspace(const ElementWorkspace&) = delete;
  ElementWorkspace& operator=(const ElementWorkspace&) = delete;
  ElementWorkspace(ElementWorkspace&&) = delete;
  ElementWorkspace& operator=(ElementWorkspace&&) = delete;

 private:
  const SolidReferenceData* const* reference_catalog_ = nullptr;  ///< Borrowed.
  /// Owned catalog-sized array; only stiffness entries own nonnull scratch.
  IncrementalElasticWorkspace** incremental_ = nullptr;
  int reference_count_ = 0;  ///< Borrowed catalog entry count.
  const SolidReferenceData* stiffness_reference_ = nullptr;  ///< Borrowed.
  const SolidReferenceData* mass_reference_ = nullptr;       ///< Borrowed.
  const SolidReferenceData* geometry_reference_ = nullptr;   ///< Borrowed.
  friend class Beam3D;
  /// Lazy owned B33 numerical scratch. Mutable only by the serial worker's
  /// joined element calls; carries no model/state ownership or result contract.
  mutable BeamKinematicWorkspace* beam_kinematics_ = nullptr;
  mutable BeamTimeKinematicWorkspace* beam_time_kinematics_ = nullptr;
  /// Owned four-point Rayleigh scratch; numerical candidate buffers only.
  mutable BeamRayleighWorkspace* beam_rayleigh_ = nullptr;
  ElementWorkspace() = default;
  Matrix* local_k_ = nullptr;      ///< Owned B3 local stiffness.
  Matrix* transform_ = nullptr;    ///< Owned B4 global-to-local transformation.
  Matrix* product_ = nullptr;      ///< Owned intermediate K_l T.
  Matrix* recovery_ = nullptr;     ///< Owned B5 16-by-12 recovery coefficients.
  Vector* local_u_ = nullptr;      ///< Owned current candidate in local axes.
  Vector* local_force_ = nullptr;  ///< Owned local generalized force.
  Vector* end_results_ = nullptr;  ///< Owned B5 candidate results.
};

}  // namespace fesa

#endif  // FESA_ELEMENTS_ELEMENT_WORKSPACE_H_
