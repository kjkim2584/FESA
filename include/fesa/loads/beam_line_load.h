/// @file
/// @brief Uniform global beam line loads: immutable input, step history and
/// work-equivalent nodal forces. Formulation: docs/Linear static/beam_line_load.md.
#ifndef FESA_LOADS_BEAM_LINE_LOAD_H_
#define FESA_LOADS_BEAM_LINE_LOAD_H_

#include "fesa/io/input_deck.h"
#include "fesa/model/domain.h"

namespace fesa {

class BeamLineLoadProgramData;
class BeamLineLoadStateData;
class BeamLineLoadWorkspaceData;
class BeamLineLoadState;

/// Immutable DLOAD definitions and block operations; owns copied source text
/// and metadata, borrows Domain until destruction. Domain definitions remain
/// unchanged. No numerical results or solver thread settings are owned here.
class BeamLineLoadProgram {
 public:
  /// Parses DLOAD only from a successfully read deck and validates element IDs,
  /// scoped ELSET targets, PX/PY/PZ, finite F/L magnitudes and per-step OP.
  /// Empty NEW blocks are retained. Other deck blocks are never modified.
  /// Copies source text; deck need not outlive this. Domain must outlive this
  /// and every produced State. Caller deletes successful *out; failure nulls
  /// *out, releases partial resources and preserves inputs. out is required,
  /// creation-only. d is optional; failure source paths borrow deck until
  /// diagnostic reporting. Unsupported options/labels return kUnsupported;
  /// invalid source/targets/history return kInputError.
  static StatusCode Create(const InputDeck* deck, const Domain* domain,
                           BeamLineLoadProgram** out, Diagnostic* d);
  ~BeamLineLoadProgram();
  BeamLineLoadProgram(const BeamLineLoadProgram&) = delete;
  BeamLineLoadProgram& operator=(const BeamLineLoadProgram&) = delete;
  BeamLineLoadProgram(BeamLineLoadProgram&&) = delete;
  BeamLineLoadProgram& operator=(BeamLineLoadProgram&&) = delete;

  /// Resolves original scope/target/label before expansion (L3). Positive step
  /// must exist in Domain. Step 1 requires null previous; later steps require
  /// this program's preceding State. New State owns independent resolved rows
  /// and sums; it never borrows previous's storage. Previous may be deleted
  /// after success. Runtime caller passes last successful step; preflight may
  /// chain resolved states without asserting solver success. Caller deletes
  /// successful *out; failure nulls it and preserves this/Domain/previous.
  /// A repeated original key cannot be MODified in a subsequent step: remove
  /// with NEW then respecify. Diagnostics borrow this program's source path.
  StatusCode ResolveStep(int step_number, const BeamLineLoadState* previous,
                        BeamLineLoadState** out, Diagnostic* d) const;
  /// Borrowed immutable model used to validate ownership attachment in Domain.
  const Domain* SourceDomain() const;

 private:
  explicit BeamLineLoadProgram(BeamLineLoadProgramData* data);
  BeamLineLoadProgramData* data_;  ///< Owned opaque definition storage.
};

/// Immutable effective line loads for one step. Owns resolved history and
/// element density vectors; borrows Program/Domain, which must outlive this.
/// Multiple readers are allowed; outputs/diagnostics are caller-exclusive.
class BeamLineLoadState {
 public:
  ~BeamLineLoadState();
  BeamLineLoadState(const BeamLineLoadState&) = delete;
  BeamLineLoadState& operator=(const BeamLineLoadState&) = delete;
  BeamLineLoadState(BeamLineLoadState&&) = delete;
  BeamLineLoadState& operator=(BeamLineLoadState&&) = delete;
  /// Copies current global [qx,qy,qz] F/L into distinct caller-owned Vector(3).
  /// Known elements without current loads receive zeros. Unknown IDs/null or
  /// wrong-sized output fail before mutation; backend failure invalidates out.
  /// Never changes this, Program or Domain. All pointers borrowed for call.
  StatusCode CopyDensity(int element_id, Vector* out, Diagnostic* d) const;
  /// Positive resolved step identity, not evidence of a completed calculation.
  int StepNumber() const;
  /// Borrowed immutable model, valid until State destruction.
  const Domain* SourceDomain() const;

 private:
  friend class BeamLineLoadProgram;
  explicit BeamLineLoadState(BeamLineLoadStateData* data);
  BeamLineLoadStateData* data_;  ///< Owned independent resolved storage.
};

/// Caller-owned reusable L1 scratch. One exclusive instance per serial worker;
/// calls overwrite scratch on success or failure. Never stored by an Element.
class BeamLineLoadWorkspace {
 public:
  /// Preallocates every buffer needed by ComputeBeamLineLoad; no per-call
  /// numerical allocation. Caller deletes success; failure nulls *out and
  /// releases partial storage. out is required and creation-only, d optional.
  static StatusCode Create(BeamLineLoadWorkspace** out, Diagnostic* d);
  ~BeamLineLoadWorkspace();
  BeamLineLoadWorkspace(const BeamLineLoadWorkspace&) = delete;
  BeamLineLoadWorkspace& operator=(const BeamLineLoadWorkspace&) = delete;
  BeamLineLoadWorkspace(BeamLineLoadWorkspace&&) = delete;
  BeamLineLoadWorkspace& operator=(BeamLineLoadWorkspace&&) = delete;

 private:
  friend StatusCode ComputeBeamLineLoad(const Element*, const Vector*,
                                       BeamLineLoadWorkspace*, Vector*, Vector*,
                                       Diagnostic*);
  explicit BeamLineLoadWorkspace(BeamLineLoadWorkspaceData* data);
  BeamLineLoadWorkspaceData* data_;  ///< Owned scratch and coefficient maps.
};

/// L1: constant global density Vector(3), F/L, over the initial element length.
/// Requires a two-node beam with immutable LocalAxes Matrix(3,3). Writes two
/// distinct caller-owned Vector(12)s in [Ua,URa,Ub,URb] order: global equivalent
/// [F,M] and local equivalent [F,M], respectively. Units F and FL; positive
/// global q follows PX/PY/PZ. Element/definitions/results remain unchanged.
/// All pointers borrowed; workspace exclusive. Invalid pointers, dimensions or
/// aliases preserve both outputs; unsupported axes return kUnsupported. Later
/// numerical/backend failure invalidates both outputs and scratch; discard
/// outputs and reuse scratch on the next call. No allocation or thread policy
/// change during evaluation. d is optional and borrowed for this call.
StatusCode ComputeBeamLineLoad(const Element* element,
                               const Vector* global_density,
                               BeamLineLoadWorkspace* workspace,
                               Vector* global_equivalent,
                               Vector* local_equivalent, Diagnostic* d);

}  // namespace fesa

#endif  // FESA_LOADS_BEAM_LINE_LOAD_H_
