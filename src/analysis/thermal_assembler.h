/// @file
/// @brief Thermal upper-COO and full residual assembly with joined element
/// work.
#ifndef FESA_SRC_ANALYSIS_THERMAL_ASSEMBLER_H_
#define FESA_SRC_ANALYSIS_THERMAL_ASSEMBLER_H_
#include <vector>

#include "fesa/analysis/first_order_integrator.h"
#include "fesa/analysis/thermal_model.h"
#include "fesa/analysis/thermal_state.h"
namespace fesa {
/// Owns fixed numerical/index storage for one bound step. Model/map/conditions
/// are borrowed and must outlive this without rebinding. No state pointer is
/// retained. Element/external flux storage retains signs; norms use absolutes.
class ThermalAssembly final {
 public:
  /// Caller owns success; failure nulls *out/releases partial buffers. Requires
  /// model's actual map/conditions. No zero-sized Matrix/Vector is created.
  /// Film/radiation requires CreateForStep with frozen step-start temperatures.
  static StatusCode Create(const ThermalModel* model, const ThermalDofMap* map,
                           const ThermalConditions* conditions,
                           ThermalAssembly** out, Diagnostic* d);
  /// As Create, additionally cloning frozen step-start temperatures for
  /// film/radiation. Borrows the full nominal Vector during call only; analysis
  /// guarantees its mapping/time. Null is valid only without film/radiation.
  /// Failure nulls out.
  static StatusCode CreateForStep(const ThermalModel* model,
                                  const ThermalDofMap* map,
                                  const ThermalConditions* conditions,
                                  const Vector* step_start,
                                  ThermalAssembly** out, Diagnostic* d);
  ~ThermalAssembly();
  ThermalAssembly(const ThermalAssembly&) = delete;
  ThermalAssembly& operator=(const ThermalAssembly&) = delete;
  ThermalAssembly(ThermalAssembly&&) = delete;
  ThermalAssembly& operator=(ThermalAssembly&&) = delete;
  bool IsValid()
      const;  ///< Successful latest assembly; false after numeric failure.
  const ThermalDofMap* DofMap() const;
  const std::vector<int>* Rows() const;
  const std::vector<int>* Columns() const;
  const Vector* Values() const;  ///< Upper COO, null when no reduced entries.
  const Vector* Residual()
      const;  ///< Q_eff-f, full nominal global order [E/T].
  const Vector* Reaction() const;  ///< f-Q on constrained rows, zero elsewhere.
  /// Borrowed unassembled internal f_e(including capacity), followed by full
  /// signed -Q per distributed original-row/element, before inactive masking.
  const Vector* ElementFluxes() const;
  const Vector* ExternalFluxes()
      const;  ///< Signed CFLUX-only per-original-row Q_eff, or null.
  const Vector* FreeResidual()
      const;  ///< Reduced equation order, null if none.
  const Vector* ParticipatingIncrement()
      const;  ///< theta-theta_n, includes constraints.
 private:
  friend class ThermalAssembler;
  ThermalAssembly() = default;
  const ThermalModel* model_ = nullptr;            ///< Borrowed.
  const ThermalDofMap* map_ = nullptr;             ///< Borrowed.
  const ThermalConditions* conditions_ = nullptr;  ///< Borrowed.
  bool valid_ = false;
  std::vector<int>* rows_ = nullptr;     ///< Owned upper-COO rows.
  std::vector<int>* columns_ = nullptr;  ///< Owned upper-COO columns.
  std::vector<int>* matrix_indices_ =
      nullptr;  ///< Owned local flattened upper entries.
  std::vector<int>* globals_ =
      nullptr;  ///< Owned element-local -> global packed IDs.
  std::vector<int>* element_offsets_ =
      nullptr;                               ///< Owned element flux prefixes.
  std::vector<int>* coo_offsets_ = nullptr;  ///< Owned upper-COO prefixes.
  /// Owned row-flux prefixes and per-element CSR routing; immutable after
  /// Create. Distributed rows retain full DOFs, including unloaded zeros.
  std::vector<int>* load_flux_offsets_ = nullptr;
  std::vector<int>* load_starts_ = nullptr;
  std::vector<int>* load_indices_ = nullptr;
  std::vector<int>* external_globals_ =
      nullptr;  ///< Owned original-row target slots.
  std::vector<int>* free_globals_ = nullptr;  ///< Owned equation-order slots.
  std::vector<int>* constrained_globals_ =
      nullptr;  ///< Owned constrained slots.
  std::vector<int>* participating_globals_ =
      nullptr;                      ///< Owned participating slots.
  Vector* values_ = nullptr;        ///< Owned optional upper values.
  Vector* residual_ = nullptr;      ///< Owned.
  Vector* reaction_ = nullptr;      ///< Owned.
  /// Owned packed internal plus signed negative distributed-row fluxes.
  Vector* element_flux_ = nullptr;
  Vector* external_flux_ =
      nullptr;  ///< Owned optional original-row external flux.
  Vector* free_residual_ = nullptr;        ///< Owned optional reduced residual.
  Vector* increment_ = nullptr;            ///< Owned participating increment.
  Vector* global_increment_ = nullptr;     ///< Owned scratch.
  Vector* constrained_scratch_ = nullptr;  ///< Owned optional reaction scratch.
  Vector* step_start_ = nullptr;  ///< Owned exchange-only frozen step-start T.
  Vector* zero_ = nullptr;                 ///< Owned immutable global zero.
};
/// No persistent resources. Borrows the integrator once per call; workers use
/// its validated coefficient. Numerical buffers are worker-local and released
/// before joining. Upper packing is disjoint; global sums use Vector
/// ScatterAdd.
class ThermalAssembler final {
 public:
  /// State candidate must already contain end-time prescribed temperatures.
  /// Borrows all inputs for the call. steady iff integrator=null/dt0; transient
  /// requires integrator and finite dt>0, matching the bound step. Calls its
  /// coefficient once. Capacity is diagonal lumped, no C_fc feature is claimed.
  /// Invalid preflight preserves out. Later failure invalidates out; all state
  /// and histories stay unchanged. Model/conditions must match out's binding.
  /// extrapolated_start is true only for an actually extrapolated iteration0
  /// evaluation (D13). Only free DOFs may then have a negative trial absolute
  /// temperature; every other call uses false. No flag/mask is retained.
  static StatusCode Assemble(const ThermalModel* model,
                             const ThermalDofMap* map,
                             const ThermalConditions* conditions,
                             const ThermalState* state, double step_time,
                             double period,
                             const FirstOrderTimeIntegrator* integrator,
                             double dt, ThermalAssembly* out, Diagnostic* d,
                             bool extrapolated_start = false);

 private:
  struct Context;
  static StatusCode ElementTask(int item, int worker, void* context,
                                Diagnostic* d);
};
}  // namespace fesa
#endif  // FESA_SRC_ANALYSIS_THERMAL_ASSEMBLER_H_
