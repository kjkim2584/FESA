/// @file
/// @brief Linear-static model ownership and immutable input definitions.
#ifndef FESA_MODEL_STEP_DEFINITION_H_
#define FESA_MODEL_STEP_DEFINITION_H_

#include <vector>

#include "fesa/model/boundary_definition.h"
#include "fesa/model/heat_transfer_step.h"
#include "fesa/model/load_definition.h"
#include "fesa/model/nonlinear_static_parameters.h"
#include "fesa/model/thermal_condition.h"

namespace fesa {

/// Immutable step identity/time with non-owning condition references.
/// Domain owns conditions and must outlive this step; append failure preserves
/// references already stored. Time units are the input model's convention.
class StepDefinition {
 public:
  /// Copies name. Caller owns *out; on failure *out is null.
  static StatusCode Create(int id, const char* name, double period,
                           StepDefinition** out, Diagnostic* d);
  /// Creates a thermal-only identity with owned explicit scalar settings.
  /// Caller owns success; failure nulls *out and cleans partial storage.
  /// Mechanical collections stay empty; detailed thermal conditions borrowed.
  static StatusCode CreateHeatTransfer(
      int id, const char* name, const HeatTransferStepParameters* parameters,
      StepDefinition** out, Diagnostic* d);
  /// True only for CreateHeatTransfer; does not imply an executable solver.
  bool IsHeatTransfer() const;
  /// Borrowed immutable owned settings, null for every mechanical factory.
  const HeatTransferStepParameters* HeatTransfer() const;
  /// Appends a borrowed same-step definition; Domain must own it before use.
  /// Failure preserves prior references and ownership; no condition evaluation.
  StatusCode AddThermalCondition(const ThermalConditionDefinition* condition,
                                  Diagnostic* d);
  /// Borrowed collection until destruction, null for a mechanical step.
  const std::vector<const ThermalConditionDefinition*>*
  ThermalConditions() const;
  /// Creates a linear dynamic identity; detailed integration settings are
  /// solver-owned. Same ownership/failure contract as Create. No execution.
  static StatusCode CreateDynamic(int id, const char* name, double period,
                                  StepDefinition** out, Diagnostic* d);
  /// Creates a finite-geometry dynamic identity with no static parameters.
  /// Same ownership/validation as Create; detailed dynamic settings remain
  /// solver-owned. Does not enable any element/CLI capability by itself.
  static StatusCode CreateNonlinearDynamic(int id, const char* name,
                                           double period, StepDefinition** out,
                                           Diagnostic* d);
  /// Procedure-independent finite geometry identity. True for nonlinear
  /// static/dynamic; false for the existing linear factories. No state owned.
  bool HasFiniteGeometry() const;
  /// True for a dynamic procedure; false for either static factory.
  bool IsDynamic() const;
  /// Owns a scalar copy of nonlinear parameters; caller may destroy input.
  /// No input mutation; failure nulls output and releases partial storage.
  static StatusCode CreateNonlinearStatic(
      int id, const char* name, const NonlinearStaticParameters* parameters,
      StepDefinition** out, Diagnostic* d);
  ~StepDefinition();
  StepDefinition(const StepDefinition&) = delete;
  StepDefinition& operator=(const StepDefinition&) = delete;
  StepDefinition(StepDefinition&&) = delete;
  StepDefinition& operator=(StepDefinition&&) = delete;
  /// Appends a borrowed definition; Domain must own it before step execution.
  /// Thermal steps reject mechanical definitions without mutation.
  StatusCode AddLoad(const LoadDefinition* load, Diagnostic* d);
  /// Appends a borrowed definition; failure preserves existing references.
  /// Thermal steps reject mechanical definitions without mutation.
  StatusCode AddBoundary(const BoundaryDefinition* boundary, Diagnostic* d);
  /// Internal positive step ID.
  int Id() const;
  /// Borrowed original name.
  const char* Name() const;
  /// Positive step period.
  double Period() const;
  /// Borrowed immutable settings until destruction; null for linear Create.
  const NonlinearStaticParameters* NonlinearParameters() const;
  /// Borrowed reference collections, valid until step destruction.
  const std::vector<const LoadDefinition*>* Loads() const;
  /// Borrowed reference collection; no conditions are owned by this step.
  const std::vector<const BoundaryDefinition*>* Boundaries() const;

  /// Records block presence/operation even for an empty NEW block.
  /// Thermal M1 uses load=true for CFLUX, false for temperature boundaries.
  /// Same-kind blocks must agree; failure preserves metadata. line is original.
  StatusCode RecordConditionBlock(bool load, ConditionOperation operation,
                                  int line, Diagnostic* d);
  /// Whether the step explicitly contains this condition kind.
  bool HasLoadBlock() const;
  /// Whether a boundary block is present, independently of row count.
  bool HasBoundaryBlock() const;
  /// Operation for the present kind; MOD is returned when absent.
  ConditionOperation LoadOperation() const;
  /// Boundary operation; MOD when absent.
  ConditionOperation BoundaryOperation() const;
  /// First block source lines; zero when absent.
  int LoadBlockLine() const;
  /// First boundary block source line; zero when absent.
  int BoundaryBlockLine() const;
  /// Records an independent thermal keyword family's block, including empty
  /// NEW. Requires a thermal step and a valid family/op/line; failure preserves
  /// metadata. No borrowed inputs survive the call.
  StatusCode RecordThermalLoadBlock(ThermalLoadFamily family,
                                    ConditionOperation operation, int line,
                                    Diagnostic* d);
  /// Absent family is MOD. Invalid family also returns MOD without mutation.
  ConditionOperation ThermalLoadOperation(ThermalLoadFamily family) const;

 private:
  bool dynamic_ = false;
  HeatTransferStepParameters* heat_transfer_ = nullptr;  ///< Owned or null.
  /// Owned collection with borrowed Domain condition entries.
  std::vector<const ThermalConditionDefinition*>* thermal_conditions_ = nullptr;
  bool finite_dynamic_ = false;  ///< Dynamic geometry; no static settings.
  int load_block_line_ = 0;
  int boundary_block_line_ = 0;
  int thermal_load_lines_[6] = {};  ///< First source line per keyword family.
  ConditionOperation thermal_load_operations_[6] = {};  ///< Scalar metadata.
  ConditionOperation load_operation_ = ConditionOperation::kMod;
  ConditionOperation boundary_operation_ = ConditionOperation::kMod;
  StepDefinition(int id, char* name, double period,
                 std::vector<const LoadDefinition*>* loads,
                 std::vector<const BoundaryDefinition*>* boundaries);
  int id_;
  char* name_;  ///< Owned.
  double period_;
  NonlinearStaticParameters* nonlinear_ = nullptr;  ///< Owned, optional.
  /// Owned list, borrowed entries.
  std::vector<const LoadDefinition*>* loads_;
  std::vector<const BoundaryDefinition*>* boundaries_;  ///< Same ownership.
};

}  // namespace fesa

#endif  // FESA_MODEL_STEP_DEFINITION_H_
