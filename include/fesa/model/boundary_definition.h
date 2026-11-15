/// @file
/// @brief Linear-static model ownership and immutable input definitions.
#ifndef FESA_MODEL_BOUNDARY_DEFINITION_H_
#define FESA_MODEL_BOUNDARY_DEFINITION_H_

#include "fesa/base/diagnostic.h"
#include "fesa/model/load_definition.h"

namespace fesa {

/// Immutable original condition. Owns target/source spelling, never a view.
/// DOF is 1..6 in global [U1,U2,U3,UR1,UR2,UR3] order.
/// Load units: force/moment; boundary units: length/radian.
class BoundaryDefinition {
 public:
  /// Copies text and preserves one original definition and source location.
  /// Target node uses node_id; set target uses name (node_id=0).
  /// Caller owns *out on success; failure nulls it. Inputs are borrowed.
  static StatusCode Create(int id, const char* scope, ConditionTarget target,
                           int node_id, const char* name, int dof, double value,
                           const char* source_file, int source_line, int step,
                           ConditionOperation operation,
                           BoundaryDefinition** out,
                           Diagnostic* d);
  ~BoundaryDefinition();
  BoundaryDefinition(const BoundaryDefinition&) = delete;
  BoundaryDefinition& operator=(const BoundaryDefinition&) = delete;
  BoundaryDefinition(BoundaryDefinition&&) = delete;
  BoundaryDefinition& operator=(BoundaryDefinition&&) = delete;
  /// Internal definition ID; never its target node ID.
  int Id() const;
  /// Borrowed interpretation scope key: PART.<name> or ASSEMBLY.<name>.
  /// The suffix retains original spelling; OriginalTarget retains the parsed token.
  const char* Scope() const;
  /// Original node/set target kind.
  ConditionTarget Target() const;
  /// Original node ID, zero for a set target.
  int NodeId() const;
  /// Borrowed original target name; empty for a direct node target.
  const char* Name() const;
  /// Global DOF 1..6.
  int Dof() const;
  /// Original value (including zero); no aggregation is performed here.
  double Value() const;
  /// Borrowed original UTF-8 file path, valid until destruction.
  const char* SourceFile() const;
  /// One-based original source line.
  int SourceLine() const;
  /// Input step index; zero denotes initial model conditions.
  int Step() const;
  /// Original OP; resolution/propagation belongs to S5/M1.
  ConditionOperation Operation() const;

  /// Copies the supplied parsed token exactly, retaining case, qualification
  /// and numeric spelling. The reader trims exterior whitespace only.
  /// Failure preserves the previously stored target. Caller retains text.
  StatusCode SetOriginalTarget(const char* text, Diagnostic* d);
  /// Borrowed parsed token until destruction; empty if not supplied by reader.
  const char* OriginalTarget() const;

 private:
  char* original_target_ = nullptr;  ///< Owned supplied parsed token.
  BoundaryDefinition(int id, char* scope, ConditionTarget target, int node_id,
                   char* name, int dof, double value, char* source_file,
                   int source_line, int step, ConditionOperation operation);
  int id_;
  char* scope_;  ///< Owned original spelling.
  ConditionTarget target_;
  int node_id_;
  char* name_;  ///< Owned original spelling.
  int dof_;
  double value_;
  char* source_file_;  ///< Owned path.
  int source_line_;
  int step_;
  ConditionOperation operation_;
};

}  // namespace fesa

#endif  // FESA_MODEL_BOUNDARY_DEFINITION_H_
