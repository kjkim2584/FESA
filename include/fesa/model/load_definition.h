/// @file
/// @brief Linear-static model ownership and immutable input definitions.
#ifndef FESA_MODEL_LOAD_DEFINITION_H_
#define FESA_MODEL_LOAD_DEFINITION_H_

#include "fesa/base/diagnostic.h"
#include "fesa/model/amplitude.h"

namespace fesa {

/// Original input target kind; never replace set identity with expanded IDs.
enum class ConditionTarget { kNode, kNodeSet };
/// Abaqus condition update operation, preserved before step resolution.
enum class ConditionOperation { kMod, kNew };

/// Immutable original condition. Owns target/source spelling, never a view.
/// DOF is 1..6 in global [U1,U2,U3,UR1,UR2,UR3] order.
/// Load units: force/moment; boundary units: length/radian.
class LoadDefinition {
 public:
  /// Copies text and preserves one original definition and source location.
  /// Target node uses node_id; set target uses name (node_id=0).
  /// Caller owns *out on success; failure nulls it. Inputs are borrowed.
  static StatusCode Create(int id, const char* scope, ConditionTarget target,
                           int node_id, const char* name, int dof, double value,
                           const char* source_file, int source_line, int step,
                           ConditionOperation operation, LoadDefinition** out,
                           Diagnostic* d);
  ~LoadDefinition();
  LoadDefinition(const LoadDefinition&) = delete;
  LoadDefinition& operator=(const LoadDefinition&) = delete;
  LoadDefinition(LoadDefinition&&) = delete;
  LoadDefinition& operator=(LoadDefinition&&) = delete;
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
  /// Binds once before Domain publication, borrowing a nonnull amplitude from
  /// the same Domain. It must outlive every use of this definition. No
  /// ownership transfer; failure preserves the binding. Only dynamic CLOAD may
  /// bind.
  StatusCode BindDynamicAmplitude(const Amplitude* amplitude, Diagnostic* d);
  /// Borrowed immutable relative amplitude; null means an unmodulated load.
  const Amplitude* DynamicAmplitude() const;

 private:
  const Amplitude* amplitude_ = nullptr;  ///< Borrowed from owning Domain.
  char* original_target_ = nullptr;  ///< Owned supplied parsed token.
  LoadDefinition(int id, char* scope, ConditionTarget target, int node_id,
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

#endif  // FESA_MODEL_LOAD_DEFINITION_H_
