/// @file
/// @brief Input-only physical constraint history before nonlinear execution.
#ifndef FESA_ANALYSIS_NONLINEAR_CONDITION_VALIDATOR_H_
#define FESA_ANALYSIS_NONLINEAR_CONDITION_VALIDATOR_H_

#include <map>
#include <utility>

#include "fesa/analysis/effective_conditions.h"
#include "fesa/math/vector.h"

namespace fesa {
/// Owns input constraint/axis history; borrows immutable Domain until
/// destruction. Consumes existing preflight snapshots, never rebuilds/evaluates
/// the model.
class NonlinearConditionValidator final {
 public:
  /// Caller owns success; failure nulls output and releases partial storage.
  static StatusCode Create(const Domain* domain,
                           NonlinearConditionValidator** out, Diagnostic* d);
  ~NonlinearConditionValidator();
  NonlinearConditionValidator(const NonlinearConditionValidator&) = delete;
  NonlinearConditionValidator& operator=(const NonlinearConditionValidator&) =
      delete;
  NonlinearConditionValidator(NonlinearConditionValidator&&) = delete;
  NonlinearConditionValidator& operator=(NonlinearConditionValidator&&) =
      delete;
  /// Borrows same-Domain nonlinear snapshots in order starting at step1.
  /// Success advances only owned input history; any failure preserves it and
  /// reports original definition context. No physical state/results mutation.
  StatusCode ValidateStep(const EffectiveConditions* conditions, Diagnostic* d);

 private:
  NonlinearConditionValidator() = default;
  const Domain* domain_ = nullptr;  ///< Borrowed, must outlive this object.
  /// Owned physical key states: missing=never, 1=active, 2=released.
  std::map<std::pair<int, int>, int>* states_ = nullptr;
  std::map<int, Vector*>* axes_ = nullptr;  ///< Owned map and unit axes.
  bool beam_ = false;  ///< Exclusive actual Beam static role permits any axis.
  int next_step_ = 1;
};
}  // namespace fesa
#endif  // FESA_ANALYSIS_NONLINEAR_CONDITION_VALIDATOR_H_
