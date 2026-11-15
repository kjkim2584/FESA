/// @file
/// @brief Internal thread-local observer for solid purpose-isolation tests.
#ifndef FESA_ELEMENTS_SOLID_ACCESS_H_
#define FESA_ELEMENTS_SOLID_ACCESS_H_
namespace fesa {
/// Logical numerical boundaries; never a production performance report.
enum class SolidEvaluationEvent {
  kGeometry,
  kJacobianFactor,
  kConstitutive,
  kStiffness,
  kMass,
  kForce,
  kRecovery,
  kRecoveryFactor,
  kHourglass,
  kMaterialTangent
};
/// Optional thread-local test callback; caller clears before its state expires.
/// Set only outside dispatch. Normal evaluation does not allocate/log here.
void SetSolidEvaluationObserver(void (*observer)(SolidEvaluationEvent));
/// Emits an internal numerical boundary to the same thread-local observer.
/// Borrows no state; callback lifetime follows SetSolidEvaluationObserver.
void ObserveSolidEvaluation(SolidEvaluationEvent event);
}  // namespace fesa
#endif  // FESA_ELEMENTS_SOLID_ACCESS_H_
