/// @file
/// @brief Internal analysis observation for adapter-failure integration tests.
#ifndef FESA_ANALYSIS_ANALYSIS_ACCESS_H_
#define FESA_ANALYSIS_ANALYSIS_ACCESS_H_
namespace fesa {
/// Observation boundaries before numerical work; no public solver options.
enum class AnalysisStage { kNodeUpdate };
/// Borrows a nonthrowing thread-local callback until reset by the caller.
/// Callbacks may arm adapter fault seams and must not mutate model objects.
void SetAnalysisObserver(void (*callback)(AnalysisStage));
}  // namespace fesa
#endif  // FESA_ANALYSIS_ANALYSIS_ACCESS_H_
