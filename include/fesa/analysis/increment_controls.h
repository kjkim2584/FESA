/// @file
/// @brief Automatic static increment factors/counts, separate from Newton.
#ifndef FESA_ANALYSIS_INCREMENT_CONTROLS_H_
#define FESA_ANALYSIS_INCREMENT_CONTROLS_H_
namespace fesa {
/// Abaqus public defaults with documented first-FESA supported coefficient
/// domain. Counts are positive except IL/IG may be0. Cutback factors(0,1),
/// success_scale(0,1], growth/maximum_growth/endpoint_stretch>=1, cutoff>=0.
/// All factors finite. Borrowed immutable within the step.
/// @see docs/MITC4/increment_execution_contract.md.
struct IncrementControls {
  int attempt_limit = 5;             ///< IA including first attempt.
  int slow_success_iterations = 10;  ///< IL; shrink if iterations>IL.
  int quick_iterations = 4;         ///< IG; two consecutive eligible successes.
  double divergence_cutback = .25;  ///< Df.
  double slow_cutback = .5;         ///< DC, also ISc cap.
  double geometry_cutback = .25;    ///< DH, also rotation branch.
  double success_scale = .75;       ///< DB after slow successful increment.
  double growth = 1.5;              ///< DD, capped by DM.
  double maximum_growth = 1.5;      ///< DM.
  double extrapolation_cutoff = .1;  ///< DE, skip predictor at ratio<=DE.
  double endpoint_stretch = 1;       ///< DT, no endpoint stretching by default.
};
}  // namespace fesa
#endif  // FESA_ANALYSIS_INCREMENT_CONTROLS_H_
