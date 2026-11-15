/// @file
/// @brief Original scalar settings for a direct HHT step.
#ifndef FESA_MODEL_DYNAMIC_STEP_PARAMETERS_H_
#define FESA_MODEL_DYNAMIC_STEP_PARAMETERS_H_
namespace fesa {
/// Call-scoped input; settings factory copies individual scalars.
/// Time quantities use model T units; alpha is dimensionless.
/// Zero min/max means unspecified to DIRECT; flags retain explicit input 0.
struct DynamicStepParameters {
  double initial_increment = 0;    ///< Positive requested h, before clipping.
  double period = 0;               ///< Positive step duration.
  double minimum = 0;              ///< Nonnegative requested lower limit.
  double maximum = 0;              ///< Nonnegative requested upper limit.
  double alpha = -.05;             ///< HHT alpha in [-1/3,0].
  int increment_limit = 100;       ///< Successful advancing increments only.
  bool minimum_specified = false;  ///< Nonempty numeric field, including 0.
  bool maximum_specified = false;  ///< Nonempty numeric field, including 0.
  bool alpha_specified = false;    ///< Explicit DYNAMIC option.
  bool increment_limit_specified = false;  ///< Explicit STEP option.
};
}  // namespace fesa
#endif  // FESA_MODEL_DYNAMIC_STEP_PARAMETERS_H_
