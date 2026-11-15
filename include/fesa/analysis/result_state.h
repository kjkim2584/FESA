/// @file
/// @brief Shared validity of the published mechanical or thermal result set.
#ifndef FESA_ANALYSIS_RESULT_STATE_H_
#define FESA_ANALYSIS_RESULT_STATE_H_

namespace fesa {
/// Common validity of the single Node/Element result set, not buffer contents.
enum class ResultState { kNoResult, kUnconfirmed, kConfirmed };
}  // namespace fesa
#endif  // FESA_ANALYSIS_RESULT_STATE_H_
