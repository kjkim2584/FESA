/// @file
/// @brief CLI syntax and preparation boundary; paths are UTF-8.
#ifndef FESA_IO_RUN_OPTIONS_H_
#define FESA_IO_RUN_OPTIONS_H_

#include <string>

#include "fesa/base/diagnostic.h"

namespace fesa {
/// Owns copied command arguments. No environment variable changes the request.
struct RunOptions {
  RunOptions() = default;
  ~RunOptions();
  RunOptions(const RunOptions&) = delete;
  RunOptions& operator=(const RunOptions&) = delete;
  RunOptions(RunOptions&&) = delete;
  RunOptions& operator=(RunOptions&&) = delete;
  std::string* input_path = nullptr;   ///< Owned; null for help.
  std::string* config_path = nullptr;  ///< Owned; null means exe/solver.cfg.
  /// Owned optional executable directory; required only when no INP CONTROLS
  /// and no explicit config_path. Direct Run never substitutes cwd.
  std::string* executable_directory = nullptr;
  int available_cpus = 0;  ///< CLI validation context; zero for direct Run.
  int threads = 1;         ///< Requested CPU count, not observation.
  bool help = false;       ///< Help has no file side effects.
  /// Opt-in pure B33 nonlinear static first-correction material Jacobian.
  bool b33_material_first = false;
};
/// Parses borrowed argv without opening files. Caller deletes *out on success;
/// failure returns null. kInvalidArgument denotes usage (exit 2); allocation
/// failure is a runtime error (exit 1). Success preserves d.
StatusCode ParseRunOptions(int argc, const char* const* argv,
                           int available_cpus, RunOptions** out, Diagnostic* d);
/// Parses/help/usage then invokes one FESASolver Run. Returns 0 for success
/// or help, 2 for usage and 1 for runtime failure; borrows all arguments.
int RunCommandLine(int argc, const char* const* argv, int available_cpus,
                   const char* executable_directory);
}  // namespace fesa
#endif  // FESA_IO_RUN_OPTIONS_H_
