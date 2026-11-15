#include <cstdio>
#include <cstring>
#include <memory>
#include <new>
#include <string>

#include "adapters/abaqus/input_syntax.h"
#include "fesa/io/analysis_log.h"
#include "fesa/io/run_options.h"
#include "fesa/solver/fesa_solver.h"

namespace fesa {
RunOptions::~RunOptions() {
  delete input_path;
  delete config_path;
  delete executable_directory;
}
namespace {
constexpr char kUsage[] =
    "FESA.exe <input.inp> [--threads N] [--config <path>] "
    "[--b33-material-first]\n"
    "FESA.exe --help\n";
StatusCode Usage(Diagnostic* d) {
  return ReportFailure(StatusCode::kInvalidArgument,
                       DiagnosticStage::kValidation, "Invalid command line", d);
}
int Failure(const Diagnostic* d) {
  bool reported = false;
  try {
    std::string text("");
    if (FormatDiagnostic(d, &text) == StatusCode::kOk) {
      std::fputs(text.c_str(), stderr);
      reported = true;
    }
  } catch (const std::bad_alloc&) {
    // Reporting allocation failure must remain safe under sustained pressure.
  }
  if (!reported)
    std::fprintf(stderr, "FESA: %s (diagnostic formatting failed)\n",
                 d->operation ? d->operation : "runtime failure");
  return 1;
}
}  // namespace
StatusCode ParseRunOptions(int argc, const char* const* argv,
                           int available_cpus, RunOptions** out,
                           Diagnostic* d) {
  if (!out) return Usage(d);
  *out = nullptr;
  if (!argv || argc < 2) return Usage(d);
  for (int i = 0; i < argc; ++i)
    if (!argv[i]) return Usage(d);
  try {
    std::unique_ptr<RunOptions> options(new RunOptions);
    if (argc == 2 && std::strcmp(argv[1], "--help") == 0) {
      options->help = true;
      *out = options.release();
      return StatusCode::kOk;
    }
    if (available_cpus < 1) return Usage(d);
    if (!*argv[1] || argv[1][0] == '-') return Usage(d);
    options->input_path = new std::string(argv[1]);
    bool threads_seen = false;
    for (int i = 2; i < argc; ++i) {
      if (std::strcmp(argv[i], "--b33-material-first") == 0) {
        if (options->b33_material_first) return Usage(d);
        options->b33_material_first = true;
        continue;
      }
      const bool threads = std::strcmp(argv[i], "--threads") == 0;
      const bool config = std::strcmp(argv[i], "--config") == 0;
      if ((!threads && !config) || (threads && threads_seen) ||
          (config && options->config_path) || i + 1 >= argc || !*argv[i + 1] ||
          std::strncmp(argv[i + 1], "--", 2) == 0)
        return Usage(d);
      ++i;
      if (threads) {
        threads_seen = true;
        if (!abaqus::Integer(argv[i], 1, available_cpus, &options->threads))
          return Usage(d);
      } else {
        options->config_path = new std::string(argv[i]);
      }
    }
    *out = options.release();
    return StatusCode::kOk;
  } catch (const std::bad_alloc&) {
    return ReportFailure(StatusCode::kAllocationFailure,
                         DiagnosticStage::kAllocation,
                         "ParseRunOptions allocation", d);
  }
}
int RunCommandLine(int argc, const char* const* argv, int available_cpus,
                   const char* executable_directory) {
  Diagnostic d;
  std::unique_ptr<RunOptions> options;
  std::unique_ptr<FESASolver> solver;
  try {
    RunOptions* parsed = nullptr;
    auto status = ParseRunOptions(argc, argv, available_cpus, &parsed, &d);
    options.reset(parsed);
    if (status == StatusCode::kInvalidArgument) {
      std::fputs(kUsage, stderr);
      std::fprintf(stderr, "--threads N must be an integer in [1, %d].\n",
                   available_cpus);
      return 2;
    }
    if (status != StatusCode::kOk) return Failure(&d);
    if (options->help) {
      std::fputs(kUsage, stdout);
      return 0;
    }
    options->available_cpus = available_cpus;
    if (executable_directory)
      options->executable_directory = new std::string(executable_directory);
    FESASolver* created = nullptr;
    status = FESASolver::Create(&created, &d);
    solver.reset(created);
    if (status != StatusCode::kOk) return Failure(&d);
    return solver->Run(options.get(), &d) == StatusCode::kOk ? 0 : 1;
  } catch (const std::bad_alloc&) {
    d.Record(StatusCode::kAllocationFailure, DiagnosticStage::kAllocation,
             "CLI execution allocation");
    return Failure(&d);
  }
}
}  // namespace fesa
