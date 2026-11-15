#define NOMINMAX
#include <Windows.h>

#include <cstdio>
#include <cwchar>
#include <memory>
#include <new>
#include <string>
#include <vector>

#include "adapters/openmp/thread_runtime.h"
#include "fesa/io/run_options.h"

namespace {
bool Utf8(const wchar_t* source, std::string* out) {
  const int size = WideCharToMultiByte(CP_UTF8, WC_ERR_INVALID_CHARS, source,
                                       -1, nullptr, 0, nullptr, nullptr);
  if (!size) return false;
  out->resize(size);
  if (!WideCharToMultiByte(CP_UTF8, WC_ERR_INVALID_CHARS, source, -1,
                           out->data(), size, nullptr, nullptr))
    return false;
  out->pop_back();
  return true;
}
}  // namespace
// Windows owns wide argv; no legacy code-page conversion touches user paths.
int wmain(int argc, wchar_t** argv) {
  try {
    std::vector<std::string> args(0);
    std::vector<const char*> pointers(0);
    args.reserve(argc);
    pointers.reserve(argc);
    for (int i = 0; i < argc; ++i) {
      std::string text("");
      if (!Utf8(argv[i], &text)) {
        std::fputs("FESA: invalid Unicode command argument\n", stderr);
        return 1;
      }
      args.push_back(text);
    }
    for (const auto& arg : args) pointers.push_back(arg.c_str());
    // Help does not depend on CPU discovery or any input/config/log path.
    if (argc == 2 && std::wcscmp(argv[1], L"--help") == 0)
      return fesa::RunCommandLine(argc, pointers.data(), 0, "");
    std::unique_ptr<wchar_t[]> filename(new wchar_t[32768]);
    const DWORD length = GetModuleFileNameW(nullptr, filename.get(), 32768);
    if (!length || length >= 32768) {
      std::fputs("FESA: executable path unavailable\n", stderr);
      return 1;
    }
    std::wstring directory(filename.get());
    const auto slash = directory.find_last_of(L"/\\");
    if (slash == std::wstring::npos) return 1;
    directory.erase(slash);
    std::string exe("");
    if (!Utf8(directory.c_str(), &exe)) return 1;
    int cpus = 0;
    fesa::Diagnostic cpu_diagnostic;
    if (fesa::QueryProcessAvailableCpus(&cpus, &cpu_diagnostic) !=
        fesa::StatusCode::kOk) {
      std::fprintf(
          stderr,
          "FESA: process CPU affinity mask unavailable or empty; native=%d\n",
          cpu_diagnostic.vendor_code);
      return 1;
    }
    return fesa::RunCommandLine(argc, pointers.data(), cpus, exe.c_str());
  } catch (const std::bad_alloc&) {
    std::fputs("FESA: command argument allocation failed\n", stderr);
    return 1;
  }
}
