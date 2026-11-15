#include "io/file_access.h"

#define NOMINMAX
#include <Windows.h>
#include <fcntl.h>
#include <io.h>

#include <cerrno>
#include <cstdint>
#include <memory>

namespace fesa {
namespace io {
namespace {
StatusCode Failure(int native, Diagnostic* d) {
  if (d)
    d->Record(StatusCode::kIoFailure, DiagnosticStage::kValidation,
              "Open UTF-8 file", native);
  return StatusCode::kIoFailure;
}
StatusCode WidePath(const char* path, wchar_t** out, Diagnostic* d) {
  const int size =
      MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, path, -1, nullptr, 0);
  if (!size) return Failure(static_cast<int>(GetLastError()), d);
  std::unique_ptr<wchar_t[]> wide(new wchar_t[size]);
  if (!MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, path, -1, wide.get(),
                           size))
    return Failure(static_cast<int>(GetLastError()), d);
  *out = wide.release();
  return StatusCode::kOk;
}
}  // namespace
StatusCode OpenFile(const char* path, FileMode mode, FILE** out,
                    Diagnostic* d) {
  *out = nullptr;
  wchar_t* converted = nullptr;
  const auto status = WidePath(path, &converted, d);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<wchar_t[]> wide(converted);
  if (mode != FileMode::kCreateNew) {
    const errno_t error =
        _wfopen_s(out, wide.get(), mode == FileMode::kRead ? L"rb" : L"wb");
    return error ? Failure(error, d) : StatusCode::kOk;
  }
  HANDLE handle = CreateFileW(wide.get(), GENERIC_WRITE, 0, nullptr, CREATE_NEW,
                              FILE_ATTRIBUTE_NORMAL, nullptr);
  if (handle == INVALID_HANDLE_VALUE)
    return Failure(static_cast<int>(GetLastError()), d);
  const int fd = _open_osfhandle(reinterpret_cast<intptr_t>(handle),
                                 _O_WRONLY | _O_BINARY);
  if (fd == -1) {
    const int error = errno;
    CloseHandle(handle);
    return Failure(error, d);
  }
  *out = _fdopen(fd, "wb");
  if (*out) return StatusCode::kOk;
  const int error = errno;
  _close(fd);
  return Failure(error, d);
}
StatusCode InspectPath(const char* path, bool* exists, bool* regular,
                       Diagnostic* d) {
  wchar_t* converted = nullptr;
  const auto status = WidePath(path, &converted, d);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<wchar_t[]> wide(converted);
  const DWORD attributes = GetFileAttributesW(wide.get());
  if (attributes == INVALID_FILE_ATTRIBUTES) {
    const auto error = GetLastError();
    if (error != ERROR_FILE_NOT_FOUND && error != ERROR_PATH_NOT_FOUND)
      return Failure(static_cast<int>(error), d);
    *exists = false;
    *regular = false;
  } else {
    *exists = true;
    *regular =
        (attributes & (FILE_ATTRIBUTE_DIRECTORY | FILE_ATTRIBUTE_DEVICE)) == 0;
  }
  return StatusCode::kOk;
}
StatusCode AbsolutePath(const char* input, std::string* out, Diagnostic* d) {
  wchar_t* converted = nullptr;
  auto status = WidePath(input, &converted, d);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<wchar_t[]> wide(converted);
  const DWORD count = GetFullPathNameW(wide.get(), 0, nullptr, nullptr);
  if (!count) return Failure(static_cast<int>(GetLastError()), d);
  std::unique_ptr<wchar_t[]> absolute(new wchar_t[count]);
  const DWORD written =
      GetFullPathNameW(wide.get(), count, absolute.get(), nullptr);
  if (!written || written >= count)
    return Failure(static_cast<int>(GetLastError()), d);
  const int bytes =
      WideCharToMultiByte(CP_UTF8, WC_ERR_INVALID_CHARS, absolute.get(), -1,
                          nullptr, 0, nullptr, nullptr);
  if (!bytes) return Failure(static_cast<int>(GetLastError()), d);
  std::unique_ptr<char[]> utf8(new char[bytes]);
  if (!WideCharToMultiByte(CP_UTF8, WC_ERR_INVALID_CHARS, absolute.get(), -1,
                           utf8.get(), bytes, nullptr, nullptr))
    return Failure(static_cast<int>(GetLastError()), d);
  out->assign(utf8.get());
  return StatusCode::kOk;
}
namespace {
StatusCode ReadIdentity(const wchar_t* path, bool* exists,
                        BY_HANDLE_FILE_INFORMATION* identity, Diagnostic* d) {
  HANDLE handle =
      CreateFileW(path, FILE_READ_ATTRIBUTES,
                  FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE,
                  nullptr, OPEN_EXISTING, FILE_FLAG_BACKUP_SEMANTICS, nullptr);
  if (handle == INVALID_HANDLE_VALUE) {
    const auto error = GetLastError();
    if (error != ERROR_FILE_NOT_FOUND && error != ERROR_PATH_NOT_FOUND)
      return Failure(static_cast<int>(error), d);
    *exists = false;
    return StatusCode::kOk;
  }
  const BOOL read = GetFileInformationByHandle(handle, identity);
  const auto error = read ? ERROR_SUCCESS : GetLastError();
  const BOOL closed = CloseHandle(handle);
  if (!read) return Failure(static_cast<int>(error), d);
  if (!closed) return Failure(static_cast<int>(GetLastError()), d);
  *exists = true;
  return StatusCode::kOk;
}
}  // namespace
StatusCode ValidateDistinctPaths(const char* source, const char* destination,
                                 Diagnostic* d) {
  std::string first(""), second("");
  auto status = AbsolutePath(source, &first, d);
  if (status != StatusCode::kOk) return status;
  status = AbsolutePath(destination, &second, d);
  if (status != StatusCode::kOk) return status;
  wchar_t* raw = nullptr;
  status = WidePath(first.c_str(), &raw, d);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<wchar_t[]> first_wide(raw);
  status = WidePath(second.c_str(), &raw, d);
  if (status != StatusCode::kOk) return status;
  std::unique_ptr<wchar_t[]> second_wide(raw);
  const int compared =
      CompareStringOrdinal(first_wide.get(), -1, second_wide.get(), -1, TRUE);
  if (!compared) return Failure(static_cast<int>(GetLastError()), d);
  bool same = compared == CSTR_EQUAL;
  if (!same) {
    bool first_exists = false, second_exists = false;
    BY_HANDLE_FILE_INFORMATION first_id = {}, second_id = {};
    status = ReadIdentity(first_wide.get(), &first_exists, &first_id, d);
    if (status != StatusCode::kOk) return status;
    status = ReadIdentity(second_wide.get(), &second_exists, &second_id, d);
    if (status != StatusCode::kOk) return status;
    same = first_exists && second_exists &&
           first_id.dwVolumeSerialNumber == second_id.dwVolumeSerialNumber &&
           first_id.nFileIndexHigh == second_id.nFileIndexHigh &&
           first_id.nFileIndexLow == second_id.nFileIndexLow;
  }
  return same ? ReportFailure(StatusCode::kIoFailure,
                              DiagnosticStage::kValidation,
                              "Source/destination path collision", d)
              : StatusCode::kOk;
}
void ResultPath(const char* input, std::string* out) {
  AnalysisLogPath(input, out);
  out->erase(out->size() - sizeof("_analysis.log") + 1);
  out->append(".h5");
}
void DefaultConfigPath(const char* directory, std::string* out) {
  out->assign(directory);
  if (!out->empty() && out->back() != '/' && out->back() != '\\')
    out->append("\\");
  out->append("solver.cfg");
}
void AnalysisLogPath(const char* input, std::string* out) {
  out->assign(input);
  const auto slash = out->find_last_of("/\\");
  const auto begin = slash == std::string::npos ? 0 : slash + 1;
  const auto dot = out->find_last_of('.');
  if (dot != std::string::npos && dot > begin) out->erase(dot);
  out->append("_analysis.log");
}
bool ReadLine(FILE* file, std::string* line) {
  line->clear();
  int c = 0;
  while ((c = std::fgetc(file)) != EOF) {
    if (c == '\n') return true;
    line->push_back(static_cast<char>(c));
  }
  return !line->empty();
}
}  // namespace io
}  // namespace fesa
