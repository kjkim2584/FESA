/// @file
/// @brief Internal UTF-8/native file boundary shared by cfg and log I/O.
#ifndef FESA_IO_FILE_ACCESS_H_
#define FESA_IO_FILE_ACCESS_H_
#include <cstdio>
#include <string>

#include "fesa/base/diagnostic.h"
namespace fesa {
namespace io {
/// Native file opening policy. CreateNew never truncates an existing file.
enum class FileMode { kRead, kReplace, kCreateNew };
/// Borrows path; returns caller-owned FILE or null with native error in d.
/// May throw bad_alloc; public call boundaries translate it to status.
StatusCode OpenFile(const char* path, FileMode mode, FILE** out, Diagnostic* d);
/// Inspects native attributes without opening or changing a file. Missing
/// files/parents yield exists=false; malformed Unicode/access errors fail.
/// Outputs change only on success; may throw bad_alloc to public boundary.
StatusCode InspectPath(const char* path, bool* exists, bool* regular,
                       Diagnostic* d);
/// Resolves UTF-8 input once using native absolute path rules. No filesystem
/// path default/move constructors; bad_alloc propagates to owning boundary.
StatusCode AbsolutePath(const char* input, std::string* out, Diagnostic* d);
/// Rejects case-insensitive absolute lexical equality or existing-file native
/// identity equality (including hard links/followed aliases). Missing paths
/// use lexical checking only. Inspection errors fail closed. Neither path is
/// modified or retained, and d's caller-owned source context is preserved.
/// Snapshot check only: no protection from concurrent replacement, nor full
/// canonicalization of missing paths through aliased parent directories.
/// May throw bad_alloc; the public owning boundary translates it to status.
StatusCode ValidateDistinctPaths(const char* source, const char* destination,
                                 Diagnostic* d);
/// Derives input stem .h5 without filesystem access; may throw bad_alloc.
void ResultPath(const char* input, std::string* out);
/// Writes directory/solver.cfg without normalization or filesystem access.
/// Uses throwing string operations so Debug allocator failure is catchable.
void DefaultConfigPath(const char* directory, std::string* out);
/// Writes input directory/base_analysis.log; preserves UTF-8 bytes.
void AnalysisLogPath(const char* input, std::string* out);
/// Reads one physical line excluding LF. False means EOF, checked via ferror.
/// Keeps NUL bytes so the caller can reject them before C-string parsing.
/// May throw bad_alloc; caller owns both borrowed parameters.
bool ReadLine(FILE* file, std::string* line);
}  // namespace io
}  // namespace fesa
#endif  // FESA_IO_FILE_ACCESS_H_
