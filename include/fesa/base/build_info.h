/// @file
/// @brief Properties of the compiled solver core; no SDK types are exposed.
#ifndef FESA_BASE_BUILD_INFO_H_
#define FESA_BASE_BUILD_INFO_H_

namespace fesa {

/// Returns the language version used to compile the core translation unit.
long CompiledLanguageVersion();

/// Returns the core pointer width in bits.
int CompiledPointerBits();

}  // namespace fesa

#endif  // FESA_BASE_BUILD_INFO_H_
