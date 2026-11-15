#include "fesa/base/build_info.h"

#include <climits>

namespace fesa {

long CompiledLanguageVersion() { return __cplusplus; }

int CompiledPointerBits() { return sizeof(void*) * CHAR_BIT; }

}  // namespace fesa
