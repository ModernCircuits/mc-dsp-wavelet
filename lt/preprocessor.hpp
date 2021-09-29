#pragma once

// #if __has_cpp_attribute(nodiscard)
// #define LT_NODISCARD [[nodiscard]] // NOLINT
#if defined(_MSC_VER)
#define LT_NODISCARD _Check_return_
#else
#define LT_NODISCARD __attribute__((warn_unused_result)) // NOLINT
#endif