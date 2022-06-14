#pragma once

#if defined(_MSC_VER)
    #define LT_NODISCARD _Check_return_
#else
    #define LT_NODISCARD __attribute__((warn_unused_result))  // NOLINT
#endif