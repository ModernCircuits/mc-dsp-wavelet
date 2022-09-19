// SPDX-License-Identifier: BSL-1.0

#pragma once

#include <mc/core/config.hpp>

#include <mc/core/exception.hpp>
#include <mc/core/format.hpp>
#include <mc/core/stdexcept.hpp>

#if defined(MC_COMPILER_GCC) || defined(MC_COMPILER_CLANG)
#define MC_NO_INLINE __attribute__((noinline))
#elif defined(_MSC_VER)
#define MC_NO_INLINE __declspec(noinline)
#else
#define MC_NO_INLINE
#endif

#if defined(MC_COMPILER_GCC) || defined(MC_COMPILER_CLANG)
#define MC_COLD __attribute__((cold))
#else
#define MC_COLD
#endif

namespace mc {

template<typename E, typename... Args>
[[noreturn]] MC_COLD MC_NO_INLINE constexpr auto raise(Args&&... args) -> void
{
    throw E{std::forward<Args>(args)...};
}

template<typename E, typename... Args>
[[noreturn]] MC_COLD MC_NO_INLINE constexpr auto
raisef(::fmt::format_string<Args...> fmtStr, Args&&... args) -> void
{
    throw E{format(fmtStr, std::forward<Args>(args)...)};
}

}  // namespace mc
