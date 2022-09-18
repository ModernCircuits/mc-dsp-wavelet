#pragma once

namespace mc::dsp {
inline constexpr auto ipow2(int n) -> int
{
    auto p = 1;
    for (auto i = 0; i < n; ++i) { p *= 2; }
    return p;
}
}  // namespace mc::dsp
