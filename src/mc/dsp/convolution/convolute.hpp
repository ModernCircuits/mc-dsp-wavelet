#pragma once

#include <mc/core/algorithm.hpp>
#include <mc/core/cstddef.hpp>

namespace mc::dsp {
template<typename Convolver>
auto convolute(
    Convolver& c,
    typename Convolver::value_type const* s,
    typename Convolver::value_type const* p,
    typename Convolver::value_type* out
) -> decltype(c.convolute(s, p, out))
{
    return c.convolute(s, p, out);
}

template<typename T>
auto convolute(
    T const* signal,
    std::size_t n,
    T const* patch,
    std::size_t l,
    T* output
) noexcept -> void
{
    auto const mm = n + l - 1;

    if (n >= l) {
        auto i = std::size_t{0};
        for (auto k = std::size_t{0}; k < l; k++) {
            output[k] = 0.0;
            for (auto m = std::size_t{0}; m <= k; m++) {
                output[k] += signal[m] * patch[k - m];
            }
        }
        for (auto k = l; k < mm; k++) {
            output[k] = 0.0;
            i++;
            auto const t1   = l + i;
            auto const tmin = static_cast<T>(std::min(t1, n));
            for (auto m = i; m < tmin; m++) { output[k] += signal[m] * patch[k - m]; }
        }
        return;
    }

    auto i = std::size_t{0};
    for (auto k = std::size_t{0}; k < n; k++) {
        output[k] = 0.0;
        for (auto m = std::size_t{0}; m <= k; m++) {
            output[k] += patch[m] * signal[k - m];
        }
    }
    for (auto k = n; k < mm; k++) {
        output[k] = 0.0;
        i++;
        auto const t1   = n + i;
        auto const tmin = static_cast<T>(std::min(t1, l));
        for (auto m = i; m < tmin; m++) { output[k] += patch[m] * signal[k - m]; }
    }
}
}  // namespace mc::dsp
