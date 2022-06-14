#pragma once

#include "mc/dsp/fft/WindowFunction.hpp"

#include "mc/algorithm.hpp"
#include "mc/cassert.hpp"
#include "mc/cmath.hpp"
#include "mc/cstddef.hpp"
#include "mc/functional.hpp"
#include "mc/numbers.hpp"
#include "mc/vector.hpp"

namespace mc::dsp
{

template<typename T>
auto ncos(std::size_t order, std::size_t i, std::size_t size) noexcept -> T
{
    return std::cos(static_cast<T>(order * i) * static_cast<T>(numbers::pi) / static_cast<T>(size - 1));
}

template<typename It>
auto fillWindow(It first, It last, WindowFunction type, bool normalise) noexcept
{
    using T = std::decay_t<decltype(*first)>;

    auto const size = static_cast<std::size_t>(std::distance(first, last));
    switch (type)
    {
        case WindowFunction::rectangular:
        {
            for (std::size_t i = 0; i < size; ++i) { first[i] = static_cast<T>(1); }
            break;
        }

        case WindowFunction::triangular:
        {
            auto const halfSlots = static_cast<T>(0.5) * static_cast<T>(size - 1);

            for (std::size_t i = 0; i < size; ++i)
            {
                first[i] = static_cast<T>(1.0) - std::abs((static_cast<T>(i) - halfSlots) / halfSlots);
            }
            break;
        }

        case WindowFunction::hann:
        {
            for (std::size_t i = 0; i < size; ++i)
            {
                auto const cos2 = ncos<T>(2, i, size);
                first[i]        = static_cast<T>(0.5 - 0.5 * cos2);
            }
            break;
        }

        case WindowFunction::hamming:
        {
            for (std::size_t i = 0; i < size; ++i)
            {
                auto const cos2 = ncos<T>(2, i, size);
                first[i]        = static_cast<T>(0.54 - 0.46 * cos2);
            }
            break;
        }

        case WindowFunction::blackman:
        {
            constexpr T alpha = 0.16F;

            for (std::size_t i = 0; i < size; ++i)
            {
                auto const cos2 = ncos<T>(2, i, size);
                auto const cos4 = ncos<T>(4, i, size);
                first[i]        = static_cast<T>(0.5 * (1 - alpha) - 0.5 * cos2 + 0.5 * alpha * cos4);
            }
            break;
        }

        case WindowFunction::blackmanHarris:
        {
            for (std::size_t i = 0; i < size; ++i)
            {
                auto const cos2 = ncos<T>(2, i, size);
                auto const cos4 = ncos<T>(4, i, size);
                auto const cos6 = ncos<T>(6, i, size);

                first[i] = static_cast<T>(0.35875 - 0.48829 * cos2 + 0.14128 * cos4 - 0.01168 * cos6);
            }
            break;
        }

        default:
        {
            MC_ASSERT(false);
            break;
        }
    }

    // DC frequency amplitude must be one
    if (normalise)
    {
        auto sum = T(0);
        for (std::size_t i = 0; i < size; ++i) { sum += first[i]; }
        auto const factor = static_cast<T>(size) / sum;
        std::transform(first, last, first, [factor](auto s) { return s * factor; });
    }
}

template<typename WindowIt, typename SignalIt>
auto multiplyWithWindow(SignalIt signalF, SignalIt signalL, WindowIt winF, WindowIt winL) noexcept
{
    auto const windowSize = static_cast<std::size_t>(std::distance(winF, winL));
    auto const signalSize = static_cast<std::size_t>(std::distance(signalF, signalL));
    auto const minSize    = std::min(signalSize, windowSize);
    std::transform(signalF, std::next(signalF, minSize), winF, signalF, std::multiplies<>{});
}

}  // namespace mc::dsp