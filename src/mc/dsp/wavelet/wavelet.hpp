// SPDX-License-Identifier: BSL-1.0

#pragma once

#include <mc/core/config.hpp>

#include <mc/core/format.hpp>
#include <mc/core/span.hpp>
#include <mc/core/string.hpp>
#include <mc/core/string_view.hpp>
#include <mc/core/vector.hpp>

namespace mc::dsp {
struct Wavelet
{
    explicit Wavelet(StringView name);

    [[nodiscard]] auto size() const noexcept -> std::size_t;
    [[nodiscard]] auto name() const noexcept -> String const&;

    [[nodiscard]] auto lpd() const noexcept -> Span<float const>;
    [[nodiscard]] auto hpd() const noexcept -> Span<float const>;
    [[nodiscard]] auto lpr() const noexcept -> Span<float const>;
    [[nodiscard]] auto hpr() const noexcept -> Span<float const>;

private:
    String _name;
    std::size_t _size;
    Vector<float> _params;
};

}  // namespace mc::dsp

template<>
struct fmt::formatter<mc::dsp::Wavelet> : formatter<string_view>
{
    template<typename FormatContext>
    auto format(mc::dsp::Wavelet const& wavelet, FormatContext& ctx) const
    {
        fmt::format_to(ctx.out(), "Wavelet: {0}\n", wavelet.name());
        fmt::format_to(ctx.out(), "  Filters length: {0}\n", wavelet.size());

        fmt::format_to(ctx.out(), "lpd: [{0}]\n", fmt::join(wavelet.lpd(), ", "));
        fmt::format_to(ctx.out(), "hpd: [{0}]\n", fmt::join(wavelet.hpd(), ", "));
        fmt::format_to(ctx.out(), "lpr: [{0}]\n", fmt::join(wavelet.lpr(), ", "));
        return fmt::format_to(ctx.out(), "hpr: [{0}]\n", fmt::join(wavelet.hpr(), ", "));
    }
};
