#pragma once

#include <mc/core/config.hpp>

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

[[nodiscard]] auto summary(Wavelet const& wavelet) -> String;

}  // namespace mc::dsp
