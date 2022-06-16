#pragma once

#include "mc/memory.hpp"
#include "mc/preprocessor.hpp"
#include "mc/span.hpp"
#include "mc/string.hpp"

namespace mc::dsp
{
struct Wavelet
{
    explicit Wavelet(char const* wname);

    [[nodiscard]] auto size() const noexcept -> std::size_t;
    [[nodiscard]] auto name() const noexcept -> std::string const&;

    [[nodiscard]] auto lpd() const noexcept -> span<float>;
    [[nodiscard]] auto hpd() const noexcept -> span<float>;
    [[nodiscard]] auto lpr() const noexcept -> span<float>;
    [[nodiscard]] auto hpr() const noexcept -> span<float>;

private:
    std::string name_;

    // When all filters are of the same length.
    // [Matlab uses zero-padding to make all filters of the same length]
    std::size_t size_;
    std::unique_ptr<float[]> params_;
};

auto summary(Wavelet const& obj) -> void;
}  // namespace mc::dsp