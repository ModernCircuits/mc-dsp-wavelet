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

    [[nodiscard]] auto size() const noexcept -> std::size_t { return size_; }
    [[nodiscard]] auto name() const noexcept -> std::string const& { return name_; }

    [[nodiscard]] auto lpd() const noexcept -> span<float> { return lpd_; }
    [[nodiscard]] auto hpd() const noexcept -> span<float> { return hpd_; }
    [[nodiscard]] auto lpr() const noexcept -> span<float> { return lpr_; }
    [[nodiscard]] auto hpr() const noexcept -> span<float> { return hpr_; }

private:
    std::string name_;

    // When all filters are of the same length.
    // [Matlab uses zero-padding to make all filters of the same length]
    std::size_t size_;
    std::unique_ptr<float[]> params_;

    span<float> lpd_;
    span<float> hpd_;
    span<float> lpr_;
    span<float> hpr_;
};

auto summary(Wavelet const& obj) -> void;
}  // namespace mc::dsp