#pragma once

#include "lt/memory.hpp"
#include "lt/preprocessor.hpp"
#include "lt/span.hpp"
#include "lt/string.hpp"

namespace lt::dsp
{
struct Wavelet
{
    explicit Wavelet(char const* wname);

    LT_NODISCARD auto size() const noexcept -> std::size_t { return size_; }
    LT_NODISCARD auto name() const noexcept -> std::string const& { return name_; }

    LT_NODISCARD auto lpd() const noexcept -> span<float> { return lpd_; }
    LT_NODISCARD auto hpd() const noexcept -> span<float> { return hpd_; }
    LT_NODISCARD auto lpr() const noexcept -> span<float> { return lpr_; }
    LT_NODISCARD auto hpr() const noexcept -> span<float> { return hpr_; }

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
}  // namespace lt::dsp