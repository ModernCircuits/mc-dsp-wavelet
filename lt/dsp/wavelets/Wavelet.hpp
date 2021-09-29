#pragma once

#include "lt/memory.hpp"
#include "lt/preprocessor.hpp"
#include "lt/span.hpp"
#include "lt/string.hpp"

struct Wavelet {
    explicit Wavelet(char const* wname);

    LT_NODISCARD auto size() const noexcept -> std::size_t { return size_; }
    LT_NODISCARD auto name() const noexcept -> std::string const& { return name_; }

    LT_NODISCARD auto lpd() const noexcept -> lt::span<float> { return lpd_; }
    LT_NODISCARD auto hpd() const noexcept -> lt::span<float> { return hpd_; }
    LT_NODISCARD auto lpr() const noexcept -> lt::span<float> { return lpr_; }
    LT_NODISCARD auto hpr() const noexcept -> lt::span<float> { return hpr_; }

private:
    std::string name_;

    // When all filters are of the same length.
    // [Matlab uses zero-padding to make all filters of the same length]
    std::size_t size_;
    std::unique_ptr<float[]> params_;

    lt::span<float> lpd_;
    lt::span<float> hpd_;
    lt::span<float> lpr_;
    lt::span<float> hpr_;
};

auto summary(Wavelet const& obj) -> void;
