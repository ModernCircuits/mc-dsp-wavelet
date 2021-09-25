#pragma once

#include "lt/memory.hpp"
#include "lt/span.hpp"
#include "lt/string.hpp"

struct Wavelet {
    explicit Wavelet(char const* wname);

    [[nodiscard]] auto size() const noexcept -> std::size_t { return size_; }
    [[nodiscard]] auto name() const noexcept -> std::string const& { return name_; }

    [[nodiscard]] auto lpd() const noexcept -> lt::span<double> { return lpd_; }
    [[nodiscard]] auto hpd() const noexcept -> lt::span<double> { return hpd_; }
    [[nodiscard]] auto lpr() const noexcept -> lt::span<double> { return lpr_; }
    [[nodiscard]] auto hpr() const noexcept -> lt::span<double> { return hpr_; }

private:
    std::string name_;

    // When all filters are of the same length.
    // [Matlab uses zero-padding to make all filters of the same length]
    std::size_t size_;
    std::unique_ptr<double[]> params_;

    lt::span<double> lpd_;
    lt::span<double> hpd_;
    lt::span<double> lpr_;
    lt::span<double> hpr_;
};

auto summary(Wavelet const& obj) -> void;
