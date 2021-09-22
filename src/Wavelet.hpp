#pragma once

#include "tcb/span.hpp"

#include <memory>
#include <string>

struct Wavelet {
    explicit Wavelet(char const* wname);

    [[nodiscard]] auto size() const noexcept -> int { return static_cast<int>(size_); }
    [[nodiscard]] auto name() const noexcept -> std::string const& { return name_; }

    [[nodiscard]] auto lpd() const noexcept -> double const* { return lpd_.data(); }
    [[nodiscard]] auto hpd() const noexcept -> double const* { return hpd_.data(); }
    [[nodiscard]] auto lpr() const noexcept -> double const* { return lpr_.data(); }
    [[nodiscard]] auto hpr() const noexcept -> double const* { return hpr_.data(); }

    [[nodiscard]] auto lpdLen() const noexcept -> int { return static_cast<int>(lpd_.size()); }
    [[nodiscard]] auto hpdLen() const noexcept -> int { return static_cast<int>(hpd_.size()); }
    [[nodiscard]] auto lprLen() const noexcept -> int { return static_cast<int>(lpr_.size()); }
    [[nodiscard]] auto hprLen() const noexcept -> int { return static_cast<int>(hpr_.size()); }

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