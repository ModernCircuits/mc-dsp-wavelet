#pragma once

#include "lt/algorithm.hpp"
#include "lt/complex.hpp"
#include "lt/memory.hpp"

constexpr auto Pi2 = 6.28318530717958647692528676655900577;

struct FFT {
    enum Direction {
        backward = -1,
        forward = 1,
    };

    FFT(int n, Direction direction);

    auto perform(Complex<double> const* inp, Complex<double>* oup) -> void;

    [[nodiscard]] auto direction() const noexcept -> Direction { return direction_; }
    auto direction(Direction newDirection) noexcept -> void { direction_ = newDirection; }

    [[nodiscard]] auto size() const noexcept -> int { return size_; }
    auto size(int newSize) noexcept -> void { size_ = newSize; }

    int factors[64] {};
    int lf;
    int lt;
    std::unique_ptr<Complex<double>[]> data;

private:
    int size_;
    Direction direction_;
};

struct RealFFT {
    RealFFT(int n, FFT::Direction direction);

    auto performRealToComplex(double const* inp, Complex<double>* oup) -> void;
    auto performComplexToReal(Complex<double> const* inp, double* oup) -> void;

private:
    std::unique_ptr<FFT> cobj_;
    std::unique_ptr<Complex<double>[]> data_;
};

auto divideby(int m, int d) -> int;
