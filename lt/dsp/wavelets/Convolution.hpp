#pragma once

#include "lt/dsp/fft/FFT.hpp"

struct Convolution {
    Convolution(int n, int l);

    auto fft(double const* inp1, double const* inp2, double* oup) const -> void;

    static auto direct(double const* inp1, int n, double const* inp2, int l, double* oup) noexcept -> void;

private:
    int ilen1_;
    int ilen2_;
    int clen_;
    std::unique_ptr<RealFFT> fobj_;
    std::unique_ptr<RealFFT> iobj_;

    std::unique_ptr<double[]> a { std::make_unique<double[]>(clen_) };
    std::unique_ptr<double[]> b { std::make_unique<double[]>(clen_) };
    std::unique_ptr<Complex<double>[]> c { std::make_unique<Complex<double>[]>(clen_) };
    std::unique_ptr<Complex<double>[]> ao { std::make_unique<Complex<double>[]>(clen_) };
    std::unique_ptr<Complex<double>[]> bo { std::make_unique<Complex<double>[]>(clen_) };
    std::unique_ptr<double[]> co { std::make_unique<double[]>(clen_) };
};
