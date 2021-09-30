#pragma once

#include "lt/algorithm.hpp"
#include "lt/complex.hpp"
#include "lt/memory.hpp"
#include "lt/preprocessor.hpp"

#include "kissfft/kissfft.hh"

namespace lt {
namespace dsp {
    constexpr auto pi2 = 6.28318530717958647692528676655900577;

    struct FFT {
        enum Direction {
            forward = 0,
            backward = 1,
        };

        FFT(int n, Direction direction);

        auto perform(Complex<float> const* inp, Complex<float>* oup) -> void;

        LT_NODISCARD auto direction() const noexcept -> Direction;
        LT_NODISCARD auto size() const noexcept -> int;
        LT_NODISCARD auto engine() -> kissfft<float>& { return fftEngine_; }

    private:
        int size_;
        Direction direction_;

        kissfft<float> fftEngine_;
    };

    struct RealFFT {
        RealFFT(int n, FFT::Direction direction);

        auto performRealToComplex(float const* inp, Complex<float>* oup) -> void;
        auto performComplexToReal(Complex<float> const* inp, float* oup) -> void;

    private:
        std::unique_ptr<FFT> cobj_;
        std::unique_ptr<Complex<float>[]> data_;
    };

    auto divideby(int m, int d) -> int;
}
}