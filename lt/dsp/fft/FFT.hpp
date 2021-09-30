#pragma once

#include "lt/dsp/fft/FFTBackend.hpp"
#include "lt/dsp/fft/FFTDirection.hpp"

#include "lt/dsp/fft/backend/KissFFT.hpp"

#include "lt/algorithm.hpp"
#include "lt/complex.hpp"
#include "lt/memory.hpp"
#include "lt/preprocessor.hpp"

namespace lt {
namespace dsp {

    template <typename T, typename BackendTag>
    struct FFT {
        using value_type = T;
        using backend_type = typename FFTBackend<T, BackendTag>::handle_type;

        FFT(int n, FFTDirection direction);

        auto perform(Complex<T> const* inp, Complex<T>* oup) -> void;

        LT_NODISCARD auto direction() const noexcept -> FFTDirection;
        LT_NODISCARD auto size() const noexcept -> int;
        LT_NODISCARD auto engine() -> backend_type& { return engine_; }

    private:
        int size_;
        FFTDirection direction_;
        backend_type engine_;
    };

    template <typename T, typename BackendTag>
    FFT<T, BackendTag>::FFT(int n, FFTDirection direction)
        : size_ { n }
        , direction_ { direction }
        , engine_ { FFTBackend<T, BackendTag>::construct(static_cast<std::size_t>(n), direction) }
    {
    }

    template <typename T, typename BackendTag>
    auto FFT<T, BackendTag>::direction() const noexcept -> FFTDirection { return direction_; }

    template <typename T, typename BackendTag>
    auto FFT<T, BackendTag>::size() const noexcept -> int { return size_; }

    template <typename T, typename BackendTag>
    auto FFT<T, BackendTag>::perform(Complex<T> const* input, Complex<T>* output) -> void
    {
        FFTBackend<T, BackendTag>::perform(engine_, input, output);
    }

    struct RFFT {
        RFFT(int n, FFTDirection direction);

        auto performRealToComplex(float const* inp, Complex<float>* oup) -> void;
        auto performComplexToReal(Complex<float> const* inp, float* oup) -> void;

    private:
        std::unique_ptr<FFT<float, KissFFT>> cobj_;
        std::unique_ptr<Complex<float>[]> data_;
    };

}
}