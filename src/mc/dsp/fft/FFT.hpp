// SPDX-License-Identifier: BSL-1.0

#pragma once

#include <mc/core/config.hpp>

#include <mc/dsp/fft/backend/KissFFT.hpp>
#include <mc/dsp/fft/backend/PFFFT.hpp>
#include <mc/dsp/fft/FFTBackend.hpp>
#include <mc/dsp/fft/FFTDirection.hpp>

#include <mc/core/algorithm.hpp>
#include <mc/core/complex.hpp>
#include <mc/core/memory.hpp>

namespace mc::dsp {

template<typename Engine, typename T>
inline auto fft(Engine& engine, Complex<T> const* input, Complex<T>* output)
    -> decltype(engine.fft(input, output))
{
    return engine.fft(input, output);
}

template<typename Engine, typename T>
inline auto ifft(Engine& engine, Complex<T> const* input, Complex<T>* output)
    -> decltype(engine.ifft(input, output))
{
    return engine.ifft(input, output);
}

template<typename FloatT>
struct FFT
{
    template<typename T>
    FFT(T&& model) : _concept{makeUnique<ModelType<T>>(std::forward<T>(model))}
    {}

    FFT(FFT const& other)                    = delete;
    auto operator=(FFT const& other) -> FFT& = delete;

    FFT(FFT&& other)                    = default;
    auto operator=(FFT&& other) -> FFT& = default;

    auto fft(Complex<FloatT> const* input, Complex<FloatT>* output)
    {
        _concept->do_fft(input, output);
    }

    auto ifft(Complex<FloatT> const* input, Complex<FloatT>* output)
    {
        _concept->do_ifft(input, output);
    }

private:
    struct ConceptType
    {
        virtual ~ConceptType() = default;
        virtual auto do_fft(Complex<FloatT> const* in, Complex<FloatT>* out) -> void  = 0;
        virtual auto do_ifft(Complex<FloatT> const* in, Complex<FloatT>* out) -> void = 0;
    };

    template<typename T>
    struct ModelType final : ConceptType
    {
        ModelType(T&& m) : model{std::forward<T>(m)} {}

        ~ModelType() override = default;

        auto do_fft(Complex<FloatT> const* input, Complex<FloatT>* output) -> void override
        {
            ::mc::dsp::fft(model, input, output);
        }

        auto do_ifft(Complex<FloatT> const* input, Complex<FloatT>* output) -> void override
        {
            ::mc::dsp::ifft(model, input, output);
        }

        T model;
    };

    UniquePtr<ConceptType> _concept{nullptr};
};

}  // namespace mc::dsp
