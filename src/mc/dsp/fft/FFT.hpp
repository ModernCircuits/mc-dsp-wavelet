// SPDX-License-Identifier: BSL-1.0

#pragma once

#include <mc/core/config.hpp>

#include <mc/dsp/fft/backend/PFFFT.hpp>

#include <mc/core/algorithm.hpp>
#include <mc/core/complex.hpp>
#include <mc/core/memory.hpp>

namespace mc::dsp {

template<typename Engine>
inline auto
fft(Engine& engine, Span<Complex<float> const> input, Span<Complex<float>> output)
    -> decltype(engine.fft(input, output))
{
    return engine.fft(input, output);
}

template<typename Engine>
inline auto
ifft(Engine& engine, Span<Complex<float> const> input, Span<Complex<float>> output)
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

    auto fft(Span<Complex<FloatT> const> input, Span<Complex<FloatT>> output)
    {
        _concept->do_fft(input, output);
    }

    auto ifft(Span<Complex<FloatT> const> input, Span<Complex<FloatT>> output)
    {
        _concept->do_ifft(input, output);
    }

private:
    struct ConceptType
    {
        virtual ~ConceptType() = default;
        virtual auto do_fft(Span<Complex<FloatT> const> in, Span<Complex<FloatT>> out)
            -> void
            = 0;
        virtual auto do_ifft(Span<Complex<FloatT> const> in, Span<Complex<FloatT>> out)
            -> void
            = 0;
    };

    template<typename T>
    struct ModelType final : ConceptType
    {
        ModelType(T&& m) : model{std::forward<T>(m)} {}

        ~ModelType() override = default;

        auto do_fft(Span<Complex<FloatT> const> input, Span<Complex<FloatT>> output)
            -> void override
        {
            ::mc::dsp::fft(model, input, output);
        }

        auto do_ifft(Span<Complex<FloatT> const> input, Span<Complex<FloatT>> output)
            -> void override
        {
            ::mc::dsp::ifft(model, input, output);
        }

        T model;
    };

    UniquePtr<ConceptType> _concept{nullptr};
};

[[nodiscard]] auto makeFFT(std::size_t size) -> FFT<float>;

}  // namespace mc::dsp
