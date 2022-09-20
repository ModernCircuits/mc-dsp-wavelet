// SPDX-License-Identifier: BSL-1.0

#pragma once

#include <mc/core/config.hpp>

#include <mc/dsp/fft/fft.hpp>

#include <mc/core/complex.hpp>
#include <mc/core/memory.hpp>
#include <mc/core/span.hpp>

#include <pffft.h>

namespace mc::dsp {

template<typename Engine>
auto rfft(Engine& engine, Span<float const> input, Span<Complex<float>> output)
    -> decltype(engine.rfft(input, output))
{
    return engine.rfft(input, output);
}

template<typename Engine>
auto irfft(Engine& engine, Span<Complex<float> const> input, Span<float> output)
    -> decltype(engine.irfft(input, output))
{
    return engine.irfft(input, output);
}

template<typename Engine>
auto rfft(Engine& engine, Span<double const> input, Span<Complex<double>> output)
    -> decltype(engine.rfft(input, output))
{
    return engine.rfft(input, output);
}

template<typename Engine>
auto irfft(Engine& engine, Span<Complex<double> const> input, Span<double> output)
    -> decltype(engine.irfft(input, output))
{
    return engine.irfft(input, output);
}

template<typename FloatT>
struct RFFT
{
    template<typename T>
    RFFT(T&& model) : _concept{makeUnique<ModelType<T>>(std::forward<T>(model))}
    {}

    RFFT(RFFT const& other)                    = delete;
    auto operator=(RFFT const& other) -> RFFT& = delete;

    RFFT(RFFT&& other)                    = default;
    auto operator=(RFFT&& other) -> RFFT& = default;

    auto rfft(Span<FloatT const> input, Span<Complex<FloatT>> output)
    {
        _concept->do_rfft(input, output);
    }

    auto irfft(Span<Complex<FloatT> const> input, Span<FloatT> output)
    {
        _concept->do_irfft(input, output);
    }

private:
    struct ConceptType
    {
        virtual ~ConceptType() = default;
        virtual auto do_rfft(Span<FloatT const> in, Span<Complex<FloatT>> out) -> void  = 0;
        virtual auto do_irfft(Span<Complex<FloatT> const> in, Span<FloatT> out) -> void = 0;
    };

    template<typename T>
    struct ModelType final : ConceptType
    {
        ModelType(T&& m) : model{std::forward<T>(m)} {}

        ~ModelType() override = default;

        auto do_rfft(Span<FloatT const> input, Span<Complex<FloatT>> output)
            -> void override
        {
            ::mc::dsp::rfft(model, input, output);
        }

        auto do_irfft(Span<Complex<FloatT> const> input, Span<FloatT> output)
            -> void override
        {
            ::mc::dsp::irfft(model, input, output);
        }

        T model;
    };

    UniquePtr<ConceptType> _concept{nullptr};
};

[[nodiscard]] auto makeRFFT(std::size_t size) -> RFFT<float>;

}  // namespace mc::dsp
