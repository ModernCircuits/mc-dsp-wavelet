// SPDX-License-Identifier: BSL-1.0

#pragma once

#include <mc/core/config.hpp>

#include <mc/dsp/fft/backend/pffft.hpp>

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

template<typename T>
struct FFT
{
    template<typename ImplT>
    FFT(ImplT&& model) : _concept{makeUnique<ModelType<ImplT>>(std::forward<ImplT>(model))}
    {}

    FFT(FFT const& other)                    = delete;
    auto operator=(FFT const& other) -> FFT& = delete;

    FFT(FFT&& other)                    = default;
    auto operator=(FFT&& other) -> FFT& = default;

    auto fft(Span<Complex<T> const> in, Span<Complex<T>> out) -> void
    {
        _concept->do_fft(in, out);
    }

    auto ifft(Span<Complex<T> const> in, Span<Complex<T>> out) -> void
    {
        _concept->do_ifft(in, out);
    }

private:
    struct ConceptType
    {
        virtual ~ConceptType() = default;
        virtual auto do_fft(Span<Complex<T> const> in, Span<Complex<T>> out) -> void  = 0;
        virtual auto do_ifft(Span<Complex<T> const> in, Span<Complex<T>> out) -> void = 0;
    };

    template<typename ImplT>
    struct ModelType final : ConceptType
    {
        ModelType(ImplT&& m) : model{std::forward<ImplT>(m)} {}

        ~ModelType() override = default;

        auto do_fft(Span<Complex<T> const> in, Span<Complex<T>> out) -> void override
        {
            ::mc::dsp::fft(model, in, out);
        }

        auto do_ifft(Span<Complex<T> const> in, Span<Complex<T>> out) -> void override
        {
            ::mc::dsp::ifft(model, in, out);
        }

        ImplT model;
    };

    UniquePtr<ConceptType> _concept{nullptr};
};

[[nodiscard]] auto makeFFT(std::size_t size) -> FFT<float>;

}  // namespace mc::dsp
