// SPDX-License-Identifier: BSL-1.0

#pragma once

#include <mc/dsp/fft/FFTBackend.hpp>
#include <mc/dsp/fft/FFTDirection.hpp>

#include <mc/core/complex.hpp>

#include <pffft.h>

namespace mc::dsp {

// struct PFFFT
// {};

// template<>
// struct FFTBackend<float, PFFFT>
// {
//     using value_type  = float;
//     using handle_type = PFFFT_Setup*;

// static auto construct(std::size_t size, FFTDirection /*direction*/)
// {
//     return pffft_new_setup(static_cast<int>(size), PFFFT_COMPLEX);
// }

// static auto destroy(handle_type handle) { pffft_destroy_setup(handle); }

// static auto perform(
//     handle_type handle,
//     Complex<float> const* in,
//     Complex<float>* out,
//     FFTDirection direction
// ) -> void
// {
//     auto const dir
//         = direction == FFTDirection::forward ? PFFFT_FORWARD : PFFFT_BACKWARD;
//     auto const* input = reinterpret_cast<float const*>(in);  // NOLINT
//     auto* output      = reinterpret_cast<float*>(out);       // NOLINT
//     pffft_transform_ordered(handle, input, output, nullptr, dir);
// }
// };

struct PFFFT
{
    PFFFT(int size) : _setup{pffft_new_setup(size, PFFFT_COMPLEX)} {}

    ~PFFFT()
    {
        if (_setup != nullptr) { pffft_destroy_setup(_setup); }
    }

    PFFFT(PFFFT const& other)                    = delete;
    auto operator=(PFFFT const& other) -> PFFFT& = delete;

    PFFFT(PFFFT&& other) noexcept : _setup{std::exchange(other._setup, nullptr)} {}

    auto operator=(PFFFT&& other) noexcept -> PFFFT&
    {
        _setup = std::exchange(other._setup, nullptr);
        return *this;
    }

    auto fft(std::complex<float> const* in, std::complex<float>* out)
    {
        auto const* input = reinterpret_cast<float const*>(in);  // NOLINT
        auto* output      = reinterpret_cast<float*>(out);       // NOLINT
        pffft_transform_ordered(_setup, input, output, nullptr, PFFFT_FORWARD);
    }

    auto ifft(std::complex<float> const* in, std::complex<float>* out)
    {
        auto const* input = reinterpret_cast<float const*>(in);  // NOLINT
        auto* output      = reinterpret_cast<float*>(out);       // NOLINT
        pffft_transform_ordered(_setup, input, output, nullptr, PFFFT_BACKWARD);
    }

private:
    PFFFT_Setup* _setup;
};

}  // namespace mc::dsp
