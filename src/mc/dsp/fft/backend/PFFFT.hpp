#pragma once

#include <mc/dsp/fft/FFTBackend.hpp>
#include <mc/dsp/fft/FFTDirection.hpp>

#include <mc/core/complex.hpp>

#include <pffft.h>

namespace mc::dsp {

struct PFFFT
{};

template<>
struct FFTBackend<float, PFFFT>
{
    using value_type  = float;
    using handle_type = PFFFT_Setup*;

    static auto construct(std::size_t size, FFTDirection /*direction*/)
    {
        return pffft_new_setup(static_cast<int>(size), PFFFT_COMPLEX);
    }

    static auto destroy(handle_type handle) { pffft_destroy_setup(handle); }

    static auto perform(
        handle_type handle,
        Complex<float> const* in,
        Complex<float>* out,
        FFTDirection direction
    ) -> void
    {
        auto const dir
            = direction == FFTDirection::forward ? PFFFT_FORWARD : PFFFT_BACKWARD;
        auto const* input = reinterpret_cast<float const*>(in);  // NOLINT
        auto* output      = reinterpret_cast<float*>(out);       // NOLINT
        pffft_transform_ordered(handle, input, output, nullptr, dir);
    }
};

}  // namespace mc::dsp
