#pragma once

#include "lt/dsp/fft/FFTBackend.hpp"

#include "lt/complex.hpp"

#include "kissfft/kissfft.hh"

namespace lt {
namespace dsp {

    struct KissFFT {
    };

    template <>
    struct FFTBackend<float, KissFFT> {
        using value_type = float;
        using handle_type = kissfft<float>;

        static auto construct(std::size_t size, FFTDirection direction)
        {
            return kissfft<float> { size, direction != FFTDirection::forward };
        }

        static auto destroy(handle_type& /*handle*/)
        {
        }

        static auto perform(handle_type& handle, Complex<float> const* in, Complex<float>* out, FFTDirection direction) -> void
        {
            (void)direction;
            handle.transform(in, out);
        }
    };

} // namespace dsp
}