#pragma once

#include <mc/core/complex.hpp>
#include <mc/core/span.hpp>

namespace mc::dsp {

/// This free function takes three complex signals a,b,c of the same size and computes the
/// complex element-wise multiplication:   a+ib * c+id = ac+iad+ibc-bd = ac-bd + i(ad+bc)
/// The computation loop isn't sent to OMP because this function itself is already expected
/// to be called by multiple threads, and it would actually slow down the process. It
/// throuws an exception if
auto spectralConvolution(
    Span<Complex<float> const> a,
    Span<Complex<float> const> b,
    Span<Complex<float>> result
) -> void;

}  // namespace mc::dsp
