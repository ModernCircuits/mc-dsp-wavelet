// SPDX-License-Identifier: BSL-1.0

#pragma once

#include <mc/core/complex.hpp>
#include <mc/core/span.hpp>

namespace mc::dsp {

/// This function behaves identically to SpectralConvolution, but computes c=a*conj(b)
/// instead of c=a*b:         a * conj(b) = a+ib * c-id = ac-iad+ibc+bd = ac+bd + i(bc-ad)
auto spectralCorrelation(
    Span<Complex<float> const> a,
    Span<Complex<float> const> b,
    Span<Complex<float>> result
) -> void;

}  // namespace mc::dsp
