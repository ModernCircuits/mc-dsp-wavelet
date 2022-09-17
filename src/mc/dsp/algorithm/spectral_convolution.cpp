#include "spectral_convolution.hpp"

#include <mc/core/cassert.hpp>

namespace mc::dsp {
auto spectralConvolution(
    Span<Complex<float> const> a,
    Span<Complex<float> const> b,
    Span<Complex<float>> result
) -> void
{
    MC_ASSERT((a.size() == result.size()) && (b.size() == result.size()));
    for (std::size_t i = 0; i < a.size(); ++i) { result[i] = a[i] * b[i]; }
}

}  // namespace mc::dsp
