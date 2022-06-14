#pragma once

#include "mc/dsp/convolution/FFTConvolver.hpp"
#include "mc/dsp/wavelets/Wavelet.hpp"

#include "mc/preprocessor.hpp"
#include "mc/span.hpp"
#include "mc/string.hpp"

namespace mc::dsp
{
struct WaveletTree
{
    WaveletTree(Wavelet* waveIn, std::size_t signalLength, std::size_t j);

    auto extension(char const* newExtension) noexcept -> void;
    [[nodiscard]] auto extension() const noexcept -> std::string const&;

    auto nodeLength(std::size_t x) -> std::size_t;
    auto coeffs(std::size_t x, std::size_t y, float* coeffs, std::size_t n) const -> void;

    Wavelet* wave;
    std::string method;
    std::size_t siglength;  // Length of the original signal.
    std::size_t outlength;  // Length of the output DWT vector
    std::size_t lenlength;  // Length of the Output Dimension Vector "length"
    std::size_t J;          // Number of decomposition Levels
    std::size_t MaxIter;    // Maximum Iterations J <= MaxIter

    std::size_t N{};  //
    std::size_t nodes;
    std::size_t cfftset;
    std::size_t zpad{};
    std::size_t length[102]{};
    float* output;
    unsigned* coeflength;
    std::unique_ptr<float[]> params;
    unsigned* nodeLength_{nullptr};

private:
    std::string ext_;
};

auto wtree(WaveletTree& wt, float const* inp) -> void;
auto summary(WaveletTree const& wt) -> void;
}  // namespace mc::dsp