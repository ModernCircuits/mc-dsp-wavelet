#pragma once

#include <mc/core/config.hpp>

#include <mc/dsp/convolution/FFTConvolver.hpp>
#include <mc/dsp/wavelets/Wavelet.hpp>

#include <mc/core/span.hpp>
#include <mc/core/string.hpp>

namespace mc::dsp {

struct WaveletPacketTransform
{
    WaveletPacketTransform(Wavelet* wave, std::size_t siglength, std::size_t j);

    [[nodiscard]] auto wave() const noexcept -> Wavelet const& { return *wave_; }

    [[nodiscard]] auto signalLength() const noexcept -> int { return signalLength_; }

    FFTConvolver* cobj{};
    int outlength{};  // Length of the output DWT vector
    int lenlength{};  // Length of the Output Dimension Vector "length"
    int J{};          // Number of decomposition Levels
    int MaxIter{};    // Maximum Iterations J <= MaxIter
    String ext;       // Type of Extension used - "per" or "sym"
    String entropy;
    float eparam{};

    int N{};  //
    int nodes{};
    int length[102]{};
    float* output{};
    float* costvalues{};
    float* basisvector{};
    int* nodeindex{};
    int* numnodeslevel{};
    int* coeflength{};
    UniquePtr<float[]> params;

private:
    Wavelet* wave_{nullptr};
    int signalLength_{};  // Length of the original signal.
};

auto dwpt(WaveletPacketTransform& wt, float const* inp) -> void;
auto idwpt(WaveletPacketTransform& wt, float* dwtop) -> void;

auto setDWPTExtension(WaveletPacketTransform& wt, char const* extension) -> void;
auto setDWPTEntropy(WaveletPacketTransform& wt, char const* entropy, float eparam) -> void;
auto getDWPTNodelength(WaveletPacketTransform& wt, int x) -> int;

auto summary(WaveletPacketTransform const& wt) -> void;

}  // namespace mc::dsp
