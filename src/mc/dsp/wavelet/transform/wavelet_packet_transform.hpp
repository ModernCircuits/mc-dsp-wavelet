// SPDX-License-Identifier: BSL-1.0

#pragma once

#include <mc/core/config.hpp>

#include <mc/dsp/algorithm/ipow2.hpp>
#include <mc/dsp/convolution.hpp>
#include <mc/dsp/wavelet/wavelet.hpp>

#include <mc/core/format.hpp>
#include <mc/core/span.hpp>
#include <mc/core/string.hpp>

namespace mc {

struct WaveletPacketTransform
{
    WaveletPacketTransform(Wavelet* wave, size_t siglength, size_t j);

    [[nodiscard]] auto wave() const noexcept -> Wavelet const& { return *_wave; }

    [[nodiscard]] auto signalLength() const noexcept -> int { return _signalLength; }

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
    Wavelet* _wave{nullptr};
    int _signalLength{};  // Length of the original signal.
};

auto dwpt(WaveletPacketTransform& wt, float const* inp) -> void;
auto idwpt(WaveletPacketTransform& wt, float* dwtop) -> void;

auto setDWPTExtension(WaveletPacketTransform& wt, char const* extension) -> void;
auto setDWPTEntropy(WaveletPacketTransform& wt, char const* entropy, float eparam) -> void;
auto getDWPTNodelength(WaveletPacketTransform& wt, int x) -> int;

}  // namespace mc

template<>
struct fmt::formatter<mc::WaveletPacketTransform> : formatter<string_view>
{
    template<typename FormatContext>
    auto format(mc::WaveletPacketTransform const& wt, FormatContext& ctx) const
    {
        fmt::format_to(ctx.out(), "{}\n", wt.wave());
        fmt::format_to(ctx.out(), "Signal Extension : {0} \n\n", wt.ext.c_str());
        fmt::format_to(ctx.out(), "Entropy : {0} \n\n", wt.entropy.c_str());
        fmt::format_to(ctx.out(), "Number of Decomposition Levels {0} \n\n", wt.J);
        fmt::format_to(ctx.out(), "Number of Active Nodes {0} \n\n", wt.nodes);
        fmt::format_to(ctx.out(), "Length of Input Signal {0} \n\n", wt.signalLength());
        fmt::format_to(ctx.out(), "Length of WT Output Vector {0} \n\n", wt.outlength);
        fmt::format_to(ctx.out(), "coefficients are contained in vector : output \n\n");
        fmt::format_to(ctx.out(), "Coefficients Access \n");

        auto const j = wt.J;

        auto it1 = 1;
        auto it2 = 0;
        for (auto i = 0; i < j; ++i) { it1 += mc::ipow2(i + 1); }
        for (auto i = j; i > 0; --i) {
            auto p2 = mc::ipow2(i);
            it1 -= p2;
            for (auto k = 0; k < p2; ++k) {
                if (wt.basisvector[it1 + k] == 1) {
                    fmt::format_to(
                        ctx.out(),
                        "Node {} {} Access : output[{}] Length : {} \n",
                        i,
                        k,
                        it2,
                        wt.length[j - i + 1]
                    );
                    it2 += wt.length[j - i + 1];
                }
            }
        }

        return fmt::format_to(ctx.out(), "\n");
    }
};
