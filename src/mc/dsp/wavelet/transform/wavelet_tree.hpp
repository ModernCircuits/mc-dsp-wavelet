#pragma once

#include <mc/core/config.hpp>

#include <mc/dsp/convolution/FFTConvolver.hpp>
#include <mc/dsp/wavelet/wavelet.hpp>

#include <mc/core/format.hpp>
#include <mc/core/span.hpp>
#include <mc/core/string.hpp>

namespace mc::dsp {
struct WaveletTree
{
    WaveletTree(Wavelet* waveIn, std::size_t signalLength, std::size_t j);

    auto extension(char const* newExtension) noexcept -> void;
    [[nodiscard]] auto extension() const noexcept -> String const&;

    auto nodeLength(std::size_t x) -> std::size_t;
    auto coeffs(std::size_t x, std::size_t y, float* coeffs, std::size_t n) const -> void;

    Wavelet* wave;
    String method;
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
    UniquePtr<float[]> params;
    unsigned* nodeLength_{nullptr};

private:
    String _ext;
};

auto wtree(WaveletTree& wt, float const* inp) -> void;

}  // namespace mc::dsp

template<>
struct fmt::formatter<mc::dsp::WaveletTree> : formatter<string_view>
{
    template<typename FormatContext>
    auto format(mc::dsp::WaveletTree const& wt, FormatContext& ctx) const
    {
        fmt::format_to(ctx.out(), "\n");
        fmt::format_to(ctx.out(), "{}\n", *wt.wave);
        fmt::format_to(ctx.out(), "Wavelet Transform : {} \n\n", wt.method.c_str());
        fmt::format_to(ctx.out(), "Signal Extension : {} \n\n", wt.extension().c_str());
        fmt::format_to(ctx.out(), "Number of Decomposition Levels {} \n\n", wt.J);
        fmt::format_to(ctx.out(), "Length of Input Signal {} \n\n", wt.siglength);
        fmt::format_to(ctx.out(), "Length of WT Output Vector {} \n\n", wt.outlength);
        fmt::format_to(ctx.out(), "coefficients are contained in vector : output \n\n");
        fmt::format_to(ctx.out(), "Coefficients Access \n");

        auto formatNode = [&](auto i, auto k, auto output, auto length) {
            fmt::format_to(
                ctx.out(),
                "Node {} {} Access : output[{}] Length : {} \n",
                i,
                k,
                output,
                length
            );
        };

        auto t  = 0;
        auto p2 = 2;
        for (std::size_t i = 0; i < wt.J; ++i) {
            for (auto k = 0; k < p2; ++k) {
                formatNode(i + 1, k, wt.nodeLength_[t], wt.length[wt.J - i]);
                t++;
            }
            p2 *= 2;
        }

        return fmt::format_to(ctx.out(), "\n");
    }
};
