// SPDX-License-Identifier: BSL-1.0

#pragma once

#include <mc/core/config.hpp>

#include <mc/algorithm/signal_extension.hpp>
#include <mc/convolution.hpp>
#include <mc/core/format.hpp>
#include <mc/core/memory.hpp>
#include <mc/core/span.hpp>
#include <mc/core/string.hpp>
#include <mc/wavelet/wavelet.hpp>

namespace mc {

struct WaveletTransform
{
    WaveletTransform(Wavelet& wave, char const* method, size_t siglength, size_t j);

    [[nodiscard]] auto wave() const noexcept -> Wavelet const&;
    [[nodiscard]] auto levels() const noexcept -> int;
    [[nodiscard]] auto signalLength() const noexcept -> size_t;
    [[nodiscard]] auto method() const noexcept -> String const&;

    auto extension(SignalExtension ext) -> void;
    [[nodiscard]] auto extension() const noexcept -> SignalExtension;

    auto convMethod(ConvolutionMethod method) -> void;
    [[nodiscard]] auto convMethod() const noexcept -> ConvolutionMethod;

    [[nodiscard]] auto output() const -> Span<float>;
    [[nodiscard]] auto approx() const -> Span<float>;
    [[nodiscard]] auto detail(size_t level) const -> Span<float>;

private:
    Wavelet* _wave;
    size_t _levels;
    size_t _signalLength;
    String _method;
    SignalExtension _ext;
    ConvolutionMethod _cmethod{ConvolutionMethod::direct};

    float* _output;

public:
    UniquePtr<FFTConvolver> convolver;
    size_t modwtsiglength;  // Modified signal length for MODWT
    size_t outlength;       // Length of the output DWT vector
    size_t lenlength;       // Length of the Output Dimension Vector "length"
    size_t MaxIter;         // Maximum Iterations J <= MaxIter

    size_t N{};  //
    size_t cfftset{0};
    size_t zpad{};
    size_t length[102]{};
    UniquePtr<float[]> params;
};

auto dwt(WaveletTransform& wt, float const* inp) -> void;
auto idwt(WaveletTransform& wt, float* dwtop) -> void;
auto swt(WaveletTransform& wt, float const* inp) -> void;
auto iswt(WaveletTransform& wt, float* swtop) -> void;
auto modwt(WaveletTransform& wt, float const* inp) -> void;
auto imodwt(WaveletTransform& wt, float* oup) -> void;

}  // namespace mc

template<>
struct fmt::formatter<mc::WaveletTransform> : formatter<string_view>
{
    template<typename FormatContext>
    auto format(mc::WaveletTransform const& wt, FormatContext& ctx) const
    {
        auto j = wt.levels();
        fmt::format_to(ctx.out(), "{}\n", wt.wave());
        fmt::format_to(ctx.out(), "Wavelet Transform : {} \n", wt.method());
        fmt::format_to(ctx.out(), "Signal Extension : {} \n", toString(wt.extension()));
        fmt::format_to(ctx.out(), "Convolution Method : {} \n", toString(wt.convMethod()));
        fmt::format_to(ctx.out(), "Number of Decomposition Levels {} \n", wt.levels());
        fmt::format_to(ctx.out(), "Length of Input Signal {} \n", wt.signalLength());
        fmt::format_to(ctx.out(), "Length of WT Output Vector {} \n", wt.outlength);
        fmt::format_to(ctx.out(), "coefficients are contained in vector : {} \n", "output");
        fmt::format_to(ctx.out(), "Approximation Coefficients \n");
        fmt::format_to(
            ctx.out(),
            "Level {} Access : output[{}] Length : {} \n",
            j,
            0,
            wt.length[0]
        );
        fmt::format_to(ctx.out(), "Detail Coefficients \n");
        auto t = wt.length[0];
        for (auto i = 0; i < j; ++i) {
            fmt::format_to(
                ctx.out(),
                "Level {} Access : output[{}] Length : {} \n",
                j - i,
                t,
                wt.length[i + 1]
            );
            t += wt.length[i + 1];
        }
        return fmt::format_to(ctx.out(), "\n");
    }
};
