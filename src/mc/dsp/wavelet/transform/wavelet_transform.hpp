#pragma once

#include <mc/core/config.hpp>

#include <mc/dsp/algorithm/signal_extension.hpp>
#include <mc/dsp/convolution/ConvolutionMethod.hpp>
#include <mc/dsp/convolution/FFTConvolver.hpp>
#include <mc/dsp/wavelet/wavelet.hpp>

#include <mc/core/format.hpp>
#include <mc/core/memory.hpp>
#include <mc/core/span.hpp>
#include <mc/core/string.hpp>

namespace mc::dsp {

struct WaveletTransform
{
    WaveletTransform(
        Wavelet& wave,
        char const* method,
        std::size_t siglength,
        std::size_t j
    );

    [[nodiscard]] auto wave() const noexcept -> Wavelet const&;
    [[nodiscard]] auto levels() const noexcept -> int;
    [[nodiscard]] auto signalLength() const noexcept -> std::size_t;
    [[nodiscard]] auto method() const noexcept -> String const&;

    auto extension(SignalExtension ext) -> void;
    [[nodiscard]] auto extension() const noexcept -> SignalExtension;

    auto convMethod(ConvolutionMethod method) -> void;
    [[nodiscard]] auto convMethod() const noexcept -> ConvolutionMethod;

    [[nodiscard]] auto output() const -> Span<float>;
    [[nodiscard]] auto approx() const -> Span<float>;
    [[nodiscard]] auto detail(std::size_t level) const -> Span<float>;

private:
    Wavelet* _wave;
    std::size_t _levels;
    std::size_t _signalLength;
    String _method;
    SignalExtension _ext;
    ConvolutionMethod _cmethod{ConvolutionMethod::direct};

    float* _output;

public:
    UniquePtr<FFTConvolver> convolver;
    std::size_t modwtsiglength;  // Modified signal length for MODWT
    std::size_t outlength;       // Length of the output DWT vector
    std::size_t lenlength;       // Length of the Output Dimension Vector "length"
    std::size_t MaxIter;         // Maximum Iterations J <= MaxIter

    std::size_t N{};  //
    std::size_t cfftset{0};
    std::size_t zpad{};
    std::size_t length[102]{};
    UniquePtr<float[]> params;
};

auto dwt(WaveletTransform& wt, float const* inp) -> void;
auto idwt(WaveletTransform& wt, float* dwtop) -> void;
auto swt(WaveletTransform& wt, float const* inp) -> void;
auto iswt(WaveletTransform& wt, float* swtop) -> void;
auto modwt(WaveletTransform& wt, float const* inp) -> void;
auto imodwt(WaveletTransform& wt, float* oup) -> void;

}  // namespace mc::dsp

template<>
struct fmt::formatter<mc::dsp::WaveletTransform> : formatter<string_view>
{
    template<typename FormatContext>
    auto format(mc::dsp::WaveletTransform const& wt, FormatContext& ctx) const
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
