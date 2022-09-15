#pragma once

#include <mc/core/config.hpp>

#include <mc/dsp/convolution/ConvolutionMethod.hpp>
#include <mc/dsp/convolution/FFTConvolver.hpp>
#include <mc/dsp/wavelets/SignalExtension.hpp>
#include <mc/dsp/wavelets/Wavelet.hpp>

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
    [[nodiscard]] auto method() const noexcept -> std::string const&;

    auto extension(SignalExtension ext) -> void;
    [[nodiscard]] auto extension() const noexcept -> SignalExtension;

    auto convMethod(ConvolutionMethod method) -> void;
    [[nodiscard]] auto convMethod() const noexcept -> ConvolutionMethod;

    [[nodiscard]] auto output() const -> Span<float>;
    [[nodiscard]] auto approx() const -> Span<float>;
    [[nodiscard]] auto detail(std::size_t level) const -> Span<float>;

private:
    Wavelet* wave_;
    std::size_t levels_;
    std::size_t signalLength_;
    std::string method_;
    SignalExtension ext_;
    ConvolutionMethod cmethod_{ConvolutionMethod::direct};

    float* output_;

public:
    std::unique_ptr<FFTConvolver> convolver;
    std::size_t modwtsiglength;  // Modified signal length for MODWT
    std::size_t outlength;       // Length of the output DWT vector
    std::size_t lenlength;       // Length of the Output Dimension Vector "length"
    std::size_t MaxIter;         // Maximum Iterations J <= MaxIter

    std::size_t N{};  //
    std::size_t cfftset{0};
    std::size_t zpad{};
    std::size_t length[102]{};
    std::unique_ptr<float[]> params;
};

auto dwt(WaveletTransform& wt, float const* inp) -> void;
auto idwt(WaveletTransform& wt, float* dwtop) -> void;
auto swt(WaveletTransform& wt, float const* inp) -> void;
auto iswt(WaveletTransform& wt, float* swtop) -> void;
auto modwt(WaveletTransform& wt, float const* inp) -> void;
auto imodwt(WaveletTransform& wt, float* oup) -> void;

auto summary(WaveletTransform const& wt) -> void;

}  // namespace mc::dsp