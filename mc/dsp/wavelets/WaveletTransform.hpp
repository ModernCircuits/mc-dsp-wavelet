#pragma once

#include "mc/dsp/convolution/ConvolutionMethod.hpp"
#include "mc/dsp/convolution/FFTConvolver.hpp"
#include "mc/dsp/wavelets/SignalExtension.hpp"
#include "mc/dsp/wavelets/Wavelet.hpp"

#include "mc/preprocessor.hpp"
#include "mc/span.hpp"
#include "mc/string.hpp"

namespace mc::dsp
{

struct WaveletTransform
{
    WaveletTransform(Wavelet& wave, char const* method, std::size_t siglength, std::size_t j);

    [[nodiscard]] auto wave() const noexcept -> Wavelet const& { return *wave_; }
    [[nodiscard]] auto levels() const noexcept -> int { return static_cast<int>(levels_); }
    [[nodiscard]] auto signalLength() const noexcept -> std::size_t { return signalLength_; }
    [[nodiscard]] auto method() const noexcept -> std::string const& { return method_; }

    auto extension(SignalExtension ext) -> void;
    [[nodiscard]] auto extension() const noexcept -> SignalExtension { return ext_; }

    auto convMethod(ConvolutionMethod method) -> void;
    [[nodiscard]] auto convMethod() const noexcept -> ConvolutionMethod { return cmethod_; }

    [[nodiscard]] auto output() const -> span<float>;
    [[nodiscard]] auto approx() const -> span<float>;
    [[nodiscard]] auto detail(std::size_t level) const -> span<float>;

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