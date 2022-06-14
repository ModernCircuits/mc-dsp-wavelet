#pragma once

#include "lt/dsp/convolution/ConvolutionMethod.hpp"
#include "lt/dsp/convolution/FFTConvolver.hpp"
#include "lt/dsp/wavelets/SignalExtension.hpp"
#include "lt/dsp/wavelets/Wavelet.hpp"

#include "lt/preprocessor.hpp"
#include "lt/span.hpp"
#include "lt/string.hpp"

namespace lt
{
namespace dsp
{

struct WaveletTransform
{
    WaveletTransform(Wavelet& wave, char const* method, std::size_t siglength, std::size_t j);

    LT_NODISCARD auto wave() const noexcept -> Wavelet const& { return *wave_; }
    LT_NODISCARD auto levels() const noexcept -> int { return static_cast<int>(levels_); }
    LT_NODISCARD auto signalLength() const noexcept -> std::size_t { return signalLength_; }
    LT_NODISCARD auto method() const noexcept -> std::string const& { return method_; }

    auto extension(SignalExtension ext) -> void;
    LT_NODISCARD auto extension() const noexcept -> SignalExtension { return ext_; }

    auto convMethod(ConvolutionMethod method) -> void;
    LT_NODISCARD auto convMethod() const noexcept -> ConvolutionMethod { return cmethod_; }

    LT_NODISCARD auto output() const -> span<float>;
    LT_NODISCARD auto approx() const -> span<float>;
    LT_NODISCARD auto detail(std::size_t level) const -> span<float>;

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

}  // namespace dsp
}  // namespace lt