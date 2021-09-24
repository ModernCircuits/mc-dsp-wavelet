#pragma once

#include "tcb/span.hpp"

#include "lt/dsp/convolution/ConvolutionMethod.hpp"
#include "lt/dsp/convolution/FFTConvolver.hpp"
#include "lt/dsp/wavelets/SignalExtension.hpp"
#include "lt/dsp/wavelets/Wavelet.hpp"

#include <string>

struct WaveletTransform {
    WaveletTransform(Wavelet& wave, char const* method, std::size_t siglength, std::size_t j);

    [[nodiscard]] auto wave() const noexcept -> Wavelet const& { return *wave_; }
    [[nodiscard]] auto levels() const noexcept -> int { return levels_; }
    [[nodiscard]] auto signalLength() const noexcept -> int { return signalLength_; }
    [[nodiscard]] auto method() const noexcept -> std::string const& { return method_; }

    auto extension(SignalExtension ext) -> void;
    [[nodiscard]] auto extension() const noexcept -> SignalExtension { return ext_; }

    auto convMethod(ConvolutionMethod method) -> void;
    [[nodiscard]] auto convMethod() const noexcept -> ConvolutionMethod { return cmethod_; }

    [[nodiscard]] auto output() const noexcept -> lt::span<double>;
    [[nodiscard]] auto approx() const noexcept -> lt::span<double>;
    [[nodiscard]] auto detail(std::size_t level) const noexcept -> lt::span<double>;

private:
    Wavelet* wave_;
    std::size_t levels_;
    std::size_t signalLength_;
    std::string method_;
    SignalExtension ext_;
    ConvolutionMethod cmethod_;

    double* output_;

public:
    std::unique_ptr<FFTConvolver> convolver;
    int modwtsiglength; // Modified signal length for MODWT
    int outlength; // Length of the output DWT vector
    int lenlength; // Length of the Output Dimension Vector "length"
    int MaxIter; // Maximum Iterations J <= MaxIter
    int even; // even = 1 if signal is of even length. even = 0 otherwise

    int N {}; //
    int cfftset;
    int zpad {};
    int length[102] {};
    std::unique_ptr<double[]> params;
};

auto dwt(WaveletTransform& wt, double const* inp) -> void;
auto idwt(WaveletTransform& wt, double* dwtop) -> void;
auto swt(WaveletTransform& wt, double const* inp) -> void;
auto iswt(WaveletTransform& wt, double* swtop) -> void;
auto modwt(WaveletTransform& wt, double const* inp) -> void;
auto imodwt(WaveletTransform& wt, double* oup) -> void;

auto summary(WaveletTransform const& wt) -> void;
