#pragma once

#include "lt/dsp/convolution/FFTConvolver.hpp"

#include "lt/preprocessor.hpp"
#include "lt/span.hpp"
#include "lt/string.hpp"

namespace lt::dsp
{

struct ContinuousWaveletTransform
{
    ContinuousWaveletTransform(char const* wave, float param, int siglength, float dt, int j);

    LT_NODISCARD auto wave() const noexcept -> std::string const& { return wave_; }
    auto scales(float s0, float dj, char const* type, int power) -> void;

    int signalLength;
    int J;             // Total Number of Scales
    float s0;          // Smallest scale. It depends on the sampling rate. s0 <= 2 * dt for most wavelets
    float dt;          // Sampling Rate
    float dj;          // Separation between scales. eg., scale = s0 * 2 ^ ( [0:N-1] *dj ) or scale = s0 *[0:N-1] * dj
    std::string type;  // Scale Type - Power or Linear
    int pow;           // Base of Power in case type = pow. Typical value is pow = 2
    int sflag;
    int pflag;
    int npad;
    int mother;
    float m;        // Wavelet parameter param
    float smean{};  // Input Signal mean

    Complex<float>* output;
    float* scale;
    float* period;
    float* coi;
    std::unique_ptr<float[]> params;

private:
    // Wavelet - morl/morlet,paul,dog/dgauss
    std::string wave_;
};

auto cwt(ContinuousWaveletTransform& wt, float const* inp) -> void;
auto icwt(ContinuousWaveletTransform& wt, float* cwtop) -> void;

auto summary(ContinuousWaveletTransform const& wt) -> void;

auto meyer(int n, float lb, float ub, float* phi, float* psi, float* tgrid) -> void;
auto gauss(int n, int p, float lb, float ub, float* psi, float* t) -> void;
auto mexhat(int n, float lb, float ub, float* psi, float* t) -> void;
auto morlet(int n, float lb, float ub, float* psi, float* t) -> void;
}  // namespace lt::dsp