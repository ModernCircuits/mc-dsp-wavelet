#pragma once

#include "tcb/span.hpp"

#include "Convolution.hpp"

#include <string>

struct ComplexWaveletTransform {
    std::string wave; // Wavelet - morl/morlet,paul,dog/dgauss
    int siglength; // Length of Input Data
    int J; // Total Number of Scales
    double s0; // Smallest scale. It depends on the sampling rate. s0 <= 2 * dt for most wavelets
    double dt; // Sampling Rate
    double dj; // Separation between scales. eg., scale = s0 * 2 ^ ( [0:N-1] *dj ) or scale = s0 *[0:N-1] * dj
    std::string type; // Scale Type - Power or Linear
    int pow; // Base of Power in case type = pow. Typical value is pow = 2
    int sflag;
    int pflag;
    int npad;
    int mother;
    double m; // Wavelet parameter param
    double smean; // Input Signal mean

    CplxData* output;
    double* scale;
    double* period;
    double* coi;
    std::unique_ptr<double[]> params;
};

auto cwtInit(char const* wave, double param, int siglength, double dt, int j) -> ComplexWaveletTransform*;

auto setCWTScales(ComplexWaveletTransform* wt, double s0, double dj, char const* type, int power) -> void;
auto cwt(ComplexWaveletTransform* wt, double const* inp) -> void;
auto icwt(ComplexWaveletTransform* wt, double* cwtop) -> void;

auto summary(ComplexWaveletTransform const& wt) -> void;
auto cwtFree(ComplexWaveletTransform* object) -> void;