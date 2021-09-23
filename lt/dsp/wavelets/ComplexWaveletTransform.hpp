#pragma once

#include "tcb/span.hpp"

#include "lt/dsp/convolution/FFTConvolver.hpp"

#include <string>

struct ComplexWaveletTransform {
    ComplexWaveletTransform(char const* wave, double param, int siglength, double dt, int j);

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

    Complex<double>* output;
    double* scale;
    double* period;
    double* coi;
    std::unique_ptr<double[]> params;
};

auto setCWTScales(ComplexWaveletTransform& wt, double s0, double dj, char const* type, int power) -> void;
auto cwt(ComplexWaveletTransform& wt, double const* inp) -> void;
auto icwt(ComplexWaveletTransform& wt, double* cwtop) -> void;

auto summary(ComplexWaveletTransform const& wt) -> void;

auto cwavelet(double const* y, int n, double dt, int mother, double param, double s0, double dj, int jtot, int npad,
    double* wave, double const* scale, double* period, double* coi) -> void;

auto psi0(int mother, double param, double* val, int* real) -> void;

auto cdelta(int mother, double param, double psi0) -> double;

auto icwavelet(double const* wave, int n, double* scale, int jtot, double dt, double dj, double cdelta, double psi0, double* oup) -> void;

auto meyer(int n, double lb, double ub, double* phi, double* psi, double* tgrid) -> void;
auto gauss(int n, int p, double lb, double ub, double* psi, double* t) -> void;
auto mexhat(int n, double lb, double ub, double* psi, double* t) -> void;
auto morlet(int n, double lb, double ub, double* psi, double* t) -> void;