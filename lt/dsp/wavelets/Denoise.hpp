/// \copyright Copyright (c) 2017, Rafat Hussain

#ifndef WAUXLIB_H_
#define WAUXLIB_H_

#include "lt/dsp/wavelets.hpp"

#include "lt/cmath.hpp"
#include <cfloat>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

struct DenoiseSet {
    DenoiseSet(int length, int j, char const* name);

    int N {}; //signal length
    int J {}; // Levels of Wavelet decomposition
    std::string wname; //Wavelet name
    std::string wmethod; //Wavelet decomposition method - dwt or swt
    std::string cmethod; //Cnvolution Method - direct or fft . Available only for modwt.
    // SWT and DWT only use direct method.
    std::string ext; // Signal Extension - sym or per
    std::string thresh; // thresholding - soft or hard
    std::string level; // Noise Estimation level - first or all
    std::string dmethod; //Denoising Method -sureshrink or visushrink
};

auto visushrink(double* signal, int n, int j, char const* wname, char const* method, char const* ext, char const* thresh, char const* level, double* denoised) -> void;

auto sureshrink(double* signal, int n, int j, char const* wname, char const* method, char const* ext, char const* thresh, char const* level, double* denoised) -> void;

auto modwtshrink(double* signal, int n, int j, char const* wname, char const* cmethod, char const* ext, char const* thresh, double* denoised) -> void;

auto denoise(DenoiseSet& obj, double* signal, double* denoised) -> void;

auto setDenoiseMethod(DenoiseSet& obj, char const* dmethod) -> void;

auto setDenoiseWTMethod(DenoiseSet& obj, char const* wmethod) -> void;

auto setDenoiseWTExtension(DenoiseSet& obj, char const* extension) -> void;

auto setDenoiseParameters(DenoiseSet& obj, char const* thresh, char const* level) -> void;

auto median(double* x, int n) -> double;

auto minindex(double const* arr, int n) -> int;

#endif /* WAUXLIB_H_ */
