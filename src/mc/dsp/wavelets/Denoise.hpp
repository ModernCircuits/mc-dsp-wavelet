#pragma once

#include "mc/dsp/wavelets.hpp"

#include <mc/core/cfloat.hpp>
#include <mc/core/cmath.hpp>
#include <mc/core/cstdlib.hpp>
#include <mc/core/cstring.hpp>
#include <mc/core/string.hpp>

namespace mc::dsp
{

struct DenoiseSet
{
    DenoiseSet(int length, int j, char const* name);

    int N{};              // signal length
    int J{};              // Levels of Wavelet decomposition
    std::string wname;    // Wavelet name
    std::string wmethod;  // Wavelet decomposition method - dwt or swt
    std::string cmethod;  // Cnvolution Method - direct or fft . Available only for modwt.
    // SWT and DWT only use direct method.
    std::string ext;      // Signal Extension - sym or per
    std::string thresh;   // thresholding - soft or hard
    std::string level;    // Noise Estimation level - first or all
    std::string dmethod;  // Denoising Method -sureshrink or visushrink
};

auto visushrink(float* signal, std::size_t n, std::size_t j, char const* wname, char const* method, char const* ext,
                char const* thresh, char const* level, float* denoised) -> void;

auto sureshrink(float* signal, std::size_t n, std::size_t j, char const* wname, char const* method, char const* ext,
                char const* thresh, char const* level, float* denoised) -> void;

auto modwtshrink(float* signal, std::size_t n, std::size_t j, char const* wname, char const* cmethod, char const* ext,
                 char const* thresh, float* denoised) -> void;

auto denoise(DenoiseSet& obj, float* signal, float* denoised) -> void;

auto setDenoiseMethod(DenoiseSet& obj, char const* dmethod) -> void;

auto setDenoiseWTMethod(DenoiseSet& obj, char const* wmethod) -> void;

auto setDenoiseWTExtension(DenoiseSet& obj, char const* extension) -> void;

auto setDenoiseParameters(DenoiseSet& obj, char const* thresh, char const* level) -> void;

auto median(float* x, int n) -> float;

auto minindex(float const* arr, int n) -> int;

}  // namespace mc::dsp