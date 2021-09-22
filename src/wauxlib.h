/// \copyright Copyright (c) 2017, Rafat Hussain

#ifndef WAUXLIB_H_
#define WAUXLIB_H_

#include "wavelib.h"

#include <string>

struct DenoiseSet {
    int N; //signal length
    int J; // Levels of Wavelet decomposition
    std::string wname; //Wavelet name
    std::string wmethod; //Wavelet decomposition method - dwt or swt
    std::string cmethod; //Cnvolution Method - direct or fft . Available only for modwt.
    // SWT and DWT only use direct method.
    std::string ext; // Signal Extension - sym or per
    std::string thresh; // thresholding - soft or hard
    std::string level; // Noise Estimation level - first or all
    std::string dmethod; //Denoising Method -sureshrink or visushrink
};

auto denoiseInit(int length, int j, char const* wname) -> DenoiseSet*;

void visushrink(double* signal, int n, int j, char const* wname, char const* method, char const* ext, char const* thresh, char const* level, double* denoised);

void sureshrink(double* signal, int n, int j, char const* wname, char const* method, char const* ext, char const* thresh, char const* level, double* denoised);

void modwtshrink(double* signal, int n, int j, char const* wname, char const* cmethod, char const* ext, char const* thresh, double* denoised);

void denoise(DenoiseSet* obj, double* signal, double* denoised);

void setDenoiseMethod(DenoiseSet* obj, char const* dmethod);

void setDenoiseWTMethod(DenoiseSet* obj, char const* wmethod);

void setDenoiseWTExtension(DenoiseSet* obj, char const* extension);

void setDenoiseParameters(DenoiseSet* obj, char const* thresh, char const* level);

void denoiseFree(DenoiseSet* object);

auto mad(double* x, int n) -> double;

#endif /* WAUXLIB_H_ */
