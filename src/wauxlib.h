/// \copyright Copyright (c) 2017, Rafat Hussain

#ifndef WAUXLIB_H_
#define WAUXLIB_H_

#include "wavelib.h"

struct denoise_set {
    int N; //signal length
    int J; // Levels of Wavelet decomposition
    char wname[10]; //Wavelet name
    char wmethod[10]; //Wavelet decomposition method - dwt or swt
    char cmethod[10]; //Cnvolution Method - direct or fft . Available only for modwt.
    // SWT and DWT only use direct method.
    char ext[10]; // Signal Extension - sym or per
    char thresh[10]; // thresholding - soft or hard
    char level[10]; // Noise Estimation level - first or all
    char dmethod[20]; //Denoising Method -sureshrink or visushrink
};

auto denoise_init(int length, int J, char const* wname) -> denoise_set*;

void visushrink(double* signal, int N, int J, char const* wname, char const* method, char const* ext, char const* thresh, char const* level, double* denoised);

void sureshrink(double* signal, int N, int J, char const* wname, char const* method, char const* ext, char const* thresh, char const* level, double* denoised);

void modwtshrink(double* signal, int N, int J, char const* wname, char const* cmethod, char const* ext, char const* thresh, double* denoised);

void denoise(denoise_set* obj, double* signal, double* denoised);

void setDenoiseMethod(denoise_set* obj, char const* dmethod);

void setDenoiseWTMethod(denoise_set* obj, char const* wmethod);

void setDenoiseWTExtension(denoise_set* obj, char const* extension);

void setDenoiseParameters(denoise_set* obj, char const* thresh, char const* level);

void denoise_free(denoise_set* object);

void getDWTRecCoeff(double const* coeff, int const* length, char const* ctype, char const* ext, int level, int J, double* lpr,
    double* hpr, int lf, int siglength, double* reccoeff);

auto mad(double* x, int N) -> double;

#endif /* WAUXLIB_H_ */
