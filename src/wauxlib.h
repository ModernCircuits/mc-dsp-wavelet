/*
Copyright (c) 2017, Rafat Hussain
*/
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
    //double params[0];
};

using denoise_object = struct denoise_set*;
auto denoise_init(int length, int J, const char* wname) -> denoise_object;

void visushrink(double* signal, int N, int J, const char* wname, const char* method, const char* ext, const char* thresh, const char* level, double* denoised);

void sureshrink(double* signal, int N, int J, const char* wname, const char* method, const char* ext, const char* thresh, const char* level, double* denoised);

void modwtshrink(double* signal, int N, int J, const char* wname, const char* cmethod, const char* ext, const char* thresh, double* denoised);

void denoise(denoise_object obj, double* signal, double* denoised);

void setDenoiseMethod(denoise_object obj, const char* dmethod);

void setDenoiseWTMethod(denoise_object obj, const char* wmethod);

void setDenoiseWTExtension(denoise_object obj, const char* extension);

void setDenoiseParameters(denoise_object obj, const char* thresh, const char* level);

void denoise_free(denoise_object object);

void getDWTRecCoeff(const double* coeff, const int* length, const char* ctype, const char* ext, int level, int J, double* lpr,
    double* hpr, int lf, int siglength, double* reccoeff);

auto mad(double* x, int N) -> double;

#endif /* WAUXLIB_H_ */
