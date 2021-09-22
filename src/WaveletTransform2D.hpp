#pragma once

#include "tcb/span.hpp"

#include "Convolution.hpp"
#include "Wavelet.hpp"

#include <string>

struct WaveletTransform2D {
    Wavelet* wave;
    std::string method;
    int rows; // Matrix Number of rows
    int cols; // Matrix Number of columns
    int outlength; // Length of the output DWT vector
    int J; // Number of decomposition Levels
    int MaxIter; // Maximum Iterations J <= MaxIter
    std::string ext; // Type of Extension used - "per" or "sym"
    int coeffaccesslength;

    int N; //
    int* dimensions;
    int* coeffaccess;
    std::unique_ptr<int[]> params;
};

auto wt2Init(Wavelet& wave, char const* method, int rows, int cols, int j) -> WaveletTransform2D*;

auto dwt2(WaveletTransform2D* wt, double* inp) -> std::unique_ptr<double[]>;
auto idwt2(WaveletTransform2D* wt, double* wavecoeff, double* oup) -> void;
auto swt2(WaveletTransform2D* wt, double* inp) -> std::unique_ptr<double[]>;
auto iswt2(WaveletTransform2D* wt, double const* wavecoeffs, double* oup) -> void;
auto modwt2(WaveletTransform2D* wt, double* inp) -> std::unique_ptr<double[]>;
auto imodwt2(WaveletTransform2D* wt, double* wavecoeff, double* oup) -> void;
auto getWT2Coeffs(WaveletTransform2D* wt, double* wcoeffs, int level, char const* type, int* rows, int* cols) -> double*;
auto setDWT2Extension(WaveletTransform2D* wt, char const* extension) -> void;
auto dispWT2Coeffs(double* a, int row, int col) -> void;

auto summary(WaveletTransform2D const& wt) -> void;
auto wt2Free(WaveletTransform2D* wt) -> void;
