#pragma once

#include "tcb/span.hpp"

#include "wavelets/Convolution.hpp"
#include "wavelets/Wavelet.hpp"

#include <string>

struct WaveletTree {
    Wavelet* wave;
    Convolution* cobj;
    std::string method;
    int siglength; // Length of the original signal.
    int outlength; // Length of the output DWT vector
    int lenlength; // Length of the Output Dimension Vector "length"
    int J; // Number of decomposition Levels
    int MaxIter; // Maximum Iterations J <= MaxIter
    int even; // even = 1 if signal is of even length. even = 0 otherwise
    std::string ext; // Type of Extension used - "per" or "sym"

    int N; //
    int nodes;
    int cfftset;
    int zpad;
    int length[102];
    double* output;
    int* nodelength;
    int* coeflength;
    std::unique_ptr<double[]> params;
};

auto wtreeInit(Wavelet* wave, int siglength, int j) -> WaveletTree*;

auto wtree(WaveletTree* wt, double const* inp) -> void;
auto setWTREEExtension(WaveletTree* wt, char const* extension) -> void;
auto getWTREENodelength(WaveletTree* wt, int x) -> int;
auto getWTREECoeffs(WaveletTree* wt, int x, int y, double* coeffs, int n) -> void;

auto summary(WaveletTree const& wt) -> void;
auto wtreeFree(WaveletTree* object) -> void;
