#pragma once

#include "tcb/span.hpp"

#include "lt/dsp/convolution/FFTConvolver.hpp"
#include "lt/dsp/wavelets/Wavelet.hpp"

#include <string>

struct WaveletPacketTransform {
    Wavelet* wave;
    FFTConvolver* cobj;
    int siglength; // Length of the original signal.
    int outlength; // Length of the output DWT vector
    int lenlength; // Length of the Output Dimension Vector "length"
    int J; // Number of decomposition Levels
    int MaxIter; // Maximum Iterations J <= MaxIter
    int even; // even = 1 if signal is of even length. even = 0 otherwise
    std::string ext; // Type of Extension used - "per" or "sym"
    std::string entropy;
    double eparam;

    int N; //
    int nodes;
    int length[102];
    double* output;
    double* costvalues;
    double* basisvector;
    int* nodeindex;
    int* numnodeslevel;
    int* coeflength;
    std::unique_ptr<double[]> params;
};

auto wptInit(Wavelet* wave, int siglength, int j) -> WaveletPacketTransform*;

auto dwpt(WaveletPacketTransform* wt, double const* inp) -> void;
auto idwpt(WaveletPacketTransform* wt, double* dwtop) -> void;
auto setDWPTExtension(WaveletPacketTransform* wt, char const* extension) -> void;
auto setDWPTEntropy(WaveletPacketTransform* wt, char const* entropy, double eparam) -> void;
auto getDWPTNodelength(WaveletPacketTransform* wt, int x) -> int;

auto summary(WaveletPacketTransform const& wt) -> void;
auto wptFree(WaveletPacketTransform* object) -> void;
