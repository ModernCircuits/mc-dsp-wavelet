#pragma once

#include "tcb/span.hpp"

#include "lt/dsp/convolution/FFTConvolver.hpp"
#include "lt/dsp/wavelets/Wavelet.hpp"

#include <string>

struct WaveletTree {
    WaveletTree(Wavelet* wave, int signalLength, int j);

    auto extension(char const* newExtension) noexcept -> void;
    [[nodiscard]] auto extension() const noexcept -> std::string const&;

    auto nodeLength(int x) -> int;
    auto coeffs(int x, int y, double* coeffs, int n) const -> void;

    Wavelet* wave;
    std::string method;
    int siglength; // Length of the original signal.
    int outlength; // Length of the output DWT vector
    int lenlength; // Length of the Output Dimension Vector "length"
    int J; // Number of decomposition Levels
    int MaxIter; // Maximum Iterations J <= MaxIter
    int even; // even = 1 if signal is of even length. even = 0 otherwise

    int N{}; //
    int nodes;
    int cfftset;
    int zpad{};
    int length[102]{};
    double* output;
    int* coeflength;
    std::unique_ptr<double[]> params;
    int* nodeLength_ { nullptr };

private:
    std::string ext_;
};

auto wtree(WaveletTree& wt, double const* inp) -> void;
auto summary(WaveletTree const& wt) -> void;
