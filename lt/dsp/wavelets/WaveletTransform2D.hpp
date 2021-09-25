#pragma once

#include "lt/dsp/convolution/FFTConvolver.hpp"
#include "lt/dsp/wavelets/Wavelet.hpp"

#include "lt/span.hpp"
#include "lt/string.hpp"

struct WaveletTransform2D {
    WaveletTransform2D(Wavelet& wave, char const* method, std::size_t rows, std::size_t cols, std::size_t j);

    [[nodiscard]] auto wave() const noexcept -> Wavelet const& { return *wave_; }
    [[nodiscard]] auto method() const noexcept -> std::string const& { return method_; }

    [[nodiscard]] auto rows() const noexcept -> int { return rows_; }
    [[nodiscard]] auto cols() const noexcept -> int { return cols_; }

    int outlength {}; // Length of the output DWT vector
    int J {}; // Number of decomposition Levels
    int MaxIter {}; // Maximum Iterations J <= MaxIter
    std::string ext; // Type of Extension used - "per" or "sym"
    int coeffaccesslength {};

    int N {}; //
    int* dimensions {};
    int* coeffaccess {};
    std::unique_ptr<int[]> params;

private:
    Wavelet* wave_ { nullptr };
    int rows_ { 0 }; // Matrix Number of rows
    int cols_ { 0 }; // Matrix Number of columns
    std::string method_;
};

auto dwt(WaveletTransform2D& wt, double* inp) -> std::unique_ptr<double[]>;
auto idwt(WaveletTransform2D& wt, double* wavecoeff, double* oup) -> void;
auto swt2(WaveletTransform2D& wt, double* inp) -> std::unique_ptr<double[]>;
auto iswt2(WaveletTransform2D& wt, double const* wavecoeffs, double* oup) -> void;
auto modwt(WaveletTransform2D& wt, double* inp) -> std::unique_ptr<double[]>;
auto imodwt(WaveletTransform2D& wt, double* wavecoeff, double* oup) -> void;
auto getWT2Coeffs(WaveletTransform2D& wt, double* wcoeffs, int level, char const* type, int* rows, int* cols) -> double*;
auto setDWT2Extension(WaveletTransform2D& wt, char const* extension) -> void;
auto dispWT2Coeffs(double* a, int row, int col) -> void;

auto summary(WaveletTransform2D const& wt) -> void;
