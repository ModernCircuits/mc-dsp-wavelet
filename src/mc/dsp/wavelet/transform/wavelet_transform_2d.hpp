#pragma once

#include <mc/core/config.hpp>

#include <mc/dsp/convolution/FFTConvolver.hpp>
#include <mc/dsp/wavelet/wavelet.hpp>

#include <mc/core/span.hpp>
#include <mc/core/string.hpp>

namespace mc::dsp {

struct WaveletTransform2D
{
    WaveletTransform2D(
        Wavelet& wave,
        char const* method,
        std::size_t rows,
        std::size_t cols,
        std::size_t j
    );

    [[nodiscard]] auto wave() const noexcept -> Wavelet const& { return *_wave; }

    [[nodiscard]] auto method() const noexcept -> String const& { return _method; }

    [[nodiscard]] auto rows() const noexcept -> int { return _rows; }

    [[nodiscard]] auto cols() const noexcept -> int { return _cols; }

    int outlength{};  // Length of the output DWT vector
    int J{};          // Number of decomposition Levels
    int MaxIter{};    // Maximum Iterations J <= MaxIter
    String ext;       // Type of Extension used - "per" or "sym"
    int coeffaccesslength{};

    int N{};  //
    int* dimensions{};
    int* coeffaccess{};
    UniquePtr<int[]> params;

private:
    Wavelet* _wave{nullptr};
    int _rows{0};  // Matrix Number of rows
    int _cols{0};  // Matrix Number of columns
    String _method;
};

auto dwt(WaveletTransform2D& wt, float const* inp) -> UniquePtr<float[]>;
auto idwt(WaveletTransform2D& wt, float* wavecoeff, float* oup) -> void;
auto swt2(WaveletTransform2D& wt, float* inp) -> UniquePtr<float[]>;
auto iswt2(WaveletTransform2D& wt, float const* wavecoeffs, float* oup) -> void;
auto modwt(WaveletTransform2D& wt, float const* inp) -> UniquePtr<float[]>;
auto imodwt(WaveletTransform2D& wt, float* wavecoeff, float* oup) -> void;
auto getWT2Coeffs(
    WaveletTransform2D& wt,
    float* wcoeffs,
    int level,
    char const* type,
    int* rows,
    int* cols
) -> float*;
auto setDWT2Extension(WaveletTransform2D& wt, char const* extension) -> void;
auto dispWT2Coeffs(float* a, int row, int col) -> void;

[[nodiscard]] auto summary(WaveletTransform2D const& wt) -> String;
}  // namespace mc::dsp
