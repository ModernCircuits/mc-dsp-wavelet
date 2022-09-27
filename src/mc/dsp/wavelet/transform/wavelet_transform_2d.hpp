// SPDX-License-Identifier: BSL-1.0

#pragma once

#include <mc/core/config.hpp>

#include <mc/dsp/convolution.hpp>
#include <mc/dsp/wavelet/wavelet.hpp>

#include <mc/core/format.hpp>
#include <mc/core/span.hpp>
#include <mc/core/string.hpp>

namespace mc {

struct WaveletTransform2D
{
    WaveletTransform2D(
        Wavelet& wave,
        char const* method,
        size_t rows,
        size_t cols,
        size_t j
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

}  // namespace mc

template<>
struct fmt::formatter<mc::WaveletTransform2D> : formatter<string_view>
{
    template<typename FormatContext>
    auto format(mc::WaveletTransform2D const& wt, FormatContext& ctx) const
    {
        fmt::format_to(ctx.out(), "{}\n", wt.wave());
        fmt::format_to(ctx.out(), "Wavelet Transform : {} \n\n", wt.method().c_str());
        fmt::format_to(ctx.out(), "Signal Extension : {} \n\n", wt.ext.c_str());
        fmt::format_to(ctx.out(), "Number of Decomposition Levels {} \n\n", wt.J);
        fmt::format_to(ctx.out(), "Input Signal Rows {} \n\n", wt.rows());
        fmt::format_to(ctx.out(), "Input Signal Cols {} \n\n", wt.cols());
        fmt::format_to(ctx.out(), "Length of Coefficients Vector {} \n\n", wt.outlength);

        auto j = wt.J;
        auto t = 0;
        for (auto i = j; i > 0; --i) {
            auto rows  = wt.dimensions[2 * (j - i)];
            auto cols  = wt.dimensions[2 * (j - i) + 1];
            auto vsize = rows * cols;
            fmt::format_to(
                ctx.out(),
                "Level {} Decomposition Rows :{} Columns:{} Vector Size (Rows*Cols):{} \n",
                i,
                rows,
                cols,
                vsize
            );
            fmt::format_to(
                ctx.out(),
                "Access Row values stored at wt.dimensions[{}]\n",
                2 * (j - i)
            );
            fmt::format_to(
                ctx.out(),
                "Access Column values stored at wt.dimensions[{}]\n\n",
                2 * (j - i) + 1
            );

            if (i == j) {
                fmt::format_to(
                    ctx.out(),
                    "Approximation Coefficients access at wt.coeffaccess[{}]={}, Vector "
                    "size:{} \n",
                    t,
                    wt.coeffaccess[t],
                    vsize
                );
            }

            t += 1;
            fmt::format_to(
                ctx.out(),
                "Horizontal Coefficients access at wt.coeffaccess[{}]={}, Vector size:{} "
                "\n",
                t,
                wt.coeffaccess[t],
                vsize
            );
            t += 1;
            fmt::format_to(
                ctx.out(),
                "Vertical Coefficients access at wt.coeffaccess[{}]={}, Vector size:{} \n",
                t,
                wt.coeffaccess[t],
                vsize
            );
            t += 1;
            fmt::format_to(
                ctx.out(),
                "Diagonal Coefficients access at wt.coeffaccess[{}]={}, Vector size:{} "
                "\n\n",
                t,
                wt.coeffaccess[t],
                vsize
            );
        }

        return fmt::format_to(ctx.out(), "\n");
    }
};
