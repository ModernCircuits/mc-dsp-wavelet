// SPDX-License-Identifier: BSL-1.0

#include "wavelet_transform_2d.hpp"

#include <mc/dsp/convolution/FFTConvolver.hpp>
#include <mc/dsp/fft/FFT.hpp>
#include <mc/dsp/wavelet/transform/common.hpp>

#include <mc/core/cassert.hpp>
#include <mc/core/cmath.hpp>
#include <mc/core/cstring.hpp>
#include <mc/core/format.hpp>
#include <mc/core/raise.hpp>
#include <mc/core/stdexcept.hpp>
#include <mc/core/string_view.hpp>
#include <mc/core/utility.hpp>

namespace mc::dsp {

namespace {

auto idwtShift(
    int shift,
    int rows,
    int cols,
    float const* lpr,
    float const* hpr,
    int lf,
    float* a,
    float* h,
    float* v,
    float* d,
    float* oup
) -> void
{
    auto const n    = rows > cols ? 2 * rows : 2 * cols;
    auto const dim1 = 2 * rows;
    auto const dim2 = 2 * cols;

    auto xLp = makeZeros<float>(n + 2 * lf - 1);
    auto cL  = makeZeros<float>(dim1 * dim2);
    auto cH  = makeZeros<float>(dim1 * dim2);

    auto ir      = rows;
    auto ic      = cols;
    auto istride = ic;
    auto ostride = 1;

    for (auto i = 0; i < ic; ++i) {
        idwtPerStride(a + i, ir, h + i, lpr, hpr, lf, xLp.get(), istride, ostride);

        for (auto k = lf / 2 - 1; k < 2 * ir + lf / 2 - 1; ++k) {
            cL[(k - lf / 2 + 1) * ic + i] = xLp[k];
        }

        idwtPerStride(v + i, ir, d + i, lpr, hpr, lf, xLp.get(), istride, ostride);

        for (auto k = lf / 2 - 1; k < 2 * ir + lf / 2 - 1; ++k) {
            cH[(k - lf / 2 + 1) * ic + i] = xLp[k];
        }
    }

    ir *= 2;
    istride = 1;
    ostride = 1;

    for (auto i = 0; i < ir; ++i) {
        idwtPerStride(
            cL.get() + i * ic,
            ic,
            cH.get() + i * ic,
            lpr,
            hpr,
            lf,
            xLp.get(),
            istride,
            ostride
        );

        for (auto k = lf / 2 - 1; k < 2 * ic + lf / 2 - 1; ++k) {
            oup[(k - lf / 2 + 1) + i * ic * 2] = xLp[k];
        }
    }

    ic *= 2;

    if (shift == -1) {
        // Save the last column
        for (auto i = 0; i < ir; ++i) { cL[i] = oup[(i + 1) * ic - 1]; }
        // Save the last row
        std::memcpy(cH.get(), oup + (ir - 1) * ic, sizeof(float) * ic);
        for (auto i = ir - 1; i > 0; --i) {
            std::memcpy(oup + i * ic + 1, oup + (i - 1) * ic, sizeof(float) * (ic - 1));
        }
        oup[0] = cL[ir - 1];
        for (auto i = 1; i < ir; ++i) { oup[i * ic] = cL[i - 1]; }

        for (auto i = 1; i < ic; ++i) { oup[i] = cH[i - 1]; }
    }
}

auto imodwtPerStride(
    int m,
    float const* cA,
    int lenCA,
    float const* cD,
    float const* filt,
    int lf,
    float* x,
    int istride,
    int ostride
) -> void
{
    int lenAvg = 0;
    int i      = 0;
    int l      = 0;
    int t      = 0;
    int is     = 0;
    int os     = 0;

    lenAvg = lf;

    for (i = 0; i < lenCA; ++i) {
        t     = i;
        os    = i * ostride;
        is    = t * istride;
        x[os] = (filt[0] * cA[is]) + (filt[lenAvg] * cD[is]);
        for (l = 1; l < lenAvg; l++) {
            t += m;
            while (t >= lenCA) { t -= lenCA; }
            while (t < 0) { t += lenCA; }
            is = t * istride;
            x[os] += (filt[l] * cA[is]) + (filt[lenAvg + l] * cD[is]);
        }
    }
}

}  // namespace

WaveletTransform2D::WaveletTransform2D(
    Wavelet& wave,
    char const* method,
    std::size_t rows,
    std::size_t cols,
    std::size_t j
)
{

    auto const maxRows = maxIterations(rows, wave.size());
    auto const maxCols = maxIterations(cols, wave.size());
    auto const maxIter = std::min(maxRows, maxCols);

    if (j > maxIter) {
        raisef<InvalidArgument>(
            "signal can only be iterated {0} times using this wavelet",
            maxIter
        );
    }

    std::size_t sumacc{0};
    if (j == 1U) {
        sumacc = 4U;
    } else if (j > 1U) {
        sumacc = j * 3U + 1U;
    } else {
        raisef<InvalidArgument>("j = {0}, should be >= 1", j);
    }

    this->params    = makeUnique<int[]>(2U * j + sumacc);
    this->outlength = 0;
    if (method == nullptr) {
        this->ext = "per";
    } else if ((method == StringView{"dwt"}) || (method == StringView{"DWT"})) {
        this->ext = "per";
    } else if ((method == StringView{"swt"}) || (method == StringView{"SWT"})) {
        if ((testSWTlength(rows, j) == 0) || (testSWTlength(cols, j) == 0)) {
            raise<InvalidArgument>(
                "\n For SWT data rows and columns must be a multiple of 2^J. \n"
            );
        }

        this->ext = "per";
    } else if ((method == StringView{"MODWT"}) || (method == StringView{"modwt"})) {
        if (strstr(wave.name().c_str(), "haar") == nullptr) {
            if (strstr(wave.name().c_str(), "db") == nullptr) {
                if (strstr(wave.name().c_str(), "sym") == nullptr) {
                    if (strstr(wave.name().c_str(), "coif") == nullptr) {
                        raise<InvalidArgument>(
                            "\n MODWT is only implemented for orthogonal wavelet "
                            "families - db, sym and coif \n"
                        );
                    }
                }
            }
        }
        this->ext = "per";
    }

    this->_wave             = &wave;
    this->_rows             = rows;
    this->_cols             = cols;
    this->J                 = j;
    this->MaxIter           = maxIter;
    this->_method           = method;
    this->coeffaccesslength = sumacc;

    this->dimensions  = &this->params[0];
    this->coeffaccess = &this->params[2 * j];
    for (std::size_t i = 0; i < (2 * j + sumacc); ++i) { this->params[i] = 0; }
}

auto setDWT2Extension(WaveletTransform2D& wt, char const* extension) -> void
{
    if (wt.method() == StringView{"dwt"}) {
        if (extension == StringView{"sym"}) {
            wt.ext = "sym";
        } else if (extension == StringView{"per"}) {
            wt.ext = "per";
        } else {
            raise<InvalidArgument>("Signal extension can be either per or sym");
        }
    } else if ((wt.method() == StringView{"swt"}) || (wt.method() == "modwt")) {
        if (extension == StringView{"per"}) {
            wt.ext = "per";
        } else {
            raise<InvalidArgument>("Signal extension can only be per");
        }
    }
}

auto dwt(WaveletTransform2D& wt, float const* inp) -> UniquePtr<float[]>
{
    int iter          = 0;
    int n             = 0;
    int rowsI         = 0;
    int colsI         = 0;
    int ir            = 0;
    int ic            = 0;
    int istride       = 0;
    int ostride       = 0;
    int aLL           = 0;
    int aLH           = 0;
    int aHL           = 0;
    int aHH           = 0;
    int cdim          = 0;
    float const* orig = nullptr;

    auto j       = wt.J;
    wt.outlength = 0;

    auto rowsN = wt.rows();
    auto colsN = wt.cols();
    auto lp    = wt.wave().lpd().size();
    auto clen  = j * 3;

    if (wt.ext == StringView{"per"}) {
        auto idx = 2 * j;
        while (idx > 0) {
            rowsN                  = (int)std::ceil((float)rowsN / 2.0F);
            colsN                  = (int)std::ceil((float)colsN / 2.0F);
            wt.dimensions[idx - 1] = colsN;
            wt.dimensions[idx - 2] = rowsN;
            wt.outlength += (rowsN * colsN) * 3;
            idx = idx - 2;
        }
        wt.outlength += (rowsN * colsN);
        n              = wt.outlength;
        auto wavecoeff = makeZeros<float>(wt.outlength);

        orig  = inp;
        ir    = wt.rows();
        ic    = wt.cols();
        colsI = wt.dimensions[2 * j - 1];

        auto lpDn1 = makeZeros<float>(ir * colsI);
        auto hpDn1 = makeZeros<float>(ir * colsI);

        for (iter = 0; iter < j; ++iter) {
            rowsI   = wt.dimensions[2 * j - 2 * iter - 2];
            colsI   = wt.dimensions[2 * j - 2 * iter - 1];
            istride = 1;
            ostride = 1;
            cdim    = rowsI * colsI;
            // Row filtering and column subsampling
            for (auto i = 0; i < ir; ++i) {
                dwtPerStride(
                    orig + i * ic,
                    ic,
                    wt.wave().lpd().data(),
                    wt.wave().hpd().data(),
                    lp,
                    lpDn1.get() + i * colsI,
                    colsI,
                    hpDn1.get() + i * colsI,
                    istride,
                    ostride
                );
            }

            // Column Filtering and Row subsampling
            aHH                      = n - cdim;
            wt.coeffaccess[clen]     = aHH;
            aHL                      = aHH - cdim;
            wt.coeffaccess[clen - 1] = aHL;
            aLH                      = aHL - cdim;
            wt.coeffaccess[clen - 2] = aLH;
            aLL                      = aLH - cdim;

            n -= 3 * cdim;
            ic      = colsI;
            istride = ic;
            ostride = ic;

            for (auto i = 0; i < ic; ++i) {
                dwtPerStride(
                    lpDn1.get() + i,
                    ir,
                    wt.wave().lpd().data(),
                    wt.wave().hpd().data(),
                    lp,
                    wavecoeff.get() + aLL + i,
                    rowsI,
                    wavecoeff.get() + aLH + i,
                    istride,
                    ostride
                );
            }

            for (auto i = 0; i < ic; ++i) {
                dwtPerStride(
                    hpDn1.get() + i,
                    ir,
                    wt.wave().lpd().data(),
                    wt.wave().hpd().data(),
                    lp,
                    wavecoeff.get() + aHL + i,
                    rowsI,
                    wavecoeff.get() + aHH + i,
                    istride,
                    ostride
                );
            }

            ir   = rowsI;
            orig = wavecoeff.get() + aLL;
            clen -= 3;
        }
        wt.coeffaccess[0] = 0;

        return wavecoeff;
    }

    MC_ASSERT(wt.ext == StringView{"sym"});

    auto idx = 2 * j;
    while (idx > 0) {
        rowsN += lp - 2;
        colsN += lp - 2;
        rowsN                  = (int)std::ceil((float)rowsN / 2.0F);
        colsN                  = (int)std::ceil((float)colsN / 2.0F);
        wt.dimensions[idx - 1] = colsN;
        wt.dimensions[idx - 2] = rowsN;
        wt.outlength += (rowsN * colsN) * 3;
        idx = idx - 2;
    }
    wt.outlength += (rowsN * colsN);
    n              = wt.outlength;
    auto wavecoeff = makeZeros<float>(wt.outlength);

    orig  = inp;
    ir    = wt.rows();
    ic    = wt.cols();
    colsI = wt.dimensions[2 * j - 1];

    auto lpDn1 = makeZeros<float>(ir * colsI);
    auto hpDn1 = makeZeros<float>(ir * colsI);

    for (iter = 0; iter < j; ++iter) {
        rowsI   = wt.dimensions[2 * j - 2 * iter - 2];
        colsI   = wt.dimensions[2 * j - 2 * iter - 1];
        istride = 1;
        ostride = 1;
        cdim    = rowsI * colsI;
        // Row filtering and column subsampling
        for (auto i = 0; i < ir; ++i) {
            dwtSymStride(
                orig + i * ic,
                ic,
                wt.wave().lpd().data(),
                wt.wave().hpd().data(),
                lp,
                lpDn1.get() + i * colsI,
                colsI,
                hpDn1.get() + i * colsI,
                istride,
                ostride
            );
        }

        // Column Filtering and Row subsampling
        aHH                      = n - cdim;
        wt.coeffaccess[clen]     = aHH;
        aHL                      = aHH - cdim;
        wt.coeffaccess[clen - 1] = aHL;
        aLH                      = aHL - cdim;
        wt.coeffaccess[clen - 2] = aLH;
        aLL                      = aLH - cdim;
        n -= 3 * cdim;
        ic      = colsI;
        istride = ic;
        ostride = ic;

        for (auto i = 0; i < ic; ++i) {
            dwtSymStride(
                lpDn1.get() + i,
                ir,
                wt.wave().lpd().data(),
                wt.wave().hpd().data(),
                lp,
                wavecoeff.get() + aLL + i,
                rowsI,
                wavecoeff.get() + aLH + i,
                istride,
                ostride
            );
        }

        for (auto i = 0; i < ic; ++i) {
            dwtSymStride(
                hpDn1.get() + i,
                ir,
                wt.wave().lpd().data(),
                wt.wave().hpd().data(),
                lp,
                wavecoeff.get() + aHL + i,
                rowsI,
                wavecoeff.get() + aHH + i,
                istride,
                ostride
            );
        }

        ir   = rowsI;
        orig = wavecoeff.get() + aLL;
        clen -= 3;
    }

    wt.coeffaccess[0] = 0;

    return wavecoeff;
}

auto idwt(WaveletTransform2D& wt, float* wavecoeff, float* oup) -> void
{

    int ir = 0;
    int ic = 0;

    int istride = 0;
    int ostride = 0;
    int iter    = 0;
    int aLL     = 0;
    int aLH     = 0;
    int aHL     = 0;
    int aHH     = 0;
    float* orig = nullptr;

    auto const rows = wt.rows();
    auto const cols = wt.cols();
    auto const j    = wt.J;

    if (wt.ext == StringView{"per"}) {
        auto const n  = rows > cols ? 2 * rows : 2 * cols;
        auto const lf = (wt.wave().lpr().size() + wt.wave().hpr().size()) / 2;

        auto idx  = j;
        auto dim1 = wt.dimensions[0];
        auto dim2 = wt.dimensions[1];
        auto k    = 0;
        while (idx > 0) {
            k += 1;
            dim1 *= 2;
            dim2 *= 2;
            idx--;
        }

        auto xLp = makeZeros<float>(n + 2 * lf - 1);
        auto cL  = makeZeros<float>(dim1 * dim2);
        auto cH  = makeZeros<float>(dim1 * dim2);
        auto out = makeZeros<float>(dim1 * dim2);

        aLL  = wt.coeffaccess[0];
        orig = wavecoeff + aLL;
        for (iter = 0; iter < j; ++iter) {
            ir      = wt.dimensions[2 * iter];
            ic      = wt.dimensions[2 * iter + 1];
            istride = ic;
            ostride = 1;
            aLH     = wt.coeffaccess[iter * 3 + 1];
            aHL     = wt.coeffaccess[iter * 3 + 2];
            aHH     = wt.coeffaccess[iter * 3 + 3];
            for (auto i = 0; i < ic; ++i) {
                idwtPerStride(
                    orig + i,
                    ir,
                    wavecoeff + aLH + i,
                    wt.wave().lpr().data(),
                    wt.wave().hpr().data(),
                    lf,
                    xLp.get(),
                    istride,
                    ostride
                );

                for (k = lf / 2 - 1; cmp_less(k, 2 * ir + lf / 2 - 1); ++k) {
                    cL[(k - lf / 2 + 1) * ic + i] = xLp[k];
                }

                idwtPerStride(
                    wavecoeff + aHL + i,
                    ir,
                    wavecoeff + aHH + i,
                    wt.wave().lpr().data(),
                    wt.wave().hpr().data(),
                    lf,
                    xLp.get(),
                    istride,
                    ostride
                );

                for (k = lf / 2 - 1; cmp_less(k, 2 * ir + lf / 2 - 1); ++k) {
                    cH[(k - lf / 2 + 1) * ic + i] = xLp[k];
                }
            }

            ir *= 2;
            istride = 1;
            ostride = 1;

            for (auto i = 0; i < ir; ++i) {
                idwtPerStride(
                    cL.get() + i * ic,
                    ic,
                    cH.get() + i * ic,
                    wt.wave().lpr().data(),
                    wt.wave().hpr().data(),
                    lf,
                    xLp.get(),
                    istride,
                    ostride
                );

                for (k = lf / 2 - 1; cmp_less(k, 2 * ic + lf / 2 - 1); ++k) {
                    out[(k - lf / 2 + 1) + i * ic * 2] = xLp[k];
                }
            }
            ic *= 2;
            if (iter == j - 1) {
                for (auto i = 0; i < wt.rows(); ++i) {
                    for (k = 0; k < wt.cols(); ++k) {
                        oup[k + i * wt.cols()] = out[k + i * ic];
                    }
                }
            } else {
                for (auto i = 0; i < wt.dimensions[2 * (iter + 1)]; ++i) {
                    for (k = 0; k < wt.dimensions[2 * (iter + 1) + 1]; ++k) {
                        oup[k + i * wt.dimensions[2 * (iter + 1) + 1]] = out[k + i * ic];
                    }
                }
            }

            orig = oup;
        }

        return;
    }
    MC_ASSERT(wt.ext == StringView{"sym"});

    auto const n  = rows > cols ? 2 * rows - 1 : 2 * cols - 1;
    auto const lf = (wt.wave().lpr().size() + wt.wave().hpr().size()) / 2;

    auto idx  = j;
    auto dim1 = wt.dimensions[0];
    auto dim2 = wt.dimensions[1];
    auto k    = 0;
    while (idx > 0) {
        k += 1;
        dim1 *= 2;
        dim2 *= 2;
        idx--;
    }

    auto xLp = makeZeros<float>(n + 2 * lf - 1);
    auto cL  = makeZeros<float>(dim1 * dim2);
    auto cH  = makeZeros<float>(dim1 * dim2);
    auto out = makeZeros<float>(dim1 * dim2);

    aLL  = wt.coeffaccess[0];
    orig = wavecoeff + aLL;
    for (iter = 0; iter < j; ++iter) {
        ir      = wt.dimensions[2 * iter];
        ic      = wt.dimensions[2 * iter + 1];
        istride = ic;
        ostride = 1;
        aLH     = wt.coeffaccess[iter * 3 + 1];
        aHL     = wt.coeffaccess[iter * 3 + 2];
        aHH     = wt.coeffaccess[iter * 3 + 3];
        for (auto i = 0; i < ic; ++i) {
            idwtSymStride(
                orig + i,
                ir,
                wavecoeff + aLH + i,
                wt.wave().lpr().data(),
                wt.wave().hpr().data(),
                lf,
                xLp.get(),
                istride,
                ostride
            );

            for (k = lf - 2; k < 2 * ir; ++k) { cL[(k - lf + 2) * ic + i] = xLp[k]; }

            idwtSymStride(
                wavecoeff + aHL + i,
                ir,
                wavecoeff + aHH + i,
                wt.wave().lpr().data(),
                wt.wave().hpr().data(),
                lf,
                xLp.get(),
                istride,
                ostride
            );

            for (k = lf - 2; k < 2 * ir; ++k) { cH[(k - lf + 2) * ic + i] = xLp[k]; }
        }

        ir *= 2;
        istride = 1;
        ostride = 1;

        for (auto i = 0; i < ir; ++i) {
            idwtSymStride(
                cL.get() + i * ic,
                ic,
                cH.get() + i * ic,
                wt.wave().lpr().data(),
                wt.wave().hpr().data(),
                lf,
                xLp.get(),
                istride,
                ostride
            );

            for (k = lf - 2; k < 2 * ic; ++k) { out[(k - lf + 2) + i * ic * 2] = xLp[k]; }
        }
        ic *= 2;
        if (iter == j - 1) {
            for (auto i = 0; i < wt.rows(); ++i) {
                for (k = 0; k < wt.cols(); ++k) {
                    oup[k + i * wt.cols()] = out[k + i * ic];
                }
            }
        } else {
            for (auto i = 0; i < wt.dimensions[2 * (iter + 1)]; ++i) {
                for (k = 0; k < wt.dimensions[2 * (iter + 1) + 1]; ++k) {
                    oup[k + i * wt.dimensions[2 * (iter + 1) + 1]] = out[k + i * ic];
                }
            }
        }

        orig = oup;
    }
}

auto swt2(WaveletTransform2D& wt, float* inp) -> UniquePtr<float[]>
{
    int j       = 0;
    int iter    = 0;
    int m       = 0;
    int n       = 0;
    int lp      = 0;
    int rowsN   = 0;
    int colsN   = 0;
    int rowsI   = 0;
    int colsI   = 0;
    int ir      = 0;
    int ic      = 0;
    int istride = 0;
    int ostride = 0;
    int aLL     = 0;
    int aLH     = 0;
    int aHL     = 0;
    int aHH     = 0;
    int cdim    = 0;
    int clen    = 0;
    float* orig = nullptr;

    j            = wt.J;
    m            = 1;
    wt.outlength = 0;

    rowsN = wt.rows();
    colsN = wt.cols();
    lp    = wt.wave().lpd().size();
    clen  = j * 3;

    auto idx = 2 * j;
    while (idx > 0) {
        wt.dimensions[idx - 1] = colsN;
        wt.dimensions[idx - 2] = rowsN;
        wt.outlength += (rowsN * colsN) * 3;
        idx = idx - 2;
    }
    wt.outlength += (rowsN * colsN);
    n              = wt.outlength;
    auto wavecoeff = makeZeros<float>(wt.outlength);

    orig  = inp;
    ir    = wt.rows();
    ic    = wt.cols();
    colsI = wt.dimensions[2 * j - 1];

    auto lpDn1 = makeUnique<float[]>(ir * colsI);
    auto hpDn1 = makeUnique<float[]>(ir * colsI);

    for (iter = 0; iter < j; ++iter) {
        if (iter > 0) { m = 2 * m; }
        rowsI   = wt.dimensions[2 * j - 2 * iter - 2];
        colsI   = wt.dimensions[2 * j - 2 * iter - 1];
        istride = 1;
        ostride = 1;
        cdim    = rowsI * colsI;
        // Row filtering and column subsampling
        for (auto i = 0; i < ir; ++i) {
            swtPerStride(
                m,
                orig + i * ic,
                ic,
                wt.wave().lpd().data(),
                wt.wave().hpd().data(),
                lp,
                lpDn1.get() + i * colsI,
                colsI,
                hpDn1.get() + i * colsI,
                istride,
                ostride
            );
        }
        // Column Filtering and Row subsampling
        aHH                      = n - cdim;
        wt.coeffaccess[clen]     = aHH;
        aHL                      = aHH - cdim;
        wt.coeffaccess[clen - 1] = aHL;
        aLH                      = aHL - cdim;
        wt.coeffaccess[clen - 2] = aLH;
        aLL                      = aLH - cdim;

        n -= 3 * cdim;
        ic      = colsI;
        istride = ic;
        ostride = ic;
        for (auto i = 0; i < ic; ++i) {
            swtPerStride(
                m,
                lpDn1.get() + i,
                ir,
                wt.wave().lpd().data(),
                wt.wave().hpd().data(),
                lp,
                wavecoeff.get() + aLL + i,
                rowsI,
                wavecoeff.get() + aLH + i,
                istride,
                ostride
            );
        }

        for (auto i = 0; i < ic; ++i) {
            swtPerStride(
                m,
                hpDn1.get() + i,
                ir,
                wt.wave().lpd().data(),
                wt.wave().hpd().data(),
                lp,
                wavecoeff.get() + aHL + i,
                rowsI,
                wavecoeff.get() + aHH + i,
                istride,
                ostride
            );
        }

        ir   = rowsI;
        orig = wavecoeff.get() + aLL;
        clen -= 3;
    }
    wt.coeffaccess[0] = 0;

    return wavecoeff;
}

auto iswt2(WaveletTransform2D& wt, float const* wavecoeffs, float* oup) -> void
{
    int k    = 0;
    int iter = 0;
    int it2  = 0;
    int it3  = 0;
    int j    = 0;
    int m    = 0;
    int rows = 0;
    int cols = 0;
    int lf   = 0;
    int ir   = 0;
    int ic   = 0;
    int k1   = 0;
    int i1   = 0;

    int aLL   = 0;
    int aLH   = 0;
    int aHL   = 0;
    int aHH   = 0;
    int shift = 0;
    j         = wt.J;
    rows      = wt.rows();
    cols      = wt.cols();
    lf        = wt.wave().lpd().size();

    auto a    = makeZeros<float>((rows + lf) * (cols + lf));
    auto h    = makeZeros<float>((rows + lf) * (cols + lf));
    auto v    = makeZeros<float>((rows + lf) * (cols + lf));
    auto d    = makeZeros<float>((rows + lf) * (cols + lf));
    auto oup1 = makeZeros<float>((rows + lf) * (cols + lf));
    auto oup2 = makeZeros<float>((rows + lf) * (cols + lf));

    aLL = wt.coeffaccess[0];

    for (auto i = 0; i < rows; ++i) {
        for (k = 0; k < cols; ++k) { oup[i * cols + k] = wavecoeffs[aLL + i * cols + k]; }
    }

    for (iter = j; iter > 0; iter--) {
        aLH = wt.coeffaccess[(j - iter) * 3 + 1];
        aHL = wt.coeffaccess[(j - iter) * 3 + 2];
        aHH = wt.coeffaccess[(j - iter) * 3 + 3];
        m   = (int)pow(2.0F, (float)iter - 1);

        for (it2 = 0; it2 < m; ++it2) {
            ir  = 0;
            ic  = 0;
            it3 = 0;
            // oup1
            for (auto i = it2; i < rows; i += 2 * m) {
                ic = 0;
                for (k = it2; k < cols; k += 2 * m) {
                    a[it3] = oup[i * cols + k];
                    h[it3] = wavecoeffs[aLH + i * cols + k];
                    v[it3] = wavecoeffs[aHL + i * cols + k];
                    d[it3] = wavecoeffs[aHH + i * cols + k];
                    it3++;
                    ic++;
                }
                ir++;
            }
            shift = 0;
            idwtShift(
                shift,
                ir,
                ic,
                wt.wave().lpr().data(),
                wt.wave().hpr().data(),
                wt.wave().lpd().size(),
                a.get(),
                h.get(),
                v.get(),
                d.get(),
                oup1.get()
            );
            // oup2
            ir  = 0;
            ic  = 0;
            it3 = 0;
            for (auto i = it2 + m; i < rows; i += 2 * m) {
                ic = 0;
                for (k = it2 + m; k < cols; k += 2 * m) {
                    a[it3] = oup[i * cols + k];
                    h[it3] = wavecoeffs[aLH + i * cols + k];
                    v[it3] = wavecoeffs[aHL + i * cols + k];
                    d[it3] = wavecoeffs[aHH + i * cols + k];
                    it3++;
                    ic++;
                }
                ir++;
            }
            shift = -1;
            idwtShift(
                shift,
                ir,
                ic,
                wt.wave().lpr().data(),
                wt.wave().hpr().data(),
                wt.wave().lpd().size(),
                a.get(),
                h.get(),
                v.get(),
                d.get(),
                oup2.get()
            );
            // Shift oup1 and oup2. Then add them to get A.
            i1 = 0;
            for (auto i = it2; i < rows; i += m) {
                k1 = 0;
                for (k = it2; k < cols; k += m) {
                    oup[i * cols + k]
                        = 0.5F * (oup1[i1 * 2 * ic + k1] + oup2[i1 * 2 * ic + k1]);
                    k1++;
                }
                i1++;
            }
        }
    }
}

auto modwt(WaveletTransform2D& wt, float const* inp) -> UniquePtr<float[]>
{
    int j             = 0;
    int iter          = 0;
    int m             = 0;
    int n             = 0;
    int lp            = 0;
    int rowsN         = 0;
    int colsN         = 0;
    int rowsI         = 0;
    int colsI         = 0;
    int ir            = 0;
    int ic            = 0;
    int istride       = 0;
    int ostride       = 0;
    int aLL           = 0;
    int aLH           = 0;
    int aHL           = 0;
    int aHH           = 0;
    int cdim          = 0;
    int clen          = 0;
    float const* orig = nullptr;

    j            = wt.J;
    m            = 1;
    wt.outlength = 0;

    rowsN = wt.rows();
    colsN = wt.cols();
    lp    = wt.wave().lpd().size();
    clen  = j * 3;

    auto idx = 2 * j;
    while (idx > 0) {
        wt.dimensions[idx - 1] = colsN;
        wt.dimensions[idx - 2] = rowsN;
        wt.outlength += (rowsN * colsN) * 3;
        idx = idx - 2;
    }
    wt.outlength += (rowsN * colsN);
    n              = wt.outlength;
    auto wavecoeff = makeZeros<float>(wt.outlength);
    auto filt      = makeUnique<float[]>(2 * lp);
    auto s         = sqrt(2.0F);
    for (auto i = 0; i < lp; ++i) {
        filt[i]      = wt.wave().lpd()[i] / s;
        filt[lp + i] = wt.wave().hpd()[i] / s;
    }

    orig  = inp;
    ir    = wt.rows();
    ic    = wt.cols();
    colsI = wt.dimensions[2 * j - 1];

    auto lpDn1 = makeUnique<float[]>(ir * colsI);
    auto hpDn1 = makeUnique<float[]>(ir * colsI);

    for (iter = 0; iter < j; ++iter) {
        if (iter > 0) { m = 2 * m; }
        rowsI   = wt.dimensions[2 * j - 2 * iter - 2];
        colsI   = wt.dimensions[2 * j - 2 * iter - 1];
        istride = 1;
        ostride = 1;
        cdim    = rowsI * colsI;
        // Row filtering and column subsampling
        for (auto i = 0; i < ir; ++i) {
            modwtPerStride(
                m,
                orig + i * ic,
                ic,
                filt.get(),
                lp,
                lpDn1.get() + i * colsI,
                colsI,
                hpDn1.get() + i * colsI,
                istride,
                ostride
            );
        }
        // Column Filtering and Row subsampling
        aHH                      = n - cdim;
        wt.coeffaccess[clen]     = aHH;
        aHL                      = aHH - cdim;
        wt.coeffaccess[clen - 1] = aHL;
        aLH                      = aHL - cdim;
        wt.coeffaccess[clen - 2] = aLH;
        aLL                      = aLH - cdim;
        n -= 3 * cdim;
        ic      = colsI;
        istride = ic;
        ostride = ic;
        for (auto i = 0; i < ic; ++i) {
            modwtPerStride(
                m,
                lpDn1.get() + i,
                ir,
                filt.get(),
                lp,
                wavecoeff.get() + aLL + i,
                rowsI,
                wavecoeff.get() + aLH + i,
                istride,
                ostride
            );
        }

        for (auto i = 0; i < ic; ++i) {
            modwtPerStride(
                m,
                hpDn1.get() + i,
                ir,
                filt.get(),
                lp,
                wavecoeff.get() + aHL + i,
                rowsI,
                wavecoeff.get() + aHH + i,
                istride,
                ostride
            );
        }

        ir   = rowsI;
        orig = wavecoeff.get() + aLL;
        clen -= 3;
    }
    wt.coeffaccess[0] = 0;

    return wavecoeff;
}

auto imodwt(WaveletTransform2D& wt, float* wavecoeff, float* oup) -> void
{
    int rows = 0;
    int cols = 0;
    int m    = 0;
    // int N;
    int ir      = 0;
    int ic      = 0;
    int lf      = 0;
    int istride = 0;
    int ostride = 0;
    int iter    = 0;
    int j       = 0;
    int aLL     = 0;
    int aLH     = 0;
    int aHL     = 0;
    int aHH     = 0;
    float* orig = nullptr;

    rows = wt.rows();
    cols = wt.cols();
    j    = wt.J;

    m = (int)pow(2.0F, (float)j - 1.0F);
    // N = rows > cols ? rows : cols;
    lf = (wt.wave().lpr().size() + wt.wave().hpr().size()) / 2;

    auto filt = makeZeros<float>(2 * lf);
    auto s    = sqrt(2.0F);
    for (auto i = 0; i < lf; ++i) {
        filt[i]      = wt.wave().lpd()[i] / s;
        filt[lf + i] = wt.wave().hpd()[i] / s;
    }

    auto cL = makeZeros<float>(rows * cols);
    auto cH = makeZeros<float>(rows * cols);
    aLL     = wt.coeffaccess[0];
    orig    = wavecoeff + aLL;
    for (iter = 0; iter < j; ++iter) {
        if (iter > 0) { m = m / 2; }
        ir      = wt.dimensions[2 * iter];
        ic      = wt.dimensions[2 * iter + 1];
        istride = ic;
        ostride = ic;
        aLH     = wt.coeffaccess[iter * 3 + 1];
        aHL     = wt.coeffaccess[iter * 3 + 2];
        aHH     = wt.coeffaccess[iter * 3 + 3];
        for (auto i = 0; i < ic; ++i) {
            imodwtPerStride(
                m,
                orig + i,
                ir,
                wavecoeff + aLH + i,
                filt.get(),
                lf,
                cL.get() + i,
                istride,
                ostride
            );
            imodwtPerStride(
                m,
                wavecoeff + aHL + i,
                ir,
                wavecoeff + aHH + i,
                filt.get(),
                lf,
                cH.get() + i,
                istride,
                ostride
            );
        }

        istride = 1;
        ostride = 1;

        for (auto i = 0; i < ir; ++i) {
            imodwtPerStride(
                m,
                cL.get() + i * ic,
                ic,
                cH.get() + i * ic,
                filt.get(),
                lf,
                oup + i * ic,
                istride,
                ostride
            );
        }

        orig = oup;
    }
}

auto getWT2Coeffs(
    WaveletTransform2D& wt,
    float* wcoeffs,
    int level,
    char const* type,
    int* rows,
    int* cols
) -> float*
{
    int j      = 0;
    int iter   = 0;
    int t      = 0;
    float* ptr = nullptr;
    j          = wt.J;
    // Error Check

    if (level > j || level < 1) {
        raisef<InvalidArgument>(
            "Error : The data is decomposed into {} levels so the acceptable values of "
            "level are between 1 and {}",
            j,
            j
        );
    }

    if ((strcmp(type, "A") == 0) && level != j) {
        raisef<InvalidArgument>(
            "Approximation Coefficients are only available for level {}",
            j
        );
    }

    if (strcmp(type, "A") == 0) {
        t    = 0;
        iter = t;
    } else if (strcmp(type, "H") == 0) {
        t    = 1;
        iter = t;
    } else if (strcmp(type, "V") == 0) {
        t    = 2;
        iter = t;
    } else if (strcmp(type, "D") == 0) {
        t    = 3;
        iter = t;
    } else {
        raise<InvalidArgument>(
            "Only four types of coefficients are accessible A, H, V and D \n"
        );
    }

    iter += (j - level) * 3;

    ptr   = wcoeffs + wt.coeffaccess[iter];
    *rows = wt.dimensions[2 * (j - level)];
    *cols = wt.dimensions[2 * (j - level) + 1];

    return ptr;
}

auto dispWT2Coeffs(float* a, int row, int col) -> void
{
    print("\n MATRIX Order : {} X {} \n \n", row, col);

    for (auto i = 0; i < row; i++) {
        print("R{}: ", i);
        for (auto j = 0; j < col; j++) { print("{} ", a[i * col + j]); }
        print(":R{} \n", i);
    }
}

}  // namespace mc::dsp
