#include "WaveletTransform2D.hpp"

#include "lt/dsp/convolution/FFTConvolver.hpp"
#include "lt/dsp/fft/FFT.hpp"
#include "lt/dsp/wavelets/common.hpp"

#include "lt/cassert.hpp"
#include "lt/cmath.hpp"
#include <cstring>
#include <string_view>

using namespace std::string_view_literals;

namespace {

auto idwtShift(int shift, int rows, int cols, double const* lpr, double const* hpr, int lf, double* a, double* h, double* v, double* d, double* oup) -> void
{
    auto const n = rows > cols ? 2 * rows : 2 * cols;
    auto const dim1 = 2 * rows;
    auto const dim2 = 2 * cols;

    auto xLp = makeZeros<double>(n + 2 * lf - 1);
    auto cL = makeZeros<double>(dim1 * dim2);
    auto cH = makeZeros<double>(dim1 * dim2);

    auto ir = rows;
    auto ic = cols;
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
        idwtPerStride(cL.get() + i * ic, ic, cH.get() + i * ic, lpr, hpr, lf, xLp.get(), istride, ostride);

        for (auto k = lf / 2 - 1; k < 2 * ic + lf / 2 - 1; ++k) {
            oup[(k - lf / 2 + 1) + i * ic * 2] = xLp[k];
        }
    }

    ic *= 2;

    if (shift == -1) {
        //Save the last column
        for (auto i = 0; i < ir; ++i) {
            cL[i] = oup[(i + 1) * ic - 1];
        }
        // Save the last row
        std::memcpy(cH.get(), oup + (ir - 1) * ic, sizeof(double) * ic);
        for (auto i = ir - 1; i > 0; --i) {
            std::memcpy(oup + i * ic + 1, oup + (i - 1) * ic, sizeof(double) * (ic - 1));
        }
        oup[0] = cL[ir - 1];
        for (auto i = 1; i < ir; ++i) {
            oup[i * ic] = cL[i - 1];
        }

        for (auto i = 1; i < ic; ++i) {
            oup[i] = cH[i - 1];
        }
    }
}

auto imodwtPerStride(int m, double const* cA, int lenCA, double const* cD, double const* filt, int lf, double* x, int istride, int ostride) -> void
{
    int lenAvg = 0;
    int i = 0;
    int l = 0;
    int t = 0;
    int is = 0;
    int os = 0;

    lenAvg = lf;

    for (i = 0; i < lenCA; ++i) {
        t = i;
        os = i * ostride;
        is = t * istride;
        x[os] = (filt[0] * cA[is]) + (filt[lenAvg] * cD[is]);
        for (l = 1; l < lenAvg; l++) {
            t += m;
            while (t >= lenCA) {
                t -= lenCA;
            }
            while (t < 0) {
                t += lenCA;
            }
            is = t * istride;
            x[os] += (filt[l] * cA[is]) + (filt[lenAvg + l] * cD[is]);
        }
    }
}

}

auto wt2Init(Wavelet& wave, char const* method, int rows, int cols, int j) -> WaveletTransform2D*
{

    auto const size = wave.size();

    auto const maxRows = wmaxiter(rows, size);
    auto const maxCols = wmaxiter(cols, size);

    auto const maxIter = (maxRows < maxCols) ? maxRows : maxCols;

    if (j > maxIter) {
        printf("\n Error - The Signal Can only be iterated %d times using this wavelet. Exiting\n", maxIter);
        exit(-1);
    }

    int sumacc { 0 };
    if (j == 1) {
        sumacc = 4;
    } else if (j > 1) {
        sumacc = j * 3 + 1;
    } else {
        printf("Error : J should be >= 1 \n");
        exit(-1);
    }

    auto obj = std::make_unique<WaveletTransform2D>();
    obj->params = std::make_unique<int[]>(2 * j + sumacc);
    obj->outlength = 0;
    if (method == nullptr) {
        obj->ext = "per";
    } else if ((method == "dwt"sv) || (method == "DWT"sv)) {
        obj->ext = "per";
    } else if ((method == "swt"sv) || (method == "SWT"sv)) {
        if ((testSWTlength(rows, j) == 0) || (testSWTlength(cols, j) == 0)) {
            printf("\n For SWT data rows and columns must be a multiple of 2^J. \n");
            exit(-1);
        }

        obj->ext = "per";
    } else if ((method == "modwt"sv) || (method == "MODWT"sv)) {
        if (strstr(wave.name().c_str(), "haar") == nullptr) {
            if (strstr(wave.name().c_str(), "db") == nullptr) {
                if (strstr(wave.name().c_str(), "sym") == nullptr) {
                    if (strstr(wave.name().c_str(), "coif") == nullptr) {
                        printf("\n MODWT is only implemented for orthogonal wavelet families - db, sym and coif \n");
                        exit(-1);
                    }
                }
            }
        }
        obj->ext = "per";
    }

    obj->wave = &wave;
    obj->rows = rows;
    obj->cols = cols;
    obj->J = j;
    obj->MaxIter = maxIter;
    obj->method = method;
    obj->coeffaccesslength = sumacc;

    obj->dimensions = &obj->params[0];
    obj->coeffaccess = &obj->params[2 * j];
    for (auto i = 0; i < (2 * j + sumacc); ++i) {
        obj->params[i] = 0;
    }

    return obj.release();
}

auto setDWT2Extension(WaveletTransform2D* wt, char const* extension) -> void
{
    if (wt->method == "dwt"sv) {
        if (extension == "sym"sv) {
            wt->ext = "sym";
        } else if (extension == "per"sv) {
            wt->ext = "per";
        } else {
            printf("Signal extension can be either per or sym");
            exit(-1);
        }
    } else if ((wt->method == "swt"sv) || (wt->method == "modwt")) {
        if (extension == "per"sv) {
            wt->ext = "per";
        } else {
            printf("Signal extension can only be per");
            exit(-1);
        }
    }
}

auto dwt(WaveletTransform2D* wt, double* inp) -> std::unique_ptr<double[]>
{
    int iter = 0;
    int n = 0;
    int rowsI = 0;
    int colsI = 0;
    int ir = 0;
    int ic = 0;
    int istride = 0;
    int ostride = 0;
    int aLL = 0;
    int aLH = 0;
    int aHL = 0;
    int aHH = 0;
    int cdim = 0;
    double* orig = nullptr;

    auto j = wt->J;
    wt->outlength = 0;

    auto rowsN = wt->rows;
    auto colsN = wt->cols;
    auto lp = wt->wave->lpdLen();
    auto clen = j * 3;

    if (wt->ext == "per"sv) {
        auto idx = 2 * j;
        while (idx > 0) {
            rowsN = (int)ceil((double)rowsN / 2.0);
            colsN = (int)ceil((double)colsN / 2.0);
            wt->dimensions[idx - 1] = colsN;
            wt->dimensions[idx - 2] = rowsN;
            wt->outlength += (rowsN * colsN) * 3;
            idx = idx - 2;
        }
        wt->outlength += (rowsN * colsN);
        n = wt->outlength;
        auto wavecoeff = makeZeros<double>(wt->outlength);

        orig = inp;
        ir = wt->rows;
        ic = wt->cols;
        colsI = wt->dimensions[2 * j - 1];

        auto lpDn1 = makeZeros<double>(ir * colsI);
        auto hpDn1 = makeZeros<double>(ir * colsI);

        for (iter = 0; iter < j; ++iter) {
            rowsI = wt->dimensions[2 * j - 2 * iter - 2];
            colsI = wt->dimensions[2 * j - 2 * iter - 1];
            istride = 1;
            ostride = 1;
            cdim = rowsI * colsI;
            // Row filtering and column subsampling
            for (auto i = 0; i < ir; ++i) {
                dwtPerStride(orig + i * ic, ic, wt->wave->lpd(), wt->wave->hpd(), lp, lpDn1.get() + i * colsI, colsI, hpDn1.get() + i * colsI, istride, ostride);
            }

            // Column Filtering and Row subsampling
            aHH = n - cdim;
            wt->coeffaccess[clen] = aHH;
            aHL = aHH - cdim;
            wt->coeffaccess[clen - 1] = aHL;
            aLH = aHL - cdim;
            wt->coeffaccess[clen - 2] = aLH;
            aLL = aLH - cdim;

            n -= 3 * cdim;
            ic = colsI;
            istride = ic;
            ostride = ic;

            for (auto i = 0; i < ic; ++i) {
                dwtPerStride(lpDn1.get() + i, ir, wt->wave->lpd(), wt->wave->hpd(), lp, wavecoeff.get() + aLL + i, rowsI, wavecoeff.get() + aLH + i, istride, ostride);
            }

            for (auto i = 0; i < ic; ++i) {
                dwtPerStride(hpDn1.get() + i, ir, wt->wave->lpd(), wt->wave->hpd(), lp, wavecoeff.get() + aHL + i, rowsI, wavecoeff.get() + aHH + i, istride, ostride);
            }

            ir = rowsI;
            orig = wavecoeff.get() + aLL;
            clen -= 3;
        }
        wt->coeffaccess[0] = 0;

        return wavecoeff;
    }

    LT_ASSERT(wt->ext == "sym"sv);

    auto idx = 2 * j;
    while (idx > 0) {
        rowsN += lp - 2;
        colsN += lp - 2;
        rowsN = (int)ceil((double)rowsN / 2.0);
        colsN = (int)ceil((double)colsN / 2.0);
        wt->dimensions[idx - 1] = colsN;
        wt->dimensions[idx - 2] = rowsN;
        wt->outlength += (rowsN * colsN) * 3;
        idx = idx - 2;
    }
    wt->outlength += (rowsN * colsN);
    n = wt->outlength;
    auto wavecoeff = makeZeros<double>(wt->outlength);

    orig = inp;
    ir = wt->rows;
    ic = wt->cols;
    colsI = wt->dimensions[2 * j - 1];

    auto lpDn1 = makeZeros<double>(ir * colsI);
    auto hpDn1 = makeZeros<double>(ir * colsI);

    for (iter = 0; iter < j; ++iter) {
        rowsI = wt->dimensions[2 * j - 2 * iter - 2];
        colsI = wt->dimensions[2 * j - 2 * iter - 1];
        istride = 1;
        ostride = 1;
        cdim = rowsI * colsI;
        // Row filtering and column subsampling
        for (auto i = 0; i < ir; ++i) {
            dwtSymStride(orig + i * ic, ic, wt->wave->lpd(), wt->wave->hpd(), lp, lpDn1.get() + i * colsI, colsI, hpDn1.get() + i * colsI, istride, ostride);
        }

        // Column Filtering and Row subsampling
        aHH = n - cdim;
        wt->coeffaccess[clen] = aHH;
        aHL = aHH - cdim;
        wt->coeffaccess[clen - 1] = aHL;
        aLH = aHL - cdim;
        wt->coeffaccess[clen - 2] = aLH;
        aLL = aLH - cdim;
        n -= 3 * cdim;
        ic = colsI;
        istride = ic;
        ostride = ic;

        for (auto i = 0; i < ic; ++i) {
            dwtSymStride(lpDn1.get() + i, ir, wt->wave->lpd(), wt->wave->hpd(), lp, wavecoeff.get() + aLL + i, rowsI, wavecoeff.get() + aLH + i, istride, ostride);
        }

        for (auto i = 0; i < ic; ++i) {
            dwtSymStride(hpDn1.get() + i, ir, wt->wave->lpd(), wt->wave->hpd(), lp, wavecoeff.get() + aHL + i, rowsI, wavecoeff.get() + aHH + i, istride, ostride);
        }

        ir = rowsI;
        orig = wavecoeff.get() + aLL;
        clen -= 3;
    }

    wt->coeffaccess[0] = 0;

    return wavecoeff;
}

auto idwt(WaveletTransform2D* wt, double* wavecoeff, double* oup) -> void
{

    int ir = 0;
    int ic = 0;

    int istride = 0;
    int ostride = 0;
    int iter = 0;
    int aLL = 0;
    int aLH = 0;
    int aHL = 0;
    int aHH = 0;
    double* orig = nullptr;

    auto const rows = wt->rows;
    auto const cols = wt->cols;
    auto const j = wt->J;

    if (wt->ext == "per"sv) {
        auto const n = rows > cols ? 2 * rows : 2 * cols;
        auto const lf = (wt->wave->lprLen() + wt->wave->hprLen()) / 2;

        auto idx = j;
        auto dim1 = wt->dimensions[0];
        auto dim2 = wt->dimensions[1];
        auto k = 0;
        while (idx > 0) {
            k += 1;
            dim1 *= 2;
            dim2 *= 2;
            idx--;
        }

        auto xLp = makeZeros<double>(n + 2 * lf - 1);
        auto cL = makeZeros<double>(dim1 * dim2);
        auto cH = makeZeros<double>(dim1 * dim2);
        auto out = makeZeros<double>(dim1 * dim2);

        aLL = wt->coeffaccess[0];
        orig = wavecoeff + aLL;
        for (iter = 0; iter < j; ++iter) {
            ir = wt->dimensions[2 * iter];
            ic = wt->dimensions[2 * iter + 1];
            istride = ic;
            ostride = 1;
            aLH = wt->coeffaccess[iter * 3 + 1];
            aHL = wt->coeffaccess[iter * 3 + 2];
            aHH = wt->coeffaccess[iter * 3 + 3];
            for (auto i = 0; i < ic; ++i) {
                idwtPerStride(orig + i, ir, wavecoeff + aLH + i, wt->wave->lpr(), wt->wave->hpr(), lf, xLp.get(), istride, ostride);

                for (k = lf / 2 - 1; k < 2 * ir + lf / 2 - 1; ++k) {
                    cL[(k - lf / 2 + 1) * ic + i] = xLp[k];
                }

                idwtPerStride(wavecoeff + aHL + i, ir, wavecoeff + aHH + i, wt->wave->lpr(), wt->wave->hpr(), lf, xLp.get(), istride, ostride);

                for (k = lf / 2 - 1; k < 2 * ir + lf / 2 - 1; ++k) {
                    cH[(k - lf / 2 + 1) * ic + i] = xLp[k];
                }
            }

            ir *= 2;
            istride = 1;
            ostride = 1;

            for (auto i = 0; i < ir; ++i) {
                idwtPerStride(cL.get() + i * ic, ic, cH.get() + i * ic, wt->wave->lpr(), wt->wave->hpr(), lf, xLp.get(), istride, ostride);

                for (k = lf / 2 - 1; k < 2 * ic + lf / 2 - 1; ++k) {
                    out[(k - lf / 2 + 1) + i * ic * 2] = xLp[k];
                }
            }
            ic *= 2;
            if (iter == j - 1) {
                for (auto i = 0; i < wt->rows; ++i) {
                    for (k = 0; k < wt->cols; ++k) {
                        oup[k + i * wt->cols] = out[k + i * ic];
                    }
                }
            } else {
                for (auto i = 0; i < wt->dimensions[2 * (iter + 1)]; ++i) {
                    for (k = 0; k < wt->dimensions[2 * (iter + 1) + 1]; ++k) {
                        oup[k + i * wt->dimensions[2 * (iter + 1) + 1]] = out[k + i * ic];
                    }
                }
            }

            orig = oup;
        }

        return;
    }
    LT_ASSERT(wt->ext == "sym"sv);

    auto const n = rows > cols ? 2 * rows - 1 : 2 * cols - 1;
    auto const lf = (wt->wave->lprLen() + wt->wave->hprLen()) / 2;

    auto idx = j;
    auto dim1 = wt->dimensions[0];
    auto dim2 = wt->dimensions[1];
    auto k = 0;
    while (idx > 0) {
        k += 1;
        dim1 *= 2;
        dim2 *= 2;
        idx--;
    }

    auto xLp = makeZeros<double>(n + 2 * lf - 1);
    auto cL = makeZeros<double>(dim1 * dim2);
    auto cH = makeZeros<double>(dim1 * dim2);
    auto out = makeZeros<double>(dim1 * dim2);

    aLL = wt->coeffaccess[0];
    orig = wavecoeff + aLL;
    for (iter = 0; iter < j; ++iter) {
        ir = wt->dimensions[2 * iter];
        ic = wt->dimensions[2 * iter + 1];
        istride = ic;
        ostride = 1;
        aLH = wt->coeffaccess[iter * 3 + 1];
        aHL = wt->coeffaccess[iter * 3 + 2];
        aHH = wt->coeffaccess[iter * 3 + 3];
        for (auto i = 0; i < ic; ++i) {
            idwtSymStride(orig + i, ir, wavecoeff + aLH + i, wt->wave->lpr(), wt->wave->hpr(), lf, xLp.get(), istride, ostride);

            for (k = lf - 2; k < 2 * ir; ++k) {
                cL[(k - lf + 2) * ic + i] = xLp[k];
            }

            idwtSymStride(wavecoeff + aHL + i, ir, wavecoeff + aHH + i, wt->wave->lpr(), wt->wave->hpr(), lf, xLp.get(), istride, ostride);

            for (k = lf - 2; k < 2 * ir; ++k) {
                cH[(k - lf + 2) * ic + i] = xLp[k];
            }
        }

        ir *= 2;
        istride = 1;
        ostride = 1;

        for (auto i = 0; i < ir; ++i) {
            idwtSymStride(cL.get() + i * ic, ic, cH.get() + i * ic, wt->wave->lpr(), wt->wave->hpr(), lf, xLp.get(), istride, ostride);

            for (k = lf - 2; k < 2 * ic; ++k) {
                out[(k - lf + 2) + i * ic * 2] = xLp[k];
            }
        }
        ic *= 2;
        if (iter == j - 1) {
            for (auto i = 0; i < wt->rows; ++i) {
                for (k = 0; k < wt->cols; ++k) {
                    oup[k + i * wt->cols] = out[k + i * ic];
                }
            }
        } else {
            for (auto i = 0; i < wt->dimensions[2 * (iter + 1)]; ++i) {
                for (k = 0; k < wt->dimensions[2 * (iter + 1) + 1]; ++k) {
                    oup[k + i * wt->dimensions[2 * (iter + 1) + 1]] = out[k + i * ic];
                }
            }
        }

        orig = oup;
    }
}

auto swt2(WaveletTransform2D* wt, double* inp) -> std::unique_ptr<double[]>
{
    int j = 0;
    int iter = 0;
    int m = 0;
    int n = 0;
    int lp = 0;
    int rowsN = 0;
    int colsN = 0;
    int rowsI = 0;
    int colsI = 0;
    int ir = 0;
    int ic = 0;
    int istride = 0;
    int ostride = 0;
    int aLL = 0;
    int aLH = 0;
    int aHL = 0;
    int aHH = 0;
    int cdim = 0;
    int clen = 0;
    double* orig = nullptr;

    j = wt->J;
    m = 1;
    wt->outlength = 0;

    rowsN = wt->rows;
    colsN = wt->cols;
    lp = wt->wave->lpdLen();
    clen = j * 3;

    auto idx = 2 * j;
    while (idx > 0) {
        wt->dimensions[idx - 1] = colsN;
        wt->dimensions[idx - 2] = rowsN;
        wt->outlength += (rowsN * colsN) * 3;
        idx = idx - 2;
    }
    wt->outlength += (rowsN * colsN);
    n = wt->outlength;
    auto wavecoeff = makeZeros<double>(wt->outlength);

    orig = inp;
    ir = wt->rows;
    ic = wt->cols;
    colsI = wt->dimensions[2 * j - 1];

    auto lpDn1 = std::make_unique<double[]>(ir * colsI);
    auto hpDn1 = std::make_unique<double[]>(ir * colsI);

    for (iter = 0; iter < j; ++iter) {
        if (iter > 0) {
            m = 2 * m;
        }
        rowsI = wt->dimensions[2 * j - 2 * iter - 2];
        colsI = wt->dimensions[2 * j - 2 * iter - 1];
        istride = 1;
        ostride = 1;
        cdim = rowsI * colsI;
        // Row filtering and column subsampling
        for (auto i = 0; i < ir; ++i) {
            swtPerStride(m, orig + i * ic, ic, wt->wave->lpd(), wt->wave->hpd(), lp, lpDn1.get() + i * colsI, colsI, hpDn1.get() + i * colsI, istride, ostride);
        }
        // Column Filtering and Row subsampling
        aHH = n - cdim;
        wt->coeffaccess[clen] = aHH;
        aHL = aHH - cdim;
        wt->coeffaccess[clen - 1] = aHL;
        aLH = aHL - cdim;
        wt->coeffaccess[clen - 2] = aLH;
        aLL = aLH - cdim;

        n -= 3 * cdim;
        ic = colsI;
        istride = ic;
        ostride = ic;
        for (auto i = 0; i < ic; ++i) {
            swtPerStride(m, lpDn1.get() + i, ir, wt->wave->lpd(), wt->wave->hpd(), lp, wavecoeff.get() + aLL + i, rowsI, wavecoeff.get() + aLH + i, istride, ostride);
        }

        for (auto i = 0; i < ic; ++i) {
            swtPerStride(m, hpDn1.get() + i, ir, wt->wave->lpd(), wt->wave->hpd(), lp, wavecoeff.get() + aHL + i, rowsI, wavecoeff.get() + aHH + i, istride, ostride);
        }

        ir = rowsI;
        orig = wavecoeff.get() + aLL;
        clen -= 3;
    }
    wt->coeffaccess[0] = 0;

    return wavecoeff;
}

auto iswt2(WaveletTransform2D* wt, double const* wavecoeffs, double* oup) -> void
{
    int k = 0;
    int iter = 0;
    int it2 = 0;
    int it3 = 0;
    int j = 0;
    int m = 0;
    int rows = 0;
    int cols = 0;
    int lf = 0;
    int ir = 0;
    int ic = 0;
    int k1 = 0;
    int i1 = 0;

    int aLL = 0;
    int aLH = 0;
    int aHL = 0;
    int aHH = 0;
    int shift = 0;
    j = wt->J;
    rows = wt->rows;
    cols = wt->cols;
    lf = wt->wave->lpdLen();

    auto a = makeZeros<double>((rows + lf) * (cols + lf));
    auto h = makeZeros<double>((rows + lf) * (cols + lf));
    auto v = makeZeros<double>((rows + lf) * (cols + lf));
    auto d = makeZeros<double>((rows + lf) * (cols + lf));
    auto oup1 = makeZeros<double>((rows + lf) * (cols + lf));
    auto oup2 = makeZeros<double>((rows + lf) * (cols + lf));

    aLL = wt->coeffaccess[0];

    for (auto i = 0; i < rows; ++i) {
        for (k = 0; k < cols; ++k) {
            oup[i * cols + k] = wavecoeffs[aLL + i * cols + k];
        }
    }

    for (iter = j; iter > 0; iter--) {
        aLH = wt->coeffaccess[(j - iter) * 3 + 1];
        aHL = wt->coeffaccess[(j - iter) * 3 + 2];
        aHH = wt->coeffaccess[(j - iter) * 3 + 3];
        m = (int)std::pow(2.0, (double)iter - 1);

        for (it2 = 0; it2 < m; ++it2) {
            ir = 0;
            ic = 0;
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
            idwtShift(shift, ir, ic, wt->wave->lpr(), wt->wave->hpr(), wt->wave->lpdLen(), a.get(), h.get(), v.get(), d.get(), oup1.get());
            //oup2
            ir = 0;
            ic = 0;
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
            idwtShift(shift, ir, ic, wt->wave->lpr(), wt->wave->hpr(), wt->wave->lpdLen(), a.get(), h.get(), v.get(), d.get(), oup2.get());
            // Shift oup1 and oup2. Then add them to get A.
            i1 = 0;
            for (auto i = it2; i < rows; i += m) {
                k1 = 0;
                for (k = it2; k < cols; k += m) {
                    oup[i * cols + k] = 0.5 * (oup1[i1 * 2 * ic + k1] + oup2[i1 * 2 * ic + k1]);
                    k1++;
                }
                i1++;
            }
        }
    }
}

auto modwt(WaveletTransform2D* wt, double* inp) -> std::unique_ptr<double[]>
{
    int j = 0;
    int iter = 0;
    int m = 0;
    int n = 0;
    int lp = 0;
    int rowsN = 0;
    int colsN = 0;
    int rowsI = 0;
    int colsI = 0;
    int ir = 0;
    int ic = 0;
    int istride = 0;
    int ostride = 0;
    int aLL = 0;
    int aLH = 0;
    int aHL = 0;
    int aHH = 0;
    int cdim = 0;
    int clen = 0;
    double* orig = nullptr;
    double s = NAN;

    j = wt->J;
    m = 1;
    wt->outlength = 0;

    rowsN = wt->rows;
    colsN = wt->cols;
    lp = wt->wave->lpdLen();
    clen = j * 3;

    auto idx = 2 * j;
    while (idx > 0) {
        wt->dimensions[idx - 1] = colsN;
        wt->dimensions[idx - 2] = rowsN;
        wt->outlength += (rowsN * colsN) * 3;
        idx = idx - 2;
    }
    wt->outlength += (rowsN * colsN);
    n = wt->outlength;
    auto wavecoeff = makeZeros<double>(wt->outlength);
    auto filt = std::make_unique<double[]>(2 * lp);
    s = std::sqrt(2.0);
    for (auto i = 0; i < lp; ++i) {
        filt[i] = wt->wave->lpd()[i] / s;
        filt[lp + i] = wt->wave->hpd()[i] / s;
    }

    orig = inp;
    ir = wt->rows;
    ic = wt->cols;
    colsI = wt->dimensions[2 * j - 1];

    auto lpDn1 = std::make_unique<double[]>(ir * colsI);
    auto hpDn1 = std::make_unique<double[]>(ir * colsI);

    for (iter = 0; iter < j; ++iter) {
        if (iter > 0) {
            m = 2 * m;
        }
        rowsI = wt->dimensions[2 * j - 2 * iter - 2];
        colsI = wt->dimensions[2 * j - 2 * iter - 1];
        istride = 1;
        ostride = 1;
        cdim = rowsI * colsI;
        // Row filtering and column subsampling
        for (auto i = 0; i < ir; ++i) {
            modwtPerStride(m, orig + i * ic, ic, filt.get(), lp, lpDn1.get() + i * colsI, colsI, hpDn1.get() + i * colsI, istride, ostride);
        }
        // Column Filtering and Row subsampling
        aHH = n - cdim;
        wt->coeffaccess[clen] = aHH;
        aHL = aHH - cdim;
        wt->coeffaccess[clen - 1] = aHL;
        aLH = aHL - cdim;
        wt->coeffaccess[clen - 2] = aLH;
        aLL = aLH - cdim;
        n -= 3 * cdim;
        ic = colsI;
        istride = ic;
        ostride = ic;
        for (auto i = 0; i < ic; ++i) {
            modwtPerStride(m, lpDn1.get() + i, ir, filt.get(), lp, wavecoeff.get() + aLL + i, rowsI, wavecoeff.get() + aLH + i, istride, ostride);
        }

        for (auto i = 0; i < ic; ++i) {
            modwtPerStride(m, hpDn1.get() + i, ir, filt.get(), lp, wavecoeff.get() + aHL + i, rowsI, wavecoeff.get() + aHH + i, istride, ostride);
        }

        ir = rowsI;
        orig = wavecoeff.get() + aLL;
        clen -= 3;
    }
    wt->coeffaccess[0] = 0;

    return wavecoeff;
}

auto imodwt(WaveletTransform2D* wt, double* wavecoeff, double* oup) -> void
{
    int rows = 0;
    int cols = 0;
    int m = 0;
    // int N;
    int ir = 0;
    int ic = 0;
    int lf = 0;
    int istride = 0;
    int ostride = 0;
    int iter = 0;
    int j = 0;
    int aLL = 0;
    int aLH = 0;
    int aHL = 0;
    int aHH = 0;
    double* orig = nullptr;
    double s = NAN;

    rows = wt->rows;
    cols = wt->cols;
    j = wt->J;

    m = (int)std::pow(2.0, (double)j - 1.0);
    // N = rows > cols ? rows : cols;
    lf = (wt->wave->lprLen() + wt->wave->hprLen()) / 2;

    auto filt = makeZeros<double>(2 * lf);
    s = std::sqrt(2.0);
    for (auto i = 0; i < lf; ++i) {
        filt[i] = wt->wave->lpd()[i] / s;
        filt[lf + i] = wt->wave->hpd()[i] / s;
    }

    auto cL = makeZeros<double>(rows * cols);
    auto cH = makeZeros<double>(rows * cols);
    aLL = wt->coeffaccess[0];
    orig = wavecoeff + aLL;
    for (iter = 0; iter < j; ++iter) {
        if (iter > 0) {
            m = m / 2;
        }
        ir = wt->dimensions[2 * iter];
        ic = wt->dimensions[2 * iter + 1];
        istride = ic;
        ostride = ic;
        aLH = wt->coeffaccess[iter * 3 + 1];
        aHL = wt->coeffaccess[iter * 3 + 2];
        aHH = wt->coeffaccess[iter * 3 + 3];
        for (auto i = 0; i < ic; ++i) {
            imodwtPerStride(m, orig + i, ir, wavecoeff + aLH + i, filt.get(), lf, cL.get() + i, istride, ostride);
            imodwtPerStride(m, wavecoeff + aHL + i, ir, wavecoeff + aHH + i, filt.get(), lf, cH.get() + i, istride, ostride);
        }

        istride = 1;
        ostride = 1;

        for (auto i = 0; i < ir; ++i) {
            imodwtPerStride(m, cL.get() + i * ic, ic, cH.get() + i * ic, filt.get(), lf, oup + i * ic, istride, ostride);
        }

        orig = oup;
    }
}

auto getWT2Coeffs(WaveletTransform2D* wt, double* wcoeffs, int level, char const* type, int* rows, int* cols) -> double*
{
    int j = 0;
    int iter = 0;
    int t = 0;
    double* ptr = nullptr;
    j = wt->J;
    // Error Check

    if (level > j || level < 1) {
        printf("Error : The data is decomposed into %d levels so the acceptable values of level are between 1 and %d", j, j);
        exit(-1);
    }

    if ((strcmp(type, "A") == 0) && level != j) {
        printf("Approximation Coefficients are only available for level %d", j);
        exit(-1);
    }

    if (strcmp(type, "A") == 0) {
        t = 0;
        iter = t;
    } else if (strcmp(type, "H") == 0) {
        t = 1;
        iter = t;
    } else if (strcmp(type, "V") == 0) {
        t = 2;
        iter = t;
    } else if (strcmp(type, "D") == 0) {
        t = 3;
        iter = t;
    } else {
        printf("Only four types of coefficients are accessible A, H, V and D \n");
        exit(-1);
    }

    iter += (j - level) * 3;

    ptr = wcoeffs + wt->coeffaccess[iter];
    *rows = wt->dimensions[2 * (j - level)];
    *cols = wt->dimensions[2 * (j - level) + 1];

    return ptr;
}

auto dispWT2Coeffs(double* a, int row, int col) -> void
{
    printf("\n MATRIX Order : %d X %d \n \n", row, col);

    for (auto i = 0; i < row; i++) {
        printf("R%d: ", i);
        for (auto j = 0; j < col; j++) {
            printf("%g ", a[i * col + j]);
        }
        printf(":R%d \n", i);
    }
}

auto summary(WaveletTransform2D const& wt) -> void
{
    int j = 0;
    int t = 0;
    int rows = 0;
    int cols = 0;
    int vsize = 0;
    j = wt.J;
    summary(*wt.wave);
    printf("\n");
    printf("Wavelet Transform : %s \n", wt.method.c_str());
    printf("\n");
    printf("Signal Extension : %s \n", wt.ext.c_str());
    printf("\n");
    printf("Number of Decomposition Levels %d \n", wt.J);
    printf("\n");
    printf("Input Signal Rows %d \n", wt.rows);
    printf("\n");
    printf("Input Signal Cols %d \n", wt.cols);
    printf("\n");
    printf("Length of Wavelet Coefficients Vector %d \n", wt.outlength);
    printf("\n");
    t = 0;
    for (auto i = j; i > 0; --i) {
        rows = wt.dimensions[2 * (j - i)];
        cols = wt.dimensions[2 * (j - i) + 1];
        vsize = rows * cols;
        printf("Level %d Decomposition Rows :%d Columns:%d Vector Size (Rows*Cols):%d \n", i, rows, cols, vsize);
        printf("Access Row values stored at wt.dimensions[%d]\n", 2 * (j - i));
        printf("Access Column values stored at wt.dimensions[%d]\n\n", 2 * (j - i) + 1);

        if (i == j) {
            printf("Approximation Coefficients access at wt.coeffaccess[%d]=%d, Vector size:%d \n", t, wt.coeffaccess[t], vsize);
        }

        t += 1;
        printf("Horizontal Coefficients access at wt.coeffaccess[%d]=%d, Vector size:%d \n", t, wt.coeffaccess[t], vsize);
        t += 1;
        printf("Vertical Coefficients access at wt.coeffaccess[%d]=%d, Vector size:%d \n", t, wt.coeffaccess[t], vsize);
        t += 1;
        printf("Diagonal Coefficients access at wt.coeffaccess[%d]=%d, Vector size:%d \n\n", t, wt.coeffaccess[t], vsize);
    }
}

auto wt2Free(WaveletTransform2D* wt) -> void
{
    delete wt;
}
