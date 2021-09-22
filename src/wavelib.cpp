/*
  Copyright (c) 2014, Rafat Hussain
*/
#include "wavelib.h"
#include "cwt.h"
#include "wtmath.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string_view>

using namespace std::string_view_literals;

WaveletTransform::WaveletTransform(Wavelet& w, char const* method, int siglength, int j)
    : wave_ { &w }
    , levels_ { j }

{
    auto const size = w.size();
    auto const maxIter = wmaxiter(siglength, w.size());

    if (levels_ > 100) {
        printf("\n The Decomposition Iterations Cannot Exceed 100. Exiting \n");
        exit(-1);
    }

    if (levels_ > maxIter) {
        printf("\n Error - The Signal Can only be iterated %d times using this wavelet. Exiting\n", maxIter);
        exit(-1);
    }

    if (method == nullptr) {
        this->params = std::make_unique<double[]>(siglength + 2 * levels_ * (size + 1));
        this->outlength = siglength + 2 * levels_ * (size + 1);
        ext_ = SignalExtension::symmetric;
    } else if ((method == "dwt"sv) || (method == "DWT"sv)) {
        this->params = std::make_unique<double[]>(siglength + 2 * levels_ * (size + 1));
        this->outlength = siglength + 2 * levels_ * (size + 1);
        ext_ = SignalExtension::symmetric;
    } else if ((method == "swt"sv) || (method == "SWT"sv)) {
        if (testSWTlength(siglength, levels_) == 0) {
            printf("\n For SWT the signal length must be a multiple of 2^levels_. \n");
            exit(-1);
        }

        this->params = std::make_unique<double[]>(siglength * (levels_ + 1));
        this->outlength = siglength * (levels_ + 1);
        ext_ = SignalExtension::periodic;
    } else if ((method == "modwt"sv) || (method == "MODWT"sv)) {

        if (strstr(w.name().c_str(), "haar") == nullptr) {
            if (strstr(w.name().c_str(), "db") == nullptr) {
                if (strstr(w.name().c_str(), "sym") == nullptr) {
                    if (strstr(w.name().c_str(), "coif") == nullptr) {
                        printf("\n MODWT is only implemented for orthogonal wavelet families - db, sym and coif \n");
                        exit(-1);
                    }
                }
            }
        }

        this->params = std::make_unique<double[]>(siglength * 2 * (levels_ + 1));
        this->outlength = siglength * (levels_ + 1);
        ext_ = SignalExtension::periodic;
    }

    this->siglength = siglength;
    this->modwtsiglength = siglength;
    this->MaxIter = maxIter;
    method_ = method;

    if (siglength % 2 == 0) {
        this->even = 1;
    } else {
        this->even = 0;
    }

    this->cobj = nullptr;

    this->cmethod_ = ConvolutionMethod::direct;
    this->cfftset = 0;
    this->lenlength = levels_ + 2;
    this->output_ = &this->params[0];
    if ((method == "dwt"sv) || (method == "DWT"sv)) {
        for (auto i = 0; i < siglength + 2 * levels() * (size + 1); ++i) {
            this->params[i] = 0.0;
        }
    } else if ((method == "swt"sv) || (method == "SWT"sv)) {
        for (auto i = 0; i < siglength * (levels() + 1); ++i) {
            this->params[i] = 0.0;
        }
    } else if ((method == "modwt"sv) || (method == "MODWT"sv)) {
        for (auto i = 0; i < siglength * 2 * (levels() + 1); ++i) {
            this->params[i] = 0.0;
        }
    }
}

auto WaveletTransform::convMethod(ConvolutionMethod method) -> void
{
    cmethod_ = method;
}

auto WaveletTransform::extension(SignalExtension ext) -> void
{
    assert((ext == SignalExtension::periodic) || (ext == SignalExtension::symmetric));
    ext_ = ext;
}

auto WaveletTransform::output() const noexcept -> lt::span<double>
{
    return lt::span<double> { output_, static_cast<std::size_t>(outlength) };
}

auto WaveletTransform::approx() const noexcept -> lt::span<double>
{
    /*
	Wavelet decomposition is stored as
	[A(J) D(J) D(J-1) ..... D(1)] in wt->output vector

	Length of A(J) , N = wt->length[0]
	*/

    return lt::span<double>(output_, static_cast<size_t>(length[0]));
}

auto WaveletTransform::detail(std::size_t level) const noexcept -> lt::span<double>
{
    /*
	returns Detail coefficents at the jth level where j = J,J-1,...,1
	and Wavelet decomposition is stored as
	[A(J) D(J) D(J-1) ..... D(1)] in wt->output vector
	Use getDWTAppx() to get A(J)
	Level 1 : Length of D(J), ie N, is stored in wt->length[1]
	Level 2 :Length of D(J-1), ie N, is stored in wt->length[2]
	....
	Level J : Length of D(1), ie N, is stored in wt->length[J]
	*/

    if (level > static_cast<std::size_t>(levels()) || level < 1U) {
        printf("The decomposition only has 1,..,%d levels", levels());
        exit(-1);
    }

    auto iter = length[0];
    for (auto i = 1U; i < levels() - level; ++i) {
        iter += length[i];
    }

    return lt::span<double>(&output_[iter], static_cast<size_t>(length[level]));
}

auto wtreeInit(Wavelet* wave, int siglength, int j) -> WaveletTree*
{
    auto const size = wave->size();
    auto const maxIter = wmaxiter(siglength, size);
    if (j > 100) {
        printf("\n The Decomposition Iterations Cannot Exceed 100. Exiting \n");
        exit(-1);
    }
    if (j > maxIter) {
        printf("\n Error - The Signal Can only be iterated %d times using this wavelet. Exiting\n", maxIter);
        exit(-1);
    }

    auto temp = 1;
    auto elength = 0;
    auto nodes = 0;
    for (auto i = 0; i < j; ++i) {
        temp *= 2;
        nodes += temp;
        auto const temp2 = (size - 2) * (temp - 1);
        elength += temp2;
    }

    auto obj = std::make_unique<WaveletTree>();
    obj->params = std::make_unique<double[]>(siglength * (j + 1) + elength + nodes + j + 1);
    obj->outlength = siglength * (j + 1) + elength;
    obj->ext = "sym";

    obj->wave = wave;
    obj->siglength = siglength;
    obj->J = j;
    obj->MaxIter = maxIter;
    obj->method = "dwt";

    if (siglength % 2 == 0) {
        obj->even = 1;
    } else {
        obj->even = 0;
    }

    obj->cobj = nullptr;
    obj->nodes = nodes;

    obj->cfftset = 0;
    obj->lenlength = j + 2;
    obj->output = &obj->params[0];
    obj->nodelength = (int*)&obj->params[siglength * (j + 1) + elength];
    obj->coeflength = (int*)&obj->params[siglength * (j + 1) + elength + nodes];

    for (auto i = 0; i < siglength * (j + 1) + elength + nodes + j + 1; ++i) {
        obj->params[i] = 0.0;
    }

    return obj.release();
}

auto wptInit(Wavelet* wave, int siglength, int j) -> WaveletPacketTransform*
{
    auto const size = wave->size();

    if (j > 100) {
        printf("\n The Decomposition Iterations Cannot Exceed 100. Exiting \n");
        exit(-1);
    }

    auto const maxIter = wmaxiter(siglength, size);
    if (j > maxIter) {
        printf("\n Error - The Signal Can only be iterated %d times using this wavelet. Exiting\n", maxIter);
        exit(-1);
    }
    auto temp = 1;
    auto nodes = 0;
    for (auto i = 0; i < j; ++i) {
        temp *= 2;
        nodes += temp;
    }

    auto idx = j;
    auto p2 = 2;
    auto n = siglength;
    auto lp = size;
    auto elength = 0;
    while (idx > 0) {
        n = n + lp - 2;
        n = (int)ceil((double)n / 2.0);
        elength = p2 * n;
        idx--;
        p2 *= 2;
    }

    auto obj = std::make_unique<WaveletPacketTransform>();
    obj->params = std::make_unique<double[]>(elength + 4 * nodes + 2 * j + 6);
    obj->outlength = siglength + 2 * (j + 1) * (size + 1);
    obj->ext = "sym";
    obj->entropy = "shannon";
    obj->eparam = 0.0;

    obj->wave = wave;
    obj->siglength = siglength;
    obj->J = j;
    obj->MaxIter = maxIter;

    if (siglength % 2 == 0) {
        obj->even = 1;
    } else {
        obj->even = 0;
    }

    obj->cobj = nullptr;
    obj->nodes = nodes;

    obj->lenlength = j + 2;
    obj->output = &obj->params[0];
    obj->costvalues = &obj->params[elength];
    obj->basisvector = &obj->params[elength + nodes + 1];
    obj->nodeindex = (int*)&obj->params[elength + 2 * nodes + 2];
    obj->numnodeslevel = (int*)&obj->params[elength + 4 * nodes + 4];
    obj->coeflength = (int*)&obj->params[elength + 4 * nodes + j + 5];

    for (auto i = 0; i < elength + 4 * nodes + 2 * j + 6; ++i) {
        obj->params[i] = 0.0;
    }

    return obj.release();
}

auto cwtInit(char const* wave, double param, int siglength, double dt, int j) -> ComplexWaveletTransform*
{
    int n;
    int nj2;
    int ibase2;
    double s0 {};
    double dj {};
    double t1;
    int m;
    int odd;
    char const* pdefault = "pow";

    m = (int)param;
    odd = 1;
    if (2 * (m / 2) == m) {
        odd = 0;
    }

    n = siglength;
    nj2 = 2 * n * j;
    auto obj = std::make_unique<ComplexWaveletTransform>();
    obj->params = std::make_unique<double[]>(nj2 + 2 * j + n);

    int mother { 0 };
    if ((wave == "morlet"sv) || (wave == "morl"sv)) {
        s0 = 2 * dt;
        dj = 0.4875;
        mother = 0;
        if (param < 0.0) {
            printf("\n Morlet Wavelet Parameter should be >= 0 \n");
            exit(-1);
        }
        if (param == 0) {
            param = 6.0;
        }
        obj->wave = "morlet";

    } else if (wave == "paul"sv) {
        s0 = 2 * dt;
        dj = 0.4875;
        mother = 1;
        if (param < 0 || param > 20) {
            printf("\n Paul Wavelet Parameter should be > 0 and <= 20 \n");
            exit(-1);
        }
        if (param == 0) {
            param = 4.0;
        }
        obj->wave = "paul";

    } else if ((wave == "dgauss"sv) || (wave == "dog"sv)) {
        s0 = 2 * dt;
        dj = 0.4875;
        mother = 2;
        if (param < 0 || odd == 1) {
            printf("\n DOG Wavelet Parameter should be > 0 and even \n");
            exit(-1);
        }
        if (param == 0) {
            param = 2.0;
        }
        obj->wave = "dog";
    }

    obj->pow = 2;
    obj->type = pdefault;

    obj->s0 = s0;
    obj->dj = dj;
    obj->dt = dt;
    obj->J = j;
    obj->siglength = siglength;
    obj->sflag = 0;
    obj->pflag = 1;
    obj->mother = mother;
    obj->m = param;

    t1 = 0.499999 + std::log((double)n) / std::log(2.0);
    ibase2 = 1 + (int)t1;

    obj->npad = (int)std::pow(2.0, (double)ibase2);

    obj->output = (CplxData*)&obj->params[0];
    obj->scale = &obj->params[nj2];
    obj->period = &obj->params[nj2 + j];
    obj->coi = &obj->params[nj2 + 2 * j];

    for (auto i = 0; i < nj2 + 2 * j + n; ++i) {
        obj->params[i] = 0.0;
    }

    return obj.release();
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

static auto wconv(WaveletTransform& wt, double* sig, int n, double const* filt, int l, double* oup) -> void
{
    if (wt.convMethod() == ConvolutionMethod::direct) {
        convDirect(sig, n, filt, l, oup);
    } else if (wt.convMethod() == ConvolutionMethod::fft) {
        if (wt.cfftset == 0) {
            wt.cobj = convInit(n, l);
            convFft(*wt.cobj, sig, filt, oup);
        } else {
            convFft(*wt.cobj, sig, filt, oup);
        }
    } else {
        printf("Convolution Only accepts two methods - direct and fft");
        exit(-1);
    }
}

static auto dwtPer(WaveletTransform& wt, double* inp, int n, double* cA, int lenCA, double* cD) -> void
{

    dwtPerStride(inp, n, wt.wave().lpd(), wt.wave().hpd(), wt.wave().lpdLen(), cA, lenCA, cD, 1, 1);
}

static auto wtreePer(WaveletTree* wt, double const* inp, int n, double* cA, int lenCA, double* cD) -> void
{
    int l;
    int l2;
    int isodd;
    int t;
    int lenAvg;

    lenAvg = wt->wave->lpdLen();
    l2 = lenAvg / 2;
    isodd = n % 2;

    for (auto i = 0; i < lenCA; ++i) {
        t = 2 * i + l2;
        cA[i] = 0.0;
        cD[i] = 0.0;
        for (l = 0; l < lenAvg; ++l) {
            if ((t - l) >= l2 && (t - l) < n) {
                cA[i] += wt->wave->lpd()[l] * inp[t - l];
                cD[i] += wt->wave->hpd()[l] * inp[t - l];
            } else if ((t - l) < l2 && (t - l) >= 0) {
                cA[i] += wt->wave->lpd()[l] * inp[t - l];
                cD[i] += wt->wave->hpd()[l] * inp[t - l];
            } else if ((t - l) < 0 && isodd == 0) {
                cA[i] += wt->wave->lpd()[l] * inp[t - l + n];
                cD[i] += wt->wave->hpd()[l] * inp[t - l + n];
            } else if ((t - l) < 0 && isodd == 1) {
                if ((t - l) != -1) {
                    cA[i] += wt->wave->lpd()[l] * inp[t - l + n + 1];
                    cD[i] += wt->wave->hpd()[l] * inp[t - l + n + 1];
                } else {
                    cA[i] += wt->wave->lpd()[l] * inp[n - 1];
                    cD[i] += wt->wave->hpd()[l] * inp[n - 1];
                }
            } else if ((t - l) >= n && isodd == 0) {
                cA[i] += wt->wave->lpd()[l] * inp[t - l - n];
                cD[i] += wt->wave->hpd()[l] * inp[t - l - n];
            } else if ((t - l) >= n && isodd == 1) {
                if (t - l != n) {
                    cA[i] += wt->wave->lpd()[l] * inp[t - l - (n + 1)];
                    cD[i] += wt->wave->hpd()[l] * inp[t - l - (n + 1)];
                } else {
                    cA[i] += wt->wave->lpd()[l] * inp[n - 1];
                    cD[i] += wt->wave->hpd()[l] * inp[n - 1];
                }
            }
        }
    }
}

static auto dwptPer(WaveletPacketTransform* wt, double const* inp, int n, double* cA, int lenCA, double* cD) -> void
{
    int l;
    int l2;
    int isodd;
    int t;
    int lenAvg;

    lenAvg = wt->wave->lpdLen();
    l2 = lenAvg / 2;
    isodd = n % 2;

    for (auto i = 0; i < lenCA; ++i) {
        t = 2 * i + l2;
        cA[i] = 0.0;
        cD[i] = 0.0;
        for (l = 0; l < lenAvg; ++l) {
            if ((t - l) >= l2 && (t - l) < n) {
                cA[i] += wt->wave->lpd()[l] * inp[t - l];
                cD[i] += wt->wave->hpd()[l] * inp[t - l];
            } else if ((t - l) < l2 && (t - l) >= 0) {
                cA[i] += wt->wave->lpd()[l] * inp[t - l];
                cD[i] += wt->wave->hpd()[l] * inp[t - l];
            } else if ((t - l) < 0 && isodd == 0) {
                cA[i] += wt->wave->lpd()[l] * inp[t - l + n];
                cD[i] += wt->wave->hpd()[l] * inp[t - l + n];
            } else if ((t - l) < 0 && isodd == 1) {
                if ((t - l) != -1) {
                    cA[i] += wt->wave->lpd()[l] * inp[t - l + n + 1];
                    cD[i] += wt->wave->hpd()[l] * inp[t - l + n + 1];
                } else {
                    cA[i] += wt->wave->lpd()[l] * inp[n - 1];
                    cD[i] += wt->wave->hpd()[l] * inp[n - 1];
                }
            } else if ((t - l) >= n && isodd == 0) {
                cA[i] += wt->wave->lpd()[l] * inp[t - l - n];
                cD[i] += wt->wave->hpd()[l] * inp[t - l - n];
            } else if ((t - l) >= n && isodd == 1) {
                if (t - l != n) {
                    cA[i] += wt->wave->lpd()[l] * inp[t - l - (n + 1)];
                    cD[i] += wt->wave->hpd()[l] * inp[t - l - (n + 1)];
                } else {
                    cA[i] += wt->wave->lpd()[l] * inp[n - 1];
                    cD[i] += wt->wave->hpd()[l] * inp[n - 1];
                }
            }
        }
    }
}

static auto dwtSym(WaveletTransform& wt, double* inp, int n, double* cA, int lenCA, double* cD) -> void
{

    dwtSymStride(inp, n, wt.wave().lpd(), wt.wave().hpd(), wt.wave().lpdLen(), cA, lenCA, cD, 1, 1);
}

static auto wtreeSym(WaveletTree* wt, double const* inp, int n, double* cA, int lenCA, double* cD) -> void
{
    int l;
    int t;
    int lenAvg;

    lenAvg = wt->wave->lpdLen();

    for (auto i = 0; i < lenCA; ++i) {
        t = 2 * i + 1;
        cA[i] = 0.0;
        cD[i] = 0.0;
        for (l = 0; l < lenAvg; ++l) {
            if ((t - l) >= 0 && (t - l) < n) {
                cA[i] += wt->wave->lpd()[l] * inp[t - l];
                cD[i] += wt->wave->hpd()[l] * inp[t - l];
            } else if ((t - l) < 0) {
                cA[i] += wt->wave->lpd()[l] * inp[-t + l - 1];
                cD[i] += wt->wave->hpd()[l] * inp[-t + l - 1];
            } else if ((t - l) >= n) {
                cA[i] += wt->wave->lpd()[l] * inp[2 * n - t + l - 1];
                cD[i] += wt->wave->hpd()[l] * inp[2 * n - t + l - 1];
            }
        }
    }
}

static auto dwptSym(WaveletPacketTransform* wt, double const* inp, int n, double* cA, int lenCA, double* cD) -> void
{
    int l;
    int t;
    int lenAvg;

    lenAvg = wt->wave->lpdLen();

    for (auto i = 0; i < lenCA; ++i) {
        t = 2 * i + 1;
        cA[i] = 0.0;
        cD[i] = 0.0;
        for (l = 0; l < lenAvg; ++l) {
            if ((t - l) >= 0 && (t - l) < n) {
                cA[i] += wt->wave->lpd()[l] * inp[t - l];
                cD[i] += wt->wave->hpd()[l] * inp[t - l];
            } else if ((t - l) < 0) {
                cA[i] += wt->wave->lpd()[l] * inp[-t + l - 1];
                cD[i] += wt->wave->hpd()[l] * inp[-t + l - 1];
            } else if ((t - l) >= n) {
                cA[i] += wt->wave->lpd()[l] * inp[2 * n - t + l - 1];
                cD[i] += wt->wave->hpd()[l] * inp[2 * n - t + l - 1];
            }
        }
    }
}

static auto dwt1(WaveletTransform& wt, double* sig, int lenSig, double* cA, double* cD) -> void
{
    constexpr auto d = 2;

    if (wt.extension() == SignalExtension::periodic) {
        auto lenAvg = (wt.wave().lpdLen() + wt.wave().hpdLen()) / 2;
        auto signal = std::make_unique<double[]>(lenSig + lenAvg + (lenSig % 2));
        lenSig = perExt(sig, lenSig, lenAvg / 2, signal.get());
        auto cAUndec = std::make_unique<double[]>(lenSig + lenAvg + wt.wave().lpdLen() - 1);

        if (wt.wave().lpdLen() == wt.wave().hpdLen() && (wt.convMethod() == ConvolutionMethod::fft)) {
            wt.cobj = convInit(lenSig + lenAvg, wt.wave().lpdLen());
            wt.cfftset = 1;
        } else if (!(wt.wave().lpdLen() == wt.wave().hpdLen())) {
            printf("Decomposition Filters must have the same length.");
            exit(-1);
        }

        wconv(wt, signal.get(), lenSig + lenAvg, wt.wave().lpd(), wt.wave().lpdLen(), cAUndec.get());
        downsamp(cAUndec.get() + lenAvg, lenSig, d, cA);
        wconv(wt, signal.get(), lenSig + lenAvg, wt.wave().hpd(), wt.wave().hpdLen(), cAUndec.get());
        downsamp(cAUndec.get() + lenAvg, lenSig, d, cD);

    } else if (wt.extension() == SignalExtension::symmetric) {
        auto lf = wt.wave().lpdLen(); // lpd and hpd have the same length
        auto signal = std::make_unique<double[]>(lenSig + 2 * (lf - 1));
        lenSig = symmExt(sig, lenSig, lf - 1, signal.get());
        auto cAUndec = std::make_unique<double[]>(lenSig + 3 * (lf - 1));

        if (wt.wave().lpdLen() == wt.wave().hpdLen() && (wt.convMethod() == ConvolutionMethod::fft)) {
            wt.cobj = convInit(lenSig + 2 * (lf - 1), lf);
            wt.cfftset = 1;
        } else if (!(wt.wave().lpdLen() == wt.wave().hpdLen())) {
            printf("Decomposition Filters must have the same length.");
            exit(-1);
        }

        wconv(wt, signal.get(), lenSig + 2 * (lf - 1), wt.wave().lpd(), wt.wave().lpdLen(), cAUndec.get());
        downsamp(cAUndec.get() + lf, lenSig + lf - 2, d, cA);
        wconv(wt, signal.get(), lenSig + 2 * (lf - 1), wt.wave().hpd(), wt.wave().hpdLen(), cAUndec.get());
        downsamp(cAUndec.get() + lf, lenSig + lf - 2, d, cD);
    } else {
        printf("Signal extension can be either per or sym");
        exit(-1);
    }

    if (wt.wave().lpdLen() == wt.wave().hpdLen() && (wt.convMethod() == ConvolutionMethod::fft)) {

        wt.cfftset = 0;
    }
}

auto dwt(WaveletTransform& wt, double const* inp) -> void
{

    auto tempLen = wt.siglength;
    auto const j = wt.levels();

    wt.length[j + 1] = tempLen;
    wt.outlength = 0;
    wt.zpad = 0;

    auto orig2 = std::make_unique<double[]>(tempLen);
    auto orig = std::make_unique<double[]>(tempLen);

    for (auto i = 0; i < wt.siglength; ++i) {
        orig[i] = inp[i];
    }

    if (wt.zpad == 1) {
        orig[tempLen - 1] = orig[tempLen - 2];
    }

    auto n = tempLen;
    auto lp = wt.wave().lpdLen();

    if (wt.extension() == SignalExtension::periodic) {
        auto idx = j;
        while (idx > 0) {
            n = (int)ceil((double)n / 2.0);
            wt.length[idx] = n;
            wt.outlength += wt.length[idx];
            idx--;
        }
        wt.length[0] = wt.length[1];
        wt.outlength += wt.length[0];
        n = wt.outlength;

        for (auto iter = 0; iter < j; ++iter) {
            auto const lenCA = wt.length[j - iter];
            n -= lenCA;
            if (wt.convMethod() == ConvolutionMethod::fft) {
                dwt1(wt, orig.get(), tempLen, orig2.get(), wt.params.get() + n);
            } else {
                dwtPer(wt, orig.get(), tempLen, orig2.get(), lenCA, wt.params.get() + n);
            }
            tempLen = wt.length[j - iter];
            if (iter == j - 1) {
                for (auto i = 0; i < lenCA; ++i) {
                    wt.params[i] = orig2[i];
                }
            } else {
                for (auto i = 0; i < lenCA; ++i) {
                    orig[i] = orig2[i];
                }
            }
        }
    } else if (wt.extension() == SignalExtension::symmetric) {
        auto idx = j;
        while (idx > 0) {
            n = n + lp - 2;
            n = (int)ceil((double)n / 2.0);
            wt.length[idx] = n;
            wt.outlength += wt.length[idx];
            idx--;
        }
        wt.length[0] = wt.length[1];
        wt.outlength += wt.length[0];
        n = wt.outlength;

        for (auto iter = 0; iter < j; ++iter) {
            auto const lenCA = wt.length[j - iter];
            n -= lenCA;
            if (wt.convMethod() == ConvolutionMethod::fft) {
                dwt1(wt, orig.get(), tempLen, orig2.get(), wt.params.get() + n);
            } else {
                dwtSym(wt, orig.get(), tempLen, orig2.get(), lenCA, wt.params.get() + n);
            }
            tempLen = wt.length[j - iter];

            if (iter == j - 1) {
                for (auto i = 0; i < lenCA; ++i) {
                    wt.params[i] = orig2[i];
                }
            } else {
                for (auto i = 0; i < lenCA; ++i) {
                    orig[i] = orig2[i];
                }
            }
        }
    } else {
        printf("Signal extension can be either per or sym");
        exit(-1);
    }
}

auto wtree(WaveletTree* wt, double const* inp) -> void
{
    int iter;
    int n;
    int lp;
    int p2;
    int k;
    int n2;
    int np;
    int lenCA;
    int t;
    int t2;
    int it1;

    auto tempLen = wt->siglength;
    auto j = wt->J;
    wt->length[j + 1] = tempLen;
    wt->outlength = 0;
    wt->zpad = 0;

    auto orig = std::make_unique<double[]>(tempLen);

    for (auto i = 0; i < wt->siglength; ++i) {
        orig[i] = inp[i];
    }

    if (wt->zpad == 1) {
        orig[tempLen - 1] = orig[tempLen - 2];
    }

    n = tempLen;
    lp = wt->wave->lpdLen();

    if (wt->ext == "per"sv) {
        auto i = j;
        p2 = 2;
        while (i > 0) {
            n = (int)ceil((double)n / 2.0);
            wt->length[i] = n;
            wt->outlength += p2 * (wt->length[i]);
            i--;
            p2 *= 2;
        }
        wt->length[0] = wt->length[1];

        n2 = wt->outlength;
        p2 = 1;
        for (iter = 0; iter < j; ++iter) {
            lenCA = wt->length[j - iter];
            n2 -= 2 * p2 * lenCA;
            n = n2;
            for (k = 0; k < p2; ++k) {
                if (iter == 0) {
                    wtreePer(wt, orig.get(), tempLen, wt->params.get() + n, lenCA, wt->params.get() + n + lenCA);
                } else {
                    wtreePer(wt, wt->params.get() + np + k * tempLen, tempLen, wt->params.get() + n, lenCA, wt->params.get() + n + lenCA);
                }
                n += 2 * lenCA;
            }

            tempLen = wt->length[j - iter];
            p2 = 2 * p2;
            np = n2;
        }
    } else if (wt->ext == "sym"sv) {
        auto i = j;
        p2 = 2;
        while (i > 0) {
            n = n + lp - 2;
            n = (int)ceil((double)n / 2.0);
            wt->length[i] = n;
            wt->outlength += p2 * (wt->length[i]);
            i--;
            p2 *= 2;
        }
        wt->length[0] = wt->length[1];

        n2 = wt->outlength;
        p2 = 1;

        for (iter = 0; iter < j; ++iter) {
            lenCA = wt->length[j - iter];
            n2 -= 2 * p2 * lenCA;
            n = n2;
            for (k = 0; k < p2; ++k) {
                if (iter == 0) {
                    wtreeSym(wt, orig.get(), tempLen, wt->params.get() + n, lenCA, wt->params.get() + n + lenCA);
                } else {
                    wtreeSym(wt, wt->params.get() + np + k * tempLen, tempLen, wt->params.get() + n, lenCA, wt->params.get() + n + lenCA);
                }
                n += 2 * lenCA;
            }

            tempLen = wt->length[j - iter];
            p2 = 2 * p2;
            np = n2;
        }

    } else {
        printf("Signal extension can be either per or sym");
        exit(-1);
    }

    j = wt->J;
    t2 = wt->outlength - 2 * wt->length[j];
    p2 = 2;
    it1 = 0;
    for (auto i = 0; i < j; ++i) {
        t = t2;
        for (k = 0; k < p2; ++k) {
            wt->nodelength[it1] = t;
            it1++;
            t += wt->length[j - i];
        }
        p2 *= 2;
        t2 = t2 - p2 * wt->length[j - i - 1];
    }

    wt->coeflength[0] = wt->siglength;

    for (auto i = 1; i < j + 1; ++i) {
        wt->coeflength[i] = wt->length[j - i + 1];
    }
}

static constexpr auto ipow2(int n) -> int
{
    auto p = 1;
    for (auto i = 0; i < n; ++i) {
        p *= 2;
    }
    return p;
}

auto dwpt(WaveletPacketTransform* wt, double const* inp) -> void
{
    int iter;
    int p2;
    int k;
    int n2;
    int np;
    int llb;
    double v1;
    double v2;
    int lenCA;
    int t;

    auto tempLen = wt->siglength;
    auto jj = wt->J;
    wt->length[jj + 1] = tempLen;
    wt->outlength = 0;
    auto temp = 1;
    auto elength = 0;
    auto size = wt->wave->size();
    auto nodes = wt->nodes;
    auto n1 = nodes + 1;
    for (auto i = 0; i < jj; ++i) {
        temp *= 2;
        auto const temp2 = (size - 2) * (temp - 1);
        elength += temp2;
    }

    auto eparam = wt->eparam;
    auto orig = std::make_unique<double[]>(tempLen);
    auto tree = std::make_unique<double[]>((tempLen * (jj + 1) + elength));
    auto nodelength = std::make_unique<int[]>(nodes);

    for (auto i = 0; i < wt->siglength; ++i) {
        orig[i] = inp[i];
    }

    for (auto i = 0; i < tempLen * (jj + 1) + elength; ++i) {
        tree[i] = 0.0;
    }

    for (auto i = 0; i < nodes + 1; ++i) {
        wt->basisvector[i] = 0.0;
        wt->costvalues[i] = 0.0;
    }

    auto n = tempLen;
    auto lp = wt->wave->lpdLen();
    // p2 = 1;

    //set eparam value here
    wt->costvalues[0] = costfunc(orig.get(), wt->siglength, wt->entropy.c_str(), eparam);
    auto it2 = 1;
    if (wt->ext == "per"sv) {
        auto i = jj;
        p2 = 2;
        while (i > 0) {
            n = (int)ceil((double)n / 2.0);
            wt->length[i] = n;
            wt->outlength += p2 * (wt->length[i]);
            i--;
            p2 *= 2;
        }
        wt->length[0] = wt->length[1];

        n2 = wt->outlength;
        p2 = 1;
        for (iter = 0; iter < jj; ++iter) {
            lenCA = wt->length[jj - iter];
            n2 -= 2 * p2 * lenCA;
            n = n2;
            for (k = 0; k < p2; ++k) {
                if (iter == 0) {
                    dwptPer(wt, orig.get(), tempLen, tree.get() + n, lenCA, tree.get() + n + lenCA);
                } else {
                    dwptPer(wt, tree.get() + np + k * tempLen, tempLen, tree.get() + n, lenCA, tree.get() + n + lenCA);
                }
                wt->costvalues[it2] = costfunc(tree.get() + n, lenCA, wt->entropy.c_str(), eparam);
                it2++;
                wt->costvalues[it2] = costfunc(tree.get() + n + lenCA, lenCA, wt->entropy.c_str(), eparam);
                it2++;
                n += 2 * lenCA;
            }

            tempLen = wt->length[jj - iter];
            p2 = 2 * p2;
            np = n2;
        }
    } else if (wt->ext == "sym"sv) {
        auto i = jj;
        p2 = 2;
        while (i > 0) {
            n = n + lp - 2;
            n = (int)ceil((double)n / 2.0);
            wt->length[i] = n;
            wt->outlength += p2 * (wt->length[i]);
            i--;
            p2 *= 2;
        }
        wt->length[0] = wt->length[1];

        n2 = wt->outlength;
        p2 = 1;

        for (iter = 0; iter < jj; ++iter) {
            lenCA = wt->length[jj - iter];
            n2 -= 2 * p2 * lenCA;
            n = n2;
            for (k = 0; k < p2; ++k) {
                if (iter == 0) {
                    dwptSym(wt, orig.get(), tempLen, tree.get() + n, lenCA, tree.get() + n + lenCA);
                } else {
                    dwptSym(wt, tree.get() + np + k * tempLen, tempLen, tree.get() + n, lenCA, tree.get() + n + lenCA);
                }
                wt->costvalues[it2] = costfunc(tree.get() + n, lenCA, wt->entropy.c_str(), eparam);
                it2++;
                wt->costvalues[it2] = costfunc(tree.get() + n + lenCA, lenCA, wt->entropy.c_str(), eparam);
                it2++;
                n += 2 * lenCA;
            }

            tempLen = wt->length[jj - iter];
            p2 = 2 * p2;
            np = n2;
        }

    } else {
        printf("Signal extension can be either per or sym");
        exit(-1);
    }

    jj = wt->J;
    auto t2 = wt->outlength - 2 * wt->length[jj];
    p2 = 2;
    auto it1 = 0;
    for (auto i = 0; i < jj; ++i) {
        t = t2;
        for (k = 0; k < p2; ++k) {
            nodelength[it1] = t;
            it1++;
            t += wt->length[jj - i];
        }
        p2 *= 2;
        t2 = t2 - p2 * wt->length[jj - i - 1];
    }

    jj = wt->J;
    llb = 1;
    for (auto i = 0; i < jj; ++i) {
        llb *= 2;
    }

    for (auto i = n1 - llb; i < n1; ++i) {
        wt->basisvector[i] = 1;
    }

    for (auto j = jj - 1; j >= 0; --j) {
        for (k = ipow2(j) - 1; k < ipow2(j + 1) - 1; ++k) {
            v1 = wt->costvalues[k];
            v2 = wt->costvalues[2 * k + 1] + wt->costvalues[2 * k + 2];
            if (v1 <= v2) {
                wt->basisvector[k] = 1;
            } else {
                wt->costvalues[k] = v2;
            }
        }
    }

    for (k = 0; k < nodes / 2; ++k) {
        if (wt->basisvector[k] == 1 || wt->basisvector[k] == 2) {
            wt->basisvector[2 * k + 1] = 2;
            wt->basisvector[2 * k + 2] = 2;
        }
    }

    for (k = 0; k < n1; ++k) {
        if (wt->basisvector[k] == 2) {
            wt->basisvector[k] = 0;
        }
    }

    // N2 = 0;
    it1 = n1;
    it2 = 0;
    wt->nodes = 0;
    wt->numnodeslevel[0] = 0;

    if (wt->basisvector[0] == 1) {
        wt->outlength = wt->siglength;
        for (auto i = 0; i < wt->siglength; ++i) {
            wt->output[i] = inp[i];
        }
        wt->nodes = 1;
        wt->nodeindex[0] = 0;
        wt->nodeindex[1] = 0;
        wt->numnodeslevel[0] = 1;
    } else {
        for (auto i = jj; i > 0; --i) {
            llb = ipow2(i);
            it1 -= llb;
            wt->numnodeslevel[i] = 0;
            for (auto j = 0; j < llb; ++j) {
                if (wt->basisvector[it1 + j] == 1) {
                    wt->nodeindex[2 * wt->nodes] = i;
                    wt->nodeindex[2 * wt->nodes + 1] = j;
                    wt->nodes += 1;
                    wt->numnodeslevel[i] += 1;
                    for (k = 0; k < wt->length[jj - i + 1]; ++k) {
                        wt->output[it2 + k] = tree[nodelength[it1 - 1 + j] + k]; // access tree
                    }
                    it2 += wt->length[jj - i + 1];
                }
            }
        }
        wt->outlength = it2;
    }

    wt->coeflength[0] = wt->siglength;

    for (auto i = 1; i < jj + 1; ++i) {
        wt->coeflength[i] = wt->length[jj - i + 1];
    }
}

/// X - Level. All Nodes at any level have the same length
auto getWTREENodelength(WaveletTree* wt, int x) -> int
{
    if (x <= 0 || x > wt->J) {
        printf("X co-ordinate must be >= 1 and <= %d", wt->J);
        exit(-1);
    }

    return wt->length[wt->J - x + 1];
}

/// X - Level. All Nodes at any level have the same length
auto getDWPTNodelength(WaveletPacketTransform* wt, int x) -> int
{
    if (x <= 0 || x > wt->J) {
        printf("X co-ordinate must be >= 1 and <= %d", wt->J);
        exit(-1);
    }

    return wt->length[wt->J - x + 1];
}

auto getWTREECoeffs(WaveletTree* wt, int x, int y, double* coeffs, int n) -> void
{
    int ymax;
    int t;
    int t2;

    if (x <= 0 || x > wt->J) {
        printf("X co-ordinate must be >= 1 and <= %d", wt->J);
        exit(-1);
    }
    ymax = 1;
    for (auto i = 0; i < x; ++i) {
        ymax *= 2;
    }

    ymax -= 1;

    if (y < 0 || y > ymax) {
        printf("Y co-ordinate must be >= 0 and <= %d", ymax);
        exit(-1);
    }

    if (x == 1) {
        t = 0;
    } else {
        t = 0;
        t2 = 1;
        for (auto i = 0; i < x - 1; ++i) {
            t2 *= 2;
            t += t2;
        }
    }

    t += y;
    t2 = wt->nodelength[t];
    for (auto i = 0; i < n; ++i) {
        coeffs[i] = wt->output[t2 + i];
    }
}

auto setCWTScales(ComplexWaveletTransform* wt, double s0, double dj, char const* type, int power) -> void
{
    wt->type = type;
    //s0*std::pow(2.0, (double)(j - 1)*dj);
    if ((wt->type == "pow"sv) || (wt->type == "power"sv)) {
        for (auto i = 0; i < wt->J; ++i) {
            wt->scale[i] = s0 * std::pow((double)power, (double)(i)*dj);
        }
        wt->sflag = 1;
        wt->pow = power;

    } else if ((wt->type == "lin"sv) || (wt->type == "linear"sv)) {
        for (auto i = 0; i < wt->J; ++i) {
            wt->scale[i] = s0 + (double)i * dj;
        }
        wt->sflag = 1;
    } else {
        printf("\n Type accepts only two values : pow and lin\n");
        exit(-1);
    }
    wt->s0 = s0;
    wt->dj = dj;
}

auto cwt(ComplexWaveletTransform* wt, double const* inp) -> void
{
    int n;
    int npad;
    int nj2;
    int j;
    int j2;
    n = wt->siglength;
    if (wt->sflag == 0) {
        for (auto i = 0; i < wt->J; ++i) {
            wt->scale[i] = wt->s0 * std::pow(2.0, (double)(i)*wt->dj);
        }
        wt->sflag = 1;
    }

    if (wt->pflag == 0) {
        npad = n;
    } else {
        npad = wt->npad;
    }

    nj2 = 2 * n * wt->J;
    j = wt->J;
    j2 = 2 * j;

    wt->smean = 0.0;

    for (auto i = 0; i < n; ++i) {
        wt->smean += inp[i];
    }
    wt->smean /= n;

    cwavelet(inp, n, wt->dt, wt->mother, wt->m, wt->s0, wt->dj, wt->J, npad, wt->params.get(), wt->params.get() + nj2, wt->params.get() + nj2 + j, wt->params.get() + nj2 + j2);
}

auto icwt(ComplexWaveletTransform* wt, double* cwtop) -> void
{
    double psi;
    double cdel;
    int real;
    int n;
    int nj2;

    n = wt->siglength;
    nj2 = n * 2 * wt->J;

    psi0(wt->mother, wt->m, &psi, &real);
    cdel = cdelta(wt->mother, wt->m, psi);

    if (((wt->type == "pow"sv) || (wt->type == "power"sv)) && wt->pow == 2) {
        icwavelet(wt->params.get(), n, wt->params.get() + nj2, wt->J, wt->dt, wt->dj, cdel, psi, cwtop);
    } else {
        printf("Inverse CWT is only available for power of 2.0 scales \n");
        exit(-1);
    }
    for (auto i = 0; i < n; ++i) {
        cwtop[i] += wt->smean;
    }
}

static auto idwt1(WaveletTransform& wt, double* temp, double* cAUp, double* cA, int lenCA, double* cD, int lenCD, double* xLp, double* xHp, double* x) -> void
{
    auto lenAvg = (wt.wave().lprLen() + wt.wave().hprLen()) / 2;
    auto n = 2 * lenCD;
    auto u = 2;

    upsamp2(cA, lenCA, u, cAUp);

    perExt(cAUp, 2 * lenCA, lenAvg / 2, temp);

    auto n2 = 2 * lenCA + lenAvg;

    if (wt.wave().lprLen() == wt.wave().hprLen() && (wt.convMethod() == ConvolutionMethod::fft)) {
        wt.cobj = convInit(n2, lenAvg);
        wt.cfftset = 1;
    } else if (!(wt.wave().lprLen() == wt.wave().hprLen())) {
        printf("Decomposition Filters must have the same length.");
        exit(-1);
    }

    wconv(wt, temp, n2, wt.wave().lpr(), lenAvg, xLp);

    upsamp2(cD, lenCD, u, cAUp);

    perExt(cAUp, 2 * lenCD, lenAvg / 2, temp);

    n2 = 2 * lenCD + lenAvg;

    wconv(wt, temp, n2, wt.wave().hpr(), lenAvg, xHp);

    for (auto i = lenAvg - 1; i < n + lenAvg - 1; ++i) {
        x[i - lenAvg + 1] = xLp[i] + xHp[i];
    }

    if (wt.wave().lprLen() == wt.wave().hprLen() && (wt.convMethod() == ConvolutionMethod::fft)) {

        wt.cfftset = 0;
    }
}

static auto idwtPer(WaveletTransform& wt, double* cA, int lenCA, double* cD, double* x) -> void
{
    idwtPerStride(cA, lenCA, cD, wt.wave().lpr(), wt.wave().hpr(), wt.wave().lprLen(), x, 1, 1);
}

static auto idwtSym(WaveletTransform& wt, double* cA, int lenCA, double* cD, double* x) -> void
{
    idwtSymStride(cA, lenCA, cD, wt.wave().lpr(), wt.wave().hpr(), wt.wave().lprLen(), x, 1, 1);
}

auto idwt(WaveletTransform& wt, double* dwtop) -> void
{

    int lf;
    int n;
    int n2;
    int iter;
    int k;
    int detLen;

    auto j = wt.levels();
    auto u = 2;
    auto appLen = wt.length[0];
    auto out = std::make_unique<double[]>(wt.siglength + 1);
    if ((wt.extension() == SignalExtension::periodic) && (wt.convMethod() == ConvolutionMethod::fft)) {
        appLen = wt.length[0];
        detLen = wt.length[1];
        n = 2 * wt.length[j];
        lf = (wt.wave().lprLen() + wt.wave().hprLen()) / 2;

        auto cAUp = std::make_unique<double[]>(n);
        auto temp = std::make_unique<double[]>((n + lf));
        auto xLp = std::make_unique<double[]>((n + 2 * lf - 1));
        auto xHp = std::make_unique<double[]>((n + 2 * lf - 1));
        iter = appLen;

        for (auto i = 0; i < appLen; ++i) {
            out[i] = wt.output()[i];
        }

        for (auto i = 0; i < j; ++i) {

            idwt1(wt, temp.get(), cAUp.get(), out.get(), detLen, wt.output().data() + iter, detLen, xLp.get(), xHp.get(), out.get());
            /*
			idwt_per(wt,out.get(), det_len, wt.output().data() + iter, det_len, X_lp);
			for (k = lf/2 - 1; k < 2 * det_len + lf/2 - 1; ++k) {
				out[k - lf/2 + 1] = X_lp[k];
			}
			*/
            iter += detLen;
            detLen = wt.length[i + 2];
        }

    } else if ((wt.extension() == SignalExtension::periodic) && (wt.convMethod() == ConvolutionMethod::direct)) {
        appLen = wt.length[0];
        detLen = wt.length[1];
        n = 2 * wt.length[j];
        lf = (wt.wave().lprLen() + wt.wave().hprLen()) / 2;

        auto xLp = std::make_unique<double[]>((n + 2 * lf - 1));
        iter = appLen;

        for (auto i = 0; i < appLen; ++i) {
            out[i] = wt.output()[i];
        }

        for (auto i = 0; i < j; ++i) {
            idwtPer(wt, out.get(), detLen, wt.output().data() + iter, xLp.get());
            for (k = lf / 2 - 1; k < 2 * detLen + lf / 2 - 1; ++k) {
                out[k - lf / 2 + 1] = xLp[k];
            }

            iter += detLen;
            detLen = wt.length[i + 2];
        }

    } else if ((wt.extension() == SignalExtension::symmetric) && (wt.convMethod() == ConvolutionMethod::direct)) {
        appLen = wt.length[0];
        detLen = wt.length[1];
        n = 2 * wt.length[j] - 1;
        lf = (wt.wave().lprLen() + wt.wave().hprLen()) / 2;

        auto xLp = std::make_unique<double[]>((n + 2 * lf - 1));
        iter = appLen;

        for (auto i = 0; i < appLen; ++i) {
            out[i] = wt.output()[i];
        }

        for (auto i = 0; i < j; ++i) {
            idwtSym(wt, out.get(), detLen, wt.output().data() + iter, xLp.get());
            for (k = lf - 2; k < 2 * detLen; ++k) {
                out[k - lf + 2] = xLp[k];
            }

            iter += detLen;
            detLen = wt.length[i + 2];
        }

    } else if ((wt.extension() == SignalExtension::symmetric) && (wt.convMethod() == ConvolutionMethod::fft)) {
        lf = wt.wave().lpdLen(); // lpd and hpd have the same length

        n = 2 * wt.length[j] - 1;
        auto cAUp = std::make_unique<double[]>(n);
        auto xLp = std::make_unique<double[]>((n + lf - 1));
        auto xHp = std::make_unique<double[]>((n + lf - 1));

        for (auto i = 0; i < appLen; ++i) {
            out[i] = wt.output()[i];
        }

        iter = appLen;

        for (auto i = 0; i < j; ++i) {
            detLen = wt.length[i + 1];
            upsamp(out.get(), detLen, u, cAUp.get());
            n2 = 2 * wt.length[i + 1] - 1;

            if (wt.wave().lprLen() == wt.wave().hprLen() && (wt.convMethod() == ConvolutionMethod::fft)) {
                wt.cobj = convInit(n2, lf);
                wt.cfftset = 1;
            } else if (!(wt.wave().lprLen() == wt.wave().hprLen())) {
                printf("Decomposition Filters must have the same length.");
                exit(-1);
            }

            wconv(wt, cAUp.get(), n2, wt.wave().lpr(), lf, xLp.get());
            upsamp(wt.output().data() + iter, detLen, u, cAUp.get());
            wconv(wt, cAUp.get(), n2, wt.wave().hpr(), lf, xHp.get());

            for (k = lf - 2; k < n2 + 1; ++k) {
                out[k - lf + 2] = xLp[k] + xHp[k];
            }
            iter += detLen;
            if (wt.wave().lprLen() == wt.wave().hprLen() && (wt.convMethod() == ConvolutionMethod::fft)) {

                wt.cfftset = 0;
            }
        }

    } else {
        printf("Signal extension can be either per or sym");
        exit(-1);
    }

    for (auto i = 0; i < wt.siglength; ++i) {
        dwtop[i] = out[i];
    }
}

static auto idwptPer(WaveletPacketTransform* wt, double const* cA, int lenCA, double const* cD, double* x) -> void
{
    int lenAvg;
    int l;
    int m;
    int n;
    int t;
    int l2;

    lenAvg = (wt->wave->lprLen() + wt->wave->hprLen()) / 2;
    l2 = lenAvg / 2;
    m = -2;
    n = -1;

    for (auto i = 0; i < lenCA + l2 - 1; ++i) {
        m += 2;
        n += 2;
        x[m] = 0.0;
        x[n] = 0.0;
        for (l = 0; l < l2; ++l) {
            t = 2 * l;
            if ((i - l) >= 0 && (i - l) < lenCA) {
                x[m] += wt->wave->lpr()[t] * cA[i - l] + wt->wave->hpr()[t] * cD[i - l];
                x[n] += wt->wave->lpr()[t + 1] * cA[i - l] + wt->wave->hpr()[t + 1] * cD[i - l];
            } else if ((i - l) >= lenCA && (i - l) < lenCA + lenAvg - 1) {
                x[m] += wt->wave->lpr()[t] * cA[i - l - lenCA] + wt->wave->hpr()[t] * cD[i - l - lenCA];
                x[n] += wt->wave->lpr()[t + 1] * cA[i - l - lenCA] + wt->wave->hpr()[t + 1] * cD[i - l - lenCA];
            } else if ((i - l) < 0 && (i - l) > -l2) {
                x[m] += wt->wave->lpr()[t] * cA[lenCA + i - l] + wt->wave->hpr()[t] * cD[lenCA + i - l];
                x[n] += wt->wave->lpr()[t + 1] * cA[lenCA + i - l] + wt->wave->hpr()[t + 1] * cD[lenCA + i - l];
            }
        }
    }
}

static auto idwptSym(WaveletPacketTransform* wt, double const* cA, int lenCA, double const* cD, double* x) -> void
{
    auto lenAvg = (wt->wave->lprLen() + wt->wave->hprLen()) / 2;
    auto m = -2;
    auto n = -1;

    for (auto v = 0; v < lenCA; ++v) {
        auto i = v;
        m += 2;
        n += 2;
        x[m] = 0.0;
        x[n] = 0.0;
        for (auto l = 0; l < lenAvg / 2; ++l) {
            auto const t = 2 * l;
            if ((i - l) >= 0 && (i - l) < lenCA) {
                x[m] += wt->wave->lpr()[t] * cA[i - l] + wt->wave->hpr()[t] * cD[i - l];
                x[n] += wt->wave->lpr()[t + 1] * cA[i - l] + wt->wave->hpr()[t + 1] * cD[i - l];
            }
        }
    }
}

auto idwpt(WaveletPacketTransform* wt, double* dwtop) -> void
{
    int k;
    int l;
    int index;

    auto j = wt->J;
    auto appLen = wt->length[0];
    auto powJ = ipow2(j);
    auto lf = (wt->wave->lprLen() + wt->wave->hprLen()) / 2;
    auto xlen = powJ * (appLen + 2 * lf);

    auto xLp = std::make_unique<double[]>(2 * (wt->length[j] + lf));
    auto x = std::make_unique<double[]>(xlen);
    auto out = std::make_unique<double[]>(wt->length[j]);
    auto out2 = std::make_unique<double[]>(wt->length[j]);
    auto prep = makeZeros<int>(powJ);
    auto ptemp = makeZeros<int>(powJ);
    auto n1 = 1;
    auto llb = 1;
    auto index2 = xlen / powJ;
    auto indexp = 0;
    if (wt->basisvector[0] == 1) {
        for (auto i = 0; i < wt->siglength; ++i) {
            dwtop[i] = wt->output[i];
        }

    } else {
        for (auto i = 0; i < j; ++i) {
            llb *= 2;
            n1 += llb;
        }

        for (auto i = 0; i < xlen; ++i) {
            x[i] = 0.0;
        }

        for (auto i = 0; i < llb; ++i) {
            prep[i] = (int)wt->basisvector[n1 - llb + i];
            ptemp[i] = 0;
        }

        if (wt->ext == "per"sv) {
            index = 0;
            for (auto i = 0; i < j; ++i) {
                auto p = ipow2(j - i - 1);
                auto detLen = wt->length[i + 1];
                index2 *= 2;
                auto index3 = 0;
                auto index4 = 0;
                n1 -= llb;
                for (l = 0; l < llb; ++l) {
                    if (ptemp[l] != 2) {
                        prep[l] = (int)wt->basisvector[n1 + l];
                    } else {
                        prep[l] = ptemp[l];
                    }
                    ptemp[l] = 0;
                }

                for (l = 0; l < p; ++l) {
                    if (prep[2 * l] == 1 && prep[2 * l + 1] == 1) {
                        for (k = 0; k < detLen; ++k) {
                            out[k] = wt->output[index + k];
                            out2[k] = wt->output[index + detLen + k];
                        }
                        idwptPer(wt, out.get(), detLen, out2.get(), xLp.get());
                        for (k = lf / 2 - 1; k < 2 * detLen + lf / 2 - 1; ++k) {
                            x[index3 + k - lf / 2 + 1] = xLp[k];
                        }
                        index += 2 * detLen;
                        index3 += index2;
                        index4 += 2 * indexp;
                        ptemp[l] = 2;
                    } else if (prep[2 * l] == 1 && prep[2 * l + 1] == 2) {
                        index4 += indexp;
                        for (k = 0; k < detLen; ++k) {
                            out[k] = wt->output[index + k];
                            out2[k] = x[index4 + k];
                        }
                        idwptPer(wt, out.get(), detLen, out2.get(), xLp.get());
                        for (k = lf / 2 - 1; k < 2 * detLen + lf / 2 - 1; ++k) {
                            x[index3 + k - lf / 2 + 1] = xLp[k];
                        }
                        index += detLen;
                        index3 += index2;
                        index4 += indexp;
                        ptemp[l] = 2;
                    } else if (prep[2 * l] == 2 && prep[2 * l + 1] == 1) {
                        for (k = 0; k < detLen; ++k) {
                            out[k] = x[index4 + k];
                            out2[k] = wt->output[index + k];
                        }
                        idwptPer(wt, out.get(), detLen, out2.get(), xLp.get());
                        for (k = lf / 2 - 1; k < 2 * detLen + lf / 2 - 1; ++k) {
                            x[index3 + k - lf / 2 + 1] = xLp[k];
                        }
                        index += detLen;
                        index3 += index2;
                        index4 += 2 * indexp;
                        ptemp[l] = 2;
                    } else if (prep[2 * l] == 2 && prep[2 * l + 1] == 2) {
                        for (k = 0; k < detLen; ++k) {
                            out[k] = x[index4 + k];
                            out2[k] = x[index4 + indexp + k];
                        }
                        idwptPer(wt, out.get(), detLen, out2.get(), xLp.get());
                        for (k = lf / 2 - 1; k < 2 * detLen + lf / 2 - 1; ++k) {
                            x[index3 + k - lf / 2 + 1] = xLp[k];
                        }
                        index4 += 2 * indexp;
                        index3 += index2;
                        ptemp[l] = 2;
                    } else {
                        index3 += index2;
                        index4 += 2 * indexp;
                    }
                }

                llb /= 2;
                indexp = index2;
            }

        } else if (wt->ext == "sym"sv) {
            index = 0;

            for (auto i = 0; i < j; ++i) {
                auto p = ipow2(j - i - 1);
                auto detLen = wt->length[i + 1];
                index2 *= 2;
                auto index3 = 0;
                auto index4 = 0;
                n1 -= llb;
                for (l = 0; l < llb; ++l) {
                    if (ptemp[l] != 2) {
                        prep[l] = (int)wt->basisvector[n1 + l];
                    } else {
                        prep[l] = ptemp[l];
                    }
                    ptemp[l] = 0;
                }

                for (l = 0; l < p; ++l) {
                    if (prep[2 * l] == 1 && prep[2 * l + 1] == 1) {
                        for (k = 0; k < detLen; ++k) {
                            out[k] = wt->output[index + k];
                            out2[k] = wt->output[index + detLen + k];
                        }
                        idwptSym(wt, out.get(), detLen, out2.get(), xLp.get());
                        for (k = lf - 2; k < 2 * detLen; ++k) {
                            x[index3 + k - lf + 2] = xLp[k];
                        }
                        index += 2 * detLen;
                        index3 += index2;
                        index4 += 2 * indexp;
                        ptemp[l] = 2;
                    } else if (prep[2 * l] == 1 && prep[2 * l + 1] == 2) {
                        index4 += indexp;
                        for (k = 0; k < detLen; ++k) {
                            out[k] = wt->output[index + k];
                            out2[k] = x[index4 + k];
                        }
                        idwptSym(wt, out.get(), detLen, out2.get(), xLp.get());
                        for (k = lf - 2; k < 2 * detLen; ++k) {
                            x[index3 + k - lf + 2] = xLp[k];
                        }
                        index += detLen;
                        index3 += index2;
                        index4 += indexp;
                        ptemp[l] = 2;
                    } else if (prep[2 * l] == 2 && prep[2 * l + 1] == 1) {
                        for (k = 0; k < detLen; ++k) {
                            out[k] = x[index4 + k];
                            out2[k] = wt->output[index + k];
                        }
                        idwptSym(wt, out.get(), detLen, out2.get(), xLp.get());
                        for (k = lf - 2; k < 2 * detLen; ++k) {
                            x[index3 + k - lf + 2] = xLp[k];
                        }
                        index += detLen;
                        index3 += index2;
                        index4 += 2 * indexp;
                        ptemp[l] = 2;
                    } else if (prep[2 * l] == 2 && prep[2 * l + 1] == 2) {
                        for (k = 0; k < detLen; ++k) {
                            out[k] = x[index4 + k];
                            out2[k] = x[index4 + indexp + k];
                        }
                        idwptSym(wt, out.get(), detLen, out2.get(), xLp.get());
                        for (k = lf - 2; k < 2 * detLen; ++k) {
                            x[index3 + k - lf + 2] = xLp[k];
                        }
                        index4 += 2 * indexp;
                        index3 += index2;
                        ptemp[l] = 2;
                    } else {
                        index3 += index2;
                        index4 += 2 * indexp;
                    }
                }

                //idwt1(wt, temp, cA_up, out, det_len, wt->output().data() + iter, det_len, X_lp.get(), X_hp, out);
                /*
				idwpt_sym(wt, out, det_len, wt->output().data() + iter, det_len, X_lp);
				for (k = lf - 2; k < 2 * det_len; ++k) {
				out[k - lf + 2] = X_lp[k];
				}

				iter += det_len;
				det_len = wt->length[i + 2];
				*/
                llb /= 2;
                indexp = index2;
            }

            //free(X_lp);

        } else {
            printf("Signal extension can be either per or sym");
            exit(-1);
        }

        for (auto i = 0; i < wt->siglength; ++i) {
            dwtop[i] = x[i];
        }
    }
}

static auto swtPer(WaveletTransform& wt, int m, double* inp, int n, double* cA, int lenCA, double* cD) -> void
{

    swtPerStride(m, inp, n, wt.wave().lpd(), wt.wave().hpd(), wt.wave().lpdLen(), cA, lenCA, cD, 1, 1);
}

static auto swtFft(WaveletTransform& wt, double const* inp) -> void
{
    int n { 0 };

    auto tempLen = wt.siglength;
    auto j = wt.levels();
    wt.length[0] = wt.length[j] = tempLen;
    wt.outlength = wt.length[j + 1] = (j + 1) * tempLen;
    auto m = 1;
    for (auto iter = 1; iter < j; ++iter) {
        m = 2 * m;
        wt.length[iter] = tempLen;
    }

    auto const lenFilt = wt.wave().size();

    auto lowPass = std::make_unique<double[]>(m * lenFilt);
    auto highPass = std::make_unique<double[]>(m * lenFilt);
    auto sig = std::make_unique<double[]>((m * lenFilt + tempLen + (tempLen % 2)));
    auto cA = std::make_unique<double[]>((2 * m * lenFilt + tempLen + (tempLen % 2)) - 1);
    auto cD = std::make_unique<double[]>((2 * m * lenFilt + tempLen + (tempLen % 2)) - 1);

    m = 1;

    for (auto i = 0; i < tempLen; ++i) {
        wt.params[i] = inp[i];
    }

    auto lenacc = wt.outlength;

    for (auto iter = 0; iter < j; ++iter) {
        lenacc -= tempLen;
        if (iter > 0) {
            m = 2 * m;
            n = m * lenFilt;
            upsamp2(wt.wave().lpd(), wt.wave().lpdLen(), m, lowPass.get());
            upsamp2(wt.wave().hpd(), wt.wave().hpdLen(), m, highPass.get());
        } else {
            n = lenFilt;
            for (auto i = 0; i < n; ++i) {
                lowPass[i] = wt.wave().lpd()[i];
                highPass[i] = wt.wave().hpd()[i];
            }
        }

        //swt_per(wt,M, wt.params.get(), temp_len, cA, temp_len, cD,temp_len);

        perExt(wt.params.get(), tempLen, n / 2, sig.get());

        if (wt.wave().lpdLen() == wt.wave().hpdLen() && (wt.convMethod() == ConvolutionMethod::fft)) {
            wt.cobj = convInit(n + tempLen + (tempLen % 2), n);
            wt.cfftset = 1;
        } else if (!(wt.wave().lpdLen() == wt.wave().hpdLen())) {
            printf("Decomposition Filters must have the same length.");
            exit(-1);
        }

        wconv(wt, sig.get(), n + tempLen + (tempLen % 2), lowPass.get(), n, cA.get());

        wconv(wt, sig.get(), n + tempLen + (tempLen % 2), highPass.get(), n, cD.get());

        if (wt.wave().lpdLen() == wt.wave().hpdLen() && (wt.convMethod() == ConvolutionMethod::fft)) {

            wt.cfftset = 0;
        }

        for (auto i = 0; i < tempLen; ++i) {
            wt.params[i] = cA[n + i];
            wt.params[lenacc + i] = cD[n + i];
        }
    }
}

static auto swtDirect(WaveletTransform& wt, double const* inp) -> void
{
    int j;
    int tempLen;
    int iter;
    int m;
    int lenacc;

    tempLen = wt.siglength;
    j = wt.levels();
    wt.length[0] = wt.length[j] = tempLen;
    wt.outlength = wt.length[j + 1] = (j + 1) * tempLen;
    m = 1;
    for (iter = 1; iter < j; ++iter) {
        m = 2 * m;
        wt.length[iter] = tempLen;
    }

    auto cA = std::make_unique<double[]>(tempLen);
    auto cD = std::make_unique<double[]>(tempLen);

    m = 1;

    for (auto i = 0; i < tempLen; ++i) {
        wt.params[i] = inp[i];
    }

    lenacc = wt.outlength;

    for (iter = 0; iter < j; ++iter) {
        lenacc -= tempLen;
        if (iter > 0) {
            m = 2 * m;
        }

        swtPer(wt, m, wt.params.get(), tempLen, cA.get(), tempLen, cD.get());

        for (auto i = 0; i < tempLen; ++i) {
            wt.params[i] = cA[i];
            wt.params[lenacc + i] = cD[i];
        }
    }
}

auto swt(WaveletTransform& wt, double const* inp) -> void
{
    if ((wt.method() == "swt"sv) && (wt.convMethod() == ConvolutionMethod::direct)) {
        swtDirect(wt, inp);
    } else if ((wt.method() == "swt"sv) && (wt.convMethod() == ConvolutionMethod::fft)) {
        swtFft(wt, inp);
    } else {
        printf("SWT Only accepts two methods - direct and fft");
        exit(-1);
    }
}

auto iswt(WaveletTransform& wt, double* swtop) -> void
{
    int n;
    int lf;
    int iter;
    int j;
    int index;
    int value;
    int count;
    int len;
    int indexShift;
    int len0;
    int u;
    int n1;
    int index2;

    n = wt.siglength;
    j = wt.levels();
    u = 2;
    lf = wt.wave().lprLen();

    auto appxSig = std::make_unique<double[]>(n);
    auto detSig = std::make_unique<double[]>(n);
    auto appx1 = std::make_unique<double[]>(n);
    auto det1 = std::make_unique<double[]>(n);
    auto appx2 = std::make_unique<double[]>(n);
    auto det2 = std::make_unique<double[]>(n);
    auto tempx = std::make_unique<double[]>(n);
    auto cL0 = std::make_unique<double[]>((n + (n % 2) + lf));
    auto cH0 = std::make_unique<double[]>((n + (n % 2) + lf));
    auto oup00L = std::make_unique<double[]>((n + 2 * lf));
    auto oup00H = std::make_unique<double[]>((n + 2 * lf));
    auto oup00 = std::make_unique<double[]>(n);
    auto oup01 = std::make_unique<double[]>(n);

    for (iter = 0; iter < j; ++iter) {
        for (auto i = 0; i < n; ++i) {
            swtop[i] = 0.0;
        }
        if (iter == 0) {
            for (auto i = 0; i < n; ++i) {
                appxSig[i] = wt.output()[i];
                detSig[i] = wt.output()[n + i];
            }
        } else {
            for (auto i = 0; i < n; ++i) {
                detSig[i] = wt.output()[(iter + 1) * n + i];
            }
        }

        value = (int)std::pow(2.0, (double)(j - 1 - iter));

        for (count = 0; count < value; count++) {
            len = 0;
            for (index = count; index < n; index += value) {
                appx1[len] = appxSig[index];
                det1[len] = detSig[index];
                len++;
            }

            //SHIFT 0
            len0 = 0;

            for (indexShift = 0; indexShift < len; indexShift += 2) {
                appx2[len0] = appx1[indexShift];
                det2[len0] = det1[indexShift];
                len0++;
            }
            upsamp2(appx2.get(), len0, u, tempx.get());
            perExt(tempx.get(), 2 * len0, lf / 2, cL0.get());

            upsamp2(det2.get(), len0, u, tempx.get());
            perExt(tempx.get(), 2 * len0, lf / 2, cH0.get());

            n1 = 2 * len0 + lf;

            if (wt.wave().lprLen() == wt.wave().hprLen() && (wt.convMethod() == ConvolutionMethod::fft)) {
                wt.cobj = convInit(n1, lf);
                wt.cfftset = 1;
            } else if (!(wt.wave().lpdLen() == wt.wave().hpdLen())) {
                printf("Decomposition Filters must have the same length.");
                exit(-1);
            }

            wconv(wt, cL0.get(), n1, wt.wave().lpr(), lf, oup00L.get());

            wconv(wt, cH0.get(), n1, wt.wave().hpr(), lf, oup00H.get());

            for (auto i = lf - 1; i < 2 * len0 + lf - 1; ++i) {
                oup00[i - lf + 1] = oup00L[i] + oup00H[i];
            }

            //SHIFT 1

            len0 = 0;

            for (indexShift = 1; indexShift < len; indexShift += 2) {
                appx2[len0] = appx1[indexShift];
                det2[len0] = det1[indexShift];
                len0++;
            }

            upsamp2(appx2.get(), len0, u, tempx.get());
            perExt(tempx.get(), 2 * len0, lf / 2, cL0.get());

            upsamp2(det2.get(), len0, u, tempx.get());
            perExt(tempx.get(), 2 * len0, lf / 2, cH0.get());

            n1 = 2 * len0 + lf;

            wconv(wt, cL0.get(), n1, wt.wave().lpr(), lf, oup00L.get());
            wconv(wt, cH0.get(), n1, wt.wave().hpr(), lf, oup00H.get());

            for (auto i = lf - 1; i < 2 * len0 + lf - 1; ++i) {
                oup01[i - lf + 1] = oup00L[i] + oup00H[i];
            }

            circshift(oup01.get(), 2 * len0, -1);

            index2 = 0;

            for (index = count; index < n; index += value) {
                swtop[index] = (oup00[index2] + oup01[index2]) / 2.0;
                index2++;
            }
        }
        for (auto i = 0; i < n; ++i) {
            appxSig[i] = swtop[i];
        }
    }
}

static auto modwtPer(WaveletTransform& wt, int m, double const* inp, double* cA, int lenCA, double* cD) -> void
{
    auto const lenAvg = wt.wave().lpdLen();
    auto filt = std::make_unique<double[]>(2 * lenAvg);
    auto s = std::sqrt(2.0);

    for (auto i = 0; i < lenAvg; ++i) {
        filt[i] = wt.wave().lpd()[i] / s;
        filt[lenAvg + i] = wt.wave().hpd()[i] / s;
    }

    for (auto i = 0; i < lenCA; ++i) {
        auto t = i;
        cA[i] = filt[0] * inp[t];
        cD[i] = filt[lenAvg] * inp[t];
        for (auto l = 1; l < lenAvg; l++) {
            t -= m;
            while (t >= lenCA) {
                t -= lenCA;
            }
            while (t < 0) {
                t += lenCA;
            }

            cA[i] += filt[l] * inp[t];
            cD[i] += filt[lenAvg + l] * inp[t];
        }
    }
}

static auto modwtDirect(WaveletTransform& wt, double const* inp) -> void
{
    if (wt.extension() != SignalExtension::periodic) {
        printf("MODWT direct method only uses periodic extension per. \n");
        printf(" Use MODWT fft method for symmetric extension sym \n");
        exit(-1);
    }

    auto tempLen = wt.siglength;
    auto j = wt.levels();
    wt.length[0] = wt.length[j] = tempLen;
    wt.outlength = wt.length[j + 1] = (j + 1) * tempLen;
    auto m = 1;
    for (auto iter = 1; iter < j; ++iter) {
        m = 2 * m;
        wt.length[iter] = tempLen;
    }

    auto cA = std::make_unique<double[]>(tempLen);
    auto cD = std::make_unique<double[]>(tempLen);

    m = 1;

    for (auto i = 0; i < tempLen; ++i) {
        wt.params[i] = inp[i];
    }

    auto lenacc = wt.outlength;

    for (auto iter = 0; iter < j; ++iter) {
        lenacc -= tempLen;
        if (iter > 0) {
            m = 2 * m;
        }

        modwtPer(wt, m, wt.params.get(), cA.get(), tempLen, cD.get());

        for (auto i = 0; i < tempLen; ++i) {
            wt.params[i] = cA[i];
            wt.params[lenacc + i] = cD[i];
        }
    }
}

static auto modwtFft(WaveletTransform& wt, double const* inp) -> void
{
    int j;
    int iter;
    int m;
    int lenacc;
    double s;
    double tmp1;
    double tmp2;

    auto tempLen = wt.siglength;
    auto lenAvg = wt.wave().lpdLen();
    int n { 0 };
    if (wt.extension() == SignalExtension::symmetric) {
        n = 2 * tempLen;
    } else if (wt.extension() == SignalExtension::periodic) {
        n = tempLen;
    }
    j = wt.levels();
    wt.modwtsiglength = n;
    wt.length[0] = wt.length[j] = n;
    wt.outlength = wt.length[j + 1] = (j + 1) * n;

    s = std::sqrt(2.0);
    for (iter = 1; iter < j; ++iter) {
        wt.length[iter] = n;
    }

    auto fftFd = fftInit(n, 1);
    auto fftBd = fftInit(n, -1);

    auto sig = std::make_unique<FftData[]>(n);
    auto cA = std::make_unique<FftData[]>(n);
    auto cD = std::make_unique<FftData[]>(n);
    auto lowPass = std::make_unique<FftData[]>(n);
    auto highPass = std::make_unique<FftData[]>(n);
    auto index = std::make_unique<int[]>(n);

    // N-point FFT of low pass and high pass filters

    // Low Pass Filter

    for (auto i = 0; i < lenAvg; ++i) {
        sig[i].re = (fft_type)wt.wave().lpd()[i] / s;
        sig[i].im = 0.0;
    }
    for (auto i = lenAvg; i < n; ++i) {
        sig[i].re = 0.0;
        sig[i].im = 0.0;
    }

    fftExec(*fftFd, sig.get(), lowPass.get());

    // High Pass Filter

    for (auto i = 0; i < lenAvg; ++i) {
        sig[i].re = (fft_type)wt.wave().hpd()[i] / s;
        sig[i].im = 0.0;
    }
    for (auto i = lenAvg; i < n; ++i) {
        sig[i].re = 0.0;
        sig[i].im = 0.0;
    }

    fftExec(*fftFd, sig.get(), highPass.get());

    // symmetric extension
    for (auto i = 0; i < tempLen; ++i) {
        sig[i].re = (fft_type)inp[i];
        sig[i].im = 0.0;
    }
    for (auto i = tempLen; i < n; ++i) {
        sig[i].re = (fft_type)inp[n - i - 1];
        sig[i].im = 0.0;
    }

    // FFT of data

    fftExec(*fftFd, sig.get(), cA.get());

    lenacc = wt.outlength;

    m = 1;

    for (iter = 0; iter < j; ++iter) {
        lenacc -= n;

        for (auto i = 0; i < n; ++i) {
            index[i] = (m * i) % n;
        }

        for (auto i = 0; i < n; ++i) {
            tmp1 = cA[i].re;
            tmp2 = cA[i].im;
            cA[i].re = lowPass[index[i]].re * tmp1 - lowPass[index[i]].im * tmp2;
            cA[i].im = lowPass[index[i]].re * tmp2 + lowPass[index[i]].im * tmp1;

            cD[i].re = highPass[index[i]].re * tmp1 - highPass[index[i]].im * tmp2;
            cD[i].im = highPass[index[i]].re * tmp2 + highPass[index[i]].im * tmp1;
        }

        fftExec(*fftBd, cD.get(), sig.get());

        for (auto i = 0; i < n; ++i) {
            wt.params[lenacc + i] = sig[i].re / n;
        }

        m *= 2;
    }

    fftExec(*fftBd, cA.get(), sig.get());

    for (auto i = 0; i < n; ++i) {
        wt.params[i] = sig[i].re / n;
    }
}

auto modwt(WaveletTransform& wt, double const* inp) -> void
{
    if (wt.convMethod() == ConvolutionMethod::direct) {
        modwtDirect(wt, inp);
        return;
    }

    modwtFft(wt, inp);
}

static auto conjComplex(FftData* x, int n) -> void
{
    for (auto i = 0; i < n; ++i) {
        x[i].im *= (-1.0);
    }
}

auto imodwtFft(WaveletTransform& wt, double* oup) -> void
{
    auto n = wt.modwtsiglength;
    auto lenAvg = wt.wave().lpdLen();
    auto j = wt.levels();

    auto s = std::sqrt(2.0);
    auto fftFd = fftInit(n, 1);
    auto fftBd = fftInit(n, -1);

    auto sig = std::make_unique<FftData[]>(n);
    auto cA = std::make_unique<FftData[]>(n);
    auto cD = std::make_unique<FftData[]>(n);
    auto lowPass = std::make_unique<FftData[]>(n);
    auto highPass = std::make_unique<FftData[]>(n);
    auto index = std::make_unique<int[]>(n);

    // N-point FFT of low pass and high pass filters

    // Low Pass Filter

    for (auto i = 0; i < lenAvg; ++i) {
        sig[i].re = (fft_type)wt.wave().lpd()[i] / s;
        sig[i].im = 0.0;
    }
    for (auto i = lenAvg; i < n; ++i) {
        sig[i].re = 0.0;
        sig[i].im = 0.0;
    }

    fftExec(*fftFd, sig.get(), lowPass.get());

    // High Pass Filter

    for (auto i = 0; i < lenAvg; ++i) {
        sig[i].re = (fft_type)wt.wave().hpd()[i] / s;
        sig[i].im = 0.0;
    }
    for (auto i = lenAvg; i < n; ++i) {
        sig[i].re = 0.0;
        sig[i].im = 0.0;
    }

    fftExec(*fftFd, sig.get(), highPass.get());

    // Complex conjugate of the two filters

    conjComplex(lowPass.get(), n);
    conjComplex(highPass.get(), n);

    auto m = (int)std::pow(2.0, (double)j - 1.0);
    auto lenacc = n;

    //
    for (auto i = 0; i < n; ++i) {
        sig[i].re = (fft_type)wt.output()[i];
        sig[i].im = 0.0;
    }

    for (auto iter = 0; iter < j; ++iter) {
        fftExec(*fftFd, sig.get(), cA.get());
        for (auto i = 0; i < n; ++i) {
            sig[i].re = wt.output()[lenacc + i];
            sig[i].im = 0.0;
        }
        fftExec(*fftFd, sig.get(), cD.get());

        for (auto i = 0; i < n; ++i) {
            index[i] = (m * i) % n;
        }

        for (auto i = 0; i < n; ++i) {
            auto const tmp1 = cA[i].re;
            auto const tmp2 = cA[i].im;
            cA[i].re = lowPass[index[i]].re * tmp1 - lowPass[index[i]].im * tmp2 + highPass[index[i]].re * cD[i].re - highPass[index[i]].im * cD[i].im;
            cA[i].im = lowPass[index[i]].re * tmp2 + lowPass[index[i]].im * tmp1 + highPass[index[i]].re * cD[i].im + highPass[index[i]].im * cD[i].re;
        }

        fftExec(*fftBd, cA.get(), sig.get());

        for (auto i = 0; i < n; ++i) {
            sig[i].re /= n;
            sig[i].im /= n;
        }
        m /= 2;
        lenacc += n;
    }

    for (auto i = 0; i < wt.siglength; ++i) {
        oup[i] = sig[i].re;
    }
}

static auto imodwtPer(WaveletTransform& wt, int m, double const* cA, int lenCA, double const* cD, double* x) -> void
{
    auto const lenAvg = wt.wave().lpdLen();
    auto filt = std::make_unique<double[]>(2 * lenAvg);
    auto s = std::sqrt(2.0);

    for (auto i = 0; i < lenAvg; ++i) {
        filt[i] = wt.wave().lpd()[i] / s;
        filt[lenAvg + i] = wt.wave().hpd()[i] / s;
    }

    for (auto i = 0; i < lenCA; ++i) {
        auto t = i;
        x[i] = (filt[0] * cA[t]) + (filt[lenAvg] * cD[t]);
        for (auto l = 1; l < lenAvg; l++) {
            t += m;
            while (t >= lenCA) {
                t -= lenCA;
            }
            while (t < 0) {
                t += lenCA;
            }

            x[i] += (filt[l] * cA[t]) + (filt[lenAvg + l] * cD[t]);
        }
    }
}

static auto imodwtDirect(WaveletTransform& wt, double* dwtop) -> void
{
    auto n = wt.siglength;
    auto lenacc = n;

    auto j = wt.levels();
    auto m = (int)std::pow(2.0, (double)j - 1.0);

    auto x = std::make_unique<double[]>(n);

    for (auto i = 0; i < n; ++i) {
        dwtop[i] = wt.output()[i];
    }

    for (auto iter = 0; iter < j; ++iter) {
        if (iter > 0) {
            m = m / 2;
        }
        imodwtPer(wt, m, dwtop, n, wt.params.get() + lenacc, x.get());
        /*
		for (auto j = lf - 1; j < N; ++j) {
			dwtop[j - lf + 1] = X[j];
		}
		for (auto j = 0; j < lf - 1; ++j) {
			dwtop[N - lf + 1 + j] = X[j];
		}
		*/
        for (auto jj = 0; jj < n; ++jj) {
            dwtop[jj] = x[jj];
        }

        lenacc += n;
    }
}

auto imodwt(WaveletTransform& wt, double* oup) -> void
{
    if (wt.convMethod() == ConvolutionMethod::direct) {
        imodwtDirect(wt, oup);
        return;
    }
    imodwtFft(wt, oup);
}

auto setWTREEExtension(WaveletTree* wt, char const* extension) -> void
{
    if (extension == "sym"sv) {
        wt->ext = "sym";
    } else if (extension == "per"sv) {
        wt->ext = "per";
    } else {
        printf("Signal extension can be either per or sym");
        exit(-1);
    }
}

auto setDWPTExtension(WaveletPacketTransform* wt, char const* extension) -> void
{
    if (extension == "sym"sv) {
        wt->ext = "sym";
    } else if (extension == "per"sv) {
        wt->ext = "per";
    } else {
        printf("Signal extension can be either per or sym");
        exit(-1);
    }
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

auto setDWPTEntropy(WaveletPacketTransform* wt, char const* entropy, double eparam) -> void
{
    if (strcmp(entropy, "shannon") == 0) {
        wt->entropy = "shannon";
    } else if (strcmp(entropy, "threshold") == 0) {
        wt->entropy = "threshold";
        wt->eparam = eparam;
    } else if (strcmp(entropy, "norm") == 0) {
        wt->entropy = "norm";
        wt->eparam = eparam;
    } else if ((strcmp(entropy, "logenergy") == 0) || (strcmp(entropy, "log energy") == 0) || (strcmp(entropy, "energy") == 0)) {
        wt->entropy = "logenergy";
    } else {
        printf("Entropy should be one of shannon, threshold, norm or logenergy");
        exit(-1);
    }
}

auto dwt2(WaveletTransform2D* wt, double* inp) -> std::unique_ptr<double[]>
{
    int iter;
    int n;
    int rowsI;
    int colsI;
    int ir;
    int ic;
    int istride;
    int ostride;
    int aLL;
    int aLH;
    int aHL;
    int aHH;
    int cdim;
    double* orig;

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

    assert(wt->ext == "sym"sv);

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

auto idwt2(WaveletTransform2D* wt, double* wavecoeff, double* oup) -> void
{

    int ir;
    int ic;

    int istride;
    int ostride;
    int iter;
    int aLL;
    int aLH;
    int aHL;
    int aHH;
    double* orig;

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
    assert(wt->ext == "sym"sv);

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
    int j;
    int iter;
    int m;
    int n;
    int lp;
    int rowsN;
    int colsN;
    int rowsI;
    int colsI;
    int ir;
    int ic;
    int istride;
    int ostride;
    int aLL;
    int aLH;
    int aHL;
    int aHH;
    int cdim;
    int clen;
    double* orig;

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
    int k;
    int iter;
    int it2;
    int it3;
    int j;
    int m;
    int rows;
    int cols;
    int lf;
    int ir;
    int ic;
    int k1;
    int i1;

    int aLL;
    int aLH;
    int aHL;
    int aHH;
    int shift;
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
            idwt2Shift(shift, ir, ic, wt->wave->lpr(), wt->wave->hpr(), wt->wave->lpdLen(), a.get(), h.get(), v.get(), d.get(), oup1.get());
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
            idwt2Shift(shift, ir, ic, wt->wave->lpr(), wt->wave->hpr(), wt->wave->lpdLen(), a.get(), h.get(), v.get(), d.get(), oup2.get());
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

auto modwt2(WaveletTransform2D* wt, double* inp) -> std::unique_ptr<double[]>
{
    int j;
    int iter;
    int m;
    int n;
    int lp;
    int rowsN;
    int colsN;
    int rowsI;
    int colsI;
    int ir;
    int ic;
    int istride;
    int ostride;
    int aLL;
    int aLH;
    int aHL;
    int aHH;
    int cdim;
    int clen;
    double* orig;
    double s;

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

auto imodwt2(WaveletTransform2D* wt, double* wavecoeff, double* oup) -> void
{
    int rows;
    int cols;
    int m;
    // int N;
    int ir;
    int ic;
    int lf;
    int istride;
    int ostride;
    int iter;
    int j;
    int aLL;
    int aLH;
    int aHL;
    int aHH;
    double* orig;
    double s;

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
    int j;
    int iter;
    int t;
    double* ptr;
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

auto summary(Wavelet const& obj) -> void
{
    auto const n = obj.size();
    printf("\n");
    printf("Wavelet Name: %s \n", obj.name().c_str());
    printf("\n");
    printf("Wavelet Filters \n");
    printf("lpd: [");
    for (auto i = 0; i < n - 1; ++i) {
        printf("%g,", obj.lpd()[i]);
    }
    printf("%g", obj.lpd()[n - 1]);
    printf("] \n");
    printf("hpd: [");
    for (auto i = 0; i < n - 1; ++i) {
        printf("%g,", obj.hpd()[i]);
    }
    printf("%g", obj.hpd()[n - 1]);
    printf("] \n");
    printf("lpr: [");
    for (auto i = 0; i < n - 1; ++i) {
        printf("%g,", obj.lpr()[i]);
    }
    printf("%g", obj.lpr()[n - 1]);
    printf("] \n");
    printf("hpr: [");
    for (auto i = 0; i < n - 1; ++i) {
        printf("%g,", obj.hpr()[i]);
    }
    printf("%g", obj.hpr()[n - 1]);
    printf("] \n");
}

auto summary(WaveletTransform const& wt) -> void
{
    int j;
    int t;
    j = wt.levels();
    summary(wt.wave());
    printf("\n");
    printf("Wavelet Transform : %s \n", wt.method().c_str());
    printf("Signal Extension : %s \n", toString(wt.extension()).c_str());
    printf("Convolutional Method : %s \n", toString(wt.convMethod()).c_str());
    printf("Number of Decomposition Levels %d \n", wt.levels());
    printf("Length of Input Signal %d \n", wt.siglength);
    printf("Length of WT Output Vector %d \n", wt.outlength);
    printf("Wavelet Coefficients are contained in vector : %s \n", "output");
    printf("Approximation Coefficients \n");
    printf("Level %d Access : output[%d] Length : %d \n", j, 0, wt.length[0]);
    printf("Detail Coefficients \n");
    t = wt.length[0];
    for (auto i = 0; i < j; ++i) {
        printf("Level %d Access : output[%d] Length : %d \n", j - i, t, wt.length[i + 1]);
        t += wt.length[i + 1];
    }
    printf("\n");
}

auto summary(WaveletTree const& wt) -> void
{
    int k;
    int p2;
    int j;
    int t;
    j = wt.J;
    summary(*wt.wave);
    printf("\n");
    printf("Wavelet Transform : %s \n", wt.method.c_str());
    printf("\n");
    printf("Signal Extension : %s \n", wt.ext.c_str());
    printf("\n");
    printf("Number of Decomposition Levels %d \n", wt.J);
    printf("\n");
    printf("Length of Input Signal %d \n", wt.siglength);
    printf("\n");
    printf("Length of WT Output Vector %d \n", wt.outlength);
    printf("\n");
    printf("Wavelet Coefficients are contained in vector : %s \n", "output");
    printf("\n");
    printf("Coefficients Access \n");
    t = 0;
    p2 = 2;
    for (auto i = 0; i < j; ++i) {
        for (k = 0; k < p2; ++k) {
            printf("Node %d %d Access : output[%d] Length : %d \n", i + 1, k, wt.nodelength[t], wt.length[j - i]);
            t++;
        }
        p2 *= 2;
    }
    printf("\n");
}

auto summary(WaveletPacketTransform const& wt) -> void
{
    int k;
    int p2;
    int j;
    int it1;
    int it2;
    j = wt.J;
    summary(*wt.wave);
    printf("\n");
    printf("Signal Extension : %s \n", wt.ext.c_str());
    printf("\n");
    printf("Entropy : %s \n", wt.entropy.c_str());
    printf("\n");
    printf("Number of Decomposition Levels %d \n", wt.J);
    printf("\n");
    printf("Number of Active Nodes %d \n", wt.nodes);
    printf("\n");
    printf("Length of Input Signal %d \n", wt.siglength);
    printf("\n");
    printf("Length of WT Output Vector %d \n", wt.outlength);
    printf("\n");
    printf("Wavelet Coefficients are contained in vector : %s \n", "output");
    printf("\n");
    printf("Coefficients Access \n");
    it1 = 1;
    it2 = 0;
    for (auto i = 0; i < j; ++i) {
        it1 += ipow2(i + 1);
    }
    for (auto i = j; i > 0; --i) {
        p2 = ipow2(i);
        it1 -= p2;
        for (k = 0; k < p2; ++k) {
            if (wt.basisvector[it1 + k] == 1) {
                printf("Node %d %d Access : output[%d] Length : %d \n", i, k, it2, wt.length[j - i + 1]);
                it2 += wt.length[j - i + 1];
            }
        }
    }

    printf("\n");
}

auto summary(ComplexWaveletTransform const& wt) -> void
{

    printf("\n");
    printf("Wavelet : %s Parameter %lf \n", wt.wave.c_str(), wt.m);
    printf("\n");
    printf("Length of Input Signal : %d \n", wt.siglength);
    printf("\n");
    printf("Sampling Rate : %g \n", wt.dt);
    printf("\n");
    printf("Total Number of Scales : %d \n", wt.J);
    printf("\n");
    printf("Smallest Scale (s0) : %lf \n", wt.s0);
    printf("\n");
    printf("Separation Between Scales (dj) %lf \n", wt.dj);
    printf("\n");
    printf("Scale Type %s \n", wt.type.c_str());
    printf("\n");
    printf("Complex CWT Output Vector is of size %d  const& %d stored in Row Major format \n", wt.J, wt.siglength);
    printf("\n");
    printf("The ith real value can be accessed using wt.output()[i].re and imaginary value by wt.output()[i].im \n");
    printf("\n");
}

auto summary(WaveletTransform2D const& wt) -> void
{
    int j;
    int t;
    int rows;
    int cols;
    int vsize;
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

auto wtreeFree(WaveletTree* object) -> void
{
    delete object;
}

auto wptFree(WaveletPacketTransform* object) -> void
{
    delete object;
}

auto cwtFree(ComplexWaveletTransform* object) -> void
{
    delete object;
}

auto wt2Free(WaveletTransform2D* wt) -> void
{
    delete wt;
}
