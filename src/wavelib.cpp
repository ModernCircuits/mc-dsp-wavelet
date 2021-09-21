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

wavelet::wavelet(char const* name)
    : name_ { name }
    , size_ { static_cast<std::size_t>(::filtlength(name)) }
    , params_ { std::make_unique<double[]>(4 * size_) }
    , lpd_ { &params_[0], size_ }
    , hpd_ { &params_[size_], size_ }
    , lpr_ { &params_[2 * size_], size_ }
    , hpr_ { &params_[3 * size_], size_ }
{
    auto* p = params_.get();
    if (name != nullptr) {
        filtcoef(name, p, p + size_, p + 2 * size_, p + 3 * size_);
    }
}

auto wt_init(wavelet& wave, char const* method, int siglength, int J) -> wt_set*
{
    int size;
    int MaxIter;
    auto obj = std::unique_ptr<wt_set>(nullptr);

    size = wave.size();

    if (J > 100) {
        printf("\n The Decomposition Iterations Cannot Exceed 100. Exiting \n");
        exit(-1);
    }

    MaxIter = wmaxiter(siglength, size);

    if (J > MaxIter) {
        printf("\n Error - The Signal Can only be iterated %d times using this wavelet. Exiting\n", MaxIter);
        exit(-1);
    }

    if (method == nullptr) {
        obj = std::make_unique<wt_set>();
        obj->params = std::make_unique<double[]>(siglength + 2 * J * (size + 1));
        obj->outlength = siglength + 2 * J * (size + 1);
        obj->ext = "sym";
    } else if ((method == "dwt"sv) || (method == "DWT"sv)) {
        obj = std::make_unique<wt_set>();
        obj->params = std::make_unique<double[]>(siglength + 2 * J * (size + 1));
        obj->outlength = siglength + 2 * J * (size + 1);
        obj->ext = "sym";
    } else if ((method == "swt"sv) || (method == "SWT"sv)) {
        if (testSWTlength(siglength, J) == 0) {
            printf("\n For SWT the signal length must be a multiple of 2^J. \n");
            exit(-1);
        }

        obj = std::make_unique<wt_set>();
        obj->params = std::make_unique<double[]>(siglength * (J + 1));
        obj->outlength = siglength * (J + 1);
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

        obj = std::make_unique<wt_set>();
        obj->params = std::make_unique<double[]>(siglength * 2 * (J + 1));
        obj->outlength = siglength * (J + 1);
        obj->ext = "per";
    }

    obj->wave = &wave;
    obj->siglength = siglength;
    obj->modwtsiglength = siglength;
    obj->J = J;
    obj->MaxIter = MaxIter;
    obj->method = method;

    if (siglength % 2 == 0) {
        obj->even = 1;
    } else {
        obj->even = 0;
    }

    obj->cobj = nullptr;

    obj->cmethod = "direct";
    obj->cfftset = 0;
    obj->lenlength = J + 2;
    obj->output = &obj->params[0];
    if ((method == "dwt"sv) || (method == "DWT"sv)) {
        for (auto i = 0; i < siglength + 2 * J * (size + 1); ++i) {
            obj->params[i] = 0.0;
        }
    } else if ((method == "swt"sv) || (method == "SWT"sv)) {
        for (auto i = 0; i < siglength * (J + 1); ++i) {
            obj->params[i] = 0.0;
        }
    } else if ((method == "modwt"sv) || (method == "MODWT"sv)) {
        for (auto i = 0; i < siglength * 2 * (J + 1); ++i) {
            obj->params[i] = 0.0;
        }
    }

    return obj.release();
}

auto wtree_init(wavelet* wave, int siglength, int J) -> wtree_set*
{
    auto const size = wave->size();
    auto const MaxIter = wmaxiter(siglength, size);
    if (J > 100) {
        printf("\n The Decomposition Iterations Cannot Exceed 100. Exiting \n");
        exit(-1);
    }
    if (J > MaxIter) {
        printf("\n Error - The Signal Can only be iterated %d times using this wavelet. Exiting\n", MaxIter);
        exit(-1);
    }

    auto temp = 1;
    auto elength = 0;
    auto nodes = 0;
    for (auto i = 0; i < J; ++i) {
        temp *= 2;
        nodes += temp;
        auto const temp2 = (size - 2) * (temp - 1);
        elength += temp2;
    }

    auto obj = std::make_unique<wtree_set>();
    obj->params = std::make_unique<double[]>(siglength * (J + 1) + elength + nodes + J + 1);
    obj->outlength = siglength * (J + 1) + elength;
    obj->ext = "sym";

    obj->wave = wave;
    obj->siglength = siglength;
    obj->J = J;
    obj->MaxIter = MaxIter;
    obj->method = "dwt";

    if (siglength % 2 == 0) {
        obj->even = 1;
    } else {
        obj->even = 0;
    }

    obj->cobj = nullptr;
    obj->nodes = nodes;

    obj->cfftset = 0;
    obj->lenlength = J + 2;
    obj->output = &obj->params[0];
    obj->nodelength = (int*)&obj->params[siglength * (J + 1) + elength];
    obj->coeflength = (int*)&obj->params[siglength * (J + 1) + elength + nodes];

    for (auto i = 0; i < siglength * (J + 1) + elength + nodes + J + 1; ++i) {
        obj->params[i] = 0.0;
    }

    return obj.release();
}

auto wpt_init(wavelet* wave, int siglength, int J) -> wpt_set*
{
    auto const size = wave->size();

    if (J > 100) {
        printf("\n The Decomposition Iterations Cannot Exceed 100. Exiting \n");
        exit(-1);
    }

    auto const MaxIter = wmaxiter(siglength, size);
    if (J > MaxIter) {
        printf("\n Error - The Signal Can only be iterated %d times using this wavelet. Exiting\n", MaxIter);
        exit(-1);
    }
    auto temp = 1;
    auto nodes = 0;
    for (auto i = 0; i < J; ++i) {
        temp *= 2;
        nodes += temp;
    }

    auto idx = J;
    auto p2 = 2;
    auto N = siglength;
    auto lp = size;
    auto elength = 0;
    while (idx > 0) {
        N = N + lp - 2;
        N = (int)ceil((double)N / 2.0);
        elength = p2 * N;
        idx--;
        p2 *= 2;
    }

    auto obj = std::make_unique<wpt_set>();
    obj->params = std::make_unique<double[]>(elength + 4 * nodes + 2 * J + 6);
    obj->outlength = siglength + 2 * (J + 1) * (size + 1);
    obj->ext = "sym";
    obj->entropy = "shannon";
    obj->eparam = 0.0;

    obj->wave = wave;
    obj->siglength = siglength;
    obj->J = J;
    obj->MaxIter = MaxIter;

    if (siglength % 2 == 0) {
        obj->even = 1;
    } else {
        obj->even = 0;
    }

    obj->cobj = nullptr;
    obj->nodes = nodes;

    obj->lenlength = J + 2;
    obj->output = &obj->params[0];
    obj->costvalues = &obj->params[elength];
    obj->basisvector = &obj->params[elength + nodes + 1];
    obj->nodeindex = (int*)&obj->params[elength + 2 * nodes + 2];
    obj->numnodeslevel = (int*)&obj->params[elength + 4 * nodes + 4];
    obj->coeflength = (int*)&obj->params[elength + 4 * nodes + J + 5];

    for (auto i = 0; i < elength + 4 * nodes + 2 * J + 6; ++i) {
        obj->params[i] = 0.0;
    }

    return obj.release();
}

auto cwt_init(char const* wave, double param, int siglength, double dt, int J) -> cwt_set*
{
    int N;
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

    N = siglength;
    nj2 = 2 * N * J;
    auto obj = std::make_unique<cwt_set>();
    obj->params = std::make_unique<double[]>(nj2 + 2 * J + N);

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
    obj->J = J;
    obj->siglength = siglength;
    obj->sflag = 0;
    obj->pflag = 1;
    obj->mother = mother;
    obj->m = param;

    t1 = 0.499999 + std::log((double)N) / std::log(2.0);
    ibase2 = 1 + (int)t1;

    obj->npad = (int)std::pow(2.0, (double)ibase2);

    obj->output = (cplx_data*)&obj->params[0];
    obj->scale = &obj->params[nj2];
    obj->period = &obj->params[nj2 + J];
    obj->coi = &obj->params[nj2 + 2 * J];

    for (auto i = 0; i < nj2 + 2 * J + N; ++i) {
        obj->params[i] = 0.0;
    }

    return obj.release();
}

auto wt2_init(wavelet& wave, char const* method, int rows, int cols, int J) -> wt2_set*
{

    auto const size = wave.size();

    auto const MaxRows = wmaxiter(rows, size);
    auto const MaxCols = wmaxiter(cols, size);

    auto const MaxIter = (MaxRows < MaxCols) ? MaxRows : MaxCols;

    if (J > MaxIter) {
        printf("\n Error - The Signal Can only be iterated %d times using this wavelet. Exiting\n", MaxIter);
        exit(-1);
    }

    int sumacc { 0 };
    if (J == 1) {
        sumacc = 4;
    } else if (J > 1) {
        sumacc = J * 3 + 1;
    } else {
        printf("Error : J should be >= 1 \n");
        exit(-1);
    }

    auto obj = std::make_unique<wt2_set>();
    obj->params = std::make_unique<int[]>(2 * J + sumacc);
    obj->outlength = 0;
    if (method == nullptr) {
        obj->ext = "per";
    } else if ((method == "dwt"sv) || (method == "DWT"sv)) {
        obj->ext = "per";
    } else if ((method == "swt"sv) || (method == "SWT"sv)) {
        if ((testSWTlength(rows, J) == 0) || (testSWTlength(cols, J) == 0)) {
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
    obj->J = J;
    obj->MaxIter = MaxIter;
    obj->method = method;
    obj->coeffaccesslength = sumacc;

    obj->dimensions = &obj->params[0];
    obj->coeffaccess = &obj->params[2 * J];
    for (auto i = 0; i < (2 * J + sumacc); ++i) {
        obj->params[i] = 0;
    }

    return obj.release();
}

static void wconv(wt_set* wt, double* sig, int N, double const* filt, int L, double* oup)
{
    if (wt->cmethod == "direct"sv) {
        conv_direct(sig, N, filt, L, oup);
    } else if ((wt->cmethod == "fft"sv) || (wt->cmethod == "FFT"sv)) {
        if (wt->cfftset == 0) {
            wt->cobj = conv_init(N, L);
            conv_fft(*wt->cobj, sig, filt, oup);
        } else {
            conv_fft(*wt->cobj, sig, filt, oup);
        }
    } else {
        printf("Convolution Only accepts two methods - direct and fft");
        exit(-1);
    }
}

static void dwt_per(wt_set* wt, double* inp, int N, double* cA, int len_cA, double* cD)
{

    dwt_per_stride(inp, N, wt->wave->lpd(), wt->wave->hpd(), wt->wave->lpd_len(), cA, len_cA, cD, 1, 1);
}

static void wtree_per(wtree_set* wt, double const* inp, int N, double* cA, int len_cA, double* cD)
{
    int l;
    int l2;
    int isodd;
    int t;
    int len_avg;

    len_avg = wt->wave->lpd_len();
    l2 = len_avg / 2;
    isodd = N % 2;

    for (auto i = 0; i < len_cA; ++i) {
        t = 2 * i + l2;
        cA[i] = 0.0;
        cD[i] = 0.0;
        for (l = 0; l < len_avg; ++l) {
            if ((t - l) >= l2 && (t - l) < N) {
                cA[i] += wt->wave->lpd()[l] * inp[t - l];
                cD[i] += wt->wave->hpd()[l] * inp[t - l];
            } else if ((t - l) < l2 && (t - l) >= 0) {
                cA[i] += wt->wave->lpd()[l] * inp[t - l];
                cD[i] += wt->wave->hpd()[l] * inp[t - l];
            } else if ((t - l) < 0 && isodd == 0) {
                cA[i] += wt->wave->lpd()[l] * inp[t - l + N];
                cD[i] += wt->wave->hpd()[l] * inp[t - l + N];
            } else if ((t - l) < 0 && isodd == 1) {
                if ((t - l) != -1) {
                    cA[i] += wt->wave->lpd()[l] * inp[t - l + N + 1];
                    cD[i] += wt->wave->hpd()[l] * inp[t - l + N + 1];
                } else {
                    cA[i] += wt->wave->lpd()[l] * inp[N - 1];
                    cD[i] += wt->wave->hpd()[l] * inp[N - 1];
                }
            } else if ((t - l) >= N && isodd == 0) {
                cA[i] += wt->wave->lpd()[l] * inp[t - l - N];
                cD[i] += wt->wave->hpd()[l] * inp[t - l - N];
            } else if ((t - l) >= N && isodd == 1) {
                if (t - l != N) {
                    cA[i] += wt->wave->lpd()[l] * inp[t - l - (N + 1)];
                    cD[i] += wt->wave->hpd()[l] * inp[t - l - (N + 1)];
                } else {
                    cA[i] += wt->wave->lpd()[l] * inp[N - 1];
                    cD[i] += wt->wave->hpd()[l] * inp[N - 1];
                }
            }
        }
    }
}

static void dwpt_per(wpt_set* wt, double const* inp, int N, double* cA, int len_cA, double* cD)
{
    int l;
    int l2;
    int isodd;
    int t;
    int len_avg;

    len_avg = wt->wave->lpd_len();
    l2 = len_avg / 2;
    isodd = N % 2;

    for (auto i = 0; i < len_cA; ++i) {
        t = 2 * i + l2;
        cA[i] = 0.0;
        cD[i] = 0.0;
        for (l = 0; l < len_avg; ++l) {
            if ((t - l) >= l2 && (t - l) < N) {
                cA[i] += wt->wave->lpd()[l] * inp[t - l];
                cD[i] += wt->wave->hpd()[l] * inp[t - l];
            } else if ((t - l) < l2 && (t - l) >= 0) {
                cA[i] += wt->wave->lpd()[l] * inp[t - l];
                cD[i] += wt->wave->hpd()[l] * inp[t - l];
            } else if ((t - l) < 0 && isodd == 0) {
                cA[i] += wt->wave->lpd()[l] * inp[t - l + N];
                cD[i] += wt->wave->hpd()[l] * inp[t - l + N];
            } else if ((t - l) < 0 && isodd == 1) {
                if ((t - l) != -1) {
                    cA[i] += wt->wave->lpd()[l] * inp[t - l + N + 1];
                    cD[i] += wt->wave->hpd()[l] * inp[t - l + N + 1];
                } else {
                    cA[i] += wt->wave->lpd()[l] * inp[N - 1];
                    cD[i] += wt->wave->hpd()[l] * inp[N - 1];
                }
            } else if ((t - l) >= N && isodd == 0) {
                cA[i] += wt->wave->lpd()[l] * inp[t - l - N];
                cD[i] += wt->wave->hpd()[l] * inp[t - l - N];
            } else if ((t - l) >= N && isodd == 1) {
                if (t - l != N) {
                    cA[i] += wt->wave->lpd()[l] * inp[t - l - (N + 1)];
                    cD[i] += wt->wave->hpd()[l] * inp[t - l - (N + 1)];
                } else {
                    cA[i] += wt->wave->lpd()[l] * inp[N - 1];
                    cD[i] += wt->wave->hpd()[l] * inp[N - 1];
                }
            }
        }
    }
}

static void dwt_sym(wt_set* wt, double* inp, int N, double* cA, int len_cA, double* cD)
{

    dwt_sym_stride(inp, N, wt->wave->lpd(), wt->wave->hpd(), wt->wave->lpd_len(), cA, len_cA, cD, 1, 1);
}

static void wtree_sym(wtree_set* wt, double const* inp, int N, double* cA, int len_cA, double* cD)
{
    int l;
    int t;
    int len_avg;

    len_avg = wt->wave->lpd_len();

    for (auto i = 0; i < len_cA; ++i) {
        t = 2 * i + 1;
        cA[i] = 0.0;
        cD[i] = 0.0;
        for (l = 0; l < len_avg; ++l) {
            if ((t - l) >= 0 && (t - l) < N) {
                cA[i] += wt->wave->lpd()[l] * inp[t - l];
                cD[i] += wt->wave->hpd()[l] * inp[t - l];
            } else if ((t - l) < 0) {
                cA[i] += wt->wave->lpd()[l] * inp[-t + l - 1];
                cD[i] += wt->wave->hpd()[l] * inp[-t + l - 1];
            } else if ((t - l) >= N) {
                cA[i] += wt->wave->lpd()[l] * inp[2 * N - t + l - 1];
                cD[i] += wt->wave->hpd()[l] * inp[2 * N - t + l - 1];
            }
        }
    }
}

static void dwpt_sym(wpt_set* wt, double const* inp, int N, double* cA, int len_cA, double* cD)
{
    int l;
    int t;
    int len_avg;

    len_avg = wt->wave->lpd_len();

    for (auto i = 0; i < len_cA; ++i) {
        t = 2 * i + 1;
        cA[i] = 0.0;
        cD[i] = 0.0;
        for (l = 0; l < len_avg; ++l) {
            if ((t - l) >= 0 && (t - l) < N) {
                cA[i] += wt->wave->lpd()[l] * inp[t - l];
                cD[i] += wt->wave->hpd()[l] * inp[t - l];
            } else if ((t - l) < 0) {
                cA[i] += wt->wave->lpd()[l] * inp[-t + l - 1];
                cD[i] += wt->wave->hpd()[l] * inp[-t + l - 1];
            } else if ((t - l) >= N) {
                cA[i] += wt->wave->lpd()[l] * inp[2 * N - t + l - 1];
                cD[i] += wt->wave->hpd()[l] * inp[2 * N - t + l - 1];
            }
        }
    }
}

static void dwt1(wt_set* wt, double* sig, int len_sig, double* cA, double* cD)
{
    constexpr auto D = 2;

    if (wt->ext == "per"sv) {
        auto len_avg = (wt->wave->lpd_len() + wt->wave->hpd_len()) / 2;
        auto signal = std::make_unique<double[]>(len_sig + len_avg + (len_sig % 2));
        len_sig = per_ext(sig, len_sig, len_avg / 2, signal.get());
        auto cA_undec = std::make_unique<double[]>(len_sig + len_avg + wt->wave->lpd_len() - 1);

        if (wt->wave->lpd_len() == wt->wave->hpd_len() && ((wt->cmethod == "fft"sv) || (wt->cmethod == "FFT"sv))) {
            wt->cobj = conv_init(len_sig + len_avg, wt->wave->lpd_len());
            wt->cfftset = 1;
        } else if (!(wt->wave->lpd_len() == wt->wave->hpd_len())) {
            printf("Decomposition Filters must have the same length.");
            exit(-1);
        }

        wconv(wt, signal.get(), len_sig + len_avg, wt->wave->lpd(), wt->wave->lpd_len(), cA_undec.get());
        downsamp(cA_undec.get() + len_avg, len_sig, D, cA);
        wconv(wt, signal.get(), len_sig + len_avg, wt->wave->hpd(), wt->wave->hpd_len(), cA_undec.get());
        downsamp(cA_undec.get() + len_avg, len_sig, D, cD);

    } else if (wt->ext == "sym"sv) {
        auto lf = wt->wave->lpd_len(); // lpd and hpd have the same length
        auto signal = std::make_unique<double[]>(len_sig + 2 * (lf - 1));
        len_sig = symm_ext(sig, len_sig, lf - 1, signal.get());
        auto cA_undec = std::make_unique<double[]>(len_sig + 3 * (lf - 1));

        if (wt->wave->lpd_len() == wt->wave->hpd_len() && ((wt->cmethod == "fft"sv) || (wt->cmethod == "FFT"sv))) {
            wt->cobj = conv_init(len_sig + 2 * (lf - 1), lf);
            wt->cfftset = 1;
        } else if (!(wt->wave->lpd_len() == wt->wave->hpd_len())) {
            printf("Decomposition Filters must have the same length.");
            exit(-1);
        }

        wconv(wt, signal.get(), len_sig + 2 * (lf - 1), wt->wave->lpd(), wt->wave->lpd_len(), cA_undec.get());
        downsamp(cA_undec.get() + lf, len_sig + lf - 2, D, cA);
        wconv(wt, signal.get(), len_sig + 2 * (lf - 1), wt->wave->hpd(), wt->wave->hpd_len(), cA_undec.get());
        downsamp(cA_undec.get() + lf, len_sig + lf - 2, D, cD);
    } else {
        printf("Signal extension can be either per or sym");
        exit(-1);
    }

    if (wt->wave->lpd_len() == wt->wave->hpd_len() && ((wt->cmethod == "fft"sv) || (wt->cmethod == "FFT"sv))) {

        wt->cfftset = 0;
    }
}

void dwt(wt_set* wt, double const* inp)
{

    auto temp_len = wt->siglength;
    auto const J = wt->J;

    wt->length[J + 1] = temp_len;
    wt->outlength = 0;
    wt->zpad = 0;

    auto orig2 = std::make_unique<double[]>(temp_len);
    auto orig = std::make_unique<double[]>(temp_len);

    for (auto i = 0; i < wt->siglength; ++i) {
        orig[i] = inp[i];
    }

    if (wt->zpad == 1) {
        orig[temp_len - 1] = orig[temp_len - 2];
    }

    auto N = temp_len;
    auto lp = wt->wave->lpd_len();

    if (wt->ext == "per"sv) {
        auto idx = J;
        while (idx > 0) {
            N = (int)ceil((double)N / 2.0);
            wt->length[idx] = N;
            wt->outlength += wt->length[idx];
            idx--;
        }
        wt->length[0] = wt->length[1];
        wt->outlength += wt->length[0];
        N = wt->outlength;

        for (auto iter = 0; iter < J; ++iter) {
            auto const len_cA = wt->length[J - iter];
            N -= len_cA;
            if ((wt->cmethod == "fft"sv) || (wt->cmethod == "FFT"sv)) {
                dwt1(wt, orig.get(), temp_len, orig2.get(), wt->params.get() + N);
            } else {
                dwt_per(wt, orig.get(), temp_len, orig2.get(), len_cA, wt->params.get() + N);
            }
            temp_len = wt->length[J - iter];
            if (iter == J - 1) {
                for (auto i = 0; i < len_cA; ++i) {
                    wt->params[i] = orig2[i];
                }
            } else {
                for (auto i = 0; i < len_cA; ++i) {
                    orig[i] = orig2[i];
                }
            }
        }
    } else if (wt->ext == "sym"sv) {
        auto idx = J;
        while (idx > 0) {
            N = N + lp - 2;
            N = (int)ceil((double)N / 2.0);
            wt->length[idx] = N;
            wt->outlength += wt->length[idx];
            idx--;
        }
        wt->length[0] = wt->length[1];
        wt->outlength += wt->length[0];
        N = wt->outlength;

        for (auto iter = 0; iter < J; ++iter) {
            auto const len_cA = wt->length[J - iter];
            N -= len_cA;
            if ((wt->cmethod == "fft"sv) || (wt->cmethod == "FFT"sv)) {
                dwt1(wt, orig.get(), temp_len, orig2.get(), wt->params.get() + N);
            } else {
                dwt_sym(wt, orig.get(), temp_len, orig2.get(), len_cA, wt->params.get() + N);
            }
            temp_len = wt->length[J - iter];

            if (iter == J - 1) {
                for (auto i = 0; i < len_cA; ++i) {
                    wt->params[i] = orig2[i];
                }
            } else {
                for (auto i = 0; i < len_cA; ++i) {
                    orig[i] = orig2[i];
                }
            }
        }
    } else {
        printf("Signal extension can be either per or sym");
        exit(-1);
    }
}

void wtree(wtree_set* wt, double const* inp)
{
    int iter;
    int N;
    int lp;
    int p2;
    int k;
    int N2;
    int Np;
    int len_cA;
    int t;
    int t2;
    int it1;

    auto temp_len = wt->siglength;
    auto J = wt->J;
    wt->length[J + 1] = temp_len;
    wt->outlength = 0;
    wt->zpad = 0;

    auto orig = std::make_unique<double[]>(temp_len);

    for (auto i = 0; i < wt->siglength; ++i) {
        orig[i] = inp[i];
    }

    if (wt->zpad == 1) {
        orig[temp_len - 1] = orig[temp_len - 2];
    }

    N = temp_len;
    lp = wt->wave->lpd_len();

    if (wt->ext == "per"sv) {
        auto i = J;
        p2 = 2;
        while (i > 0) {
            N = (int)ceil((double)N / 2.0);
            wt->length[i] = N;
            wt->outlength += p2 * (wt->length[i]);
            i--;
            p2 *= 2;
        }
        wt->length[0] = wt->length[1];

        N2 = wt->outlength;
        p2 = 1;
        for (iter = 0; iter < J; ++iter) {
            len_cA = wt->length[J - iter];
            N2 -= 2 * p2 * len_cA;
            N = N2;
            for (k = 0; k < p2; ++k) {
                if (iter == 0) {
                    wtree_per(wt, orig.get(), temp_len, wt->params.get() + N, len_cA, wt->params.get() + N + len_cA);
                } else {
                    wtree_per(wt, wt->params.get() + Np + k * temp_len, temp_len, wt->params.get() + N, len_cA, wt->params.get() + N + len_cA);
                }
                N += 2 * len_cA;
            }

            temp_len = wt->length[J - iter];
            p2 = 2 * p2;
            Np = N2;
        }
    } else if (wt->ext == "sym"sv) {
        auto i = J;
        p2 = 2;
        while (i > 0) {
            N = N + lp - 2;
            N = (int)ceil((double)N / 2.0);
            wt->length[i] = N;
            wt->outlength += p2 * (wt->length[i]);
            i--;
            p2 *= 2;
        }
        wt->length[0] = wt->length[1];

        N2 = wt->outlength;
        p2 = 1;

        for (iter = 0; iter < J; ++iter) {
            len_cA = wt->length[J - iter];
            N2 -= 2 * p2 * len_cA;
            N = N2;
            for (k = 0; k < p2; ++k) {
                if (iter == 0) {
                    wtree_sym(wt, orig.get(), temp_len, wt->params.get() + N, len_cA, wt->params.get() + N + len_cA);
                } else {
                    wtree_sym(wt, wt->params.get() + Np + k * temp_len, temp_len, wt->params.get() + N, len_cA, wt->params.get() + N + len_cA);
                }
                N += 2 * len_cA;
            }

            temp_len = wt->length[J - iter];
            p2 = 2 * p2;
            Np = N2;
        }

    } else {
        printf("Signal extension can be either per or sym");
        exit(-1);
    }

    J = wt->J;
    t2 = wt->outlength - 2 * wt->length[J];
    p2 = 2;
    it1 = 0;
    for (auto i = 0; i < J; ++i) {
        t = t2;
        for (k = 0; k < p2; ++k) {
            wt->nodelength[it1] = t;
            it1++;
            t += wt->length[J - i];
        }
        p2 *= 2;
        t2 = t2 - p2 * wt->length[J - i - 1];
    }

    wt->coeflength[0] = wt->siglength;

    for (auto i = 1; i < J + 1; ++i) {
        wt->coeflength[i] = wt->length[J - i + 1];
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

void dwpt(wpt_set* wt, double const* inp)
{
    int iter;
    int p2;
    int k;
    int N2;
    int Np;
    int llb;
    double v1;
    double v2;
    int len_cA;
    int t;

    auto temp_len = wt->siglength;
    auto J = wt->J;
    wt->length[J + 1] = temp_len;
    wt->outlength = 0;
    auto temp = 1;
    auto elength = 0;
    auto size = wt->wave->size();
    auto nodes = wt->nodes;
    auto n1 = nodes + 1;
    for (auto i = 0; i < J; ++i) {
        temp *= 2;
        auto const temp2 = (size - 2) * (temp - 1);
        elength += temp2;
    }

    auto eparam = wt->eparam;
    auto orig = std::make_unique<double[]>(temp_len);
    auto tree = std::make_unique<double[]>((temp_len * (J + 1) + elength));
    auto nodelength = std::make_unique<int[]>(nodes);

    for (auto i = 0; i < wt->siglength; ++i) {
        orig[i] = inp[i];
    }

    for (auto i = 0; i < temp_len * (J + 1) + elength; ++i) {
        tree[i] = 0.0;
    }

    for (auto i = 0; i < nodes + 1; ++i) {
        wt->basisvector[i] = 0.0;
        wt->costvalues[i] = 0.0;
    }

    auto N = temp_len;
    auto lp = wt->wave->lpd_len();
    // p2 = 1;

    //set eparam value here
    wt->costvalues[0] = costfunc(orig.get(), wt->siglength, wt->entropy.c_str(), eparam);
    auto it2 = 1;
    if (wt->ext == "per"sv) {
        auto i = J;
        p2 = 2;
        while (i > 0) {
            N = (int)ceil((double)N / 2.0);
            wt->length[i] = N;
            wt->outlength += p2 * (wt->length[i]);
            i--;
            p2 *= 2;
        }
        wt->length[0] = wt->length[1];

        N2 = wt->outlength;
        p2 = 1;
        for (iter = 0; iter < J; ++iter) {
            len_cA = wt->length[J - iter];
            N2 -= 2 * p2 * len_cA;
            N = N2;
            for (k = 0; k < p2; ++k) {
                if (iter == 0) {
                    dwpt_per(wt, orig.get(), temp_len, tree.get() + N, len_cA, tree.get() + N + len_cA);
                } else {
                    dwpt_per(wt, tree.get() + Np + k * temp_len, temp_len, tree.get() + N, len_cA, tree.get() + N + len_cA);
                }
                wt->costvalues[it2] = costfunc(tree.get() + N, len_cA, wt->entropy.c_str(), eparam);
                it2++;
                wt->costvalues[it2] = costfunc(tree.get() + N + len_cA, len_cA, wt->entropy.c_str(), eparam);
                it2++;
                N += 2 * len_cA;
            }

            temp_len = wt->length[J - iter];
            p2 = 2 * p2;
            Np = N2;
        }
    } else if (wt->ext == "sym"sv) {
        auto i = J;
        p2 = 2;
        while (i > 0) {
            N = N + lp - 2;
            N = (int)ceil((double)N / 2.0);
            wt->length[i] = N;
            wt->outlength += p2 * (wt->length[i]);
            i--;
            p2 *= 2;
        }
        wt->length[0] = wt->length[1];

        N2 = wt->outlength;
        p2 = 1;

        for (iter = 0; iter < J; ++iter) {
            len_cA = wt->length[J - iter];
            N2 -= 2 * p2 * len_cA;
            N = N2;
            for (k = 0; k < p2; ++k) {
                if (iter == 0) {
                    dwpt_sym(wt, orig.get(), temp_len, tree.get() + N, len_cA, tree.get() + N + len_cA);
                } else {
                    dwpt_sym(wt, tree.get() + Np + k * temp_len, temp_len, tree.get() + N, len_cA, tree.get() + N + len_cA);
                }
                wt->costvalues[it2] = costfunc(tree.get() + N, len_cA, wt->entropy.c_str(), eparam);
                it2++;
                wt->costvalues[it2] = costfunc(tree.get() + N + len_cA, len_cA, wt->entropy.c_str(), eparam);
                it2++;
                N += 2 * len_cA;
            }

            temp_len = wt->length[J - iter];
            p2 = 2 * p2;
            Np = N2;
        }

    } else {
        printf("Signal extension can be either per or sym");
        exit(-1);
    }

    J = wt->J;
    auto t2 = wt->outlength - 2 * wt->length[J];
    p2 = 2;
    auto it1 = 0;
    for (auto i = 0; i < J; ++i) {
        t = t2;
        for (k = 0; k < p2; ++k) {
            nodelength[it1] = t;
            it1++;
            t += wt->length[J - i];
        }
        p2 *= 2;
        t2 = t2 - p2 * wt->length[J - i - 1];
    }

    J = wt->J;
    llb = 1;
    for (auto i = 0; i < J; ++i) {
        llb *= 2;
    }

    for (auto i = n1 - llb; i < n1; ++i) {
        wt->basisvector[i] = 1;
    }

    for (auto j = J - 1; j >= 0; --j) {
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
        for (auto i = J; i > 0; --i) {
            llb = ipow2(i);
            it1 -= llb;
            wt->numnodeslevel[i] = 0;
            for (auto j = 0; j < llb; ++j) {
                if (wt->basisvector[it1 + j] == 1) {
                    wt->nodeindex[2 * wt->nodes] = i;
                    wt->nodeindex[2 * wt->nodes + 1] = j;
                    wt->nodes += 1;
                    wt->numnodeslevel[i] += 1;
                    for (k = 0; k < wt->length[J - i + 1]; ++k) {
                        wt->output[it2 + k] = tree[nodelength[it1 - 1 + j] + k]; // access tree
                    }
                    it2 += wt->length[J - i + 1];
                }
            }
        }
        wt->outlength = it2;
    }

    wt->coeflength[0] = wt->siglength;

    for (auto i = 1; i < J + 1; ++i) {
        wt->coeflength[i] = wt->length[J - i + 1];
    }
}

/// X - Level. All Nodes at any level have the same length
auto getWTREENodelength(wtree_set* wt, int X) -> int
{
    if (X <= 0 || X > wt->J) {
        printf("X co-ordinate must be >= 1 and <= %d", wt->J);
        exit(-1);
    }

    return wt->length[wt->J - X + 1];
}

/// X - Level. All Nodes at any level have the same length
auto getDWPTNodelength(wpt_set* wt, int X) -> int
{
    if (X <= 0 || X > wt->J) {
        printf("X co-ordinate must be >= 1 and <= %d", wt->J);
        exit(-1);
    }

    return wt->length[wt->J - X + 1];
}

void getWTREECoeffs(wtree_set* wt, int X, int Y, double* coeffs, int N)
{
    int ymax;
    int t;
    int t2;

    if (X <= 0 || X > wt->J) {
        printf("X co-ordinate must be >= 1 and <= %d", wt->J);
        exit(-1);
    }
    ymax = 1;
    for (auto i = 0; i < X; ++i) {
        ymax *= 2;
    }

    ymax -= 1;

    if (Y < 0 || Y > ymax) {
        printf("Y co-ordinate must be >= 0 and <= %d", ymax);
        exit(-1);
    }

    if (X == 1) {
        t = 0;
    } else {
        t = 0;
        t2 = 1;
        for (auto i = 0; i < X - 1; ++i) {
            t2 *= 2;
            t += t2;
        }
    }

    t += Y;
    t2 = wt->nodelength[t];
    for (auto i = 0; i < N; ++i) {
        coeffs[i] = wt->output[t2 + i];
    }
}

void setCWTScales(cwt_set* wt, double s0, double dj, char const* type, int power)
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

void cwt(cwt_set* wt, double const* inp)
{
    int N;
    int npad;
    int nj2;
    int j;
    int j2;
    N = wt->siglength;
    if (wt->sflag == 0) {
        for (auto i = 0; i < wt->J; ++i) {
            wt->scale[i] = wt->s0 * std::pow(2.0, (double)(i)*wt->dj);
        }
        wt->sflag = 1;
    }

    if (wt->pflag == 0) {
        npad = N;
    } else {
        npad = wt->npad;
    }

    nj2 = 2 * N * wt->J;
    j = wt->J;
    j2 = 2 * j;

    wt->smean = 0.0;

    for (auto i = 0; i < N; ++i) {
        wt->smean += inp[i];
    }
    wt->smean /= N;

    cwavelet(inp, N, wt->dt, wt->mother, wt->m, wt->s0, wt->dj, wt->J, npad, wt->params.get(), wt->params.get() + nj2, wt->params.get() + nj2 + j, wt->params.get() + nj2 + j2);
}

void icwt(cwt_set* wt, double* cwtop)
{
    double psi;
    double cdel;
    int real;
    int N;
    int nj2;

    N = wt->siglength;
    nj2 = N * 2 * wt->J;

    psi0(wt->mother, wt->m, &psi, &real);
    cdel = cdelta(wt->mother, wt->m, psi);

    if (((wt->type == "pow"sv) || (wt->type == "power"sv)) && wt->pow == 2) {
        icwavelet(wt->params.get(), N, wt->params.get() + nj2, wt->J, wt->dt, wt->dj, cdel, psi, cwtop);
    } else {
        printf("Inverse CWT is only available for power of 2.0 scales \n");
        exit(-1);
    }
    for (auto i = 0; i < N; ++i) {
        cwtop[i] += wt->smean;
    }
}

static void idwt1(wt_set* wt, double* temp, double* cA_up, double* cA, int len_cA, double* cD, int len_cD, double* X_lp, double* X_hp, double* X)
{
    auto len_avg = (wt->wave->lpr_len() + wt->wave->hpr_len()) / 2;
    auto N = 2 * len_cD;
    auto U = 2;

    upsamp2(cA, len_cA, U, cA_up);

    per_ext(cA_up, 2 * len_cA, len_avg / 2, temp);

    auto N2 = 2 * len_cA + len_avg;

    if (wt->wave->lpr_len() == wt->wave->hpr_len() && ((wt->cmethod == "fft"sv) || (wt->cmethod == "FFT"sv))) {
        wt->cobj = conv_init(N2, len_avg);
        wt->cfftset = 1;
    } else if (!(wt->wave->lpr_len() == wt->wave->hpr_len())) {
        printf("Decomposition Filters must have the same length.");
        exit(-1);
    }

    wconv(wt, temp, N2, wt->wave->lpr(), len_avg, X_lp);

    upsamp2(cD, len_cD, U, cA_up);

    per_ext(cA_up, 2 * len_cD, len_avg / 2, temp);

    N2 = 2 * len_cD + len_avg;

    wconv(wt, temp, N2, wt->wave->hpr(), len_avg, X_hp);

    for (auto i = len_avg - 1; i < N + len_avg - 1; ++i) {
        X[i - len_avg + 1] = X_lp[i] + X_hp[i];
    }

    if (wt->wave->lpr_len() == wt->wave->hpr_len() && ((wt->cmethod == "fft"sv) || (wt->cmethod == "FFT"sv))) {

        wt->cfftset = 0;
    }
}

static void idwt_per(wt_set* wt, double* cA, int len_cA, double* cD, double* X)
{
    idwt_per_stride(cA, len_cA, cD, wt->wave->lpr(), wt->wave->hpr(), wt->wave->lpr_len(), X, 1, 1);
}

static void idwt_sym(wt_set* wt, double* cA, int len_cA, double* cD, double* X)
{
    idwt_sym_stride(cA, len_cA, cD, wt->wave->lpr(), wt->wave->hpr(), wt->wave->lpr_len(), X, 1, 1);
}

void idwt(wt_set* wt, double* dwtop)
{

    int lf;
    int N;
    int N2;
    int iter;
    int k;
    int det_len;

    auto J = wt->J;
    auto U = 2;
    auto app_len = wt->length[0];
    auto out = std::make_unique<double[]>(wt->siglength + 1);
    if ((wt->ext == "per"sv) && ((wt->cmethod == "fft"sv) || (wt->cmethod == "FFT"sv))) {
        app_len = wt->length[0];
        det_len = wt->length[1];
        N = 2 * wt->length[J];
        lf = (wt->wave->lpr_len() + wt->wave->hpr_len()) / 2;

        auto cA_up = std::make_unique<double[]>(N);
        auto temp = std::make_unique<double[]>((N + lf));
        auto X_lp = std::make_unique<double[]>((N + 2 * lf - 1));
        auto X_hp = std::make_unique<double[]>((N + 2 * lf - 1));
        iter = app_len;

        for (auto i = 0; i < app_len; ++i) {
            out[i] = wt->output[i];
        }

        for (auto i = 0; i < J; ++i) {

            idwt1(wt, temp.get(), cA_up.get(), out.get(), det_len, wt->output + iter, det_len, X_lp.get(), X_hp.get(), out.get());
            /*
			idwt_per(wt,out.get(), det_len, wt->output + iter, det_len, X_lp);
			for (k = lf/2 - 1; k < 2 * det_len + lf/2 - 1; ++k) {
				out[k - lf/2 + 1] = X_lp[k];
			}
			*/
            iter += det_len;
            det_len = wt->length[i + 2];
        }

    } else if ((wt->ext == "per"sv) && (wt->cmethod == "direct"sv)) {
        app_len = wt->length[0];
        det_len = wt->length[1];
        N = 2 * wt->length[J];
        lf = (wt->wave->lpr_len() + wt->wave->hpr_len()) / 2;

        auto X_lp = std::make_unique<double[]>((N + 2 * lf - 1));
        iter = app_len;

        for (auto i = 0; i < app_len; ++i) {
            out[i] = wt->output[i];
        }

        for (auto i = 0; i < J; ++i) {
            idwt_per(wt, out.get(), det_len, wt->output + iter, X_lp.get());
            for (k = lf / 2 - 1; k < 2 * det_len + lf / 2 - 1; ++k) {
                out[k - lf / 2 + 1] = X_lp[k];
            }

            iter += det_len;
            det_len = wt->length[i + 2];
        }

    } else if ((wt->ext == "sym"sv) && (wt->cmethod == "direct"sv)) {
        app_len = wt->length[0];
        det_len = wt->length[1];
        N = 2 * wt->length[J] - 1;
        lf = (wt->wave->lpr_len() + wt->wave->hpr_len()) / 2;

        auto X_lp = std::make_unique<double[]>((N + 2 * lf - 1));
        iter = app_len;

        for (auto i = 0; i < app_len; ++i) {
            out[i] = wt->output[i];
        }

        for (auto i = 0; i < J; ++i) {
            idwt_sym(wt, out.get(), det_len, wt->output + iter, X_lp.get());
            for (k = lf - 2; k < 2 * det_len; ++k) {
                out[k - lf + 2] = X_lp[k];
            }

            iter += det_len;
            det_len = wt->length[i + 2];
        }

    } else if ((wt->ext == "sym"sv) && ((wt->cmethod == "fft"sv) || (wt->cmethod == "FFT"sv))) {
        lf = wt->wave->lpd_len(); // lpd and hpd have the same length

        N = 2 * wt->length[J] - 1;
        auto cA_up = std::make_unique<double[]>(N);
        auto X_lp = std::make_unique<double[]>((N + lf - 1));
        auto X_hp = std::make_unique<double[]>((N + lf - 1));

        for (auto i = 0; i < app_len; ++i) {
            out[i] = wt->output[i];
        }

        iter = app_len;

        for (auto i = 0; i < J; ++i) {
            det_len = wt->length[i + 1];
            upsamp(out.get(), det_len, U, cA_up.get());
            N2 = 2 * wt->length[i + 1] - 1;

            if (wt->wave->lpr_len() == wt->wave->hpr_len() && ((wt->cmethod == "fft"sv) || (wt->cmethod == "FFT"sv))) {
                wt->cobj = conv_init(N2, lf);
                wt->cfftset = 1;
            } else if (!(wt->wave->lpr_len() == wt->wave->hpr_len())) {
                printf("Decomposition Filters must have the same length.");
                exit(-1);
            }

            wconv(wt, cA_up.get(), N2, wt->wave->lpr(), lf, X_lp.get());
            upsamp(wt->output + iter, det_len, U, cA_up.get());
            wconv(wt, cA_up.get(), N2, wt->wave->hpr(), lf, X_hp.get());

            for (k = lf - 2; k < N2 + 1; ++k) {
                out[k - lf + 2] = X_lp[k] + X_hp[k];
            }
            iter += det_len;
            if (wt->wave->lpr_len() == wt->wave->hpr_len() && ((wt->cmethod == "fft"sv) || (wt->cmethod == "FFT"sv))) {

                wt->cfftset = 0;
            }
        }

    } else {
        printf("Signal extension can be either per or sym");
        exit(-1);
    }

    for (auto i = 0; i < wt->siglength; ++i) {
        dwtop[i] = out[i];
    }
}

static void idwpt_per(wpt_set* wt, double const* cA, int len_cA, double const* cD, double* X)
{
    int len_avg;
    int l;
    int m;
    int n;
    int t;
    int l2;

    len_avg = (wt->wave->lpr_len() + wt->wave->hpr_len()) / 2;
    l2 = len_avg / 2;
    m = -2;
    n = -1;

    for (auto i = 0; i < len_cA + l2 - 1; ++i) {
        m += 2;
        n += 2;
        X[m] = 0.0;
        X[n] = 0.0;
        for (l = 0; l < l2; ++l) {
            t = 2 * l;
            if ((i - l) >= 0 && (i - l) < len_cA) {
                X[m] += wt->wave->lpr()[t] * cA[i - l] + wt->wave->hpr()[t] * cD[i - l];
                X[n] += wt->wave->lpr()[t + 1] * cA[i - l] + wt->wave->hpr()[t + 1] * cD[i - l];
            } else if ((i - l) >= len_cA && (i - l) < len_cA + len_avg - 1) {
                X[m] += wt->wave->lpr()[t] * cA[i - l - len_cA] + wt->wave->hpr()[t] * cD[i - l - len_cA];
                X[n] += wt->wave->lpr()[t + 1] * cA[i - l - len_cA] + wt->wave->hpr()[t + 1] * cD[i - l - len_cA];
            } else if ((i - l) < 0 && (i - l) > -l2) {
                X[m] += wt->wave->lpr()[t] * cA[len_cA + i - l] + wt->wave->hpr()[t] * cD[len_cA + i - l];
                X[n] += wt->wave->lpr()[t + 1] * cA[len_cA + i - l] + wt->wave->hpr()[t + 1] * cD[len_cA + i - l];
            }
        }
    }
}

static void idwpt_sym(wpt_set* wt, double const* cA, int len_cA, double const* cD, double* X)
{
    auto len_avg = (wt->wave->lpr_len() + wt->wave->hpr_len()) / 2;
    auto m = -2;
    auto n = -1;

    for (auto v = 0; v < len_cA; ++v) {
        auto i = v;
        m += 2;
        n += 2;
        X[m] = 0.0;
        X[n] = 0.0;
        for (auto l = 0; l < len_avg / 2; ++l) {
            auto const t = 2 * l;
            if ((i - l) >= 0 && (i - l) < len_cA) {
                X[m] += wt->wave->lpr()[t] * cA[i - l] + wt->wave->hpr()[t] * cD[i - l];
                X[n] += wt->wave->lpr()[t + 1] * cA[i - l] + wt->wave->hpr()[t + 1] * cD[i - l];
            }
        }
    }
}

void idwpt(wpt_set* wt, double* dwtop)
{
    int k;
    int l;
    int index;

    auto J = wt->J;
    auto app_len = wt->length[0];
    auto powJ = ipow2(J);
    auto lf = (wt->wave->lpr_len() + wt->wave->hpr_len()) / 2;
    auto xlen = powJ * (app_len + 2 * lf);

    auto X_lp = std::make_unique<double[]>(2 * (wt->length[J] + lf));
    auto X = std::make_unique<double[]>(xlen);
    auto out = std::make_unique<double[]>(wt->length[J]);
    auto out2 = std::make_unique<double[]>(wt->length[J]);
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
        for (auto i = 0; i < J; ++i) {
            llb *= 2;
            n1 += llb;
        }

        for (auto i = 0; i < xlen; ++i) {
            X[i] = 0.0;
        }

        for (auto i = 0; i < llb; ++i) {
            prep[i] = (int)wt->basisvector[n1 - llb + i];
            ptemp[i] = 0;
        }

        if (wt->ext == "per"sv) {
            index = 0;
            for (auto i = 0; i < J; ++i) {
                auto p = ipow2(J - i - 1);
                auto det_len = wt->length[i + 1];
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
                        for (k = 0; k < det_len; ++k) {
                            out[k] = wt->output[index + k];
                            out2[k] = wt->output[index + det_len + k];
                        }
                        idwpt_per(wt, out.get(), det_len, out2.get(), X_lp.get());
                        for (k = lf / 2 - 1; k < 2 * det_len + lf / 2 - 1; ++k) {
                            X[index3 + k - lf / 2 + 1] = X_lp[k];
                        }
                        index += 2 * det_len;
                        index3 += index2;
                        index4 += 2 * indexp;
                        ptemp[l] = 2;
                    } else if (prep[2 * l] == 1 && prep[2 * l + 1] == 2) {
                        index4 += indexp;
                        for (k = 0; k < det_len; ++k) {
                            out[k] = wt->output[index + k];
                            out2[k] = X[index4 + k];
                        }
                        idwpt_per(wt, out.get(), det_len, out2.get(), X_lp.get());
                        for (k = lf / 2 - 1; k < 2 * det_len + lf / 2 - 1; ++k) {
                            X[index3 + k - lf / 2 + 1] = X_lp[k];
                        }
                        index += det_len;
                        index3 += index2;
                        index4 += indexp;
                        ptemp[l] = 2;
                    } else if (prep[2 * l] == 2 && prep[2 * l + 1] == 1) {
                        for (k = 0; k < det_len; ++k) {
                            out[k] = X[index4 + k];
                            out2[k] = wt->output[index + k];
                        }
                        idwpt_per(wt, out.get(), det_len, out2.get(), X_lp.get());
                        for (k = lf / 2 - 1; k < 2 * det_len + lf / 2 - 1; ++k) {
                            X[index3 + k - lf / 2 + 1] = X_lp[k];
                        }
                        index += det_len;
                        index3 += index2;
                        index4 += 2 * indexp;
                        ptemp[l] = 2;
                    } else if (prep[2 * l] == 2 && prep[2 * l + 1] == 2) {
                        for (k = 0; k < det_len; ++k) {
                            out[k] = X[index4 + k];
                            out2[k] = X[index4 + indexp + k];
                        }
                        idwpt_per(wt, out.get(), det_len, out2.get(), X_lp.get());
                        for (k = lf / 2 - 1; k < 2 * det_len + lf / 2 - 1; ++k) {
                            X[index3 + k - lf / 2 + 1] = X_lp[k];
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

            for (auto i = 0; i < J; ++i) {
                auto p = ipow2(J - i - 1);
                auto det_len = wt->length[i + 1];
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
                        for (k = 0; k < det_len; ++k) {
                            out[k] = wt->output[index + k];
                            out2[k] = wt->output[index + det_len + k];
                        }
                        idwpt_sym(wt, out.get(), det_len, out2.get(), X_lp.get());
                        for (k = lf - 2; k < 2 * det_len; ++k) {
                            X[index3 + k - lf + 2] = X_lp[k];
                        }
                        index += 2 * det_len;
                        index3 += index2;
                        index4 += 2 * indexp;
                        ptemp[l] = 2;
                    } else if (prep[2 * l] == 1 && prep[2 * l + 1] == 2) {
                        index4 += indexp;
                        for (k = 0; k < det_len; ++k) {
                            out[k] = wt->output[index + k];
                            out2[k] = X[index4 + k];
                        }
                        idwpt_sym(wt, out.get(), det_len, out2.get(), X_lp.get());
                        for (k = lf - 2; k < 2 * det_len; ++k) {
                            X[index3 + k - lf + 2] = X_lp[k];
                        }
                        index += det_len;
                        index3 += index2;
                        index4 += indexp;
                        ptemp[l] = 2;
                    } else if (prep[2 * l] == 2 && prep[2 * l + 1] == 1) {
                        for (k = 0; k < det_len; ++k) {
                            out[k] = X[index4 + k];
                            out2[k] = wt->output[index + k];
                        }
                        idwpt_sym(wt, out.get(), det_len, out2.get(), X_lp.get());
                        for (k = lf - 2; k < 2 * det_len; ++k) {
                            X[index3 + k - lf + 2] = X_lp[k];
                        }
                        index += det_len;
                        index3 += index2;
                        index4 += 2 * indexp;
                        ptemp[l] = 2;
                    } else if (prep[2 * l] == 2 && prep[2 * l + 1] == 2) {
                        for (k = 0; k < det_len; ++k) {
                            out[k] = X[index4 + k];
                            out2[k] = X[index4 + indexp + k];
                        }
                        idwpt_sym(wt, out.get(), det_len, out2.get(), X_lp.get());
                        for (k = lf - 2; k < 2 * det_len; ++k) {
                            X[index3 + k - lf + 2] = X_lp[k];
                        }
                        index4 += 2 * indexp;
                        index3 += index2;
                        ptemp[l] = 2;
                    } else {
                        index3 += index2;
                        index4 += 2 * indexp;
                    }
                }

                //idwt1(wt, temp, cA_up, out, det_len, wt->output + iter, det_len, X_lp.get(), X_hp, out);
                /*
				idwpt_sym(wt, out, det_len, wt->output + iter, det_len, X_lp);
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
            dwtop[i] = X[i];
        }
    }
}

static void swt_per(wt_set* wt, int M, double* inp, int N, double* cA, int len_cA, double* cD)
{

    swt_per_stride(M, inp, N, wt->wave->lpd(), wt->wave->hpd(), wt->wave->lpd_len(), cA, len_cA, cD, 1, 1);
}

static void swt_fft(wt_set* wt, double const* inp)
{
    int N { 0 };

    auto temp_len = wt->siglength;
    auto J = wt->J;
    wt->length[0] = wt->length[J] = temp_len;
    wt->outlength = wt->length[J + 1] = (J + 1) * temp_len;
    auto M = 1;
    for (auto iter = 1; iter < J; ++iter) {
        M = 2 * M;
        wt->length[iter] = temp_len;
    }

    auto const len_filt = wt->wave->size();

    auto low_pass = std::make_unique<double[]>(M * len_filt);
    auto high_pass = std::make_unique<double[]>(M * len_filt);
    auto sig = std::make_unique<double[]>((M * len_filt + temp_len + (temp_len % 2)));
    auto cA = std::make_unique<double[]>((2 * M * len_filt + temp_len + (temp_len % 2)) - 1);
    auto cD = std::make_unique<double[]>((2 * M * len_filt + temp_len + (temp_len % 2)) - 1);

    M = 1;

    for (auto i = 0; i < temp_len; ++i) {
        wt->params[i] = inp[i];
    }

    auto lenacc = wt->outlength;

    for (auto iter = 0; iter < J; ++iter) {
        lenacc -= temp_len;
        if (iter > 0) {
            M = 2 * M;
            N = M * len_filt;
            upsamp2(wt->wave->lpd(), wt->wave->lpd_len(), M, low_pass.get());
            upsamp2(wt->wave->hpd(), wt->wave->hpd_len(), M, high_pass.get());
        } else {
            N = len_filt;
            for (auto i = 0; i < N; ++i) {
                low_pass[i] = wt->wave->lpd()[i];
                high_pass[i] = wt->wave->hpd()[i];
            }
        }

        //swt_per(wt,M, wt->params.get(), temp_len, cA, temp_len, cD,temp_len);

        per_ext(wt->params.get(), temp_len, N / 2, sig.get());

        if (wt->wave->lpd_len() == wt->wave->hpd_len() && ((wt->cmethod == "fft"sv) || (wt->cmethod == "FFT"sv))) {
            wt->cobj = conv_init(N + temp_len + (temp_len % 2), N);
            wt->cfftset = 1;
        } else if (!(wt->wave->lpd_len() == wt->wave->hpd_len())) {
            printf("Decomposition Filters must have the same length.");
            exit(-1);
        }

        wconv(wt, sig.get(), N + temp_len + (temp_len % 2), low_pass.get(), N, cA.get());

        wconv(wt, sig.get(), N + temp_len + (temp_len % 2), high_pass.get(), N, cD.get());

        if (wt->wave->lpd_len() == wt->wave->hpd_len() && ((wt->cmethod == "fft"sv) || (wt->cmethod == "FFT"sv))) {

            wt->cfftset = 0;
        }

        for (auto i = 0; i < temp_len; ++i) {
            wt->params[i] = cA[N + i];
            wt->params[lenacc + i] = cD[N + i];
        }
    }
}

static void swt_direct(wt_set* wt, double const* inp)
{
    int J;
    int temp_len;
    int iter;
    int M;
    int lenacc;

    temp_len = wt->siglength;
    J = wt->J;
    wt->length[0] = wt->length[J] = temp_len;
    wt->outlength = wt->length[J + 1] = (J + 1) * temp_len;
    M = 1;
    for (iter = 1; iter < J; ++iter) {
        M = 2 * M;
        wt->length[iter] = temp_len;
    }

    auto cA = std::make_unique<double[]>(temp_len);
    auto cD = std::make_unique<double[]>(temp_len);

    M = 1;

    for (auto i = 0; i < temp_len; ++i) {
        wt->params[i] = inp[i];
    }

    lenacc = wt->outlength;

    for (iter = 0; iter < J; ++iter) {
        lenacc -= temp_len;
        if (iter > 0) {
            M = 2 * M;
        }

        swt_per(wt, M, wt->params.get(), temp_len, cA.get(), temp_len, cD.get());

        for (auto i = 0; i < temp_len; ++i) {
            wt->params[i] = cA[i];
            wt->params[lenacc + i] = cD[i];
        }
    }
}

void swt(wt_set* wt, double const* inp)
{
    if ((wt->method == "swt"sv) && (wt->cmethod == "direct"sv)) {
        swt_direct(wt, inp);
    } else if ((wt->method == "swt"sv) && ((wt->cmethod == "fft"sv) || (wt->cmethod == "FFT"sv))) {
        swt_fft(wt, inp);
    } else {
        printf("SWT Only accepts two methods - direct and fft");
        exit(-1);
    }
}

void iswt(wt_set* wt, double* swtop)
{
    int N;
    int lf;
    int iter;
    int J;
    int index;
    int value;
    int count;
    int len;
    int index_shift;
    int len0;
    int U;
    int N1;
    int index2;

    N = wt->siglength;
    J = wt->J;
    U = 2;
    lf = wt->wave->lpr_len();

    auto appx_sig = std::make_unique<double[]>(N);
    auto det_sig = std::make_unique<double[]>(N);
    auto appx1 = std::make_unique<double[]>(N);
    auto det1 = std::make_unique<double[]>(N);
    auto appx2 = std::make_unique<double[]>(N);
    auto det2 = std::make_unique<double[]>(N);
    auto tempx = std::make_unique<double[]>(N);
    auto cL0 = std::make_unique<double[]>((N + (N % 2) + lf));
    auto cH0 = std::make_unique<double[]>((N + (N % 2) + lf));
    auto oup00L = std::make_unique<double[]>((N + 2 * lf));
    auto oup00H = std::make_unique<double[]>((N + 2 * lf));
    auto oup00 = std::make_unique<double[]>(N);
    auto oup01 = std::make_unique<double[]>(N);

    for (iter = 0; iter < J; ++iter) {
        for (auto i = 0; i < N; ++i) {
            swtop[i] = 0.0;
        }
        if (iter == 0) {
            for (auto i = 0; i < N; ++i) {
                appx_sig[i] = wt->output[i];
                det_sig[i] = wt->output[N + i];
            }
        } else {
            for (auto i = 0; i < N; ++i) {
                det_sig[i] = wt->output[(iter + 1) * N + i];
            }
        }

        value = (int)std::pow(2.0, (double)(J - 1 - iter));

        for (count = 0; count < value; count++) {
            len = 0;
            for (index = count; index < N; index += value) {
                appx1[len] = appx_sig[index];
                det1[len] = det_sig[index];
                len++;
            }

            //SHIFT 0
            len0 = 0;

            for (index_shift = 0; index_shift < len; index_shift += 2) {
                appx2[len0] = appx1[index_shift];
                det2[len0] = det1[index_shift];
                len0++;
            }
            upsamp2(appx2.get(), len0, U, tempx.get());
            per_ext(tempx.get(), 2 * len0, lf / 2, cL0.get());

            upsamp2(det2.get(), len0, U, tempx.get());
            per_ext(tempx.get(), 2 * len0, lf / 2, cH0.get());

            N1 = 2 * len0 + lf;

            if (wt->wave->lpr_len() == wt->wave->hpr_len() && ((wt->cmethod == "fft"sv) || (wt->cmethod == "FFT"sv))) {
                wt->cobj = conv_init(N1, lf);
                wt->cfftset = 1;
            } else if (!(wt->wave->lpd_len() == wt->wave->hpd_len())) {
                printf("Decomposition Filters must have the same length.");
                exit(-1);
            }

            wconv(wt, cL0.get(), N1, wt->wave->lpr(), lf, oup00L.get());

            wconv(wt, cH0.get(), N1, wt->wave->hpr(), lf, oup00H.get());

            for (auto i = lf - 1; i < 2 * len0 + lf - 1; ++i) {
                oup00[i - lf + 1] = oup00L[i] + oup00H[i];
            }

            //SHIFT 1

            len0 = 0;

            for (index_shift = 1; index_shift < len; index_shift += 2) {
                appx2[len0] = appx1[index_shift];
                det2[len0] = det1[index_shift];
                len0++;
            }

            upsamp2(appx2.get(), len0, U, tempx.get());
            per_ext(tempx.get(), 2 * len0, lf / 2, cL0.get());

            upsamp2(det2.get(), len0, U, tempx.get());
            per_ext(tempx.get(), 2 * len0, lf / 2, cH0.get());

            N1 = 2 * len0 + lf;

            wconv(wt, cL0.get(), N1, wt->wave->lpr(), lf, oup00L.get());
            wconv(wt, cH0.get(), N1, wt->wave->hpr(), lf, oup00H.get());

            for (auto i = lf - 1; i < 2 * len0 + lf - 1; ++i) {
                oup01[i - lf + 1] = oup00L[i] + oup00H[i];
            }

            circshift(oup01.get(), 2 * len0, -1);

            index2 = 0;

            for (index = count; index < N; index += value) {
                swtop[index] = (oup00[index2] + oup01[index2]) / 2.0;
                index2++;
            }
        }
        for (auto i = 0; i < N; ++i) {
            appx_sig[i] = swtop[i];
        }
    }
}

static void modwt_per(wt_set* wt, int M, double const* inp, double* cA, int len_cA, double* cD)
{
    auto const len_avg = wt->wave->lpd_len();
    auto filt = std::make_unique<double[]>(2 * len_avg);
    auto s = std::sqrt(2.0);

    for (auto i = 0; i < len_avg; ++i) {
        filt[i] = wt->wave->lpd()[i] / s;
        filt[len_avg + i] = wt->wave->hpd()[i] / s;
    }

    for (auto i = 0; i < len_cA; ++i) {
        auto t = i;
        cA[i] = filt[0] * inp[t];
        cD[i] = filt[len_avg] * inp[t];
        for (auto l = 1; l < len_avg; l++) {
            t -= M;
            while (t >= len_cA) {
                t -= len_cA;
            }
            while (t < 0) {
                t += len_cA;
            }

            cA[i] += filt[l] * inp[t];
            cD[i] += filt[len_avg + l] * inp[t];
        }
    }
}

static void modwt_direct(wt_set* wt, double const* inp)
{
    if (wt->ext != "per"sv) {
        printf("MODWT direct method only uses periodic extension per. \n");
        printf(" Use MODWT fft method for symmetric extension sym \n");
        exit(-1);
    }

    auto temp_len = wt->siglength;
    auto J = wt->J;
    wt->length[0] = wt->length[J] = temp_len;
    wt->outlength = wt->length[J + 1] = (J + 1) * temp_len;
    auto M = 1;
    for (auto iter = 1; iter < J; ++iter) {
        M = 2 * M;
        wt->length[iter] = temp_len;
    }

    auto cA = std::make_unique<double[]>(temp_len);
    auto cD = std::make_unique<double[]>(temp_len);

    M = 1;

    for (auto i = 0; i < temp_len; ++i) {
        wt->params[i] = inp[i];
    }

    auto lenacc = wt->outlength;

    for (auto iter = 0; iter < J; ++iter) {
        lenacc -= temp_len;
        if (iter > 0) {
            M = 2 * M;
        }

        modwt_per(wt, M, wt->params.get(), cA.get(), temp_len, cD.get());

        for (auto i = 0; i < temp_len; ++i) {
            wt->params[i] = cA[i];
            wt->params[lenacc + i] = cD[i];
        }
    }
}

static void modwt_fft(wt_set* wt, double const* inp)
{
    int J;
    int iter;
    int M;
    int lenacc;
    double s;
    double tmp1;
    double tmp2;

    auto temp_len = wt->siglength;
    auto len_avg = wt->wave->lpd_len();
    int N { 0 };
    if (wt->ext == "sym"sv) {
        N = 2 * temp_len;
    } else if (wt->ext == "per"sv) {
        N = temp_len;
    }
    J = wt->J;
    wt->modwtsiglength = N;
    wt->length[0] = wt->length[J] = N;
    wt->outlength = wt->length[J + 1] = (J + 1) * N;

    s = std::sqrt(2.0);
    for (iter = 1; iter < J; ++iter) {
        wt->length[iter] = N;
    }

    auto fft_fd = fft_init(N, 1);
    auto fft_bd = fft_init(N, -1);

    auto sig = std::make_unique<fft_data[]>(N);
    auto cA = std::make_unique<fft_data[]>(N);
    auto cD = std::make_unique<fft_data[]>(N);
    auto low_pass = std::make_unique<fft_data[]>(N);
    auto high_pass = std::make_unique<fft_data[]>(N);
    auto index = std::make_unique<int[]>(N);

    // N-point FFT of low pass and high pass filters

    // Low Pass Filter

    for (auto i = 0; i < len_avg; ++i) {
        sig[i].re = (fft_type)wt->wave->lpd()[i] / s;
        sig[i].im = 0.0;
    }
    for (auto i = len_avg; i < N; ++i) {
        sig[i].re = 0.0;
        sig[i].im = 0.0;
    }

    fft_exec(*fft_fd, sig.get(), low_pass.get());

    // High Pass Filter

    for (auto i = 0; i < len_avg; ++i) {
        sig[i].re = (fft_type)wt->wave->hpd()[i] / s;
        sig[i].im = 0.0;
    }
    for (auto i = len_avg; i < N; ++i) {
        sig[i].re = 0.0;
        sig[i].im = 0.0;
    }

    fft_exec(*fft_fd, sig.get(), high_pass.get());

    // symmetric extension
    for (auto i = 0; i < temp_len; ++i) {
        sig[i].re = (fft_type)inp[i];
        sig[i].im = 0.0;
    }
    for (auto i = temp_len; i < N; ++i) {
        sig[i].re = (fft_type)inp[N - i - 1];
        sig[i].im = 0.0;
    }

    // FFT of data

    fft_exec(*fft_fd, sig.get(), cA.get());

    lenacc = wt->outlength;

    M = 1;

    for (iter = 0; iter < J; ++iter) {
        lenacc -= N;

        for (auto i = 0; i < N; ++i) {
            index[i] = (M * i) % N;
        }

        for (auto i = 0; i < N; ++i) {
            tmp1 = cA[i].re;
            tmp2 = cA[i].im;
            cA[i].re = low_pass[index[i]].re * tmp1 - low_pass[index[i]].im * tmp2;
            cA[i].im = low_pass[index[i]].re * tmp2 + low_pass[index[i]].im * tmp1;

            cD[i].re = high_pass[index[i]].re * tmp1 - high_pass[index[i]].im * tmp2;
            cD[i].im = high_pass[index[i]].re * tmp2 + high_pass[index[i]].im * tmp1;
        }

        fft_exec(*fft_bd, cD.get(), sig.get());

        for (auto i = 0; i < N; ++i) {
            wt->params[lenacc + i] = sig[i].re / N;
        }

        M *= 2;
    }

    fft_exec(*fft_bd, cA.get(), sig.get());

    for (auto i = 0; i < N; ++i) {
        wt->params[i] = sig[i].re / N;
    }
}

void modwt(wt_set* wt, double const* inp)
{
    if (wt->cmethod == "direct"sv) {
        modwt_direct(wt, inp);
    } else if (wt->cmethod == "fft"sv) {
        modwt_fft(wt, inp);
    } else {
        printf("Error- Available Choices for this method are - direct and fft \n");
        exit(-1);
    }
}

static void conj_complex(fft_data* x, int N)
{
    for (auto i = 0; i < N; ++i) {
        x[i].im *= (-1.0);
    }
}

void imodwt_fft(wt_set* wt, double* oup)
{
    auto N = wt->modwtsiglength;
    auto len_avg = wt->wave->lpd_len();
    auto J = wt->J;

    auto s = std::sqrt(2.0);
    auto fft_fd = fft_init(N, 1);
    auto fft_bd = fft_init(N, -1);

    auto sig = std::make_unique<fft_data[]>(N);
    auto cA = std::make_unique<fft_data[]>(N);
    auto cD = std::make_unique<fft_data[]>(N);
    auto low_pass = std::make_unique<fft_data[]>(N);
    auto high_pass = std::make_unique<fft_data[]>(N);
    auto index = std::make_unique<int[]>(N);

    // N-point FFT of low pass and high pass filters

    // Low Pass Filter

    for (auto i = 0; i < len_avg; ++i) {
        sig[i].re = (fft_type)wt->wave->lpd()[i] / s;
        sig[i].im = 0.0;
    }
    for (auto i = len_avg; i < N; ++i) {
        sig[i].re = 0.0;
        sig[i].im = 0.0;
    }

    fft_exec(*fft_fd, sig.get(), low_pass.get());

    // High Pass Filter

    for (auto i = 0; i < len_avg; ++i) {
        sig[i].re = (fft_type)wt->wave->hpd()[i] / s;
        sig[i].im = 0.0;
    }
    for (auto i = len_avg; i < N; ++i) {
        sig[i].re = 0.0;
        sig[i].im = 0.0;
    }

    fft_exec(*fft_fd, sig.get(), high_pass.get());

    // Complex conjugate of the two filters

    conj_complex(low_pass.get(), N);
    conj_complex(high_pass.get(), N);

    auto M = (int)std::pow(2.0, (double)J - 1.0);
    auto lenacc = N;

    //
    for (auto i = 0; i < N; ++i) {
        sig[i].re = (fft_type)wt->output[i];
        sig[i].im = 0.0;
    }

    for (auto iter = 0; iter < J; ++iter) {
        fft_exec(*fft_fd, sig.get(), cA.get());
        for (auto i = 0; i < N; ++i) {
            sig[i].re = wt->output[lenacc + i];
            sig[i].im = 0.0;
        }
        fft_exec(*fft_fd, sig.get(), cD.get());

        for (auto i = 0; i < N; ++i) {
            index[i] = (M * i) % N;
        }

        for (auto i = 0; i < N; ++i) {
            auto const tmp1 = cA[i].re;
            auto const tmp2 = cA[i].im;
            cA[i].re = low_pass[index[i]].re * tmp1 - low_pass[index[i]].im * tmp2 + high_pass[index[i]].re * cD[i].re - high_pass[index[i]].im * cD[i].im;
            cA[i].im = low_pass[index[i]].re * tmp2 + low_pass[index[i]].im * tmp1 + high_pass[index[i]].re * cD[i].im + high_pass[index[i]].im * cD[i].re;
        }

        fft_exec(*fft_bd, cA.get(), sig.get());

        for (auto i = 0; i < N; ++i) {
            sig[i].re /= N;
            sig[i].im /= N;
        }
        M /= 2;
        lenacc += N;
    }

    for (auto i = 0; i < wt->siglength; ++i) {
        oup[i] = sig[i].re;
    }
}

static void imodwt_per(wt_set* wt, int M, double const* cA, int len_cA, double const* cD, double* X)
{
    auto const len_avg = wt->wave->lpd_len();
    auto filt = std::make_unique<double[]>(2 * len_avg);
    auto s = std::sqrt(2.0);

    for (auto i = 0; i < len_avg; ++i) {
        filt[i] = wt->wave->lpd()[i] / s;
        filt[len_avg + i] = wt->wave->hpd()[i] / s;
    }

    for (auto i = 0; i < len_cA; ++i) {
        auto t = i;
        X[i] = (filt[0] * cA[t]) + (filt[len_avg] * cD[t]);
        for (auto l = 1; l < len_avg; l++) {
            t += M;
            while (t >= len_cA) {
                t -= len_cA;
            }
            while (t < 0) {
                t += len_cA;
            }

            X[i] += (filt[l] * cA[t]) + (filt[len_avg + l] * cD[t]);
        }
    }
}

static void imodwt_direct(wt_set* wt, double* dwtop)
{
    auto N = wt->siglength;
    auto lenacc = N;

    auto J = wt->J;
    auto M = (int)std::pow(2.0, (double)J - 1.0);

    auto X = std::make_unique<double[]>(N);

    for (auto i = 0; i < N; ++i) {
        dwtop[i] = wt->output[i];
    }

    for (auto iter = 0; iter < J; ++iter) {
        if (iter > 0) {
            M = M / 2;
        }
        imodwt_per(wt, M, dwtop, N, wt->params.get() + lenacc, X.get());
        /*
		for (auto j = lf - 1; j < N; ++j) {
			dwtop[j - lf + 1] = X[j];
		}
		for (auto j = 0; j < lf - 1; ++j) {
			dwtop[N - lf + 1 + j] = X[j];
		}
		*/
        for (auto j = 0; j < N; ++j) {
            dwtop[j] = X[j];
        }

        lenacc += N;
    }
}

void imodwt(wt_set* wt, double* oup)
{
    if (wt->cmethod == "direct"sv) {
        imodwt_direct(wt, oup);
    } else if (wt->cmethod == "fft"sv) {
        imodwt_fft(wt, oup);
    } else {
        printf("Error- Available Choices for this method are - direct and fft \n");
        exit(-1);
    }
}

void setDWTExtension(wt_set* wt, char const* extension)
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

void setWTREEExtension(wtree_set* wt, char const* extension)
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

void setDWPTExtension(wpt_set* wt, char const* extension)
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

void setDWT2Extension(wt2_set* wt, char const* extension)
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

void setDWPTEntropy(wpt_set* wt, char const* entropy, double eparam)
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

void setWTConv(wt_set* wt, char const* cmethod)
{
    if ((cmethod == "fft"sv) || (cmethod == "FFT"sv)) {
        wt->cmethod = "fft";
    } else if (cmethod == "direct"sv) {
        wt->cmethod = "direct";
    } else {
        printf("Convolution Only accepts two methods - direct and fft");
        exit(-1);
    }
}

auto dwt2(wt2_set* wt, double* inp) -> std::unique_ptr<double[]>
{
    int iter;
    int N;
    int rows_i;
    int cols_i;
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

    auto J = wt->J;
    wt->outlength = 0;

    auto rows_n = wt->rows;
    auto cols_n = wt->cols;
    auto lp = wt->wave->lpd_len();
    auto clen = J * 3;

    if (wt->ext == "per"sv) {
        auto idx = 2 * J;
        while (idx > 0) {
            rows_n = (int)ceil((double)rows_n / 2.0);
            cols_n = (int)ceil((double)cols_n / 2.0);
            wt->dimensions[idx - 1] = cols_n;
            wt->dimensions[idx - 2] = rows_n;
            wt->outlength += (rows_n * cols_n) * 3;
            idx = idx - 2;
        }
        wt->outlength += (rows_n * cols_n);
        N = wt->outlength;
        auto wavecoeff = makeZeros<double>(wt->outlength);

        orig = inp;
        ir = wt->rows;
        ic = wt->cols;
        cols_i = wt->dimensions[2 * J - 1];

        auto lp_dn1 = makeZeros<double>(ir * cols_i);
        auto hp_dn1 = makeZeros<double>(ir * cols_i);

        for (iter = 0; iter < J; ++iter) {
            rows_i = wt->dimensions[2 * J - 2 * iter - 2];
            cols_i = wt->dimensions[2 * J - 2 * iter - 1];
            istride = 1;
            ostride = 1;
            cdim = rows_i * cols_i;
            // Row filtering and column subsampling
            for (auto i = 0; i < ir; ++i) {
                dwt_per_stride(orig + i * ic, ic, wt->wave->lpd(), wt->wave->hpd(), lp, lp_dn1.get() + i * cols_i, cols_i, hp_dn1.get() + i * cols_i, istride, ostride);
            }

            // Column Filtering and Row subsampling
            aHH = N - cdim;
            wt->coeffaccess[clen] = aHH;
            aHL = aHH - cdim;
            wt->coeffaccess[clen - 1] = aHL;
            aLH = aHL - cdim;
            wt->coeffaccess[clen - 2] = aLH;
            aLL = aLH - cdim;

            N -= 3 * cdim;
            ic = cols_i;
            istride = ic;
            ostride = ic;

            for (auto i = 0; i < ic; ++i) {
                dwt_per_stride(lp_dn1.get() + i, ir, wt->wave->lpd(), wt->wave->hpd(), lp, wavecoeff.get() + aLL + i, rows_i, wavecoeff.get() + aLH + i, istride, ostride);
            }

            for (auto i = 0; i < ic; ++i) {
                dwt_per_stride(hp_dn1.get() + i, ir, wt->wave->lpd(), wt->wave->hpd(), lp, wavecoeff.get() + aHL + i, rows_i, wavecoeff.get() + aHH + i, istride, ostride);
            }

            ir = rows_i;
            orig = wavecoeff.get() + aLL;
            clen -= 3;
        }
        wt->coeffaccess[0] = 0;

        return wavecoeff;
    }

    assert(wt->ext == "sym"sv);

    auto idx = 2 * J;
    while (idx > 0) {
        rows_n += lp - 2;
        cols_n += lp - 2;
        rows_n = (int)ceil((double)rows_n / 2.0);
        cols_n = (int)ceil((double)cols_n / 2.0);
        wt->dimensions[idx - 1] = cols_n;
        wt->dimensions[idx - 2] = rows_n;
        wt->outlength += (rows_n * cols_n) * 3;
        idx = idx - 2;
    }
    wt->outlength += (rows_n * cols_n);
    N = wt->outlength;
    auto wavecoeff = makeZeros<double>(wt->outlength);

    orig = inp;
    ir = wt->rows;
    ic = wt->cols;
    cols_i = wt->dimensions[2 * J - 1];

    auto lp_dn1 = makeZeros<double>(ir * cols_i);
    auto hp_dn1 = makeZeros<double>(ir * cols_i);

    for (iter = 0; iter < J; ++iter) {
        rows_i = wt->dimensions[2 * J - 2 * iter - 2];
        cols_i = wt->dimensions[2 * J - 2 * iter - 1];
        istride = 1;
        ostride = 1;
        cdim = rows_i * cols_i;
        // Row filtering and column subsampling
        for (auto i = 0; i < ir; ++i) {
            dwt_sym_stride(orig + i * ic, ic, wt->wave->lpd(), wt->wave->hpd(), lp, lp_dn1.get() + i * cols_i, cols_i, hp_dn1.get() + i * cols_i, istride, ostride);
        }

        // Column Filtering and Row subsampling
        aHH = N - cdim;
        wt->coeffaccess[clen] = aHH;
        aHL = aHH - cdim;
        wt->coeffaccess[clen - 1] = aHL;
        aLH = aHL - cdim;
        wt->coeffaccess[clen - 2] = aLH;
        aLL = aLH - cdim;
        N -= 3 * cdim;
        ic = cols_i;
        istride = ic;
        ostride = ic;

        for (auto i = 0; i < ic; ++i) {
            dwt_sym_stride(lp_dn1.get() + i, ir, wt->wave->lpd(), wt->wave->hpd(), lp, wavecoeff.get() + aLL + i, rows_i, wavecoeff.get() + aLH + i, istride, ostride);
        }

        for (auto i = 0; i < ic; ++i) {
            dwt_sym_stride(hp_dn1.get() + i, ir, wt->wave->lpd(), wt->wave->hpd(), lp, wavecoeff.get() + aHL + i, rows_i, wavecoeff.get() + aHH + i, istride, ostride);
        }

        ir = rows_i;
        orig = wavecoeff.get() + aLL;
        clen -= 3;
    }

    wt->coeffaccess[0] = 0;

    return wavecoeff;
}

void idwt2(wt2_set* wt, double* wavecoeff, double* oup)
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
    auto const J = wt->J;

    if (wt->ext == "per"sv) {
        auto const N = rows > cols ? 2 * rows : 2 * cols;
        auto const lf = (wt->wave->lpr_len() + wt->wave->hpr_len()) / 2;

        auto idx = J;
        auto dim1 = wt->dimensions[0];
        auto dim2 = wt->dimensions[1];
        auto k = 0;
        while (idx > 0) {
            k += 1;
            dim1 *= 2;
            dim2 *= 2;
            idx--;
        }

        auto X_lp = makeZeros<double>(N + 2 * lf - 1);
        auto cL = makeZeros<double>(dim1 * dim2);
        auto cH = makeZeros<double>(dim1 * dim2);
        auto out = makeZeros<double>(dim1 * dim2);

        aLL = wt->coeffaccess[0];
        orig = wavecoeff + aLL;
        for (iter = 0; iter < J; ++iter) {
            ir = wt->dimensions[2 * iter];
            ic = wt->dimensions[2 * iter + 1];
            istride = ic;
            ostride = 1;
            aLH = wt->coeffaccess[iter * 3 + 1];
            aHL = wt->coeffaccess[iter * 3 + 2];
            aHH = wt->coeffaccess[iter * 3 + 3];
            for (auto i = 0; i < ic; ++i) {
                idwt_per_stride(orig + i, ir, wavecoeff + aLH + i, wt->wave->lpr(), wt->wave->hpr(), lf, X_lp.get(), istride, ostride);

                for (k = lf / 2 - 1; k < 2 * ir + lf / 2 - 1; ++k) {
                    cL[(k - lf / 2 + 1) * ic + i] = X_lp[k];
                }

                idwt_per_stride(wavecoeff + aHL + i, ir, wavecoeff + aHH + i, wt->wave->lpr(), wt->wave->hpr(), lf, X_lp.get(), istride, ostride);

                for (k = lf / 2 - 1; k < 2 * ir + lf / 2 - 1; ++k) {
                    cH[(k - lf / 2 + 1) * ic + i] = X_lp[k];
                }
            }

            ir *= 2;
            istride = 1;
            ostride = 1;

            for (auto i = 0; i < ir; ++i) {
                idwt_per_stride(cL.get() + i * ic, ic, cH.get() + i * ic, wt->wave->lpr(), wt->wave->hpr(), lf, X_lp.get(), istride, ostride);

                for (k = lf / 2 - 1; k < 2 * ic + lf / 2 - 1; ++k) {
                    out[(k - lf / 2 + 1) + i * ic * 2] = X_lp[k];
                }
            }
            ic *= 2;
            if (iter == J - 1) {
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

    auto const N = rows > cols ? 2 * rows - 1 : 2 * cols - 1;
    auto const lf = (wt->wave->lpr_len() + wt->wave->hpr_len()) / 2;

    auto idx = J;
    auto dim1 = wt->dimensions[0];
    auto dim2 = wt->dimensions[1];
    auto k = 0;
    while (idx > 0) {
        k += 1;
        dim1 *= 2;
        dim2 *= 2;
        idx--;
    }

    auto X_lp = makeZeros<double>(N + 2 * lf - 1);
    auto cL = makeZeros<double>(dim1 * dim2);
    auto cH = makeZeros<double>(dim1 * dim2);
    auto out = makeZeros<double>(dim1 * dim2);

    aLL = wt->coeffaccess[0];
    orig = wavecoeff + aLL;
    for (iter = 0; iter < J; ++iter) {
        ir = wt->dimensions[2 * iter];
        ic = wt->dimensions[2 * iter + 1];
        istride = ic;
        ostride = 1;
        aLH = wt->coeffaccess[iter * 3 + 1];
        aHL = wt->coeffaccess[iter * 3 + 2];
        aHH = wt->coeffaccess[iter * 3 + 3];
        for (auto i = 0; i < ic; ++i) {
            idwt_sym_stride(orig + i, ir, wavecoeff + aLH + i, wt->wave->lpr(), wt->wave->hpr(), lf, X_lp.get(), istride, ostride);

            for (k = lf - 2; k < 2 * ir; ++k) {
                cL[(k - lf + 2) * ic + i] = X_lp[k];
            }

            idwt_sym_stride(wavecoeff + aHL + i, ir, wavecoeff + aHH + i, wt->wave->lpr(), wt->wave->hpr(), lf, X_lp.get(), istride, ostride);

            for (k = lf - 2; k < 2 * ir; ++k) {
                cH[(k - lf + 2) * ic + i] = X_lp[k];
            }
        }

        ir *= 2;
        istride = 1;
        ostride = 1;

        for (auto i = 0; i < ir; ++i) {
            idwt_sym_stride(cL.get() + i * ic, ic, cH.get() + i * ic, wt->wave->lpr(), wt->wave->hpr(), lf, X_lp.get(), istride, ostride);

            for (k = lf - 2; k < 2 * ic; ++k) {
                out[(k - lf + 2) + i * ic * 2] = X_lp[k];
            }
        }
        ic *= 2;
        if (iter == J - 1) {
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

auto swt2(wt2_set* wt, double* inp) -> std::unique_ptr<double[]>
{
    int J;
    int iter;
    int M;
    int N;
    int lp;
    int rows_n;
    int cols_n;
    int rows_i;
    int cols_i;
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

    J = wt->J;
    M = 1;
    wt->outlength = 0;

    rows_n = wt->rows;
    cols_n = wt->cols;
    lp = wt->wave->lpd_len();
    clen = J * 3;

    auto idx = 2 * J;
    while (idx > 0) {
        wt->dimensions[idx - 1] = cols_n;
        wt->dimensions[idx - 2] = rows_n;
        wt->outlength += (rows_n * cols_n) * 3;
        idx = idx - 2;
    }
    wt->outlength += (rows_n * cols_n);
    N = wt->outlength;
    auto wavecoeff = makeZeros<double>(wt->outlength);

    orig = inp;
    ir = wt->rows;
    ic = wt->cols;
    cols_i = wt->dimensions[2 * J - 1];

    auto lp_dn1 = std::make_unique<double[]>(ir * cols_i);
    auto hp_dn1 = std::make_unique<double[]>(ir * cols_i);

    for (iter = 0; iter < J; ++iter) {
        if (iter > 0) {
            M = 2 * M;
        }
        rows_i = wt->dimensions[2 * J - 2 * iter - 2];
        cols_i = wt->dimensions[2 * J - 2 * iter - 1];
        istride = 1;
        ostride = 1;
        cdim = rows_i * cols_i;
        // Row filtering and column subsampling
        for (auto i = 0; i < ir; ++i) {
            swt_per_stride(M, orig + i * ic, ic, wt->wave->lpd(), wt->wave->hpd(), lp, lp_dn1.get() + i * cols_i, cols_i, hp_dn1.get() + i * cols_i, istride, ostride);
        }
        // Column Filtering and Row subsampling
        aHH = N - cdim;
        wt->coeffaccess[clen] = aHH;
        aHL = aHH - cdim;
        wt->coeffaccess[clen - 1] = aHL;
        aLH = aHL - cdim;
        wt->coeffaccess[clen - 2] = aLH;
        aLL = aLH - cdim;

        N -= 3 * cdim;
        ic = cols_i;
        istride = ic;
        ostride = ic;
        for (auto i = 0; i < ic; ++i) {
            swt_per_stride(M, lp_dn1.get() + i, ir, wt->wave->lpd(), wt->wave->hpd(), lp, wavecoeff.get() + aLL + i, rows_i, wavecoeff.get() + aLH + i, istride, ostride);
        }

        for (auto i = 0; i < ic; ++i) {
            swt_per_stride(M, hp_dn1.get() + i, ir, wt->wave->lpd(), wt->wave->hpd(), lp, wavecoeff.get() + aHL + i, rows_i, wavecoeff.get() + aHH + i, istride, ostride);
        }

        ir = rows_i;
        orig = wavecoeff.get() + aLL;
        clen -= 3;
    }
    wt->coeffaccess[0] = 0;

    return wavecoeff;
}

void iswt2(wt2_set* wt, double const* wavecoeffs, double* oup)
{
    int k;
    int iter;
    int it2;
    int it3;
    int J;
    int M;
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
    J = wt->J;
    rows = wt->rows;
    cols = wt->cols;
    lf = wt->wave->lpd_len();

    auto A = makeZeros<double>((rows + lf) * (cols + lf));
    auto H = makeZeros<double>((rows + lf) * (cols + lf));
    auto V = makeZeros<double>((rows + lf) * (cols + lf));
    auto D = makeZeros<double>((rows + lf) * (cols + lf));
    auto oup1 = makeZeros<double>((rows + lf) * (cols + lf));
    auto oup2 = makeZeros<double>((rows + lf) * (cols + lf));

    aLL = wt->coeffaccess[0];

    for (auto i = 0; i < rows; ++i) {
        for (k = 0; k < cols; ++k) {
            oup[i * cols + k] = wavecoeffs[aLL + i * cols + k];
        }
    }

    for (iter = J; iter > 0; iter--) {
        aLH = wt->coeffaccess[(J - iter) * 3 + 1];
        aHL = wt->coeffaccess[(J - iter) * 3 + 2];
        aHH = wt->coeffaccess[(J - iter) * 3 + 3];
        M = (int)std::pow(2.0, (double)iter - 1);

        for (it2 = 0; it2 < M; ++it2) {
            ir = 0;
            ic = 0;
            it3 = 0;
            // oup1
            for (auto i = it2; i < rows; i += 2 * M) {
                ic = 0;
                for (k = it2; k < cols; k += 2 * M) {
                    A[it3] = oup[i * cols + k];
                    H[it3] = wavecoeffs[aLH + i * cols + k];
                    V[it3] = wavecoeffs[aHL + i * cols + k];
                    D[it3] = wavecoeffs[aHH + i * cols + k];
                    it3++;
                    ic++;
                }
                ir++;
            }
            shift = 0;
            idwt2_shift(shift, ir, ic, wt->wave->lpr(), wt->wave->hpr(), wt->wave->lpd_len(), A.get(), H.get(), V.get(), D.get(), oup1.get());
            //oup2
            ir = 0;
            ic = 0;
            it3 = 0;
            for (auto i = it2 + M; i < rows; i += 2 * M) {
                ic = 0;
                for (k = it2 + M; k < cols; k += 2 * M) {
                    A[it3] = oup[i * cols + k];
                    H[it3] = wavecoeffs[aLH + i * cols + k];
                    V[it3] = wavecoeffs[aHL + i * cols + k];
                    D[it3] = wavecoeffs[aHH + i * cols + k];
                    it3++;
                    ic++;
                }
                ir++;
            }
            shift = -1;
            idwt2_shift(shift, ir, ic, wt->wave->lpr(), wt->wave->hpr(), wt->wave->lpd_len(), A.get(), H.get(), V.get(), D.get(), oup2.get());
            // Shift oup1 and oup2. Then add them to get A.
            i1 = 0;
            for (auto i = it2; i < rows; i += M) {
                k1 = 0;
                for (k = it2; k < cols; k += M) {
                    oup[i * cols + k] = 0.5 * (oup1[i1 * 2 * ic + k1] + oup2[i1 * 2 * ic + k1]);
                    k1++;
                }
                i1++;
            }
        }
    }
}

auto modwt2(wt2_set* wt, double* inp) -> std::unique_ptr<double[]>
{
    int J;
    int iter;
    int M;
    int N;
    int lp;
    int rows_n;
    int cols_n;
    int rows_i;
    int cols_i;
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

    J = wt->J;
    M = 1;
    wt->outlength = 0;

    rows_n = wt->rows;
    cols_n = wt->cols;
    lp = wt->wave->lpd_len();
    clen = J * 3;

    auto idx = 2 * J;
    while (idx > 0) {
        wt->dimensions[idx - 1] = cols_n;
        wt->dimensions[idx - 2] = rows_n;
        wt->outlength += (rows_n * cols_n) * 3;
        idx = idx - 2;
    }
    wt->outlength += (rows_n * cols_n);
    N = wt->outlength;
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
    cols_i = wt->dimensions[2 * J - 1];

    auto lp_dn1 = std::make_unique<double[]>(ir * cols_i);
    auto hp_dn1 = std::make_unique<double[]>(ir * cols_i);

    for (iter = 0; iter < J; ++iter) {
        if (iter > 0) {
            M = 2 * M;
        }
        rows_i = wt->dimensions[2 * J - 2 * iter - 2];
        cols_i = wt->dimensions[2 * J - 2 * iter - 1];
        istride = 1;
        ostride = 1;
        cdim = rows_i * cols_i;
        // Row filtering and column subsampling
        for (auto i = 0; i < ir; ++i) {
            modwt_per_stride(M, orig + i * ic, ic, filt.get(), lp, lp_dn1.get() + i * cols_i, cols_i, hp_dn1.get() + i * cols_i, istride, ostride);
        }
        // Column Filtering and Row subsampling
        aHH = N - cdim;
        wt->coeffaccess[clen] = aHH;
        aHL = aHH - cdim;
        wt->coeffaccess[clen - 1] = aHL;
        aLH = aHL - cdim;
        wt->coeffaccess[clen - 2] = aLH;
        aLL = aLH - cdim;
        N -= 3 * cdim;
        ic = cols_i;
        istride = ic;
        ostride = ic;
        for (auto i = 0; i < ic; ++i) {
            modwt_per_stride(M, lp_dn1.get() + i, ir, filt.get(), lp, wavecoeff.get() + aLL + i, rows_i, wavecoeff.get() + aLH + i, istride, ostride);
        }

        for (auto i = 0; i < ic; ++i) {
            modwt_per_stride(M, hp_dn1.get() + i, ir, filt.get(), lp, wavecoeff.get() + aHL + i, rows_i, wavecoeff.get() + aHH + i, istride, ostride);
        }

        ir = rows_i;
        orig = wavecoeff.get() + aLL;
        clen -= 3;
    }
    wt->coeffaccess[0] = 0;

    return wavecoeff;
}

void imodwt2(wt2_set* wt, double* wavecoeff, double* oup)
{
    int rows;
    int cols;
    int M;
    // int N;
    int ir;
    int ic;
    int lf;
    int istride;
    int ostride;
    int iter;
    int J;
    int aLL;
    int aLH;
    int aHL;
    int aHH;
    double* orig;
    double s;

    rows = wt->rows;
    cols = wt->cols;
    J = wt->J;

    M = (int)std::pow(2.0, (double)J - 1.0);
    // N = rows > cols ? rows : cols;
    lf = (wt->wave->lpr_len() + wt->wave->hpr_len()) / 2;

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
    for (iter = 0; iter < J; ++iter) {
        if (iter > 0) {
            M = M / 2;
        }
        ir = wt->dimensions[2 * iter];
        ic = wt->dimensions[2 * iter + 1];
        istride = ic;
        ostride = ic;
        aLH = wt->coeffaccess[iter * 3 + 1];
        aHL = wt->coeffaccess[iter * 3 + 2];
        aHH = wt->coeffaccess[iter * 3 + 3];
        for (auto i = 0; i < ic; ++i) {
            imodwt_per_stride(M, orig + i, ir, wavecoeff + aLH + i, filt.get(), lf, cL.get() + i, istride, ostride);
            imodwt_per_stride(M, wavecoeff + aHL + i, ir, wavecoeff + aHH + i, filt.get(), lf, cH.get() + i, istride, ostride);
        }

        istride = 1;
        ostride = 1;

        for (auto i = 0; i < ir; ++i) {
            imodwt_per_stride(M, cL.get() + i * ic, ic, cH.get() + i * ic, filt.get(), lf, oup + i * ic, istride, ostride);
        }

        orig = oup;
    }
}

auto getWT2Coeffs(wt2_set* wt, double* wcoeffs, int level, char const* type, int* rows, int* cols) -> double*
{
    int J;
    int iter;
    int t;
    double* ptr;
    J = wt->J;
    // Error Check

    if (level > J || level < 1) {
        printf("Error : The data is decomposed into %d levels so the acceptable values of level are between 1 and %d", J, J);
        exit(-1);
    }

    if ((strcmp(type, "A") == 0) && level != J) {
        printf("Approximation Coefficients are only available for level %d", J);
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

    iter += (J - level) * 3;

    ptr = wcoeffs + wt->coeffaccess[iter];
    *rows = wt->dimensions[2 * (J - level)];
    *cols = wt->dimensions[2 * (J - level) + 1];

    return ptr;
}

void dispWT2Coeffs(double* A, int row, int col)
{
    printf("\n MATRIX Order : %d X %d \n \n", row, col);

    for (auto i = 0; i < row; i++) {
        printf("R%d: ", i);
        for (auto j = 0; j < col; j++) {
            printf("%g ", A[i * col + j]);
        }
        printf(":R%d \n", i);
    }
}

void wave_summary(wavelet const& obj)
{
    auto const N = obj.size();
    printf("\n");
    printf("Wavelet Name: %s \n", obj.name().c_str());
    printf("\n");
    printf("Wavelet Filters \n");
    printf("lpd: [");
    for (auto i = 0; i < N - 1; ++i) {
        printf("%g,", obj.lpd()[i]);
    }
    printf("%g", obj.lpd()[N - 1]);
    printf("] \n");
    printf("hpd: [");
    for (auto i = 0; i < N - 1; ++i) {
        printf("%g,", obj.hpd()[i]);
    }
    printf("%g", obj.hpd()[N - 1]);
    printf("] \n");
    printf("lpr: [");
    for (auto i = 0; i < N - 1; ++i) {
        printf("%g,", obj.lpr()[i]);
    }
    printf("%g", obj.lpr()[N - 1]);
    printf("] \n");
    printf("hpr: [");
    for (auto i = 0; i < N - 1; ++i) {
        printf("%g,", obj.hpr()[i]);
    }
    printf("%g", obj.hpr()[N - 1]);
    printf("] \n");
}

void wt_summary(wt_set* wt)
{
    int J;
    int t;
    J = wt->J;
    wave_summary(*wt->wave);
    printf("\n");
    printf("Wavelet Transform : %s \n", wt->method.c_str());
    printf("Signal Extension : %s \n", wt->ext.c_str());
    printf("Convolutional Method : %s \n", wt->cmethod.c_str());
    printf("Number of Decomposition Levels %d \n", wt->J);
    printf("Length of Input Signal %d \n", wt->siglength);
    printf("Length of WT Output Vector %d \n", wt->outlength);
    printf("Wavelet Coefficients are contained in vector : %s \n", "output");
    printf("Approximation Coefficients \n");
    printf("Level %d Access : output[%d] Length : %d \n", J, 0, wt->length[0]);
    printf("Detail Coefficients \n");
    t = wt->length[0];
    for (auto i = 0; i < J; ++i) {
        printf("Level %d Access : output[%d] Length : %d \n", J - i, t, wt->length[i + 1]);
        t += wt->length[i + 1];
    }
    printf("\n");
}

void wtree_summary(wtree_set* wt)
{
    int k;
    int p2;
    int J;
    int t;
    J = wt->J;
    wave_summary(*wt->wave);
    printf("\n");
    printf("Wavelet Transform : %s \n", wt->method.c_str());
    printf("\n");
    printf("Signal Extension : %s \n", wt->ext.c_str());
    printf("\n");
    printf("Number of Decomposition Levels %d \n", wt->J);
    printf("\n");
    printf("Length of Input Signal %d \n", wt->siglength);
    printf("\n");
    printf("Length of WT Output Vector %d \n", wt->outlength);
    printf("\n");
    printf("Wavelet Coefficients are contained in vector : %s \n", "output");
    printf("\n");
    printf("Coefficients Access \n");
    t = 0;
    p2 = 2;
    for (auto i = 0; i < J; ++i) {
        for (k = 0; k < p2; ++k) {
            printf("Node %d %d Access : output[%d] Length : %d \n", i + 1, k, wt->nodelength[t], wt->length[J - i]);
            t++;
        }
        p2 *= 2;
    }
    printf("\n");
}

void wpt_summary(wpt_set* wt)
{
    int k;
    int p2;
    int J;
    int it1;
    int it2;
    J = wt->J;
    wave_summary(*wt->wave);
    printf("\n");
    printf("Signal Extension : %s \n", wt->ext.c_str());
    printf("\n");
    printf("Entropy : %s \n", wt->entropy.c_str());
    printf("\n");
    printf("Number of Decomposition Levels %d \n", wt->J);
    printf("\n");
    printf("Number of Active Nodes %d \n", wt->nodes);
    printf("\n");
    printf("Length of Input Signal %d \n", wt->siglength);
    printf("\n");
    printf("Length of WT Output Vector %d \n", wt->outlength);
    printf("\n");
    printf("Wavelet Coefficients are contained in vector : %s \n", "output");
    printf("\n");
    printf("Coefficients Access \n");
    it1 = 1;
    it2 = 0;
    for (auto i = 0; i < J; ++i) {
        it1 += ipow2(i + 1);
    }
    for (auto i = J; i > 0; --i) {
        p2 = ipow2(i);
        it1 -= p2;
        for (k = 0; k < p2; ++k) {
            if (wt->basisvector[it1 + k] == 1) {
                printf("Node %d %d Access : output[%d] Length : %d \n", i, k, it2, wt->length[J - i + 1]);
                it2 += wt->length[J - i + 1];
            }
        }
    }

    printf("\n");
}

void cwt_summary(cwt_set* wt)
{

    printf("\n");
    printf("Wavelet : %s Parameter %lf \n", wt->wave.c_str(), wt->m);
    printf("\n");
    printf("Length of Input Signal : %d \n", wt->siglength);
    printf("\n");
    printf("Sampling Rate : %g \n", wt->dt);
    printf("\n");
    printf("Total Number of Scales : %d \n", wt->J);
    printf("\n");
    printf("Smallest Scale (s0) : %lf \n", wt->s0);
    printf("\n");
    printf("Separation Between Scales (dj) %lf \n", wt->dj);
    printf("\n");
    printf("Scale Type %s \n", wt->type.c_str());
    printf("\n");
    printf("Complex CWT Output Vector is of size %d * %d stored in Row Major format \n", wt->J, wt->siglength);
    printf("\n");
    printf("The ith real value can be accessed using wt->output[i].re and imaginary value by wt->output[i].im \n");
    printf("\n");
}

void wt2_summary(wt2_set* wt)
{
    int J;
    int t;
    int rows;
    int cols;
    int vsize;
    J = wt->J;
    wave_summary(*wt->wave);
    printf("\n");
    printf("Wavelet Transform : %s \n", wt->method.c_str());
    printf("\n");
    printf("Signal Extension : %s \n", wt->ext.c_str());
    printf("\n");
    printf("Number of Decomposition Levels %d \n", wt->J);
    printf("\n");
    printf("Input Signal Rows %d \n", wt->rows);
    printf("\n");
    printf("Input Signal Cols %d \n", wt->cols);
    printf("\n");
    printf("Length of Wavelet Coefficients Vector %d \n", wt->outlength);
    printf("\n");
    t = 0;
    for (auto i = J; i > 0; --i) {
        rows = wt->dimensions[2 * (J - i)];
        cols = wt->dimensions[2 * (J - i) + 1];
        vsize = rows * cols;
        printf("Level %d Decomposition Rows :%d Columns:%d Vector Size (Rows*Cols):%d \n", i, rows, cols, vsize);
        printf("Access Row values stored at wt->dimensions[%d]\n", 2 * (J - i));
        printf("Access Column values stored at wt->dimensions[%d]\n\n", 2 * (J - i) + 1);

        if (i == J) {
            printf("Approximation Coefficients access at wt->coeffaccess[%d]=%d, Vector size:%d \n", t, wt->coeffaccess[t], vsize);
        }

        t += 1;
        printf("Horizontal Coefficients access at wt->coeffaccess[%d]=%d, Vector size:%d \n", t, wt->coeffaccess[t], vsize);
        t += 1;
        printf("Vertical Coefficients access at wt->coeffaccess[%d]=%d, Vector size:%d \n", t, wt->coeffaccess[t], vsize);
        t += 1;
        printf("Diagonal Coefficients access at wt->coeffaccess[%d]=%d, Vector size:%d \n\n", t, wt->coeffaccess[t], vsize);
    }
}

void wt_free(wt_set* object)
{
    delete object;
}

void wtree_free(wtree_set* object)
{
    delete object;
}

void wpt_free(wpt_set* object)
{
    delete object;
}

void cwt_free(cwt_set* object)
{
    delete object;
}

void wt2_free(wt2_set* wt)
{
    delete wt;
}
