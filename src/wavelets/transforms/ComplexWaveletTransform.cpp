#include "ComplexWaveletTransform.hpp"

#include "wavelets/conv.h"
#include "wavelets/cwt.h"
#include "wavelets/hsfft.h"
#include "wavelets/wtmath.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string_view>

using namespace std::string_view_literals;

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

auto cwtFree(ComplexWaveletTransform* object) -> void
{
    delete object;
}
