#include "ComplexWaveletTransform.hpp"

#include "lt/dsp/convolution/Convolution.hpp"
#include "lt/dsp/fft/FFT.hpp"
#include "lt/dsp/wavelets/transforms/common.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
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

    obj->output = (Complex<double>*)&obj->params[0];
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
    printf("The ith real value can be accessed using wt.output()[i].real() and imaginary value by wt.output()[i].imag() \n");
    printf("\n");
}

auto cwtFree(ComplexWaveletTransform* object) -> void
{
    delete object;
}

static auto factorial(int n) -> double
{
    static const double fact[41] = {
        1.0,
        1.0,
        2.0,
        6.0,
        24.0,
        120.0,
        720.0,
        5040.0,
        40320.0,
        362880.0,
        3628800.0,
        39916800.0,
        479001600.0,
        6227020800.0,
        87178291200.0,
        1307674368000.0,
        20922789888000.0,
        355687428096000.0,
        6402373705728000.0,
        121645100408832000.0,
        2432902008176640000.0,
        51090942171709440000.0,
        1124000727777607680000.0,
        25852016738884976640000.0,
        620448401733239439360000.0,
        15511210043330985984000000.0,
        403291461126605635584000000.0,
        10888869450418352160768000000.0,
        304888344611713860501504000000.0,
        8841761993739701954543616000000.0,
        265252859812191058636308480000000.0,
        8222838654177922817725562880000000.0,
        263130836933693530167218012160000000.0,
        8683317618811886495518194401280000000.0,
        295232799039604140847618609643520000000.0,
        10333147966386144929666651337523200000000.0,
        371993326789901217467999448150835200000000.0,
        13763753091226345046315979581580902400000000.0,
        523022617466601111760007224100074291200000000.0,
        20397882081197443358640281739902897356800000000.0,
        815915283247897734345611269596115894272000000000.0,
    };

    if (n > 40 || n < 0) {
        printf("This program is only valid for 0 <= N <= 40 \n");
        return -1.0;
    }

    return fact[n];
}

static auto waveFunction(int nk, double dt, int mother, double param, double scale1, double const* kwave, double pi, double* period1,
    double* coi1, Complex<double>* daughter) -> void
{

    double norm;
    double expnt;
    double fourierFactor;
    int k;
    int m;
    double temp;
    int sign;
    int re;

    if (mother == 0) {
        //MORLET
        if (param < 0.0) {
            param = 6.0;
        }
        norm = std::sqrt(2.0 * pi * scale1 / dt) * std::pow(pi, -0.25);

        for (k = 1; k <= nk / 2 + 1; ++k) {
            temp = (scale1 * kwave[k - 1] - param);
            expnt = -0.5 * temp * temp;
            daughter[k - 1].real(norm * std::exp(expnt));
            daughter[k - 1].imag(0.0);
        }
        for (k = nk / 2 + 2; k <= nk; ++k) {
            daughter[k - 1].real(0.0);
            daughter[k - 1].imag(0.0);
        }
        fourierFactor = (4.0 * pi) / (param + std::sqrt(2.0 + param * param));
        *period1 = scale1 * fourierFactor;
        *coi1 = fourierFactor / std::sqrt(2.0);
    } else if (mother == 1) {
        // PAUL
        if (param < 0.0) {
            param = 4.0;
        }
        m = (int)param;
        norm = std::sqrt(2.0 * pi * scale1 / dt) * (std::pow(2.0, (double)m) / std::sqrt((double)(m * factorial(2 * m - 1))));
        for (k = 1; k <= nk / 2 + 1; ++k) {
            temp = scale1 * kwave[k - 1];
            expnt = -temp;
            daughter[k - 1].real(norm * std::pow(temp, (double)m) * std::exp(expnt));
            daughter[k - 1].imag(0.0);
        }
        for (k = nk / 2 + 2; k <= nk; ++k) {
            daughter[k - 1].real(0.0);
            daughter[k - 1].imag(0.0);
        }
        fourierFactor = (4.0 * pi) / (2.0 * m + 1.0);
        *period1 = scale1 * fourierFactor;
        *coi1 = fourierFactor * std::sqrt(2.0);
    } else if (mother == 2) {
        if (param < 0.0) {
            param = 2.0;
        }
        m = (int)param;

        if (m % 2 == 0) {
            re = 1;
        } else {
            re = 0;
        }

        if (m % 4 == 0 || m % 4 == 1) {
            sign = -1;
        } else {
            sign = 1;
        }

        norm = std::sqrt(2.0 * pi * scale1 / dt) * std::sqrt(1.0 / cwtGamma(m + 0.50));
        norm *= sign;

        if (re == 1) {
            for (k = 1; k <= nk; ++k) {
                temp = scale1 * kwave[k - 1];
                daughter[k - 1].real(norm * std::pow(temp, (double)m) * std::exp(-0.50 * std::pow(temp, 2.0)));
                daughter[k - 1].imag(0.0);
            }
        } else if (re == 0) {
            for (k = 1; k <= nk; ++k) {
                temp = scale1 * kwave[k - 1];
                daughter[k - 1].real(0.0);
                daughter[k - 1].imag(norm * std::pow(temp, (double)m) * std::exp(-0.50 * std::pow(temp, 2.0)));
            }
        }
        fourierFactor = (2.0 * pi) * std::sqrt(2.0 / (2.0 * m + 1.0));
        *period1 = scale1 * fourierFactor;
        *coi1 = fourierFactor / std::sqrt(2.0);
    }
}

auto cwavelet(double const* y, int n, double dt, int mother, double param, double s0, double dj, int jtot, int npad,
    double* wave, double const* scale, double* period, double* coi) -> void
{

    int j;
    int k;
    int iter;
    double ymean;
    double freq1;
    double pi;
    double period1;
    double coi1;
    double tmp1;
    double tmp2;
    double scale1;

    (void)s0;
    (void)dj; /* yes, we need these parameters unused */

    pi = 4.0 * atan(1.0);

    if (npad < n) {
        printf("npad must be >= N \n");
        exit(-1);
    }

    auto obj = std::make_unique<FFT>(npad, 1);
    auto iobj = std::make_unique<FFT>(npad, -1);

    auto ypad = std::make_unique<Complex<double>[]>(npad);
    auto yfft = std::make_unique<Complex<double>[]>(npad);
    auto daughter = std::make_unique<Complex<double>[]>(npad);
    auto kwave = std::make_unique<double[]>(npad);
    ymean = 0.0;

    for (auto i = 0; i < n; ++i) {
        ymean += y[i];
    }

    ymean /= n;

    for (auto i = 0; i < n; ++i) {
        ypad[i].real(y[i] - ymean);
        ypad[i].imag(0.0);
    }

    for (auto i = n; i < npad; ++i) {
        ypad[i].real(0.0);
        ypad[i].imag(0.0);
    }

    // Find FFT of the input y (ypad)

    obj->perform(ypad.get(), yfft.get());

    for (auto i = 0; i < npad; ++i) {
        yfft[i].real(yfft[i].real() / (double)npad);
        yfft[i].imag(yfft[i].imag() / (double)npad);
    }

    //Construct the wavenumber array

    freq1 = 2.0 * pi / ((double)npad * dt);
    kwave[0] = 0.0;

    for (auto i = 1; i < npad / 2 + 1; ++i) {
        kwave[i] = i * freq1;
    }

    for (auto i = npad / 2 + 1; i < npad; ++i) {
        kwave[i] = -kwave[npad - i];
    }

    // Main loop

    for (j = 1; j <= jtot; ++j) {
        scale1 = scale[j - 1]; // = s0*std::pow(2.0, (double)(j - 1)*dj);
        waveFunction(npad, dt, mother, param, scale1, kwave.get(), pi, &period1, &coi1, daughter.get());
        period[j - 1] = period1;
        for (k = 0; k < npad; ++k) {
            tmp1 = daughter[k].real() * yfft[k].real() - daughter[k].imag() * yfft[k].imag();
            tmp2 = daughter[k].real() * yfft[k].imag() + daughter[k].imag() * yfft[k].real();
            daughter[k].real(tmp1);
            daughter[k].imag(tmp2);
        }
        iobj->perform(daughter.get(), ypad.get());
        iter = 2 * (j - 1) * n;
        for (auto i = 0; i < n; ++i) {
            wave[iter + 2 * i] = ypad[i].real();
            wave[iter + 2 * i + 1] = ypad[i].imag();
        }
    }

    for (auto i = 1; i <= (n + 1) / 2; ++i) {
        coi[i - 1] = coi1 * dt * ((double)i - 1.0);
        coi[n - i] = coi[i - 1];
    }
}

auto psi0(int mother, double param, double* val, int* real) -> void
{
    double pi;
    double coeff;
    int m;
    int sign;

    m = (int)param;
    pi = 4.0 * atan(1.0);

    if (mother == 0) {
        // Morlet
        *val = 1.0 / std::sqrt(std::sqrt(pi));
        *real = 1;
    } else if (mother == 1) {
        //Paul
        if (m % 2 == 0) {
            *real = 1;
        } else {
            *real = 0;
        }

        if (m % 4 == 0 || m % 4 == 1) {
            sign = 1;
        } else {
            sign = -1;
        }
        *val = sign * std::pow(2.0, (double)m) * factorial(m) / (std::sqrt(pi * factorial(2 * m)));

    } else if (mother == 2) {
        // D.O.G
        *real = 1;

        if (m % 2 == 0) {
            if (m % 4 == 0) {
                sign = -1;
            } else {
                sign = 1;
            }
            coeff = sign * std::pow(2.0, (double)m / 2) / cwtGamma(0.5);
            *val = coeff * cwtGamma(((double)m + 1.0) / 2.0) / std::sqrt(cwtGamma(m + 0.50));
        } else {
            *val = 0;
        }
    }
}

static auto maxabs(double* array, int n) -> int
{
    double maxval;
    double temp;

    int index;
    maxval = 0.0;
    index = -1;

    for (auto i = 0; i < n; ++i) {
        temp = fabs(array[i]);
        if (temp >= maxval) {
            maxval = temp;
            index = i;
        }
    }

    return index;
}

auto cdelta(int mother, double param, double psi0) -> double
{
    auto s0 { 0.0 };
    auto n { 1 };
    auto subscale = 8.0;
    auto dt = 0.25;

    if (mother == 0) {
        n = 16;
        s0 = dt / 4;
    } else if (mother == 1) {
        n = 16;
        s0 = dt / 4.0;
    } else if (mother == 2) {
        s0 = dt / 8.0;
        n = 256;
        if (param == 2.0) {
            subscale = 16.0;
            s0 = dt / 16.0;
            n = 2048;
        }
    }

    auto const dj = 1.0 / subscale;
    auto const jtot = 16 * (int)subscale;

    auto delta = std::make_unique<double[]>(n);
    auto wave = std::make_unique<double[]>(2 * n * jtot);
    auto coi = std::make_unique<double[]>(n);
    auto scale = std::make_unique<double[]>(jtot);
    auto period = std::make_unique<double[]>(jtot);
    auto mval = std::make_unique<double[]>(n);

    delta[0] = 1;

    for (auto i = 1; i < n; ++i) {
        delta[i] = 0;
    }

    for (auto i = 0; i < jtot; ++i) {
        scale[i] = s0 * std::pow(2.0, (double)(i)*dj);
    }

    cwavelet(delta.get(), n, dt, mother, param, s0, dj, jtot, n, wave.get(), scale.get(), period.get(), coi.get());

    for (auto i = 0; i < n; ++i) {
        mval[i] = 0;
    }

    for (auto j = 0; j < jtot; ++j) {
        auto const iter = 2 * j * n;
        auto const den = std::sqrt(scale[j]);
        for (auto i = 0; i < n; ++i) {
            mval[i] += wave[iter + 2 * i] / den;
        }
    }

    auto const maxarr = maxabs(mval.get(), n);
    auto const cdel = std::sqrt(dt) * dj * mval[maxarr] / psi0;
    return cdel;
}

auto icwavelet(double const* wave, int n, double* scale, int jtot, double dt, double dj, double cdelta, double psi0, double* oup) -> void
{

    int j;
    int iter;
    double den;
    double coeff;

    coeff = std::sqrt(dt) * dj / (cdelta * psi0);

    for (auto i = 0; i < n; ++i) {
        oup[i] = 0.0;
    }

    for (j = 0; j < jtot; ++j) {
        iter = 2 * j * n;
        den = std::sqrt(scale[j]);
        for (auto i = 0; i < n; ++i) {
            oup[i] += wave[iter + 2 * i] / den;
        }
    }

    for (auto i = 0; i < n; ++i) {
        oup[i] *= coeff;
    }
}
