#include "ComplexWaveletTransform.hpp"

#include "lt/dsp/convolution/FFTConvolver.hpp"
#include "lt/dsp/fft/FFT.hpp"
#include "lt/dsp/wavelets/common.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string_view>

using namespace std::string_view_literals;

namespace {
#include "lt/dsp/fft/FFT.hpp"

#include <cmath>
#include <cstdio>
#include <memory>

namespace {

    auto factorial(int n) -> double
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

    auto nsfftFd(FFT* obj, Complex<double>* inp, Complex<double>* oup, double lb, double ub, double* w) -> void
    {
        auto const n = obj->size();
        auto const l = n / 2;

        auto temp1 = std::make_unique<double[]>(l);
        auto temp2 = std::make_unique<double[]>(l);

        auto const delta = (ub - lb) / n;
        auto const den = 2 * (ub - lb);
        auto j = -n;

        for (auto i = 0; i < n; ++i) {
            w[i] = (double)j / den;
            j += 2;
        }

        obj->perform(inp, oup);

        for (auto i = 0; i < l; ++i) {
            temp1[i] = oup[i].real();
            temp2[i] = oup[i].imag();
        }

        for (auto i = 0; i < n - l; ++i) {
            oup[i].real(oup[i + l].real());
            oup[i].imag(oup[i + l].imag());
        }

        for (auto i = 0; i < l; ++i) {
            oup[n - l + i].real(temp1[i]);
            oup[n - l + i].imag(temp2[i]);
        }

        auto const plb = PI2 * lb;

        for (auto i = 0; i < n; ++i) {
            auto const tempr = oup[i].real();
            auto const tempi = oup[i].imag();
            auto const theta = w[i] * plb;

            oup[i].real(delta * (tempr * cos(theta) + tempi * sin(theta)));
            oup[i].imag(delta * (tempi * cos(theta) - tempr * sin(theta)));
        }
    }

    auto nsfftBk(FFT* obj, Complex<double>* inp, Complex<double>* oup, double lb, double ub, double* t) -> void
    {

        auto const n = obj->size();
        auto const l = n / 2;

        auto const m = divideby(n, 2);

        if (m == 0) {
            printf("The Non-Standard FFT Length must be a power of 2");
            exit(1);
        }

        auto temp1 = std::make_unique<double[]>(l);
        auto temp2 = std::make_unique<double[]>(l);
        auto w = std::make_unique<double[]>(n);
        auto inpt = std::make_unique<Complex<double>[]>(n);

        auto const delta = (ub - lb) / n;
        auto const den = 2 * (ub - lb);
        auto j = -n;

        for (auto i = 0; i < n; ++i) {
            w[i] = (double)j / den;
            j += 2;
        }

        auto const plb = PI2 * lb;

        for (auto i = 0; i < n; ++i) {
            auto const theta = w[i] * plb;
            inpt[i].real((inp[i].real() * cos(theta) - inp[i].imag() * sin(theta)) / delta);
            inpt[i].imag((inp[i].imag() * cos(theta) + inp[i].real() * sin(theta)) / delta);
        }

        for (auto i = 0; i < l; ++i) {
            temp1[i] = inpt[i].real();
            temp2[i] = inpt[i].imag();
        }

        for (auto i = 0; i < n - l; ++i) {
            inpt[i].real(inpt[i + l].real());
            inpt[i].imag(inpt[i + l].imag());
        }

        for (auto i = 0; i < l; ++i) {
            inpt[n - l + i].real(temp1[i]);
            inpt[n - l + i].imag(temp2[i]);
        }

        obj->perform(inpt.get(), oup);

        for (auto i = 0; i < n; ++i) {
            t[i] = lb + i * delta;
        }
    }

    auto nsfftExec(FFT* fft, Complex<double>* inp, Complex<double>* oup, double lb, double ub, double* w) -> void
    {
        if (fft->direction() == FFT::forward) {
            nsfftFd(fft, inp, oup, lb, ub, w);
        } else {
            nsfftBk(fft, inp, oup, lb, ub, w);
        }
    }

}

auto roundTowardsZero(double x) -> double
{
    return x >= 0.0 ? std::floor(x) : std::ceil(x);
}

auto cwtGamma(double x) -> double
{
    /*
	 * This C program code is based on  W J Cody's fortran code.
	 * http://www.netlib.org/specfun/gamma
	 *
	 * References:
   "An Overview of Software Development for Special Functions",
	W. J. Cody, Lecture Notes in Mathematics, 506,
	Numerical Analysis Dundee, 1975, G. A. Watson (ed.),
	Springer Verlag, Berlin, 1976.

   Computer Approximations, Hart, Et. Al., Wiley and sons, New York, 1968.
   */

    // numerator and denominator coefficients for 1 <= x <= 2

    double y;
    double oup;
    double fact;
    double sum;
    double y2;
    double yi;
    double z;
    double nsum;
    double dsum;
    int swi;
    int n;

    constexpr double spi = 0.9189385332046727417803297;
    constexpr double pi = 3.1415926535897932384626434;
    constexpr double xmax = 171.624e+0;
    constexpr double xinf = 1.79e308;
    constexpr double eps = 2.22e-16;
    constexpr double xninf = 1.79e-308;

    double num[8] = {
        -1.71618513886549492533811e+0,
        2.47656508055759199108314e+1,
        -3.79804256470945635097577e+2,
        6.29331155312818442661052e+2,
        8.66966202790413211295064e+2,
        -3.14512729688483675254357e+4,
        -3.61444134186911729807069e+4,
        6.64561438202405440627855e+4,
    };

    double den[8] = {
        -3.08402300119738975254353e+1,
        3.15350626979604161529144e+2,
        -1.01515636749021914166146e+3,
        -3.10777167157231109440444e+3,
        2.25381184209801510330112e+4,
        4.75584627752788110767815e+3,
        -1.34659959864969306392456e+5,
        -1.15132259675553483497211e+5,
    };

    // Coefficients for Hart's Minimax approximation x >= 12

    double c[7] = {
        -1.910444077728e-03,
        8.4171387781295e-04,
        -5.952379913043012e-04,
        7.93650793500350248e-04,
        -2.777777777777681622553e-03,
        8.333333333333333331554247e-02,
        5.7083835261e-03,
    };

    y = x;
    swi = 0;
    fact = 1.0;
    n = 0;

    if (y < 0.) {
        // Negative x
        y = -x;
        yi = roundTowardsZero(y);
        oup = y - yi;

        if (oup != 0.0) {
            if (yi != roundTowardsZero(yi * .5) * 2.) {
                swi = 1;
            }
            fact = -pi / sin(pi * oup);
            y += 1.;
        } else {
            return xinf;
        }
    }

    if (y < eps) {
        if (y >= xninf) {
            oup = 1.0 / y;
        } else {
            return xinf;
        }

    } else if (y < 12.) {
        yi = y;
        if (y < 1.) {
            z = y;
            y += 1.;
        } else {
            n = (int)y - 1;
            y -= (double)n;
            z = y - 1.0;
        }
        nsum = 0.;
        dsum = 1.;
        for (auto i = 0; i < 8; ++i) {
            nsum = (nsum + num[i]) * z;
            dsum = dsum * z + den[i];
        }
        oup = nsum / dsum + 1.;

        if (yi < y) {

            oup /= yi;
        } else if (yi > y) {

            for (auto i = 0; i < n; ++i) {
                oup *= y;
                y += 1.;
            }
        }

    } else {
        if (y <= xmax) {
            y2 = y * y;
            sum = c[6];
            for (auto i = 0; i < 6; ++i) {
                sum = sum / y2 + c[i];
            }
            sum = sum / y - y + spi;
            sum += (y - .5) * std::log(y);
            oup = std::exp(sum);
        } else {
            return (xinf);
        }
    }

    if (swi != 0) {
        oup = -oup;
    }
    if (fact != 1.) {
        oup = fact / oup;
    }

    return oup;
}

auto waveFunction(int nk, double dt, int mother, double param, double scale1, double const* kwave, double pi, double* period1,
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

    pi = 4.0 * std::atan(1.0);

    if (npad < n) {
        printf("npad must be >= N \n");
        exit(-1);
    }

    auto obj = std::make_unique<FFT>(npad, FFT::forward);
    auto iobj = std::make_unique<FFT>(npad, FFT::backward);

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
    pi = 4.0 * std::atan(1.0);

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

ComplexWaveletTransform::ComplexWaveletTransform(char const* wave, double param, int siglength, double dtIn, int j)
{

    double s0Tmp {};
    double djTmp {};
    char const* pdefault = "pow";

    auto const mm = (int)param;
    auto odd = 1;
    if (2 * (mm / 2) == mm) {
        odd = 0;
    }

    auto n = siglength;
    auto nj2 = 2 * n * j;
    params = std::make_unique<double[]>(nj2 + 2 * j + n);

    int motherTmp { 0 };
    if ((wave == "morlet"sv) || (wave == "morl"sv)) {
        s0Tmp = 2 * dtIn;
        djTmp = 0.4875;
        motherTmp = 0;
        if (param < 0.0) {
            printf("\n Morlet Wavelet Parameter should be >= 0 \n");
            exit(-1);
        }
        if (param == 0) {
            param = 6.0;
        }
        wave_ = "morlet";

    } else if (wave == "paul"sv) {
        s0Tmp = 2 * dtIn;
        djTmp = 0.4875;
        motherTmp = 1;
        if (param < 0 || param > 20) {
            printf("\n Paul Wavelet Parameter should be > 0 and <= 20 \n");
            exit(-1);
        }
        if (param == 0) {
            param = 4.0;
        }
        wave_ = "paul";

    } else if ((wave == "dgauss"sv) || (wave == "dog"sv)) {
        s0Tmp = 2 * dtIn;
        djTmp = 0.4875;
        motherTmp = 2;
        if (param < 0 || odd == 1) {
            printf("\n DOG Wavelet Parameter should be > 0 and even \n");
            exit(-1);
        }
        if (param == 0) {
            param = 2.0;
        }
        wave_ = "dog";
    }

    pow = 2;
    type = pdefault;

    s0 = s0Tmp;
    dj = djTmp;
    dt = dtIn;
    J = j;
    signalLength = siglength;
    sflag = 0;
    pflag = 1;
    mother = motherTmp;
    m = param;

    auto const t1 = 0.499999 + std::log((double)n) / std::log(2.0);
    auto const ibase2 = 1 + (int)t1;

    npad = (int)std::pow(2.0, (double)ibase2);

    output = (Complex<double>*)&params[0];
    scale = &params[nj2];
    period = &params[nj2 + j];
    coi = &params[nj2 + 2 * j];

    for (auto i = 0; i < nj2 + 2 * j + n; ++i) {
        params[i] = 0.0;
    }
}

auto ComplexWaveletTransform::scales(double newS0, double newDj, char const* newType, int power) -> void
{
    type = newType;
    //newS0*std::pow(2.0, (double)(j - 1)*newDj);
    if ((type == "pow"sv) || (type == "power"sv)) {
        for (auto i = 0; i < J; ++i) {
            scale[i] = newS0 * std::pow((double)power, (double)(i)*newDj);
        }
        sflag = 1;
        pow = power;

    } else if ((type == "lin"sv) || (type == "linear"sv)) {
        for (auto i = 0; i < J; ++i) {
            scale[i] = newS0 + (double)i * newDj;
        }
        sflag = 1;
    } else {
        printf("\n Type accepts only two values : pow and lin\n");
        exit(-1);
    }
    s0 = newS0;
    dj = newDj;
}

auto cwt(ComplexWaveletTransform& wt, double const* inp) -> void
{
    int n;
    int npad;
    int nj2;
    int j;
    int j2;
    n = wt.signalLength;
    if (wt.sflag == 0) {
        for (auto i = 0; i < wt.J; ++i) {
            wt.scale[i] = wt.s0 * std::pow(2.0, (double)(i)*wt.dj);
        }
        wt.sflag = 1;
    }

    if (wt.pflag == 0) {
        npad = n;
    } else {
        npad = wt.npad;
    }

    nj2 = 2 * n * wt.J;
    j = wt.J;
    j2 = 2 * j;

    wt.smean = 0.0;

    for (auto i = 0; i < n; ++i) {
        wt.smean += inp[i];
    }
    wt.smean /= n;

    cwavelet(inp, n, wt.dt, wt.mother, wt.m, wt.s0, wt.dj, wt.J, npad, wt.params.get(), wt.params.get() + nj2, wt.params.get() + nj2 + j, wt.params.get() + nj2 + j2);
}

auto icwt(ComplexWaveletTransform& wt, double* cwtop) -> void
{
    double psi;
    double cdel;
    int real;
    int n;
    int nj2;

    n = wt.signalLength;
    nj2 = n * 2 * wt.J;

    psi0(wt.mother, wt.m, &psi, &real);
    cdel = cdelta(wt.mother, wt.m, psi);

    if (((wt.type == "pow"sv) || (wt.type == "power"sv)) && wt.pow == 2) {
        icwavelet(wt.params.get(), n, wt.params.get() + nj2, wt.J, wt.dt, wt.dj, cdel, psi, cwtop);
    } else {
        printf("Inverse CWT is only available for power of 2.0 scales \n");
        exit(-1);
    }
    for (auto i = 0; i < n; ++i) {
        cwtop[i] += wt.smean;
    }
}

auto summary(ComplexWaveletTransform const& wt) -> void
{
    printf("\nWavelet : %s Parameter %lf \n", wt.wave().c_str(), wt.m);
    printf("\nLength of Input Signal : %d \n", wt.signalLength);
    printf("\nSampling Rate : %g \n", wt.dt);
    printf("\nTotal Number of Scales : %d \n", wt.J);
    printf("\nSmallest Scale (s0) : %lf \n", wt.s0);
    printf("\nSeparation Between Scales (dj) %lf \n", wt.dj);
    printf("\nScale Type %s \n", wt.type.c_str());
    printf("\nComplex CWT Output Vector is of size %d  const& %d stored in Row Major format \n", wt.J, wt.signalLength);
    printf("\nThe ith real value can be accessed using wt.output()[i].real() and imaginary value by wt.output()[i].imag() \n");
    printf("\n");
}

auto meyer(int n, double lb, double ub, double* phi, double* psi, double* tgrid) -> void
{
    auto const m = divideby(n, 2);
    if (m == 0) {
        printf("Size of Wavelet must be a power of 2");
        exit(1);
    }
    if (lb >= ub) {
        printf("upper bound must be greater than lower bound");
        exit(1);
    }

    auto obj = std::make_unique<FFT>(n, FFT::backward);
    auto w = std::make_unique<double[]>(n);
    auto phiw = std::make_unique<Complex<double>[]>(n);
    auto psiw = std::make_unique<Complex<double>[]>(n);
    auto oup = std::make_unique<Complex<double>[]>(n);

    auto const delta = 2 * (ub - lb) / PI2;

    auto j = (double)n;
    j *= -1.0;

    for (auto i = 0; i < n; ++i) {
        w[i] = j / delta;
        j += 2.0;
        psiw[i].real(0.0);
        psiw[i].imag(0.0);
        phiw[i].real(0.0);
        phiw[i].imag(0.0);
    }

    for (auto i = 0; i < n; ++i) {
        auto const wf = fabs(w[i]);
        if (wf <= PI2 / 3.0) {
            phiw[i].real(1.0);
        }
        if (wf > PI2 / 3.0 && wf <= 2 * PI2 / 3.0) {
            auto const x = (3 * wf / PI2) - 1.0;
            auto const x2 = x * x;
            auto const x3 = x2 * x;
            auto const x4 = x3 * x;
            auto const v = x4 * (35 - 84 * x + 70 * x2 - 20 * x3);
            auto const theta = v * PI2 / 4.0;
            auto const cs = std::cos(theta);
            auto const sn = std::sin(theta);

            phiw[i].real(cs);
            psiw[i].real(std::cos(w[i] / 2.0) * sn);
            psiw[i].imag(std::sin(w[i] / 2.0) * sn);
        }
        if (wf > 2.0 * PI2 / 3.0 && wf <= 4 * PI2 / 3.0) {
            auto const x = (1.5 * wf / PI2) - 1.0;
            auto const x2 = x * x;
            auto const x3 = x2 * x;
            auto const x4 = x3 * x;
            auto const v = x4 * (35 - 84 * x + 70 * x2 - 20 * x3);
            auto const theta = v * PI2 / 4.0;
            auto const cs = std::cos(theta);

            psiw[i].real(std::cos(w[i] / 2.0) * cs);
            psiw[i].imag(std::sin(w[i] / 2.0) * cs);
        }
    }

    nsfftExec(obj.get(), phiw.get(), oup.get(), lb, ub, tgrid);

    for (auto i = 0; i < n; ++i) {
        phi[i] = oup[i].real() / n;
    }

    nsfftExec(obj.get(), psiw.get(), oup.get(), lb, ub, tgrid);

    for (auto i = 0; i < n; ++i) {
        psi[i] = oup[i].real() / n;
    }
}

auto gauss(int n, int p, double lb, double ub, double* psi, double* t) -> void
{
    double delta;
    double num;
    double den;
    double t2;
    double t4;
    int i;

    if (lb >= ub) {
        printf("upper bound must be greater than lower bound");
        exit(1);
    }

    t[0] = lb;
    t[n - 1] = ub;
    delta = (ub - lb) / (n - 1);
    for (i = 1; i < n - 1; ++i) {
        t[i] = lb + delta * i;
    }

    den = std::sqrt(cwtGamma(p + 0.5));

    if ((p + 1) % 2 == 0) {
        num = 1.0;
    } else {
        num = -1.0;
    }

    num /= den;

    //printf("\n%g\n",num);

    if (p == 1) {
        for (i = 0; i < n; ++i) {
            psi[i] = -t[i] * std::exp(-t[i] * t[i] / 2.0) * num;
        }
    } else if (p == 2) {
        for (i = 0; i < n; ++i) {
            t2 = t[i] * t[i];
            psi[i] = (-1.0 + t2) * std::exp(-t2 / 2.0) * num;
        }
    } else if (p == 3) {
        for (i = 0; i < n; ++i) {
            t2 = t[i] * t[i];
            psi[i] = t[i] * (3.0 - t2) * std::exp(-t2 / 2.0) * num;
        }
    } else if (p == 4) {
        for (i = 0; i < n; ++i) {
            t2 = t[i] * t[i];
            psi[i] = (t2 * t2 - 6.0 * t2 + 3.0) * std::exp(-t2 / 2.0) * num;
        }
    } else if (p == 5) {
        for (i = 0; i < n; ++i) {
            t2 = t[i] * t[i];
            psi[i] = t[i] * (-t2 * t2 + 10.0 * t2 - 15.0) * std::exp(-t2 / 2.0) * num;
        }
    } else if (p == 6) {
        for (i = 0; i < n; ++i) {
            t2 = t[i] * t[i];
            psi[i] = (t2 * t2 * t2 - 15.0 * t2 * t2 + 45.0 * t2 - 15.0) * std::exp(-t2 / 2.0) * num;
        }
    } else if (p == 7) {
        for (i = 0; i < n; ++i) {
            t2 = t[i] * t[i];
            psi[i] = t[i] * (-t2 * t2 * t2 + 21.0 * t2 * t2 - 105.0 * t2 + 105.0) * std::exp(-t2 / 2.0) * num;
        }
    } else if (p == 8) {
        for (i = 0; i < n; ++i) {
            t2 = t[i] * t[i];
            t4 = t2 * t2;
            psi[i] = (t4 * t4 - 28.0 * t4 * t2 + 210.0 * t4 - 420.0 * t2 + 105.0) * std::exp(-t2 / 2.0) * num;
        }
    } else if (p == 9) {
        for (i = 0; i < n; ++i) {
            t2 = t[i] * t[i];
            t4 = t2 * t2;
            psi[i] = t[i] * (-t4 * t4 + 36.0 * t4 * t2 - 378.0 * t4 + 1260.0 * t2 - 945.0) * std::exp(-t2 / 2.0) * num;
        }
    } else if (p == 10) {
        for (i = 0; i < n; ++i) {
            t2 = t[i] * t[i];
            t4 = t2 * t2;
            psi[i] = (t4 * t4 * t2 - 45.0 * t4 * t4 + 630.0 * t4 * t2 - 3150.0 * t4 + 4725.0 * t2 - 945.0) * std::exp(-t2 / 2.0) * num;
        }
    } else {
        printf("\n The Gaussian Derivative Wavelet is only available for Derivatives 1 to 10");
        exit(1);
    }
}

auto mexhat(int n, double lb, double ub, double* psi, double* t) -> void
{
    gauss(n, 2, lb, ub, psi, t);
}

auto morlet(int n, double lb, double ub, double* psi, double* t) -> void
{
    int i;
    double delta;

    if (lb >= ub) {
        printf("upper bound must be greater than lower bound");
        exit(1);
    }

    t[0] = lb;
    t[n - 1] = ub;
    delta = (ub - lb) / (n - 1);
    for (i = 1; i < n - 1; ++i) {
        t[i] = lb + delta * i;
    }

    for (i = 0; i < n; ++i) {
        psi[i] = std::exp(-t[i] * t[i] / 2.0) * std::cos(5 * t[i]);
    }
}