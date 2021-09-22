#include "cwtmath.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <memory>


static auto nsfftFd(FftSet* obj, FftData* inp, FftData* oup, double lb, double ub, double* w) -> void
{
    auto const n = obj->N;
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

    fftExec(*obj, inp, oup);

    for (auto i = 0; i < l; ++i) {
        temp1[i] = oup[i].re;
        temp2[i] = oup[i].im;
    }

    for (auto i = 0; i < n - l; ++i) {
        oup[i].re = oup[i + l].re;
        oup[i].im = oup[i + l].im;
    }

    for (auto i = 0; i < l; ++i) {
        oup[n - l + i].re = temp1[i];
        oup[n - l + i].im = temp2[i];
    }

    auto const plb = PI2 * lb;

    for (auto i = 0; i < n; ++i) {
        auto const tempr = oup[i].re;
        auto const tempi = oup[i].im;
        auto const theta = w[i] * plb;

        oup[i].re = delta * (tempr * cos(theta) + tempi * sin(theta));
        oup[i].im = delta * (tempi * cos(theta) - tempr * sin(theta));
    }
}

static auto nsfftBk(FftSet* obj, FftData* inp, FftData* oup, double lb, double ub, double* t) -> void
{

    auto const n = obj->N;
    auto const l = n / 2;

    auto const m = divideby(n, 2);

    if (m == 0) {
        printf("The Non-Standard FFT Length must be a power of 2");
        exit(1);
    }

    auto temp1 = std::make_unique<double[]>(l);
    auto temp2 = std::make_unique<double[]>(l);
    auto w = std::make_unique<double[]>(n);
    auto inpt = std::make_unique<FftData[]>(n);

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
        inpt[i].re = (inp[i].re * cos(theta) - inp[i].im * sin(theta)) / delta;
        inpt[i].im = (inp[i].im * cos(theta) + inp[i].re * sin(theta)) / delta;
    }

    for (auto i = 0; i < l; ++i) {
        temp1[i] = inpt[i].re;
        temp2[i] = inpt[i].im;
    }

    for (auto i = 0; i < n - l; ++i) {
        inpt[i].re = inpt[i + l].re;
        inpt[i].im = inpt[i + l].im;
    }

    for (auto i = 0; i < l; ++i) {
        inpt[n - l + i].re = temp1[i];
        inpt[n - l + i].im = temp2[i];
    }

    fftExec(*obj, inpt.get(), oup);

    for (auto i = 0; i < n; ++i) {
        t[i] = lb + i * delta;
    }
}

auto nsfftExec(FftSet* obj, FftData* inp, FftData* oup, double lb, double ub, double* w) -> void
{
    if (obj->sgn == 1) {
        nsfftFd(obj, inp, oup, lb, ub, w);
    } else if (obj->sgn == -1) {
        nsfftBk(obj, inp, oup, lb, ub, w);
    }
}

static auto roundTowardsZero(double x) -> double
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
