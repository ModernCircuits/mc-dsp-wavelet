#include "wavefunc.h"

#include "lt/dsp/fft/FFT.hpp"

#include <cmath>
#include <cstdio>
#include <memory>

namespace {

auto nsfftFd(FFT* obj, Complex<double>* inp, Complex<double>* oup, double lb, double ub, double* w) -> void
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

auto nsfftExec(FFT* obj, Complex<double>* inp, Complex<double>* oup, double lb, double ub, double* w) -> void
{
    if (obj->sgn == 1) {
        nsfftFd(obj, inp, oup, lb, ub, w);
    } else if (obj->sgn == -1) {
        nsfftBk(obj, inp, oup, lb, ub, w);
    }
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

    auto obj = std::make_unique<FFT>(n, -1);
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
