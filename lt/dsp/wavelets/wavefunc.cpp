#include "wavefunc.h"

#include <cmath>
#include <cstdio>
#include <memory>

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

    auto obj = fftInit(n, -1);
    auto w = std::make_unique<double[]>(n);
    auto phiw = std::make_unique<FftData[]>(n);
    auto psiw = std::make_unique<FftData[]>(n);
    auto oup = std::make_unique<FftData[]>(n);

    auto const delta = 2 * (ub - lb) / PI2;

    auto j = (double)n;
    j *= -1.0;

    for (auto i = 0; i < n; ++i) {
        w[i] = j / delta;
        j += 2.0;
        psiw[i].re = psiw[i].im = 0.0;
        phiw[i].re = phiw[i].im = 0.0;
    }

    for (auto i = 0; i < n; ++i) {
        auto const wf = fabs(w[i]);
        if (wf <= PI2 / 3.0) {
            phiw[i].re = 1.0;
        }
        if (wf > PI2 / 3.0 && wf <= 2 * PI2 / 3.0) {
            auto const x = (3 * wf / PI2) - 1.0;
            auto const x2 = x * x;
            auto const x3 = x2 * x;
            auto const x4 = x3 * x;
            auto const v = x4 * (35 - 84 * x + 70 * x2 - 20 * x3);
            auto const theta = v * PI2 / 4.0;
            auto const cs = cos(theta);
            auto const sn = sin(theta);

            phiw[i].re = cs;
            psiw[i].re = cos(w[i] / 2.0) * sn;
            psiw[i].im = sin(w[i] / 2.0) * sn;
        }
        if (wf > 2.0 * PI2 / 3.0 && wf <= 4 * PI2 / 3.0) {
            auto const x = (1.5 * wf / PI2) - 1.0;
            auto const x2 = x * x;
            auto const x3 = x2 * x;
            auto const x4 = x3 * x;
            auto const v = x4 * (35 - 84 * x + 70 * x2 - 20 * x3);
            auto const theta = v * PI2 / 4.0;
            auto const cs = cos(theta);

            psiw[i].re = cos(w[i] / 2.0) * cs;
            psiw[i].im = sin(w[i] / 2.0) * cs;
        }
    }

    nsfftExec(obj.get(), phiw.get(), oup.get(), lb, ub, tgrid);

    for (auto i = 0; i < n; ++i) {
        phi[i] = oup[i].re / n;
    }

    nsfftExec(obj.get(), psiw.get(), oup.get(), lb, ub, tgrid);

    for (auto i = 0; i < n; ++i) {
        psi[i] = oup[i].re / n;
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
        psi[i] = std::exp(-t[i] * t[i] / 2.0) * cos(5 * t[i]);
    }
}
