// SPDX-License-Identifier: BSL-1.0

#include "common.hpp"

#include <mc/dsp/fft/fft.hpp>

#include <mc/core/cmath.hpp>
#include <mc/core/cstdlib.hpp>
#include <mc/core/cstring.hpp>

namespace mc::dsp {

auto dwtPerStride(
    float const* inp,
    int n,
    float const* lpd,
    float const* hpd,
    int lpdLen,
    float* cA,
    int lenCA,
    float* cD,
    int istride,
    int ostride
) -> void
{

    auto const l2    = lpdLen / 2;
    auto const isodd = n % 2;

    for (auto i = 0; i < lenCA; ++i) {
        auto const t  = 2 * i + l2;
        auto const os = i * ostride;

        cA[os] = 0.0F;
        cD[os] = 0.0F;

        for (auto l = 0; l < lpdLen; ++l) {
            if ((t - l) >= l2 && (t - l) < n) {
                auto const is = (t - l) * istride;
                cA[os] += lpd[l] * inp[is];
                cD[os] += hpd[l] * inp[is];
            } else if ((t - l) < l2 && (t - l) >= 0) {
                auto const is = (t - l) * istride;
                cA[os] += lpd[l] * inp[is];
                cD[os] += hpd[l] * inp[is];
            } else if ((t - l) < 0 && isodd == 0) {
                auto const is = (t - l + n) * istride;
                cA[os] += lpd[l] * inp[is];
                cD[os] += hpd[l] * inp[is];
            } else if ((t - l) < 0 && isodd == 1) {
                if ((t - l) != -1) {
                    auto const is = (t - l + n + 1) * istride;
                    cA[os] += lpd[l] * inp[is];
                    cD[os] += hpd[l] * inp[is];
                } else {
                    auto const is = (n - 1) * istride;
                    cA[os] += lpd[l] * inp[is];
                    cD[os] += hpd[l] * inp[is];
                }
            } else if ((t - l) >= n && isodd == 0) {
                auto const is = (t - l - n) * istride;
                cA[os] += lpd[l] * inp[is];
                cD[os] += hpd[l] * inp[is];
            } else if ((t - l) >= n && isodd == 1) {
                if (t - l != n) {
                    auto const is = (t - l - (n + 1)) * istride;
                    cA[os] += lpd[l] * inp[is];
                    cD[os] += hpd[l] * inp[is];
                } else {
                    auto const is = (n - 1) * istride;
                    cA[os] += lpd[l] * inp[is];
                    cD[os] += hpd[l] * inp[is];
                }
            }
        }
    }
}

auto dwtSymStride(
    float const* inp,
    int n,
    float const* lpd,
    float const* hpd,
    int lpdLen,
    float* cA,
    int lenCA,
    float* cD,
    int istride,
    int ostride
) -> void
{
    for (auto i = 0; i < lenCA; ++i) {
        auto const t  = 2 * i + 1;
        auto const os = i * ostride;

        cA[os] = 0.0F;
        cD[os] = 0.0F;

        for (auto l = 0; l < lpdLen; ++l) {
            auto const is = [&] {
                if ((t - l) >= 0 && (t - l) < n) { return (t - l) * istride; }
                if ((t - l) < 0) { return (-t + l - 1) * istride; }
                return (2 * n - t + l - 1) * istride;
            }();

            cA[os] += lpd[l] * inp[is];
            cD[os] += hpd[l] * inp[is];
        }
    }
}

auto modwtPerStride(
    int m,
    float const* inp,
    int /*N*/,
    float const* filt,
    int lpdLen,
    float* cA,
    int lenCA,
    float* cD,
    int istride,
    int ostride
) -> void
{
    int l      = 0;
    int i      = 0;
    int t      = 0;
    int lenAvg = 0;
    int is     = 0;
    int os     = 0;
    lenAvg     = lpdLen;

    for (i = 0; i < lenCA; ++i) {
        t      = i;
        os     = i * ostride;
        is     = t * istride;
        cA[os] = filt[0] * inp[is];
        cD[os] = filt[lenAvg] * inp[is];
        for (l = 1; l < lenAvg; l++) {
            t -= m;
            while (t >= lenCA) { t -= lenCA; }
            while (t < 0) { t += lenCA; }
            os = i * ostride;
            is = t * istride;
            cA[os] += filt[l] * inp[is];
            cD[os] += filt[lenAvg + l] * inp[is];
        }
    }
}

auto swtPerStride(
    int m,
    float const* inp,
    int n,
    float const* lpd,
    float const* hpd,
    int lpdLen,
    float* cA,
    int lenCA,
    float* cD,
    int istride,
    int ostride
) -> void
{
    int l      = 0;
    int l2     = 0;
    int isodd  = 0;
    int i      = 0;
    int t      = 0;
    int lenAvg = 0;
    int j      = 0;
    int is     = 0;
    int os     = 0;
    lenAvg     = m * lpdLen;
    l2         = lenAvg / 2;
    isodd      = n % 2;

    for (i = 0; i < lenCA; ++i) {
        t      = i + l2;
        os     = i * ostride;
        cA[os] = 0.0F;
        cD[os] = 0.0F;
        l      = -1;
        for (j = 0; j < lenAvg; j += m) {
            l++;
            while (j >= lenCA) { j -= lenCA; }
            if ((t - j) >= l2 && (t - j) < n) {
                is = (t - j) * istride;
                cA[os] += lpd[l] * inp[is];
                cD[os] += hpd[l] * inp[is];
            } else if ((t - j) < l2 && (t - j) >= 0) {
                is = (t - j) * istride;
                cA[os] += lpd[l] * inp[is];
                cD[os] += hpd[l] * inp[is];
            } else if ((t - j) < 0) {
                is = (t - j + n) * istride;
                cA[os] += lpd[l] * inp[is];
                cD[os] += hpd[l] * inp[is];
            } else if ((t - j) >= n && isodd == 0) {
                is = (t - j - n) * istride;
                cA[os] += lpd[l] * inp[is];
                cD[os] += hpd[l] * inp[is];
            } else if ((t - j) >= n && isodd == 1) {
                if (t - l != n) {
                    is = (t - j - (n + 1)) * istride;
                    cA[os] += lpd[l] * inp[is];
                    cD[os] += hpd[l] * inp[is];
                } else {
                    is = (n - 1) * istride;
                    cA[os] += lpd[l] * inp[is];
                    cD[os] += hpd[l] * inp[n - 1];
                }
            }
        }
    }
}

auto idwtPerStride(
    float const* cA,
    int lenCA,
    float const* cD,
    float const* lpr,
    float const* hpr,
    int lprLen,
    float* x,
    int istride,
    int ostride
) -> void
{
    int lenAvg = 0;
    int i      = 0;
    int l      = 0;
    int m      = 0;
    int n      = 0;
    int t      = 0;
    int l2     = 0;
    int is     = 0;
    int ms     = 0;
    int ns     = 0;

    lenAvg = lprLen;
    l2     = lenAvg / 2;
    m      = -2;
    n      = -1;

    for (i = 0; i < lenCA + l2 - 1; ++i) {
        m += 2;
        n += 2;
        ms    = m * ostride;
        ns    = n * ostride;
        x[ms] = 0.0F;
        x[ns] = 0.0F;
        for (l = 0; l < l2; ++l) {
            t = 2 * l;
            if ((i - l) >= 0 && (i - l) < lenCA) {
                is = (i - l) * istride;
                x[ms] += lpr[t] * cA[is] + hpr[t] * cD[is];
                x[ns] += lpr[t + 1] * cA[is] + hpr[t + 1] * cD[is];
            } else if ((i - l) >= lenCA && (i - l) < lenCA + lenAvg - 1) {
                is = (i - l - lenCA) * istride;
                x[ms] += lpr[t] * cA[is] + hpr[t] * cD[is];
                x[ns] += lpr[t + 1] * cA[is] + hpr[t + 1] * cD[is];
            } else if ((i - l) < 0 && (i - l) > -l2) {
                is = (lenCA + i - l) * istride;
                x[ms] += lpr[t] * cA[is] + hpr[t] * cD[is];
                x[ns] += lpr[t + 1] * cA[is] + hpr[t + 1] * cD[is];
            }
        }
    }
}

auto idwtSymStride(
    float const* cA,
    int lenCA,
    float const* cD,
    float const* lpr,
    float const* hpr,
    int lprLen,
    float* x,
    int istride,
    int ostride
) -> void
{
    int lenAvg = 0;
    int i      = 0;
    int l      = 0;
    int m      = 0;
    int n      = 0;
    int t      = 0;
    int v      = 0;
    int ms     = 0;
    int ns     = 0;
    int is     = 0;
    lenAvg     = lprLen;
    m          = -2;
    n          = -1;

    for (v = 0; v < lenCA; ++v) {
        i = v;
        m += 2;
        n += 2;
        ms    = m * ostride;
        ns    = n * ostride;
        x[ms] = 0.0F;
        x[ns] = 0.0F;
        for (l = 0; l < lenAvg / 2; ++l) {
            t = 2 * l;
            if ((i - l) >= 0 && (i - l) < lenCA) {
                is = (i - l) * istride;
                x[ms] += lpr[t] * cA[is] + hpr[t] * cD[is];
                x[ns] += lpr[t + 1] * cA[is] + hpr[t + 1] * cD[is];
            }
        }
    }
}

auto testSWTlength(int n, int j) -> int
{
    int ret = 0;
    int div = 0;
    int i   = 0;
    ret     = 1;

    div = 1;
    for (i = 0; i < j; ++i) { div *= 2; }

    if ((n % div) != 0) { ret = 0; }

    return ret;
}

auto maxIterations(std::size_t sigLen, std::size_t filtLen) -> std::size_t
{
    auto const temp
        = std::log(static_cast<float>(sigLen) / (static_cast<float>(filtLen) - 1.0F))
        / std::log(2.0F);
    return static_cast<std::size_t>(temp);
}
}  // namespace mc::dsp
