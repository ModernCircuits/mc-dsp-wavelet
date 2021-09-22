/*
 * conv.c
 *
 *  Created on: May 1, 2013
 *      Author: Rafat Hussain
 */

#include "conv.h"

#include <algorithm>
#include <memory>

namespace {
[[nodiscard]] auto factorf(int m) -> int
{
    int n;
    n = m;
    while (n % 7 == 0) {
        n = n / 7;
    }
    while (n % 3 == 0) {
        n = n / 3;
    }
    while (n % 5 == 0) {
        n = n / 5;
    }
    while (n % 2 == 0) {
        n = n / 2;
    }

    return n;
}

[[nodiscard]] auto findnexte(int m) -> int
{
    int n;
    n = m;

    while (factorf(n) != 1 || n % 2 != 0) {
        ++n;
    }

    return n;
}
}

auto convInit(int n, int l) -> std::unique_ptr<Convolution>
{
    auto obj = std::make_unique<Convolution>();

    auto const convLen = n + l - 1;
    obj->clen = findnexte(convLen);
    obj->ilen1 = n;
    obj->ilen2 = l;

    obj->fobj = fftRealInit(obj->clen, 1);
    obj->iobj = fftRealInit(obj->clen, -1);

    return obj;
}

auto convDirect(fft_type const* inp1, int n, fft_type const* inp2, int l, fft_type* oup) -> void
{

    int k;
    int m;
    fft_type t1;
    fft_type tmin;

    auto const mm = n + l - 1;
    auto i = 0;

    if (n >= l) {

        for (k = 0; k < l; k++) {
            oup[k] = 0.0;
            for (m = 0; m <= k; m++) {
                oup[k] += inp1[m] * inp2[k - m];
            }
        }

        for (k = l; k < mm; k++) {
            oup[k] = 0.0;
            i++;
            t1 = l + i;
            tmin = std::min<double>(t1, n);
            for (m = i; m < tmin; m++) {
                oup[k] += inp1[m] * inp2[k - m];
            }
        }

    } else {
        for (k = 0; k < n; k++) {
            oup[k] = 0.0;
            for (m = 0; m <= k; m++) {
                oup[k] += inp2[m] * inp1[k - m];
            }
        }

        for (k = n; k < mm; k++) {
            oup[k] = 0.0;
            i++;
            t1 = n + i;
            tmin = std::min<double>(t1, l);
            for (m = i; m < tmin; m++) {
                oup[k] += inp2[m] * inp1[k - m];
            }
        }
    }
}

auto convFft(Convolution const& obj, fft_type const* inp1, fft_type const* inp2, fft_type* oup) -> void
{

    auto n = obj.clen;
    auto l1 = obj.ilen1;
    auto l2 = obj.ilen2;
    auto ls = l1 + l2 - 1;

    auto a = std::make_unique<fft_type[]>(n);
    auto b = std::make_unique<fft_type[]>(n);
    auto c = std::make_unique<FftData[]>(n);
    auto ao = std::make_unique<FftData[]>(n);
    auto bo = std::make_unique<FftData[]>(n);
    auto co = std::make_unique<fft_type[]>(n);

    for (auto i = 0; i < n; i++) {
        if (i < l1) {
            a[i] = inp1[i];
        } else {
            a[i] = 0.0;
        }

        if (i < l2) {
            b[i] = inp2[i];
        } else {
            b[i] = 0.0;
        }
    }

    fftR2cExec(obj.fobj.get(), a.get(), ao.get());
    fftR2cExec(obj.fobj.get(), b.get(), bo.get());

    for (auto i = 0; i < n; i++) {
        c[i].re = ao[i].re * bo[i].re - ao[i].im * bo[i].im;
        c[i].im = ao[i].im * bo[i].re + ao[i].re * bo[i].im;
    }

    fftC2rExec(obj.iobj.get(), c.get(), co.get());

    for (auto i = 0; i < ls; i++) {
        oup[i] = co[i] / n;
    }
}
