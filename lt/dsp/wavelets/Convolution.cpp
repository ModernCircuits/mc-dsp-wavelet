#include "Convolution.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
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

auto fftRealInit(int n, int sgn) -> std::unique_ptr<FftRealSet>
{
    auto obj = std::make_unique<FftRealSet>();
    obj->data = std::make_unique<FftData[]>(n / 2);
    obj->cobj = fftInit(n / 2, sgn);

    for (auto k = 0; k < n / 2; ++k) {
        auto const theta = PI2 * k / n;
        obj->data[k].re = cos(theta);
        obj->data[k].im = sin(theta);
    }

    return obj;
}

auto fftR2cExec(FftRealSet* obj, fft_type const* inp, FftData* oup) -> void
{
    int i;
    int n2;
    int n;
    fft_type temp1;
    fft_type temp2;
    n2 = obj->cobj->N;
    n = n2 * 2;

    auto cinp = std::make_unique<FftData[]>(n2);
    auto coup = std::make_unique<FftData[]>(n2);

    for (i = 0; i < n2; ++i) {
        cinp[i].re = inp[2 * i];
        cinp[i].im = inp[2 * i + 1];
    }

    fftExec(*(obj->cobj), cinp.get(), coup.get());

    oup[0].re = coup[0].re + coup[0].im;
    oup[0].im = 0.0;

    for (i = 1; i < n2; ++i) {
        temp1 = coup[i].im + coup[n2 - i].im;
        temp2 = coup[n2 - i].re - coup[i].re;
        oup[i].re = (coup[i].re + coup[n2 - i].re + (temp1 * obj->data[i].re) + (temp2 * obj->data[i].im)) / 2.0;
        oup[i].im = (coup[i].im - coup[n2 - i].im + (temp2 * obj->data[i].re) - (temp1 * obj->data[i].im)) / 2.0;
    }

    oup[n2].re = coup[0].re - coup[0].im;
    oup[n2].im = 0.0;

    for (i = 1; i < n2; ++i) {
        oup[n - i].re = oup[i].re;
        oup[n - i].im = -oup[i].im;
    }
}

auto fftC2rExec(FftRealSet* obj, FftData* inp, fft_type* oup) -> void
{
    int i;
    int n2;
    fft_type temp1;
    fft_type temp2;
    n2 = obj->cobj->N;

    auto cinp = std::make_unique<FftData[]>(n2);
    auto coup = std::make_unique<FftData[]>(n2);

    for (i = 0; i < n2; ++i) {
        temp1 = -inp[i].im - inp[n2 - i].im;
        temp2 = -inp[n2 - i].re + inp[i].re;
        cinp[i].re = inp[i].re + inp[n2 - i].re + (temp1 * obj->data[i].re) - (temp2 * obj->data[i].im);
        cinp[i].im = inp[i].im - inp[n2 - i].im + (temp2 * obj->data[i].re) + (temp1 * obj->data[i].im);
    }

    fftExec(*(obj->cobj), cinp.get(), coup.get());
    for (i = 0; i < n2; ++i) {
        oup[2 * i] = coup[i].re;
        oup[2 * i + 1] = coup[i].im;
    }
}
