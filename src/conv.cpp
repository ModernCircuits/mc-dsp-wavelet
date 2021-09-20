/*
 * conv.c
 *
 *  Created on: May 1, 2013
 *      Author: Rafat Hussain
 */

#include "conv.h"

#include <algorithm>
#include <memory>

auto factorf(int M) -> int
{
    int N;
    N = M;
    while (N % 7 == 0) {
        N = N / 7;
    }
    while (N % 3 == 0) {
        N = N / 3;
    }
    while (N % 5 == 0) {
        N = N / 5;
    }
    while (N % 2 == 0) {
        N = N / 2;
    }

    return N;
}

auto findnext(int M) -> int
{
    int N;
    N = M;

    while (factorf(N) != 1) {
        ++N;
    }

    return N;
}

auto findnexte(int M) -> int
{
    int N;
    N = M;

    while (factorf(N) != 1 || N % 2 != 0) {
        ++N;
    }

    return N;
}

auto conv_init(int N, int L) -> conv_set*
{

    conv_set* obj = nullptr;
    int conv_len;
    conv_len = N + L - 1;

    obj = (conv_set*)malloc(sizeof(struct conv_set));

    //obj->clen = npow2(conv_len);
    //obj->clen = conv_len;
    obj->clen = findnexte(conv_len);
    obj->ilen1 = N;
    obj->ilen2 = L;

    obj->fobj = fft_real_init(obj->clen, 1);
    obj->iobj = fft_real_init(obj->clen, -1);

    return obj;
}

void conv_direct(fft_type const* inp1, int N, fft_type const* inp2, int L, fft_type* oup)
{

    int M;
    int k;
    int m;
    fft_type t1;
    fft_type tmin;

    M = N + L - 1;
    auto i = 0;

    if (N >= L) {

        for (k = 0; k < L; k++) {
            oup[k] = 0.0;
            for (m = 0; m <= k; m++) {
                oup[k] += inp1[m] * inp2[k - m];
            }
        }

        for (k = L; k < M; k++) {
            oup[k] = 0.0;
            i++;
            t1 = L + i;
            tmin = std::min<double>(t1, N);
            for (m = i; m < tmin; m++) {
                oup[k] += inp1[m] * inp2[k - m];
            }
        }

    } else {
        for (k = 0; k < N; k++) {
            oup[k] = 0.0;
            for (m = 0; m <= k; m++) {
                oup[k] += inp2[m] * inp1[k - m];
            }
        }

        for (k = N; k < M; k++) {
            oup[k] = 0.0;
            i++;
            t1 = N + i;
            tmin = std::min<double>(t1, L);
            for (m = i; m < tmin; m++) {
                oup[k] += inp2[m] * inp1[k - m];
            }
        }
    }
}

void conv_fft(const conv_set* obj, fft_type const* inp1, fft_type const* inp2, fft_type* oup)
{
    int i;
    int N;
    int L1;
    int L2;
    int ls;
    fft_type* a;
    fft_type* b;
    fft_type* co;

    N = obj->clen;
    L1 = obj->ilen1;
    L2 = obj->ilen2;
    ls = L1 + L2 - 1;

    a = (fft_type*)malloc(sizeof(fft_data) * N);
    b = (fft_type*)malloc(sizeof(fft_data) * N);
    auto c = std::make_unique<fft_data[]>(N);
    auto ao = std::make_unique<fft_data[]>(N);
    auto bo = std::make_unique<fft_data[]>(N);
    co = (fft_type*)malloc(sizeof(fft_data) * N);

    for (i = 0; i < N; i++) {
        if (i < L1) {
            a[i] = inp1[i];
        } else {
            a[i] = 0.0;
        }

        if (i < L2) {
            b[i] = inp2[i];
        } else {
            b[i] = 0.0;
        }
    }

    fft_r2c_exec(obj->fobj, a, ao.get());
    fft_r2c_exec(obj->fobj, b, bo.get());

    for (i = 0; i < N; i++) {
        c[i].re = ao[i].re * bo[i].re - ao[i].im * bo[i].im;
        c[i].im = ao[i].im * bo[i].re + ao[i].re * bo[i].im;
    }

    fft_c2r_exec(obj->iobj, c.get(), co);

    for (i = 0; i < ls; i++) {
        oup[i] = co[i] / N;
    }

    free(a);
    free(b);
    free(co);
}

void free_conv(conv_set* object)
{
    free_real_fft(object->fobj);
    free_real_fft(object->iobj);
    free(object);
}
