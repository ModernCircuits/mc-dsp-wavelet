/*
 * real.c
 *
 *  Created on: Apr 20, 2013
 *      Author: Rafat Hussain
 */
#include "real.h"

#include <cstdio>
#include <memory>

auto fft_real_init(int N, int sgn) -> fft_real_object
{
    fft_real_object obj = nullptr;
    fft_type theta;
    int k;

    obj = (fft_real_object)malloc(sizeof(struct fft_real_set) + sizeof(fft_data) * (N / 2));

    obj->cobj = fft_init(N / 2, sgn);

    for (k = 0; k < N / 2; ++k) {
        theta = PI2 * k / N;
        obj->twiddle2[k].re = cos(theta);
        obj->twiddle2[k].im = sin(theta);
    }

    return obj;
}

void fft_r2c_exec(fft_real_object obj, fft_type const* inp, fft_data* oup)
{
    int i;
    int N2;
    int N;
    fft_type temp1;
    fft_type temp2;
    N2 = obj->cobj->N;
    N = N2 * 2;

    auto cinp = std::make_unique<fft_data[]>(N2);
    auto coup = std::make_unique<fft_data[]>(N2);

    for (i = 0; i < N2; ++i) {
        cinp[i].re = inp[2 * i];
        cinp[i].im = inp[2 * i + 1];
    }

    fft_exec(obj->cobj, cinp.get(), coup.get());

    oup[0].re = coup[0].re + coup[0].im;
    oup[0].im = 0.0;

    for (i = 1; i < N2; ++i) {
        temp1 = coup[i].im + coup[N2 - i].im;
        temp2 = coup[N2 - i].re - coup[i].re;
        oup[i].re = (coup[i].re + coup[N2 - i].re + (temp1 * obj->twiddle2[i].re) + (temp2 * obj->twiddle2[i].im)) / 2.0;
        oup[i].im = (coup[i].im - coup[N2 - i].im + (temp2 * obj->twiddle2[i].re) - (temp1 * obj->twiddle2[i].im)) / 2.0;
    }

    oup[N2].re = coup[0].re - coup[0].im;
    oup[N2].im = 0.0;

    for (i = 1; i < N2; ++i) {
        oup[N - i].re = oup[i].re;
        oup[N - i].im = -oup[i].im;
    }
}

void fft_c2r_exec(fft_real_object obj, fft_data* inp, fft_type* oup)
{
    int i;
    int N2;
    fft_type temp1;
    fft_type temp2;
    N2 = obj->cobj->N;

    auto cinp = std::make_unique<fft_data[]>(N2);
    auto coup = std::make_unique<fft_data[]>(N2);

    for (i = 0; i < N2; ++i) {
        temp1 = -inp[i].im - inp[N2 - i].im;
        temp2 = -inp[N2 - i].re + inp[i].re;
        cinp[i].re = inp[i].re + inp[N2 - i].re + (temp1 * obj->twiddle2[i].re) - (temp2 * obj->twiddle2[i].im);
        cinp[i].im = inp[i].im - inp[N2 - i].im + (temp2 * obj->twiddle2[i].re) + (temp1 * obj->twiddle2[i].im);
    }

    fft_exec(obj->cobj, cinp.get(), coup.get());
    for (i = 0; i < N2; ++i) {
        oup[2 * i] = coup[i].re;
        oup[2 * i + 1] = coup[i].im;
    }
}

void free_real_fft(fft_real_object object)
{
    free_fft(object->cobj);
    free(object);
}
