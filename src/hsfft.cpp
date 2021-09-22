/*
 * hsfft.c
 *
 *  Created on: Apr 14, 2013
 *      Author: Rafat Hussain
 */

#include "hsfft.h"

#include <cstdio>
#include <memory>

auto fftInit(int n, int sgn) -> std::unique_ptr<FftSet>
{
    auto obj = std::unique_ptr<FftSet>();
    int twiLen { 0 };
    auto const out = dividebyN(n);

    if (out == 1) {
        obj = std::make_unique<FftSet>();
        obj->data = makeZeros<FftData>(n);
        obj->lf = factors(n, obj->factors);
        longvectorN(obj->data.get(), obj->factors, obj->lf);
        twiLen = n;
        obj->lt = 0;
    } else {
        int k;
        int m;
        k = (int)std::pow(2.0, ceil(log10(n) / log10(2.0)));

        if (k < 2 * n - 2) {
            m = k * 2;
        } else {
            m = k;
        }
        obj = std::make_unique<FftSet>();
        obj->data = makeZeros<FftData>(m);
        obj->lf = factors(m, obj->factors);
        longvectorN(obj->data.get(), obj->factors, obj->lf);
        obj->lt = 1;
        twiLen = m;
    }

    obj->N = n;
    obj->sgn = sgn;

    if (sgn == -1) {
        for (auto ct = 0; ct < twiLen; ct++) {
            (obj->data.get() + ct)->im = -(obj->data.get() + ct)->im;
        }
    }

    return obj;
}

static void mixedRadixDitRec(FftData* op, FftData* ip, const FftSet* obj, int sgn, int n, int l, int inc)
{

    auto const radix = n > 1 ? obj->factors[inc] : 0;

    if (n == 1) {

        op[0].re = ip[0].re;
        op[0].im = ip[0].im;

    } else if (n == 2) {
        fft_type tau1r;
        fft_type tau1i;
        op[0].re = ip[0].re;
        op[0].im = ip[0].im;

        op[1].re = ip[l].re;
        op[1].im = ip[l].im;

        tau1r = op[0].re;
        tau1i = op[0].im;

        op[0].re = tau1r + op[1].re;
        op[0].im = tau1i + op[1].im;

        op[1].re = tau1r - op[1].re;
        op[1].im = tau1i - op[1].im;

    } else if (n == 3) {
        fft_type tau0r;
        fft_type tau0i;
        fft_type tau1r;
        fft_type tau1i;
        fft_type tau2r;
        fft_type tau2i;
        op[0].re = ip[0].re;
        op[0].im = ip[0].im;

        op[1].re = ip[l].re;
        op[1].im = ip[l].im;

        op[2].re = ip[2 * l].re;
        op[2].im = ip[2 * l].im;

        tau0r = op[1].re + op[2].re;
        tau0i = op[1].im + op[2].im;

        tau1r = sgn * 0.86602540378 * (op[1].re - op[2].re);
        tau1i = sgn * 0.86602540378 * (op[1].im - op[2].im);

        tau2r = op[0].re - tau0r * 0.5000000000;
        tau2i = op[0].im - tau0i * 0.5000000000;

        op[0].re = tau0r + op[0].re;
        op[0].im = tau0i + op[0].im;

        op[1].re = tau2r + tau1i;
        op[1].im = tau2i - tau1r;

        op[2].re = tau2r - tau1i;
        op[2].im = tau2i + tau1r;

        return;

    } else if (n == 4) {
        fft_type tau0r;
        fft_type tau0i;
        fft_type tau1r;
        fft_type tau1i;
        fft_type tau2r;
        fft_type tau2i;
        fft_type tau3r;
        fft_type tau3i;
        op[0].re = ip[0].re;
        op[0].im = ip[0].im;

        op[1].re = ip[l].re;
        op[1].im = ip[l].im;

        op[2].re = ip[2 * l].re;
        op[2].im = ip[2 * l].im;

        op[3].re = ip[3 * l].re;
        op[3].im = ip[3 * l].im;

        tau0r = op[0].re + op[2].re;
        tau0i = op[0].im + op[2].im;

        tau1r = op[0].re - op[2].re;
        tau1i = op[0].im - op[2].im;

        tau2r = op[1].re + op[3].re;
        tau2i = op[1].im + op[3].im;

        tau3r = sgn * (op[1].re - op[3].re);
        tau3i = sgn * (op[1].im - op[3].im);

        op[0].re = tau0r + tau2r;
        op[0].im = tau0i + tau2i;

        op[1].re = tau1r + tau3i;
        op[1].im = tau1i - tau3r;

        op[2].re = tau0r - tau2r;
        op[2].im = tau0i - tau2i;

        op[3].re = tau1r - tau3i;
        op[3].im = tau1i + tau3r;

    } else if (n == 5) {
        fft_type tau0r;
        fft_type tau0i;
        fft_type tau1r;
        fft_type tau1i;
        fft_type tau2r;
        fft_type tau2i;
        fft_type tau3r;
        fft_type tau3i;
        fft_type tau4r;
        fft_type tau4i;
        fft_type tau5r;
        fft_type tau5i;
        fft_type tau6r;
        fft_type tau6i;
        fft_type c1;
        fft_type c2;
        fft_type s1;
        fft_type s2;
        op[0].re = ip[0].re;
        op[0].im = ip[0].im;

        op[1].re = ip[l].re;
        op[1].im = ip[l].im;

        op[2].re = ip[2 * l].re;
        op[2].im = ip[2 * l].im;

        op[3].re = ip[3 * l].re;
        op[3].im = ip[3 * l].im;

        op[4].re = ip[4 * l].re;
        op[4].im = ip[4 * l].im;

        c1 = 0.30901699437;
        c2 = -0.80901699437;
        s1 = 0.95105651629;
        s2 = 0.58778525229;

        tau0r = op[1].re + op[4].re;
        tau2r = op[1].re - op[4].re;
        tau0i = op[1].im + op[4].im;
        tau2i = op[1].im - op[4].im;

        tau1r = op[2].re + op[3].re;
        tau3r = op[2].re - op[3].re;
        tau1i = op[2].im + op[3].im;
        tau3i = op[2].im - op[3].im;

        tau4r = c1 * tau0r + c2 * tau1r;
        tau4i = c1 * tau0i + c2 * tau1i;

        //tau5r = sgn * ( s1 * tau2r + s2 * tau3r);
        //tau5i = sgn * ( s1 * tau2i + s2 * tau3i);

        if (sgn == 1) {
            tau5r = s1 * tau2r + s2 * tau3r;
            tau5i = s1 * tau2i + s2 * tau3i;

        } else {
            tau5r = -s1 * tau2r - s2 * tau3r;
            tau5i = -s1 * tau2i - s2 * tau3i;
        }

        tau6r = op[0].re + tau4r;
        tau6i = op[0].im + tau4i;

        op[1].re = tau6r + tau5i;
        op[1].im = tau6i - tau5r;

        op[4].re = tau6r - tau5i;
        op[4].im = tau6i + tau5r;

        tau4r = c2 * tau0r + c1 * tau1r;
        tau4i = c2 * tau0i + c1 * tau1i;

        //tau5r = sgn * ( s2 * tau2r - s1 * tau3r);
        //tau5i = sgn * ( s2 * tau2i - s1 * tau3i);

        if (sgn == 1) {
            tau5r = s2 * tau2r - s1 * tau3r;
            tau5i = s2 * tau2i - s1 * tau3i;

        } else {
            tau5r = -s2 * tau2r + s1 * tau3r;
            tau5i = -s2 * tau2i + s1 * tau3i;
        }

        tau6r = op[0].re + tau4r;
        tau6i = op[0].im + tau4i;

        op[2].re = tau6r + tau5i;
        op[2].im = tau6i - tau5r;

        op[3].re = tau6r - tau5i;
        op[3].im = tau6i + tau5r;

        op[0].re += tau0r + tau1r;
        op[0].im += tau0i + tau1i;

    } else if (n == 7) {
        fft_type tau0r;
        fft_type tau0i;
        fft_type tau1r;
        fft_type tau1i;
        fft_type tau2r;
        fft_type tau2i;
        fft_type tau3r;
        fft_type tau3i;
        fft_type tau4r;
        fft_type tau4i;
        fft_type tau5r;
        fft_type tau5i;
        fft_type tau6r;
        fft_type tau6i;
        fft_type tau7r;
        fft_type tau7i;
        fft_type c1;
        fft_type c2;
        fft_type c3;
        fft_type s1;
        fft_type s2;
        fft_type s3;
        op[0].re = ip[0].re;
        op[0].im = ip[0].im;

        op[1].re = ip[l].re;
        op[1].im = ip[l].im;

        op[2].re = ip[2 * l].re;
        op[2].im = ip[2 * l].im;

        op[3].re = ip[3 * l].re;
        op[3].im = ip[3 * l].im;

        op[4].re = ip[4 * l].re;
        op[4].im = ip[4 * l].im;

        op[5].re = ip[5 * l].re;
        op[5].im = ip[5 * l].im;

        op[6].re = ip[6 * l].re;
        op[6].im = ip[6 * l].im;

        c1 = 0.62348980185;
        c2 = -0.22252093395;
        c3 = -0.9009688679;
        s1 = 0.78183148246;
        s2 = 0.97492791218;
        s3 = 0.43388373911;

        tau0r = op[1].re + op[6].re;
        tau3r = op[1].re - op[6].re;

        tau0i = op[1].im + op[6].im;
        tau3i = op[1].im - op[6].im;

        tau1r = op[2].re + op[5].re;
        tau4r = op[2].re - op[5].re;

        tau1i = op[2].im + op[5].im;
        tau4i = op[2].im - op[5].im;

        tau2r = op[3].re + op[4].re;
        tau5r = op[3].re - op[4].re;

        tau2i = op[3].im + op[4].im;
        tau5i = op[3].im - op[4].im;

        tau6r = op[0].re + c1 * tau0r + c2 * tau1r + c3 * tau2r;
        tau6i = op[0].im + c1 * tau0i + c2 * tau1i + c3 * tau2i;

        //tau7r = sgn * ( -s1 * tau3r - s2 * tau4r - s3 * tau5r);
        //tau7i = sgn * ( -s1 * tau3i - s2 * tau4i - s3 * tau5i);

        if (sgn == 1) {
            tau7r = -s1 * tau3r - s2 * tau4r - s3 * tau5r;
            tau7i = -s1 * tau3i - s2 * tau4i - s3 * tau5i;

        } else {
            tau7r = s1 * tau3r + s2 * tau4r + s3 * tau5r;
            tau7i = s1 * tau3i + s2 * tau4i + s3 * tau5i;
        }

        op[1].re = tau6r - tau7i;
        op[6].re = tau6r + tau7i;

        op[1].im = tau6i + tau7r;
        op[6].im = tau6i - tau7r;

        tau6r = op[0].re + c2 * tau0r + c3 * tau1r + c1 * tau2r;
        tau6i = op[0].im + c2 * tau0i + c3 * tau1i + c1 * tau2i;

        //tau7r = sgn * ( -s2 * tau3r + s3 * tau4r + s1 * tau5r);
        //tau7i = sgn * ( -s2 * tau3i + s3 * tau4i + s1 * tau5i);

        if (sgn == 1) {
            tau7r = -s2 * tau3r + s3 * tau4r + s1 * tau5r;
            tau7i = -s2 * tau3i + s3 * tau4i + s1 * tau5i;
        } else {
            tau7r = s2 * tau3r - s3 * tau4r - s1 * tau5r;
            tau7i = s2 * tau3i - s3 * tau4i - s1 * tau5i;
        }

        op[2].re = tau6r - tau7i;
        op[5].re = tau6r + tau7i;
        op[2].im = tau6i + tau7r;
        op[5].im = tau6i - tau7r;

        tau6r = op[0].re + c3 * tau0r + c1 * tau1r + c2 * tau2r;
        tau6i = op[0].im + c3 * tau0i + c1 * tau1i + c2 * tau2i;

        //tau7r = sgn * ( -s3 * tau3r + s1 * tau4r - s2 * tau5r);
        //tau7i = sgn * ( -s3 * tau3i + s1 * tau4i - s2 * tau5i);

        if (sgn == 1) {
            tau7r = -s3 * tau3r + s1 * tau4r - s2 * tau5r;
            tau7i = -s3 * tau3i + s1 * tau4i - s2 * tau5i;

        } else {
            tau7r = s3 * tau3r - s1 * tau4r + s2 * tau5r;
            tau7i = s3 * tau3i - s1 * tau4i + s2 * tau5i;
        }

        op[3].re = tau6r - tau7i;
        op[4].re = tau6r + tau7i;
        op[3].im = tau6i + tau7r;
        op[4].im = tau6i - tau7r;

        op[0].re += tau0r + tau1r + tau2r;
        op[0].im += tau0i + tau1i + tau2i;

    } else if (n == 8) {
        fft_type tau0r;
        fft_type tau0i;
        fft_type tau1r;
        fft_type tau1i;
        fft_type tau2r;
        fft_type tau2i;
        fft_type tau3r;
        fft_type tau3i;
        fft_type tau4r;
        fft_type tau4i;
        fft_type tau5r;
        fft_type tau5i;
        fft_type tau6r;
        fft_type tau6i;
        fft_type tau7r;
        fft_type tau7i;
        fft_type tau8r;
        fft_type tau8i;
        fft_type tau9r;
        fft_type tau9i;
        fft_type c1;
        fft_type s1;
        fft_type temp1r;
        fft_type temp1i;
        fft_type temp2r;
        fft_type temp2i;
        op[0].re = ip[0].re;
        op[0].im = ip[0].im;

        op[1].re = ip[l].re;
        op[1].im = ip[l].im;

        op[2].re = ip[2 * l].re;
        op[2].im = ip[2 * l].im;

        op[3].re = ip[3 * l].re;
        op[3].im = ip[3 * l].im;

        op[4].re = ip[4 * l].re;
        op[4].im = ip[4 * l].im;

        op[5].re = ip[5 * l].re;
        op[5].im = ip[5 * l].im;

        op[6].re = ip[6 * l].re;
        op[6].im = ip[6 * l].im;

        op[7].re = ip[7 * l].re;
        op[7].im = ip[7 * l].im;

        c1 = 0.70710678118654752440084436210485;
        s1 = 0.70710678118654752440084436210485;

        tau0r = op[0].re + op[4].re;
        tau4r = op[0].re - op[4].re;

        tau0i = op[0].im + op[4].im;
        tau4i = op[0].im - op[4].im;

        tau1r = op[1].re + op[7].re;
        tau5r = op[1].re - op[7].re;

        tau1i = op[1].im + op[7].im;
        tau5i = op[1].im - op[7].im;

        tau2r = op[3].re + op[5].re;
        tau6r = op[3].re - op[5].re;

        tau2i = op[3].im + op[5].im;
        tau6i = op[3].im - op[5].im;

        tau3r = op[2].re + op[6].re;
        tau7r = op[2].re - op[6].re;

        tau3i = op[2].im + op[6].im;
        tau7i = op[2].im - op[6].im;

        op[0].re = tau0r + tau1r + tau2r + tau3r;
        op[0].im = tau0i + tau1i + tau2i + tau3i;

        op[4].re = tau0r - tau1r - tau2r + tau3r;
        op[4].im = tau0i - tau1i - tau2i + tau3i;

        temp1r = tau1r - tau2r;
        temp1i = tau1i - tau2i;

        temp2r = tau5r + tau6r;
        temp2i = tau5i + tau6i;

        tau8r = tau4r + c1 * temp1r;
        tau8i = tau4i + c1 * temp1i;

        //tau9r = sgn * ( -s1 * temp2r - tau7r);
        //tau9i = sgn * ( -s1 * temp2i - tau7i);

        if (sgn == 1) {
            tau9r = -s1 * temp2r - tau7r;
            tau9i = -s1 * temp2i - tau7i;

        } else {
            tau9r = s1 * temp2r + tau7r;
            tau9i = s1 * temp2i + tau7i;
        }

        op[1].re = tau8r - tau9i;
        op[1].im = tau8i + tau9r;

        op[7].re = tau8r + tau9i;
        op[7].im = tau8i - tau9r;

        tau8r = tau0r - tau3r;
        tau8i = tau0i - tau3i;

        //tau9r = sgn * ( -tau5r + tau6r);
        //tau9i = sgn * ( -tau5i + tau6i);

        if (sgn == 1) {
            tau9r = -tau5r + tau6r;
            tau9i = -tau5i + tau6i;

        } else {
            tau9r = tau5r - tau6r;
            tau9i = tau5i - tau6i;
        }

        op[2].re = tau8r - tau9i;
        op[2].im = tau8i + tau9r;

        op[6].re = tau8r + tau9i;
        op[6].im = tau8i - tau9r;

        tau8r = tau4r - c1 * temp1r;
        tau8i = tau4i - c1 * temp1i;

        //tau9r = sgn * ( -s1 * temp2r + tau7r);
        //tau9i = sgn * ( -s1 * temp2i + tau7i);

        if (sgn == 1) {
            tau9r = -s1 * temp2r + tau7r;
            tau9i = -s1 * temp2i + tau7i;

        } else {
            tau9r = s1 * temp2r - tau7r;
            tau9i = s1 * temp2i - tau7i;
        }

        op[3].re = tau8r - tau9i;
        op[3].im = tau8i + tau9r;

        op[5].re = tau8r + tau9i;
        op[5].im = tau8i - tau9r;

    } else if (radix == 2) {
        int k;
        int tkm1;
        int ind;
        fft_type wlr;
        fft_type wli;
        fft_type tau1r;
        fft_type tau1i;
        fft_type tau2r;
        fft_type tau2i;
        auto const m = n / 2;
        auto const ll = 2 * l;
        mixedRadixDitRec(op, ip, obj, sgn, m, ll, inc + 1);
        mixedRadixDitRec(op + m, ip + l, obj, sgn, m, ll, inc + 1);

        for (k = 0; k < m; k++) {
            ind = m - 1 + k;
            wlr = (obj->data.get() + ind)->re;
            wli = (obj->data.get() + ind)->im;

            tkm1 = k + m;

            tau1r = op[k].re;
            tau1i = op[k].im;

            tau2r = op[tkm1].re * wlr - op[tkm1].im * wli;
            tau2i = op[tkm1].im * wlr + op[tkm1].re * wli;

            op[k].re = tau1r + tau2r;
            op[k].im = tau1i + tau2i;

            op[tkm1].re = tau1r - tau2r;
            op[tkm1].im = tau1i - tau2i;
        }

    } else if (radix == 3) {
        int k;
        int tkm1;
        int tkm2;
        int ind;
        fft_type wlr;
        fft_type wli;
        fft_type wl2r;
        fft_type wl2i;
        fft_type tau0r;
        fft_type tau0i;
        fft_type tau1r;
        fft_type tau1i;
        fft_type tau2r;
        fft_type tau2i;
        fft_type ar;
        fft_type ai;
        fft_type br;
        fft_type bi;
        fft_type cr;
        fft_type ci;
        auto const m = n / 3;
        auto const ll = 3 * l;
        mixedRadixDitRec(op, ip, obj, sgn, m, ll, inc + 1);
        mixedRadixDitRec(op + m, ip + l, obj, sgn, m, ll, inc + 1);
        mixedRadixDitRec(op + 2 * m, ip + 2 * l, obj, sgn, m, ll, inc + 1);
        //mixed_radix3_dit_rec(op,ip,obj,sgn,ll,m);

        for (k = 0; k < m; ++k) {
            ind = m - 1 + 2 * k;
            wlr = (obj->data.get() + ind)->re;
            wli = (obj->data.get() + ind)->im;
            ind++;
            wl2r = (obj->data.get() + ind)->re;
            wl2i = (obj->data.get() + ind)->im;
            tkm1 = k + m;
            tkm2 = tkm1 + m;

            ar = op[k].re;
            ai = op[k].im;

            br = op[tkm1].re * wlr - op[tkm1].im * wli;
            bi = op[tkm1].im * wlr + op[tkm1].re * wli;

            cr = op[tkm2].re * wl2r - op[tkm2].im * wl2i;
            ci = op[tkm2].im * wl2r + op[tkm2].re * wl2i;

            tau0r = br + cr;
            tau0i = bi + ci;

            tau1r = sgn * 0.86602540378 * (br - cr);
            tau1i = sgn * 0.86602540378 * (bi - ci);

            tau2r = ar - tau0r * 0.5000000000;
            tau2i = ai - tau0i * 0.5000000000;

            op[k].re = ar + tau0r;
            op[k].im = ai + tau0i;

            op[tkm1].re = tau2r + tau1i;
            op[tkm1].im = tau2i - tau1r;

            op[tkm2].re = tau2r - tau1i;
            op[tkm2].im = tau2i + tau1r;
        }

    } else if (radix == 4) {
        int k;
        int tkm1;
        int tkm2;
        int tkm3;
        int ind;
        fft_type wlr;
        fft_type wli;
        fft_type wl2r;
        fft_type wl2i;
        fft_type wl3r;
        fft_type wl3i;
        fft_type tau0r;
        fft_type tau0i;
        fft_type tau1r;
        fft_type tau1i;
        fft_type tau2r;
        fft_type tau2i;
        fft_type tau3r;
        fft_type tau3i;
        fft_type ar;
        fft_type ai;
        fft_type br;
        fft_type bi;
        fft_type cr;
        fft_type ci;
        fft_type dr;
        fft_type di;
        auto const m = n / 4;
        auto const ll = 4 * l;
        mixedRadixDitRec(op, ip, obj, sgn, m, ll, inc + 1);
        mixedRadixDitRec(op + m, ip + l, obj, sgn, m, ll, inc + 1);
        mixedRadixDitRec(op + 2 * m, ip + 2 * l, obj, sgn, m, ll, inc + 1);
        mixedRadixDitRec(op + 3 * m, ip + 3 * l, obj, sgn, m, ll, inc + 1);

        //mixed_radix4_dit_rec(op,ip,obj,sgn,ll,m);

        tkm1 = m;
        tkm2 = tkm1 + m;
        tkm3 = tkm2 + m;

        ar = op[0].re;
        ai = op[0].im;

        br = op[tkm1].re;
        bi = op[tkm1].im;

        cr = op[tkm2].re;
        ci = op[tkm2].im;

        dr = op[tkm3].re;
        di = op[tkm3].im;

        tau0r = ar + cr;
        tau0i = ai + ci;

        tau1r = ar - cr;
        tau1i = ai - ci;

        tau2r = br + dr;
        tau2i = bi + di;

        tau3r = sgn * (br - dr);
        tau3i = sgn * (bi - di);

        op[0].re = tau0r + tau2r;
        op[0].im = tau0i + tau2i;

        op[tkm1].re = tau1r + tau3i;
        op[tkm1].im = tau1i - tau3r;

        op[tkm2].re = tau0r - tau2r;
        op[tkm2].im = tau0i - tau2i;

        op[tkm3].re = tau1r - tau3i;
        op[tkm3].im = tau1i + tau3r;

        for (k = 1; k < m; k++) {
            ind = m - 1 + 3 * k;
            wlr = (obj->data.get() + ind)->re;
            wli = (obj->data.get() + ind)->im;
            ind++;
            wl2r = (obj->data.get() + ind)->re;
            wl2i = (obj->data.get() + ind)->im;
            ind++;
            wl3r = (obj->data.get() + ind)->re;
            wl3i = (obj->data.get() + ind)->im;

            tkm1 = k + m;
            tkm2 = tkm1 + m;
            tkm3 = tkm2 + m;

            ar = op[k].re;
            ai = op[k].im;

            br = op[tkm1].re * wlr - op[tkm1].im * wli;
            bi = op[tkm1].im * wlr + op[tkm1].re * wli;

            cr = op[tkm2].re * wl2r - op[tkm2].im * wl2i;
            ci = op[tkm2].im * wl2r + op[tkm2].re * wl2i;

            dr = op[tkm3].re * wl3r - op[tkm3].im * wl3i;
            di = op[tkm3].im * wl3r + op[tkm3].re * wl3i;

            tau0r = ar + cr;
            tau0i = ai + ci;

            tau1r = ar - cr;
            tau1i = ai - ci;

            tau2r = br + dr;
            tau2i = bi + di;

            tau3r = sgn * (br - dr);
            tau3i = sgn * (bi - di);

            op[k].re = tau0r + tau2r;
            op[k].im = tau0i + tau2i;

            op[tkm1].re = tau1r + tau3i;
            op[tkm1].im = tau1i - tau3r;

            op[tkm2].re = tau0r - tau2r;
            op[tkm2].im = tau0i - tau2i;

            op[tkm3].re = tau1r - tau3i;
            op[tkm3].im = tau1i + tau3r;
        }

    } else if (radix == 5) {
        int k;
        int tkm1;
        int tkm2;
        int tkm3;
        int tkm4;
        int ind;
        fft_type wlr;
        fft_type wli;
        fft_type wl2r;
        fft_type wl2i;
        fft_type wl3r;
        fft_type wl3i;
        fft_type wl4r;
        fft_type wl4i;
        fft_type tau0r;
        fft_type tau0i;
        fft_type tau1r;
        fft_type tau1i;
        fft_type tau2r;
        fft_type tau2i;
        fft_type tau3r;
        fft_type tau3i;
        fft_type ar;
        fft_type ai;
        fft_type br;
        fft_type bi;
        fft_type cr;
        fft_type ci;
        fft_type dr;
        fft_type di;
        fft_type er;
        fft_type ei;
        fft_type tau4r;
        fft_type tau4i;
        fft_type tau5r;
        fft_type tau5i;
        fft_type tau6r;
        fft_type tau6i;
        fft_type c1;
        fft_type c2;
        fft_type s1;
        fft_type s2;
        auto const m = n / 5;
        auto const ll = 5 * l;
        mixedRadixDitRec(op, ip, obj, sgn, m, ll, inc + 1);
        mixedRadixDitRec(op + m, ip + l, obj, sgn, m, ll, inc + 1);
        mixedRadixDitRec(op + 2 * m, ip + 2 * l, obj, sgn, m, ll, inc + 1);
        mixedRadixDitRec(op + 3 * m, ip + 3 * l, obj, sgn, m, ll, inc + 1);
        mixedRadixDitRec(op + 4 * m, ip + 4 * l, obj, sgn, m, ll, inc + 1);
        //mixed_radix3_dit_rec(op,ip,obj,sgn,ll,m);

        c1 = 0.30901699437;
        c2 = -0.80901699437;
        s1 = 0.95105651629;
        s2 = 0.58778525229;

        tkm1 = m;
        tkm2 = tkm1 + m;
        tkm3 = tkm2 + m;
        tkm4 = tkm3 + m;

        ar = op[0].re;
        ai = op[0].im;

        br = op[tkm1].re;
        bi = op[tkm1].im;

        cr = op[tkm2].re;
        ci = op[tkm2].im;

        dr = op[tkm3].re;
        di = op[tkm3].im;

        er = op[tkm4].re;
        ei = op[tkm4].im;

        tau0r = br + er;
        tau0i = bi + ei;

        tau1r = cr + dr;
        tau1i = ci + di;

        tau2r = br - er;
        tau2i = bi - ei;

        tau3r = cr - dr;
        tau3i = ci - di;

        op[0].re = ar + tau0r + tau1r;
        op[0].im = ai + tau0i + tau1i;

        tau4r = c1 * tau0r + c2 * tau1r;
        tau4i = c1 * tau0i + c2 * tau1i;

        tau5r = sgn * (s1 * tau2r + s2 * tau3r);
        tau5i = sgn * (s1 * tau2i + s2 * tau3i);

        tau6r = ar + tau4r;
        tau6i = ai + tau4i;

        op[tkm1].re = tau6r + tau5i;
        op[tkm1].im = tau6i - tau5r;

        op[tkm4].re = tau6r - tau5i;
        op[tkm4].im = tau6i + tau5r;

        tau4r = c2 * tau0r + c1 * tau1r;
        tau4i = c2 * tau0i + c1 * tau1i;

        tau5r = sgn * (s2 * tau2r - s1 * tau3r);
        tau5i = sgn * (s2 * tau2i - s1 * tau3i);

        tau6r = ar + tau4r;
        tau6i = ai + tau4i;

        op[tkm2].re = tau6r + tau5i;
        op[tkm2].im = tau6i - tau5r;

        op[tkm3].re = tau6r - tau5i;
        op[tkm3].im = tau6i + tau5r;

        for (k = 1; k < m; k++) {
            ind = m - 1 + 4 * k;
            wlr = (obj->data.get() + ind)->re;
            wli = (obj->data.get() + ind)->im;
            ind++;
            wl2r = (obj->data.get() + ind)->re;
            wl2i = (obj->data.get() + ind)->im;
            ind++;
            wl3r = (obj->data.get() + ind)->re;
            wl3i = (obj->data.get() + ind)->im;
            ind++;
            wl4r = (obj->data.get() + ind)->re;
            wl4i = (obj->data.get() + ind)->im;

            tkm1 = k + m;
            tkm2 = tkm1 + m;
            tkm3 = tkm2 + m;
            tkm4 = tkm3 + m;

            ar = op[k].re;
            ai = op[k].im;

            br = op[tkm1].re * wlr - op[tkm1].im * wli;
            bi = op[tkm1].im * wlr + op[tkm1].re * wli;

            cr = op[tkm2].re * wl2r - op[tkm2].im * wl2i;
            ci = op[tkm2].im * wl2r + op[tkm2].re * wl2i;

            dr = op[tkm3].re * wl3r - op[tkm3].im * wl3i;
            di = op[tkm3].im * wl3r + op[tkm3].re * wl3i;

            er = op[tkm4].re * wl4r - op[tkm4].im * wl4i;
            ei = op[tkm4].im * wl4r + op[tkm4].re * wl4i;

            tau0r = br + er;
            tau0i = bi + ei;

            tau1r = cr + dr;
            tau1i = ci + di;

            tau2r = br - er;
            tau2i = bi - ei;

            tau3r = cr - dr;
            tau3i = ci - di;

            op[k].re = ar + tau0r + tau1r;
            op[k].im = ai + tau0i + tau1i;

            tau4r = c1 * tau0r + c2 * tau1r;
            tau4i = c1 * tau0i + c2 * tau1i;

            //tau5r = sgn * ( s1 * tau2r + s2 * tau3r);
            //tau5i = sgn * ( s1 * tau2i + s2 * tau3i);

            if (sgn == 1) {
                tau5r = s1 * tau2r + s2 * tau3r;
                tau5i = s1 * tau2i + s2 * tau3i;

            } else {
                tau5r = -s1 * tau2r - s2 * tau3r;
                tau5i = -s1 * tau2i - s2 * tau3i;
            }

            tau6r = ar + tau4r;
            tau6i = ai + tau4i;

            op[tkm1].re = tau6r + tau5i;
            op[tkm1].im = tau6i - tau5r;

            op[tkm4].re = tau6r - tau5i;
            op[tkm4].im = tau6i + tau5r;

            tau4r = c2 * tau0r + c1 * tau1r;
            tau4i = c2 * tau0i + c1 * tau1i;

            //tau5r = sgn * ( s2 * tau2r - s1 * tau3r);
            //tau5i = sgn * ( s2 * tau2i - s1 * tau3i);

            if (sgn == 1) {
                tau5r = s2 * tau2r - s1 * tau3r;
                tau5i = s2 * tau2i - s1 * tau3i;

            } else {
                tau5r = -s2 * tau2r + s1 * tau3r;
                tau5i = -s2 * tau2i + s1 * tau3i;
            }

            tau6r = ar + tau4r;
            tau6i = ai + tau4i;

            op[tkm2].re = tau6r + tau5i;
            op[tkm2].im = tau6i - tau5r;

            op[tkm3].re = tau6r - tau5i;
            op[tkm3].im = tau6i + tau5r;
        }

    } else if (radix == 7) {
        int k;
        int tkm1;
        int tkm2;
        int tkm3;
        int tkm4;
        int tkm5;
        int tkm6;
        int ind;
        fft_type wlr;
        fft_type wli;
        fft_type wl2r;
        fft_type wl2i;
        fft_type wl3r;
        fft_type wl3i;
        fft_type wl4r;
        fft_type wl4i;
        fft_type wl5r;
        fft_type wl5i;
        fft_type wl6r;
        fft_type wl6i;
        fft_type tau0r;
        fft_type tau0i;
        fft_type tau1r;
        fft_type tau1i;
        fft_type tau2r;
        fft_type tau2i;
        fft_type tau3r;
        fft_type tau3i;
        fft_type ar;
        fft_type ai;
        fft_type br;
        fft_type bi;
        fft_type cr;
        fft_type ci;
        fft_type dr;
        fft_type di;
        fft_type er;
        fft_type ei;
        fft_type fr;
        fft_type fi;
        fft_type gr;
        fft_type gi;
        fft_type tau4r;
        fft_type tau4i;
        fft_type tau5r;
        fft_type tau5i;
        fft_type tau6r;
        fft_type tau6i;
        fft_type tau7r;
        fft_type tau7i;
        fft_type c1;
        fft_type c2;
        fft_type c3;
        fft_type s1;
        fft_type s2;
        fft_type s3;
        auto const m = n / 7;
        auto const ll = 7 * l;
        mixedRadixDitRec(op, ip, obj, sgn, m, ll, inc + 1);
        mixedRadixDitRec(op + m, ip + l, obj, sgn, m, ll, inc + 1);
        mixedRadixDitRec(op + 2 * m, ip + 2 * l, obj, sgn, m, ll, inc + 1);
        mixedRadixDitRec(op + 3 * m, ip + 3 * l, obj, sgn, m, ll, inc + 1);
        mixedRadixDitRec(op + 4 * m, ip + 4 * l, obj, sgn, m, ll, inc + 1);
        mixedRadixDitRec(op + 5 * m, ip + 5 * l, obj, sgn, m, ll, inc + 1);
        mixedRadixDitRec(op + 6 * m, ip + 6 * l, obj, sgn, m, ll, inc + 1);
        //mixed_radix3_dit_rec(op,ip,obj,sgn,ll,m);

        c1 = 0.62348980185;
        c2 = -0.22252093395;
        c3 = -0.9009688679;
        s1 = 0.78183148246;
        s2 = 0.97492791218;
        s3 = 0.43388373911;

        tkm1 = m;
        tkm2 = tkm1 + m;
        tkm3 = tkm2 + m;
        tkm4 = tkm3 + m;
        tkm5 = tkm4 + m;
        tkm6 = tkm5 + m;

        ar = op[0].re;
        ai = op[0].im;

        br = op[tkm1].re;
        bi = op[tkm1].im;

        cr = op[tkm2].re;
        ci = op[tkm2].im;

        dr = op[tkm3].re;
        di = op[tkm3].im;

        er = op[tkm4].re;
        ei = op[tkm4].im;

        fr = op[tkm5].re;
        fi = op[tkm5].im;

        gr = op[tkm6].re;
        gi = op[tkm6].im;

        tau0r = br + gr;
        tau3r = br - gr;
        tau0i = bi + gi;
        tau3i = bi - gi;

        tau1r = cr + fr;
        tau4r = cr - fr;
        tau1i = ci + fi;
        tau4i = ci - fi;

        tau2r = dr + er;
        tau5r = dr - er;
        tau2i = di + ei;
        tau5i = di - ei;

        op[0].re = ar + tau0r + tau1r + tau2r;
        op[0].im = ai + tau0i + tau1i + tau2i;

        tau6r = ar + c1 * tau0r + c2 * tau1r + c3 * tau2r;
        tau6i = ai + c1 * tau0i + c2 * tau1i + c3 * tau2i;

        //tau7r = sgn * ( -s1 * tau3r - s2 * tau4r - s3 * tau5r);
        //tau7i = sgn * ( -s1 * tau3i - s2 * tau4i - s3 * tau5i);

        if (sgn == 1) {
            tau7r = -s1 * tau3r - s2 * tau4r - s3 * tau5r;
            tau7i = -s1 * tau3i - s2 * tau4i - s3 * tau5i;

        } else {
            tau7r = s1 * tau3r + s2 * tau4r + s3 * tau5r;
            tau7i = s1 * tau3i + s2 * tau4i + s3 * tau5i;
        }

        op[tkm1].re = tau6r - tau7i;
        op[tkm1].im = tau6i + tau7r;

        op[tkm6].re = tau6r + tau7i;
        op[tkm6].im = tau6i - tau7r;

        tau6r = ar + c2 * tau0r + c3 * tau1r + c1 * tau2r;
        tau6i = ai + c2 * tau0i + c3 * tau1i + c1 * tau2i;

        //tau7r = sgn * ( -s2 * tau3r + s3 * tau4r + s1 * tau5r);
        //tau7i = sgn * ( -s2 * tau3i + s3 * tau4i + s1 * tau5i);

        if (sgn == 1) {
            tau7r = -s2 * tau3r + s3 * tau4r + s1 * tau5r;
            tau7i = -s2 * tau3i + s3 * tau4i + s1 * tau5i;

        } else {
            tau7r = s2 * tau3r - s3 * tau4r - s1 * tau5r;
            tau7i = s2 * tau3i - s3 * tau4i - s1 * tau5i;
        }

        op[tkm2].re = tau6r - tau7i;
        op[tkm2].im = tau6i + tau7r;

        op[tkm5].re = tau6r + tau7i;
        op[tkm5].im = tau6i - tau7r;

        tau6r = ar + c3 * tau0r + c1 * tau1r + c2 * tau2r;
        tau6i = ai + c3 * tau0i + c1 * tau1i + c2 * tau2i;

        //tau7r = sgn * ( -s3 * tau3r + s1 * tau4r - s2 * tau5r);
        //tau7i = sgn * ( -s3 * tau3i + s1 * tau4i - s2 * tau5i);

        if (sgn == 1) {
            tau7r = -s3 * tau3r + s1 * tau4r - s2 * tau5r;
            tau7i = -s3 * tau3i + s1 * tau4i - s2 * tau5i;

        } else {
            tau7r = s3 * tau3r - s1 * tau4r + s2 * tau5r;
            tau7i = s3 * tau3i - s1 * tau4i + s2 * tau5i;
        }

        op[tkm3].re = tau6r - tau7i;
        op[tkm3].im = tau6i + tau7r;

        op[tkm4].re = tau6r + tau7i;
        op[tkm4].im = tau6i - tau7r;

        for (k = 1; k < m; k++) {
            ind = m - 1 + 6 * k;
            wlr = (obj->data.get() + ind)->re;
            wli = (obj->data.get() + ind)->im;
            ind++;
            wl2r = (obj->data.get() + ind)->re;
            wl2i = (obj->data.get() + ind)->im;
            ind++;
            wl3r = (obj->data.get() + ind)->re;
            wl3i = (obj->data.get() + ind)->im;
            ind++;
            wl4r = (obj->data.get() + ind)->re;
            wl4i = (obj->data.get() + ind)->im;
            ind++;
            wl5r = (obj->data.get() + ind)->re;
            wl5i = (obj->data.get() + ind)->im;
            ind++;
            wl6r = (obj->data.get() + ind)->re;
            wl6i = (obj->data.get() + ind)->im;

            tkm1 = k + m;
            tkm2 = tkm1 + m;
            tkm3 = tkm2 + m;
            tkm4 = tkm3 + m;
            tkm5 = tkm4 + m;
            tkm6 = tkm5 + m;

            ar = op[k].re;
            ai = op[k].im;

            br = op[tkm1].re * wlr - op[tkm1].im * wli;
            bi = op[tkm1].im * wlr + op[tkm1].re * wli;

            cr = op[tkm2].re * wl2r - op[tkm2].im * wl2i;
            ci = op[tkm2].im * wl2r + op[tkm2].re * wl2i;

            dr = op[tkm3].re * wl3r - op[tkm3].im * wl3i;
            di = op[tkm3].im * wl3r + op[tkm3].re * wl3i;

            er = op[tkm4].re * wl4r - op[tkm4].im * wl4i;
            ei = op[tkm4].im * wl4r + op[tkm4].re * wl4i;

            fr = op[tkm5].re * wl5r - op[tkm5].im * wl5i;
            fi = op[tkm5].im * wl5r + op[tkm5].re * wl5i;

            gr = op[tkm6].re * wl6r - op[tkm6].im * wl6i;
            gi = op[tkm6].im * wl6r + op[tkm6].re * wl6i;

            tau0r = br + gr;
            tau3r = br - gr;
            tau0i = bi + gi;
            tau3i = bi - gi;

            tau1r = cr + fr;
            tau4r = cr - fr;
            tau1i = ci + fi;
            tau4i = ci - fi;

            tau2r = dr + er;
            tau5r = dr - er;
            tau2i = di + ei;
            tau5i = di - ei;

            op[k].re = ar + tau0r + tau1r + tau2r;
            op[k].im = ai + tau0i + tau1i + tau2i;

            tau6r = ar + c1 * tau0r + c2 * tau1r + c3 * tau2r;
            tau6i = ai + c1 * tau0i + c2 * tau1i + c3 * tau2i;

            //tau7r = sgn * ( -s1 * tau3r - s2 * tau4r - s3 * tau5r);
            //tau7i = sgn * ( -s1 * tau3i - s2 * tau4i - s3 * tau5i);

            if (sgn == 1) {
                tau7r = -s1 * tau3r - s2 * tau4r - s3 * tau5r;
                tau7i = -s1 * tau3i - s2 * tau4i - s3 * tau5i;

            } else {
                tau7r = s1 * tau3r + s2 * tau4r + s3 * tau5r;
                tau7i = s1 * tau3i + s2 * tau4i + s3 * tau5i;
            }

            op[tkm1].re = tau6r - tau7i;
            op[tkm1].im = tau6i + tau7r;

            op[tkm6].re = tau6r + tau7i;
            op[tkm6].im = tau6i - tau7r;

            tau6r = ar + c2 * tau0r + c3 * tau1r + c1 * tau2r;
            tau6i = ai + c2 * tau0i + c3 * tau1i + c1 * tau2i;

            //tau7r = sgn * ( -s2 * tau3r + s3 * tau4r + s1 * tau5r);
            //tau7i = sgn * ( -s2 * tau3i + s3 * tau4i + s1 * tau5i);

            if (sgn == 1) {
                tau7r = -s2 * tau3r + s3 * tau4r + s1 * tau5r;
                tau7i = -s2 * tau3i + s3 * tau4i + s1 * tau5i;

            } else {
                tau7r = s2 * tau3r - s3 * tau4r - s1 * tau5r;
                tau7i = s2 * tau3i - s3 * tau4i - s1 * tau5i;
            }

            op[tkm2].re = tau6r - tau7i;
            op[tkm2].im = tau6i + tau7r;

            op[tkm5].re = tau6r + tau7i;
            op[tkm5].im = tau6i - tau7r;

            tau6r = ar + c3 * tau0r + c1 * tau1r + c2 * tau2r;
            tau6i = ai + c3 * tau0i + c1 * tau1i + c2 * tau2i;

            //tau7r = sgn * ( -s3 * tau3r + s1 * tau4r - s2 * tau5r);
            //tau7i = sgn * ( -s3 * tau3i + s1 * tau4i - s2 * tau5i);

            if (sgn == 1) {
                tau7r = -s3 * tau3r + s1 * tau4r - s2 * tau5r;
                tau7i = -s3 * tau3i + s1 * tau4i - s2 * tau5i;

            } else {
                tau7r = s3 * tau3r - s1 * tau4r + s2 * tau5r;
                tau7i = s3 * tau3i - s1 * tau4i + s2 * tau5i;
            }

            op[tkm3].re = tau6r - tau7i;
            op[tkm3].im = tau6i + tau7r;

            op[tkm4].re = tau6r + tau7i;
            op[tkm4].im = tau6i - tau7r;
        }

    } else if (radix == 8) {
        int k;
        int tkm1;
        int tkm2;
        int tkm3;
        int tkm4;
        int tkm5;
        int tkm6;
        int tkm7;
        int ind;
        fft_type wlr;
        fft_type wli;
        fft_type wl2r;
        fft_type wl2i;
        fft_type wl3r;
        fft_type wl3i;
        fft_type wl4r;
        fft_type wl4i;
        fft_type wl5r;
        fft_type wl5i;
        fft_type wl6r;
        fft_type wl6i;
        fft_type wl7r;
        fft_type wl7i;
        fft_type tau0r;
        fft_type tau0i;
        fft_type tau1r;
        fft_type tau1i;
        fft_type tau2r;
        fft_type tau2i;
        fft_type tau3r;
        fft_type tau3i;
        fft_type ar;
        fft_type ai;
        fft_type br;
        fft_type bi;
        fft_type cr;
        fft_type ci;
        fft_type dr;
        fft_type di;
        fft_type er;
        fft_type ei;
        fft_type fr;
        fft_type fi;
        fft_type gr;
        fft_type gi;
        fft_type hr;
        fft_type hi;
        fft_type tau4r;
        fft_type tau4i;
        fft_type tau5r;
        fft_type tau5i;
        fft_type tau6r;
        fft_type tau6i;
        fft_type tau7r;
        fft_type tau7i;
        fft_type tau8r;
        fft_type tau8i;
        fft_type tau9r;
        fft_type tau9i;
        fft_type c1;
        fft_type s1;
        fft_type temp1r;
        fft_type temp1i;
        fft_type temp2r;
        fft_type temp2i;

        auto const m = n / 8;
        auto const ll = 8 * l;
        mixedRadixDitRec(op, ip, obj, sgn, m, ll, inc + 1);
        mixedRadixDitRec(op + m, ip + l, obj, sgn, m, ll, inc + 1);
        mixedRadixDitRec(op + 2 * m, ip + 2 * l, obj, sgn, m, ll, inc + 1);
        mixedRadixDitRec(op + 3 * m, ip + 3 * l, obj, sgn, m, ll, inc + 1);
        mixedRadixDitRec(op + 4 * m, ip + 4 * l, obj, sgn, m, ll, inc + 1);
        mixedRadixDitRec(op + 5 * m, ip + 5 * l, obj, sgn, m, ll, inc + 1);
        mixedRadixDitRec(op + 6 * m, ip + 6 * l, obj, sgn, m, ll, inc + 1);
        mixedRadixDitRec(op + 7 * m, ip + 7 * l, obj, sgn, m, ll, inc + 1);
        //mixed_radix3_dit_rec(op,ip,obj,sgn,ll,m);

        c1 = 0.70710678118654752440084436210485;
        s1 = 0.70710678118654752440084436210485;

        for (k = 0; k < m; k++) {
            ind = m - 1 + 7 * k;
            wlr = (obj->data.get() + ind)->re;
            wli = (obj->data.get() + ind)->im;
            ind++;
            wl2r = (obj->data.get() + ind)->re;
            wl2i = (obj->data.get() + ind)->im;
            ind++;
            wl3r = (obj->data.get() + ind)->re;
            wl3i = (obj->data.get() + ind)->im;
            ind++;
            wl4r = (obj->data.get() + ind)->re;
            wl4i = (obj->data.get() + ind)->im;
            ind++;
            wl5r = (obj->data.get() + ind)->re;
            wl5i = (obj->data.get() + ind)->im;
            ind++;
            wl6r = (obj->data.get() + ind)->re;
            wl6i = (obj->data.get() + ind)->im;
            ind++;
            wl7r = (obj->data.get() + ind)->re;
            wl7i = (obj->data.get() + ind)->im;

            tkm1 = k + m;
            tkm2 = tkm1 + m;
            tkm3 = tkm2 + m;
            tkm4 = tkm3 + m;
            tkm5 = tkm4 + m;
            tkm6 = tkm5 + m;
            tkm7 = tkm6 + m;

            ar = op[k].re;
            ai = op[k].im;

            br = op[tkm1].re * wlr - op[tkm1].im * wli;
            bi = op[tkm1].im * wlr + op[tkm1].re * wli;

            cr = op[tkm2].re * wl2r - op[tkm2].im * wl2i;
            ci = op[tkm2].im * wl2r + op[tkm2].re * wl2i;

            dr = op[tkm3].re * wl3r - op[tkm3].im * wl3i;
            di = op[tkm3].im * wl3r + op[tkm3].re * wl3i;

            er = op[tkm4].re * wl4r - op[tkm4].im * wl4i;
            ei = op[tkm4].im * wl4r + op[tkm4].re * wl4i;

            fr = op[tkm5].re * wl5r - op[tkm5].im * wl5i;
            fi = op[tkm5].im * wl5r + op[tkm5].re * wl5i;

            gr = op[tkm6].re * wl6r - op[tkm6].im * wl6i;
            gi = op[tkm6].im * wl6r + op[tkm6].re * wl6i;

            hr = op[tkm7].re * wl7r - op[tkm7].im * wl7i;
            hi = op[tkm7].im * wl7r + op[tkm7].re * wl7i;

            tau0r = ar + er;
            tau4r = ar - er;
            tau0i = ai + ei;
            tau4i = ai - ei;

            tau1r = br + hr;
            tau5r = br - hr;
            tau1i = bi + hi;
            tau5i = bi - hi;

            tau2r = dr + fr;
            tau6r = dr - fr;
            tau6i = di - fi;
            tau2i = di + fi;

            tau3r = cr + gr;
            tau7r = cr - gr;
            tau7i = ci - gi;
            tau3i = ci + gi;

            op[k].re = tau0r + tau1r + tau2r + tau3r;
            op[k].im = tau0i + tau1i + tau2i + tau3i;

            op[tkm4].re = tau0r - tau1r - tau2r + tau3r;
            op[tkm4].im = tau0i - tau1i - tau2i + tau3i;

            temp1r = tau1r - tau2r;
            temp1i = tau1i - tau2i;

            temp2r = tau5r + tau6r;
            temp2i = tau5i + tau6i;

            tau8r = tau4r + c1 * temp1r;
            tau8i = tau4i + c1 * temp1i;

            //tau9r = sgn * ( -s1 * temp2r - tau7r);
            //tau9i = sgn * ( -s1 * temp2i - tau7i);

            if (sgn == 1) {
                tau9r = -s1 * temp2r - tau7r;
                tau9i = -s1 * temp2i - tau7i;

            } else {
                tau9r = s1 * temp2r + tau7r;
                tau9i = s1 * temp2i + tau7i;
            }

            op[tkm1].re = tau8r - tau9i;
            op[tkm1].im = tau8i + tau9r;

            op[tkm7].re = tau8r + tau9i;
            op[tkm7].im = tau8i - tau9r;

            tau8r = tau0r - tau3r;
            tau8i = tau0i - tau3i;

            //tau9r = sgn * ( -tau5r + tau6r);
            //tau9i = sgn * ( -tau5i + tau6i);

            if (sgn == 1) {
                tau9r = -tau5r + tau6r;
                tau9i = -tau5i + tau6i;

            } else {
                tau9r = tau5r - tau6r;
                tau9i = tau5i - tau6i;
            }

            op[tkm2].re = tau8r - tau9i;
            op[tkm2].im = tau8i + tau9r;

            op[tkm6].re = tau8r + tau9i;
            op[tkm6].im = tau8i - tau9r;

            tau8r = tau4r - c1 * temp1r;
            tau8i = tau4i - c1 * temp1i;

            //tau9r = sgn * ( -s1 * temp2r + tau7r);
            //tau9i = sgn * ( -s1 * temp2i + tau7i);

            if (sgn == 1) {
                tau9r = -s1 * temp2r + tau7r;
                tau9i = -s1 * temp2i + tau7i;

            } else {
                tau9r = s1 * temp2r - tau7r;
                tau9i = s1 * temp2i - tau7i;
            }

            op[tkm3].re = tau8r - tau9i;
            op[tkm3].im = tau8i + tau9r;

            op[tkm5].re = tau8r + tau9i;
            op[tkm5].im = tau8i - tau9r;
        }

    } else {
        std::puts("Should never be reached");
        std::exit(EXIT_FAILURE);
    }

    // else {
    //     int k;
    //     int i;
    //     int ind;
    //     int M;
    //     int tkm;
    //     int u;
    //     int v;
    //     int t;
    //     int tt;
    //     fft_type temp1r;
    //     fft_type temp1i;
    //     fft_type temp2r;
    //     fft_type temp2i;
    //     auto wlr = std::make_unique<fft_type[]>(radix - 1);
    //     auto wli = std::make_unique<fft_type[]>(radix - 1);
    //     auto taur = std::make_unique<fft_type[]>(radix - 1);
    //     auto taui = std::make_unique<fft_type[]>(radix - 1);
    //     auto c1 = std::make_unique<fft_type[]>(radix - 1);
    //     auto s1 = std::make_unique<fft_type[]>(radix - 1);
    //     auto yr = std::make_unique<fft_type[]>(radix);
    //     auto yi = std::make_unique<fft_type[]>(radix);

    //     auto const m = N / radix;
    //     auto const ll = radix * l;

    //     for (i = 0; i < radix; ++i) {
    //         mixed_radix_dit_rec(op + i * m, ip + i * l, obj, sgn, m, ll, inc + 1);
    //     }

    //     M = (radix - 1) / 2;

    //     for (i = 1; i < M + 1; ++i) {
    //         c1[i - 1] = cos(i * PI2 / radix);
    //         s1[i - 1] = sin(i * PI2 / radix);
    //     }

    //     for (i = 0; i < M; ++i) {
    //         s1[i + M] = -s1[M - 1 - i];
    //         c1[i + M] = c1[M - 1 - i];
    //     }

    //     for (k = 0; k < m; ++k) {
    //         ind = m - 1 + (radix - 1) * k;
    //         yr[0] = op[k].re;
    //         yi[0] = op[k].im;
    //         for (i = 0; i < radix - 1; ++i) {
    //             wlr[i] = (obj->data.get() + ind)->re;
    //             wli[i] = (obj->data.get() + ind)->im;
    //             tkm = k + (i + 1) * m;
    //             yr[i + 1] = op[tkm].re * wlr[i] - op[tkm].im * wli[i];
    //             yi[i + 1] = op[tkm].im * wlr[i] + op[tkm].re * wli[i];
    //             ind++;
    //         }

    //         for (i = 0; i < M; ++i) {
    //             taur[i] = yr[i + 1] + yr[radix - 1 - i];
    //             taui[i + M] = yi[i + 1] - yi[radix - 1 - i];
    //             taui[i] = yi[i + 1] + yi[radix - 1 - i];
    //             taur[i + M] = yr[i + 1] - yr[radix - 1 - i];
    //         }

    //         temp1r = yr[0];
    //         temp1i = yi[0];

    //         for (i = 0; i < M; ++i) {
    //             temp1r += taur[i];
    //             temp1i += taui[i];
    //         }

    //         op[k].re = temp1r;
    //         op[k].im = temp1i;

    //         for (u = 0; u < M; u++) {
    //             temp1r = yr[0];
    //             temp1i = yi[0];
    //             temp2r = 0.0;
    //             temp2i = 0.0;
    //             for (v = 0; v < M; v++) {
    //                 //int ind2 = (u+v)%M;
    //                 t = (u + 1) * (v + 1);
    //                 while (t >= radix) {
    //                     t -= radix;
    //                 }
    //                 tt = t - 1;

    //                 temp1r += c1[tt] * taur[v];
    //                 temp1i += c1[tt] * taui[v];
    //                 temp2r -= s1[tt] * taur[v + M];
    //                 temp2i -= s1[tt] * taui[v + M];
    //             }
    //             temp2r = sgn * temp2r;
    //             temp2i = sgn * temp2i;

    //             op[k + (u + 1) * m].re = temp1r - temp2i;
    //             op[k + (u + 1) * m].im = temp1i + temp2r;

    //             op[k + (radix - u - 1) * m].re = temp1r + temp2i;
    //             op[k + (radix - u - 1) * m].im = temp1i - temp2r;
    //         }
    //     }
    // }
}

static void bluesteinExp(FftData* hl, FftData* hlt, int len, int m)
{
    fft_type pi;
    fft_type theta;
    fft_type angle;
    int l2;
    int len2;
    int i;
    pi = 3.1415926535897932384626433832795;
    theta = pi / len;
    l2 = 0;
    len2 = 2 * len;

    for (i = 0; i < len; ++i) {
        angle = theta * l2;
        hlt[i].re = cos(angle);
        hlt[i].im = sin(angle);
        hl[i].re = hlt[i].re;
        hl[i].im = hlt[i].im;
        l2 += 2 * i + 1;
        while (l2 > len2) {
            l2 -= len2;
        }
    }

    for (i = len; i < m - len + 1; i++) {
        hl[i].re = 0.0;
        hl[i].im = 0.0;
    }

    for (i = m - len + 1; i < m; i++) {
        hl[i].re = hlt[m - i].re;
        hl[i].im = hlt[m - i].im;
    }
}

static void bluesteinFft(FftData* data, FftData* oup, FftSet* obj, int sgn, int n)
{

    int m;
    int ii;
    int i;
    fft_type scale;
    fft_type temp;

    obj->lt = 0;
    auto k = (int)std::pow(2.0, ceil((double)log10((double)n) / log10((double)2.0)));
    auto defLt = 1;
    auto defSgn = obj->sgn;
    auto defN = obj->N;

    if (k < 2 * n - 2) {
        m = k * 2;
    } else {
        m = k;
    }
    obj->N = m;

    auto yn = std::make_unique<FftData[]>(m);
    auto hk = std::make_unique<FftData[]>(m);
    auto tempop = std::make_unique<FftData[]>(m);
    auto yno = std::make_unique<FftData[]>(m);
    auto hlt = std::make_unique<FftData[]>(n);

    bluesteinExp(tempop.get(), hlt.get(), n, m);
    scale = 1.0 / m;

    for (ii = 0; ii < m; ++ii) {
        tempop[ii].im *= scale;
        tempop[ii].re *= scale;
    }

    //fft_set* obj = initialize_fft2(M,1);
    fftExec(*obj, tempop.get(), hk.get());

    if (sgn == 1) {
        for (i = 0; i < n; i++) {
            tempop[i].re = data[i].re * hlt[i].re + data[i].im * hlt[i].im;
            tempop[i].im = -data[i].re * hlt[i].im + data[i].im * hlt[i].re;
        }
    } else {
        for (i = 0; i < n; i++) {
            tempop[i].re = data[i].re * hlt[i].re - data[i].im * hlt[i].im;
            tempop[i].im = data[i].re * hlt[i].im + data[i].im * hlt[i].re;
        }
    }

    for (i = n; i < m; i++) {
        tempop[i].re = 0.0;
        tempop[i].im = 0.0;
    }

    fftExec(*obj, tempop.get(), yn.get());

    if (sgn == 1) {
        for (i = 0; i < m; i++) {
            temp = yn[i].re * hk[i].re - yn[i].im * hk[i].im;
            yn[i].im = yn[i].re * hk[i].im + yn[i].im * hk[i].re;
            yn[i].re = temp;
        }
    } else {
        for (i = 0; i < m; i++) {
            temp = yn[i].re * hk[i].re + yn[i].im * hk[i].im;
            yn[i].im = -yn[i].re * hk[i].im + yn[i].im * hk[i].re;
            yn[i].re = temp;
        }
    }

    //IFFT

    for (ii = 0; ii < m; ++ii) {
        (obj->data.get() + ii)->im = -(obj->data.get() + ii)->im;
    }

    obj->sgn = -1 * sgn;

    fftExec(*obj, yn.get(), yno.get());

    if (sgn == 1) {
        for (i = 0; i < n; i++) {
            oup[i].re = yno[i].re * hlt[i].re + yno[i].im * hlt[i].im;
            oup[i].im = -yno[i].re * hlt[i].im + yno[i].im * hlt[i].re;
        }
    } else {
        for (i = 0; i < n; i++) {
            oup[i].re = yno[i].re * hlt[i].re - yno[i].im * hlt[i].im;
            oup[i].im = yno[i].re * hlt[i].im + yno[i].im * hlt[i].re;
        }
    }

    obj->sgn = defSgn;
    obj->N = defN;
    obj->lt = defLt;
    for (ii = 0; ii < m; ++ii) {
        (obj->data.get() + ii)->im = -(obj->data.get() + ii)->im;
    }
}

void fftExec(FftSet& obj, FftData* inp, FftData* oup)
{
    if (obj.lt == 0) {
        //fftct_radix3_dit_rec(inp,oup,obj, obj.sgn, obj.N);
        //fftct_mixed_rec(inp,oup,obj, obj.sgn, obj.N);
        int l;
        int inc;
        int nn;
        int sgn1;
        nn = obj.N;
        sgn1 = obj.sgn;
        l = 1;
        inc = 0;
        //radix3_dit_rec(oup,inp,obj,sgn1,nn,l);
        mixedRadixDitRec(oup, inp, &obj, sgn1, nn, l, inc);
    } else if (obj.lt == 1) {
        int nn;
        int sgn1;
        nn = obj.N;
        sgn1 = obj.sgn;
        bluesteinFft(inp, oup, &obj, sgn1, nn);
    }
}

auto divideby(int m, int d) -> int
{
    while (m % d == 0) {
        m = m / d;
    }
    if (m == 1) {
        return 1;
    }
    return 0;
}

auto dividebyN(int n) -> int
{
    while (n % 53 == 0) {
        n = n / 53;
    }
    while (n % 47 == 0) {
        n = n / 47;
    }
    while (n % 43 == 0) {
        n = n / 43;
    }
    while (n % 41 == 0) {
        n = n / 41;
    }
    while (n % 37 == 0) {
        n = n / 37;
    }
    while (n % 31 == 0) {
        n = n / 31;
    }
    while (n % 29 == 0) {
        n = n / 29;
    }
    while (n % 23 == 0) {
        n = n / 23;
    }
    while (n % 17 == 0) {
        n = n / 17;
    }
    while (n % 13 == 0) {
        n = n / 13;
    }
    while (n % 11 == 0) {
        n = n / 11;
    }
    while (n % 8 == 0) {
        n = n / 8;
    }
    while (n % 7 == 0) {
        n = n / 7;
    }
    while (n % 5 == 0) {
        n = n / 5;
    }
    while (n % 4 == 0) {
        n = n / 4;
    }
    while (n % 3 == 0) {
        n = n / 3;
    }
    while (n % 2 == 0) {
        n = n / 2;
    }
    if (n == 1) {
        return 1;
    }
    return 0;
}

auto factors(int m, int* arr) -> int
{
    int i;
    int n;
    int num;
    int mult;
    int m1;
    int m2;
    i = 0;
    n = m;
    while (n % 53 == 0) {
        n = n / 53;
        arr[i] = 53;
        i++;
    }
    while (n % 47 == 0) {
        n = n / 47;
        arr[i] = 47;
        i++;
    }
    while (n % 43 == 0) {
        n = n / 43;
        arr[i] = 43;
        i++;
    }
    while (n % 41 == 0) {
        n = n / 41;
        arr[i] = 41;
        i++;
    }
    while (n % 37 == 0) {
        n = n / 37;
        arr[i] = 37;
        i++;
    }
    while (n % 31 == 0) {
        n = n / 31;
        arr[i] = 31;
        i++;
    }
    while (n % 29 == 0) {
        n = n / 29;
        arr[i] = 29;
        i++;
    }
    while (n % 23 == 0) {
        n = n / 23;
        arr[i] = 23;
        i++;
    }
    while (n % 19 == 0) {
        n = n / 19;
        arr[i] = 19;
        i++;
    }
    while (n % 17 == 0) {
        n = n / 17;
        arr[i] = 17;
        i++;
    }
    while (n % 13 == 0) {
        n = n / 13;
        arr[i] = 13;
        i++;
    }
    while (n % 11 == 0) {
        n = n / 11;
        arr[i] = 11;
        i++;
    }
    while (n % 8 == 0) {
        n = n / 8;
        arr[i] = 8;
        i++;
    }
    while (n % 7 == 0) {
        n = n / 7;
        arr[i] = 7;
        i++;
    }
    while (n % 5 == 0) {
        n = n / 5;
        arr[i] = 5;
        i++;
    }
    while (n % 4 == 0) {
        n = n / 4;
        arr[i] = 4;
        i++;
    }
    while (n % 3 == 0) {
        n = n / 3;
        arr[i] = 3;
        i++;
    }
    while (n % 2 == 0) {
        n = n / 2;
        arr[i] = 2;
        i++;
    }
    if (n > 31) {
        num = 2;

        while (n > 1) {
            mult = num * 6;
            m1 = mult - 1;
            m2 = mult + 1;
            while (n % m1 == 0) {
                arr[i] = m1;
                i++;
                n = n / m1;
            }
            while (n % m2 == 0) {
                arr[i] = m2;
                i++;
                n = n / m2;
            }
            num += 1;
        }
    }
    return i;
}

void twiddle(FftData* vec, int n, int radix)
{
    int k;
    int kl;
    fft_type theta;
    fft_type theta2;
    theta = PI2 / n;
    kl = n / radix;
    vec[0].re = 1.0;
    vec[0].im = 0.0;

    for (k = 1; k < kl; k++) {
        theta2 = theta * k;
        vec[k].re = cos(theta2);
        vec[k].im = -sin(theta2);
    }
}

void longvectorN(FftData* sig, int const* array, int tx)
{
    int l;
    int i;
    int ls;
    int ct;
    int j;
    int k;
    fft_type theta;
    l = 1;
    ct = 0;
    for (i = 0; i < tx; i++) {
        l = l * array[tx - 1 - i];
        ls = l / array[tx - 1 - i];
        theta = -1.0 * PI2 / l;
        for (j = 0; j < ls; j++) {
            for (k = 0; k < array[tx - 1 - i] - 1; k++) {
                sig[ct].re = cos((k + 1) * j * theta);
                sig[ct].im = sin((k + 1) * j * theta);
                ct++;
            }
        }
    }
}
