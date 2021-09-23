#include "FFT.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <memory>

namespace {
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

auto factor(int m, int* arr) -> int
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

auto longvectorN(Complex<double>* sig, int const* array, int tx) -> void
{
    auto l = 1;
    auto ct = 0;
    for (auto i = 0; i < tx; i++) {
        l = l * array[tx - 1 - i];
        auto const ls = l / array[tx - 1 - i];
        auto const theta = -1.0 * PI2 / l;
        for (auto j = 0; j < ls; j++) {
            for (auto k = 0; k < array[tx - 1 - i] - 1; k++) {
                sig[ct].real(cos((k + 1) * j * theta));
                sig[ct].imag(sin((k + 1) * j * theta));
                ct++;
            }
        }
    }
}
}

FFT::FFT(int n, int sgn)
    : N { n }
    , sgn { sgn }
{
    int twiLen { 0 };
    auto const out = dividebyN(n);

    if (out == 1) {
        this->data = makeZeros<Complex<double>>(n);
        this->lf = factor(n, this->factors);
        longvectorN(this->data.get(), this->factors, this->lf);
        twiLen = n;
        this->lt = 0;
    } else {
        int k;
        int m;
        k = (int)std::pow(2.0, ceil(log10(n) / log10(2.0)));

        if (k < 2 * n - 2) {
            m = k * 2;
        } else {
            m = k;
        }
        this->data = makeZeros<Complex<double>>(m);
        this->lf = factor(m, this->factors);
        longvectorN(this->data.get(), this->factors, this->lf);
        this->lt = 1;
        twiLen = m;
    }

    if (sgn == -1) {
        for (auto ct = 0; ct < twiLen; ct++) {
            (this->data.get() + ct)->imag(-(this->data.get() + ct)->imag());
        }
    }
}

static auto mixedRadixDitRec(Complex<double>* op, Complex<double>* ip, const FFT* obj, int sgn, int n, int l, int inc) -> void
{

    auto const radix = n > 1 ? obj->factors[inc] : 0;

    if (n == 1) {

        op[0].real(ip[0].real());
        op[0].imag(ip[0].imag());

    } else if (n == 2) {
        double tau1r;
        double tau1i;
        op[0].real(ip[0].real());
        op[0].imag(ip[0].imag());

        op[1].real(ip[l].real());
        op[1].imag(ip[l].imag());

        tau1r = op[0].real();
        tau1i = op[0].imag();

        op[0].real(tau1r + op[1].real());
        op[0].imag(tau1i + op[1].imag());

        op[1].real(tau1r - op[1].real());
        op[1].imag(tau1i - op[1].imag());

    } else if (n == 3) {
        double tau0r;
        double tau0i;
        double tau1r;
        double tau1i;
        double tau2r;
        double tau2i;
        op[0].real(ip[0].real());
        op[0].imag(ip[0].imag());

        op[1].real(ip[l].real());
        op[1].imag(ip[l].imag());

        op[2].real(ip[2 * l].real());
        op[2].imag(ip[2 * l].imag());

        tau0r = op[1].real() + op[2].real();
        tau0i = op[1].imag() + op[2].imag();

        tau1r = sgn * 0.86602540378 * (op[1].real() - op[2].real());
        tau1i = sgn * 0.86602540378 * (op[1].imag() - op[2].imag());

        tau2r = op[0].real() - tau0r * 0.5000000000;
        tau2i = op[0].imag() - tau0i * 0.5000000000;

        op[0].real(tau0r + op[0].real());
        op[0].imag(tau0i + op[0].imag());

        op[1].real(tau2r + tau1i);
        op[1].imag(tau2i - tau1r);

        op[2].real(tau2r - tau1i);
        op[2].imag(tau2i + tau1r);

        return;

    } else if (n == 4) {
        double tau0r;
        double tau0i;
        double tau1r;
        double tau1i;
        double tau2r;
        double tau2i;
        double tau3r;
        double tau3i;
        op[0].real(ip[0].real());
        op[0].imag(ip[0].imag());

        op[1].real(ip[l].real());
        op[1].imag(ip[l].imag());

        op[2].real(ip[2 * l].real());
        op[2].imag(ip[2 * l].imag());

        op[3].real(ip[3 * l].real());
        op[3].imag(ip[3 * l].imag());

        tau0r = op[0].real() + op[2].real();
        tau0i = op[0].imag() + op[2].imag();

        tau1r = op[0].real() - op[2].real();
        tau1i = op[0].imag() - op[2].imag();

        tau2r = op[1].real() + op[3].real();
        tau2i = op[1].imag() + op[3].imag();

        tau3r = sgn * (op[1].real() - op[3].real());
        tau3i = sgn * (op[1].imag() - op[3].imag());

        op[0].real(tau0r + tau2r);
        op[0].imag(tau0i + tau2i);

        op[1].real(tau1r + tau3i);
        op[1].imag(tau1i - tau3r);

        op[2].real(tau0r - tau2r);
        op[2].imag(tau0i - tau2i);

        op[3].real(tau1r - tau3i);
        op[3].imag(tau1i + tau3r);

    } else if (n == 5) {
        double tau0r;
        double tau0i;
        double tau1r;
        double tau1i;
        double tau2r;
        double tau2i;
        double tau3r;
        double tau3i;
        double tau4r;
        double tau4i;
        double tau5r;
        double tau5i;
        double tau6r;
        double tau6i;
        double c1;
        double c2;
        double s1;
        double s2;
        op[0].real(ip[0].real());
        op[0].imag(ip[0].imag());

        op[1].real(ip[l].real());
        op[1].imag(ip[l].imag());

        op[2].real(ip[2 * l].real());
        op[2].imag(ip[2 * l].imag());

        op[3].real(ip[3 * l].real());
        op[3].imag(ip[3 * l].imag());

        op[4].real(ip[4 * l].real());
        op[4].imag(ip[4 * l].imag());

        c1 = 0.30901699437;
        c2 = -0.80901699437;
        s1 = 0.95105651629;
        s2 = 0.58778525229;

        tau0r = op[1].real() + op[4].real();
        tau2r = op[1].real() - op[4].real();
        tau0i = op[1].imag() + op[4].imag();
        tau2i = op[1].imag() - op[4].imag();

        tau1r = op[2].real() + op[3].real();
        tau3r = op[2].real() - op[3].real();
        tau1i = op[2].imag() + op[3].imag();
        tau3i = op[2].imag() - op[3].imag();

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

        tau6r = op[0].real() + tau4r;
        tau6i = op[0].imag() + tau4i;

        op[1].real(tau6r + tau5i);
        op[1].imag(tau6i - tau5r);

        op[4].real(tau6r - tau5i);
        op[4].imag(tau6i + tau5r);

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

        tau6r = op[0].real() + tau4r;
        tau6i = op[0].imag() + tau4i;

        op[2].real(tau6r + tau5i);
        op[2].imag(tau6i - tau5r);

        op[3].real(tau6r - tau5i);
        op[3].imag(tau6i + tau5r);

        op[0].real(op[0].real() + tau0r + tau1r);
        op[0].imag(op[0].imag() + tau0i + tau1i);

    } else if (n == 7) {
        double tau0r;
        double tau0i;
        double tau1r;
        double tau1i;
        double tau2r;
        double tau2i;
        double tau3r;
        double tau3i;
        double tau4r;
        double tau4i;
        double tau5r;
        double tau5i;
        double tau6r;
        double tau6i;
        double tau7r;
        double tau7i;
        double c1;
        double c2;
        double c3;
        double s1;
        double s2;
        double s3;
        op[0].real(ip[0].real());
        op[0].imag(ip[0].imag());

        op[1].real(ip[l].real());
        op[1].imag(ip[l].imag());

        op[2].real(ip[2 * l].real());
        op[2].imag(ip[2 * l].imag());

        op[3].real(ip[3 * l].real());
        op[3].imag(ip[3 * l].imag());

        op[4].real(ip[4 * l].real());
        op[4].imag(ip[4 * l].imag());

        op[5].real(ip[5 * l].real());
        op[5].imag(ip[5 * l].imag());

        op[6].real(ip[6 * l].real());
        op[6].imag(ip[6 * l].imag());

        c1 = 0.62348980185;
        c2 = -0.22252093395;
        c3 = -0.9009688679;
        s1 = 0.78183148246;
        s2 = 0.97492791218;
        s3 = 0.43388373911;

        tau0r = op[1].real() + op[6].real();
        tau3r = op[1].real() - op[6].real();

        tau0i = op[1].imag() + op[6].imag();
        tau3i = op[1].imag() - op[6].imag();

        tau1r = op[2].real() + op[5].real();
        tau4r = op[2].real() - op[5].real();

        tau1i = op[2].imag() + op[5].imag();
        tau4i = op[2].imag() - op[5].imag();

        tau2r = op[3].real() + op[4].real();
        tau5r = op[3].real() - op[4].real();

        tau2i = op[3].imag() + op[4].imag();
        tau5i = op[3].imag() - op[4].imag();

        tau6r = op[0].real() + c1 * tau0r + c2 * tau1r + c3 * tau2r;
        tau6i = op[0].imag() + c1 * tau0i + c2 * tau1i + c3 * tau2i;

        //tau7r = sgn * ( -s1 * tau3r - s2 * tau4r - s3 * tau5r);
        //tau7i = sgn * ( -s1 * tau3i - s2 * tau4i - s3 * tau5i);

        if (sgn == 1) {
            tau7r = -s1 * tau3r - s2 * tau4r - s3 * tau5r;
            tau7i = -s1 * tau3i - s2 * tau4i - s3 * tau5i;

        } else {
            tau7r = s1 * tau3r + s2 * tau4r + s3 * tau5r;
            tau7i = s1 * tau3i + s2 * tau4i + s3 * tau5i;
        }

        op[1].real(tau6r - tau7i);
        op[6].real(tau6r + tau7i);

        op[1].imag(tau6i + tau7r);
        op[6].imag(tau6i - tau7r);

        tau6r = op[0].real() + c2 * tau0r + c3 * tau1r + c1 * tau2r;
        tau6i = op[0].imag() + c2 * tau0i + c3 * tau1i + c1 * tau2i;

        //tau7r = sgn * ( -s2 * tau3r + s3 * tau4r + s1 * tau5r);
        //tau7i = sgn * ( -s2 * tau3i + s3 * tau4i + s1 * tau5i);

        if (sgn == 1) {
            tau7r = -s2 * tau3r + s3 * tau4r + s1 * tau5r;
            tau7i = -s2 * tau3i + s3 * tau4i + s1 * tau5i;
        } else {
            tau7r = s2 * tau3r - s3 * tau4r - s1 * tau5r;
            tau7i = s2 * tau3i - s3 * tau4i - s1 * tau5i;
        }

        op[2].real(tau6r - tau7i);
        op[5].real(tau6r + tau7i);
        op[2].imag(tau6i + tau7r);
        op[5].imag(tau6i - tau7r);

        tau6r = op[0].real() + c3 * tau0r + c1 * tau1r + c2 * tau2r;
        tau6i = op[0].imag() + c3 * tau0i + c1 * tau1i + c2 * tau2i;

        //tau7r = sgn * ( -s3 * tau3r + s1 * tau4r - s2 * tau5r);
        //tau7i = sgn * ( -s3 * tau3i + s1 * tau4i - s2 * tau5i);

        if (sgn == 1) {
            tau7r = -s3 * tau3r + s1 * tau4r - s2 * tau5r;
            tau7i = -s3 * tau3i + s1 * tau4i - s2 * tau5i;

        } else {
            tau7r = s3 * tau3r - s1 * tau4r + s2 * tau5r;
            tau7i = s3 * tau3i - s1 * tau4i + s2 * tau5i;
        }

        op[3].real(tau6r - tau7i);
        op[4].real(tau6r + tau7i);
        op[3].imag(tau6i + tau7r);
        op[4].imag(tau6i - tau7r);

        op[0].real(op[0].real() + tau0r + tau1r + tau2r);
        op[0].imag(op[0].imag() + tau0i + tau1i + tau2i);

    } else if (n == 8) {
        double tau0r;
        double tau0i;
        double tau1r;
        double tau1i;
        double tau2r;
        double tau2i;
        double tau3r;
        double tau3i;
        double tau4r;
        double tau4i;
        double tau5r;
        double tau5i;
        double tau6r;
        double tau6i;
        double tau7r;
        double tau7i;
        double tau8r;
        double tau8i;
        double tau9r;
        double tau9i;
        double c1;
        double s1;
        double temp1r;
        double temp1i;
        double temp2r;
        double temp2i;
        op[0].real(ip[0].real());
        op[0].imag(ip[0].imag());

        op[1].real(ip[l].real());
        op[1].imag(ip[l].imag());

        op[2].real(ip[2 * l].real());
        op[2].imag(ip[2 * l].imag());

        op[3].real(ip[3 * l].real());
        op[3].imag(ip[3 * l].imag());

        op[4].real(ip[4 * l].real());
        op[4].imag(ip[4 * l].imag());

        op[5].real(ip[5 * l].real());
        op[5].imag(ip[5 * l].imag());

        op[6].real(ip[6 * l].real());
        op[6].imag(ip[6 * l].imag());

        op[7].real(ip[7 * l].real());
        op[7].imag(ip[7 * l].imag());

        c1 = 0.70710678118654752440084436210485;
        s1 = 0.70710678118654752440084436210485;

        tau0r = op[0].real() + op[4].real();
        tau4r = op[0].real() - op[4].real();

        tau0i = op[0].imag() + op[4].imag();
        tau4i = op[0].imag() - op[4].imag();

        tau1r = op[1].real() + op[7].real();
        tau5r = op[1].real() - op[7].real();

        tau1i = op[1].imag() + op[7].imag();
        tau5i = op[1].imag() - op[7].imag();

        tau2r = op[3].real() + op[5].real();
        tau6r = op[3].real() - op[5].real();

        tau2i = op[3].imag() + op[5].imag();
        tau6i = op[3].imag() - op[5].imag();

        tau3r = op[2].real() + op[6].real();
        tau7r = op[2].real() - op[6].real();

        tau3i = op[2].imag() + op[6].imag();
        tau7i = op[2].imag() - op[6].imag();

        op[0].real(tau0r + tau1r + tau2r + tau3r);
        op[0].imag(tau0i + tau1i + tau2i + tau3i);

        op[4].real(tau0r - tau1r - tau2r + tau3r);
        op[4].imag(tau0i - tau1i - tau2i + tau3i);

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

        op[1].real(tau8r - tau9i);
        op[1].imag(tau8i + tau9r);

        op[7].real(tau8r + tau9i);
        op[7].imag(tau8i - tau9r);

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

        op[2].real(tau8r - tau9i);
        op[2].imag(tau8i + tau9r);

        op[6].real(tau8r + tau9i);
        op[6].imag(tau8i - tau9r);

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

        op[3].real(tau8r - tau9i);
        op[3].imag(tau8i + tau9r);

        op[5].real(tau8r + tau9i);
        op[5].imag(tau8i - tau9r);

    } else if (radix == 2) {
        int k;
        int tkm1;
        int ind;
        double wlr;
        double wli;
        double tau1r;
        double tau1i;
        double tau2r;
        double tau2i;
        auto const m = n / 2;
        auto const ll = 2 * l;
        mixedRadixDitRec(op, ip, obj, sgn, m, ll, inc + 1);
        mixedRadixDitRec(op + m, ip + l, obj, sgn, m, ll, inc + 1);

        for (k = 0; k < m; k++) {
            ind = m - 1 + k;
            wlr = (obj->data.get() + ind)->real();
            wli = (obj->data.get() + ind)->imag();

            tkm1 = k + m;

            tau1r = op[k].real();
            tau1i = op[k].imag();

            tau2r = op[tkm1].real() * wlr - op[tkm1].imag() * wli;
            tau2i = op[tkm1].imag() * wlr + op[tkm1].real() * wli;

            op[k].real(tau1r + tau2r);
            op[k].imag(tau1i + tau2i);

            op[tkm1].real(tau1r - tau2r);
            op[tkm1].imag(tau1i - tau2i);
        }

    } else if (radix == 3) {
        int k;
        int tkm1;
        int tkm2;
        int ind;
        double wlr;
        double wli;
        double wl2r;
        double wl2i;
        double tau0r;
        double tau0i;
        double tau1r;
        double tau1i;
        double tau2r;
        double tau2i;
        double ar;
        double ai;
        double br;
        double bi;
        double cr;
        double ci;
        auto const m = n / 3;
        auto const ll = 3 * l;
        mixedRadixDitRec(op, ip, obj, sgn, m, ll, inc + 1);
        mixedRadixDitRec(op + m, ip + l, obj, sgn, m, ll, inc + 1);
        mixedRadixDitRec(op + 2 * m, ip + 2 * l, obj, sgn, m, ll, inc + 1);
        //mixed_radix3_dit_rec(op,ip,obj,sgn,ll,m);

        for (k = 0; k < m; ++k) {
            ind = m - 1 + 2 * k;
            wlr = (obj->data.get() + ind)->real();
            wli = (obj->data.get() + ind)->imag();
            ind++;
            wl2r = (obj->data.get() + ind)->real();
            wl2i = (obj->data.get() + ind)->imag();
            tkm1 = k + m;
            tkm2 = tkm1 + m;

            ar = op[k].real();
            ai = op[k].imag();

            br = op[tkm1].real() * wlr - op[tkm1].imag() * wli;
            bi = op[tkm1].imag() * wlr + op[tkm1].real() * wli;

            cr = op[tkm2].real() * wl2r - op[tkm2].imag() * wl2i;
            ci = op[tkm2].imag() * wl2r + op[tkm2].real() * wl2i;

            tau0r = br + cr;
            tau0i = bi + ci;

            tau1r = sgn * 0.86602540378 * (br - cr);
            tau1i = sgn * 0.86602540378 * (bi - ci);

            tau2r = ar - tau0r * 0.5000000000;
            tau2i = ai - tau0i * 0.5000000000;

            op[k].real(ar + tau0r);
            op[k].imag(ai + tau0i);

            op[tkm1].real(tau2r + tau1i);
            op[tkm1].imag(tau2i - tau1r);

            op[tkm2].real(tau2r - tau1i);
            op[tkm2].imag(tau2i + tau1r);
        }

    } else if (radix == 4) {
        int k;
        int tkm1;
        int tkm2;
        int tkm3;
        int ind;
        double wlr;
        double wli;
        double wl2r;
        double wl2i;
        double wl3r;
        double wl3i;
        double tau0r;
        double tau0i;
        double tau1r;
        double tau1i;
        double tau2r;
        double tau2i;
        double tau3r;
        double tau3i;
        double ar;
        double ai;
        double br;
        double bi;
        double cr;
        double ci;
        double dr;
        double di;
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

        ar = op[0].real();
        ai = op[0].imag();

        br = op[tkm1].real();
        bi = op[tkm1].imag();

        cr = op[tkm2].real();
        ci = op[tkm2].imag();

        dr = op[tkm3].real();
        di = op[tkm3].imag();

        tau0r = ar + cr;
        tau0i = ai + ci;

        tau1r = ar - cr;
        tau1i = ai - ci;

        tau2r = br + dr;
        tau2i = bi + di;

        tau3r = sgn * (br - dr);
        tau3i = sgn * (bi - di);

        op[0].real(tau0r + tau2r);
        op[0].imag(tau0i + tau2i);

        op[tkm1].real(tau1r + tau3i);
        op[tkm1].imag(tau1i - tau3r);

        op[tkm2].real(tau0r - tau2r);
        op[tkm2].imag(tau0i - tau2i);

        op[tkm3].real(tau1r - tau3i);
        op[tkm3].imag(tau1i + tau3r);

        for (k = 1; k < m; k++) {
            ind = m - 1 + 3 * k;
            wlr = (obj->data.get() + ind)->real();
            wli = (obj->data.get() + ind)->imag();
            ind++;
            wl2r = (obj->data.get() + ind)->real();
            wl2i = (obj->data.get() + ind)->imag();
            ind++;
            wl3r = (obj->data.get() + ind)->real();
            wl3i = (obj->data.get() + ind)->imag();

            tkm1 = k + m;
            tkm2 = tkm1 + m;
            tkm3 = tkm2 + m;

            ar = op[k].real();
            ai = op[k].imag();

            br = op[tkm1].real() * wlr - op[tkm1].imag() * wli;
            bi = op[tkm1].imag() * wlr + op[tkm1].real() * wli;

            cr = op[tkm2].real() * wl2r - op[tkm2].imag() * wl2i;
            ci = op[tkm2].imag() * wl2r + op[tkm2].real() * wl2i;

            dr = op[tkm3].real() * wl3r - op[tkm3].imag() * wl3i;
            di = op[tkm3].imag() * wl3r + op[tkm3].real() * wl3i;

            tau0r = ar + cr;
            tau0i = ai + ci;

            tau1r = ar - cr;
            tau1i = ai - ci;

            tau2r = br + dr;
            tau2i = bi + di;

            tau3r = sgn * (br - dr);
            tau3i = sgn * (bi - di);

            op[k].real(tau0r + tau2r);
            op[k].imag(tau0i + tau2i);

            op[tkm1].real(tau1r + tau3i);
            op[tkm1].imag(tau1i - tau3r);

            op[tkm2].real(tau0r - tau2r);
            op[tkm2].imag(tau0i - tau2i);

            op[tkm3].real(tau1r - tau3i);
            op[tkm3].imag(tau1i + tau3r);
        }

    } else if (radix == 5) {
        int k;
        int tkm1;
        int tkm2;
        int tkm3;
        int tkm4;
        int ind;
        double wlr;
        double wli;
        double wl2r;
        double wl2i;
        double wl3r;
        double wl3i;
        double wl4r;
        double wl4i;
        double tau0r;
        double tau0i;
        double tau1r;
        double tau1i;
        double tau2r;
        double tau2i;
        double tau3r;
        double tau3i;
        double ar;
        double ai;
        double br;
        double bi;
        double cr;
        double ci;
        double dr;
        double di;
        double er;
        double ei;
        double tau4r;
        double tau4i;
        double tau5r;
        double tau5i;
        double tau6r;
        double tau6i;
        double c1;
        double c2;
        double s1;
        double s2;
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

        ar = op[0].real();
        ai = op[0].imag();

        br = op[tkm1].real();
        bi = op[tkm1].imag();

        cr = op[tkm2].real();
        ci = op[tkm2].imag();

        dr = op[tkm3].real();
        di = op[tkm3].imag();

        er = op[tkm4].real();
        ei = op[tkm4].imag();

        tau0r = br + er;
        tau0i = bi + ei;

        tau1r = cr + dr;
        tau1i = ci + di;

        tau2r = br - er;
        tau2i = bi - ei;

        tau3r = cr - dr;
        tau3i = ci - di;

        op[0].real(ar + tau0r + tau1r);
        op[0].imag(ai + tau0i + tau1i);

        tau4r = c1 * tau0r + c2 * tau1r;
        tau4i = c1 * tau0i + c2 * tau1i;

        tau5r = sgn * (s1 * tau2r + s2 * tau3r);
        tau5i = sgn * (s1 * tau2i + s2 * tau3i);

        tau6r = ar + tau4r;
        tau6i = ai + tau4i;

        op[tkm1].real(tau6r + tau5i);
        op[tkm1].imag(tau6i - tau5r);

        op[tkm4].real(tau6r - tau5i);
        op[tkm4].imag(tau6i + tau5r);

        tau4r = c2 * tau0r + c1 * tau1r;
        tau4i = c2 * tau0i + c1 * tau1i;

        tau5r = sgn * (s2 * tau2r - s1 * tau3r);
        tau5i = sgn * (s2 * tau2i - s1 * tau3i);

        tau6r = ar + tau4r;
        tau6i = ai + tau4i;

        op[tkm2].real(tau6r + tau5i);
        op[tkm2].imag(tau6i - tau5r);

        op[tkm3].real(tau6r - tau5i);
        op[tkm3].imag(tau6i + tau5r);

        for (k = 1; k < m; k++) {
            ind = m - 1 + 4 * k;
            wlr = (obj->data.get() + ind)->real();
            wli = (obj->data.get() + ind)->imag();
            ind++;
            wl2r = (obj->data.get() + ind)->real();
            wl2i = (obj->data.get() + ind)->imag();
            ind++;
            wl3r = (obj->data.get() + ind)->real();
            wl3i = (obj->data.get() + ind)->imag();
            ind++;
            wl4r = (obj->data.get() + ind)->real();
            wl4i = (obj->data.get() + ind)->imag();

            tkm1 = k + m;
            tkm2 = tkm1 + m;
            tkm3 = tkm2 + m;
            tkm4 = tkm3 + m;

            ar = op[k].real();
            ai = op[k].imag();

            br = op[tkm1].real() * wlr - op[tkm1].imag() * wli;
            bi = op[tkm1].imag() * wlr + op[tkm1].real() * wli;

            cr = op[tkm2].real() * wl2r - op[tkm2].imag() * wl2i;
            ci = op[tkm2].imag() * wl2r + op[tkm2].real() * wl2i;

            dr = op[tkm3].real() * wl3r - op[tkm3].imag() * wl3i;
            di = op[tkm3].imag() * wl3r + op[tkm3].real() * wl3i;

            er = op[tkm4].real() * wl4r - op[tkm4].imag() * wl4i;
            ei = op[tkm4].imag() * wl4r + op[tkm4].real() * wl4i;

            tau0r = br + er;
            tau0i = bi + ei;

            tau1r = cr + dr;
            tau1i = ci + di;

            tau2r = br - er;
            tau2i = bi - ei;

            tau3r = cr - dr;
            tau3i = ci - di;

            op[k].real(ar + tau0r + tau1r);
            op[k].imag(ai + tau0i + tau1i);

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

            op[tkm1].real(tau6r + tau5i);
            op[tkm1].imag(tau6i - tau5r);

            op[tkm4].real(tau6r - tau5i);
            op[tkm4].imag(tau6i + tau5r);

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

            op[tkm2].real(tau6r + tau5i);
            op[tkm2].imag(tau6i - tau5r);

            op[tkm3].real(tau6r - tau5i);
            op[tkm3].imag(tau6i + tau5r);
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
        double wlr;
        double wli;
        double wl2r;
        double wl2i;
        double wl3r;
        double wl3i;
        double wl4r;
        double wl4i;
        double wl5r;
        double wl5i;
        double wl6r;
        double wl6i;
        double tau0r;
        double tau0i;
        double tau1r;
        double tau1i;
        double tau2r;
        double tau2i;
        double tau3r;
        double tau3i;
        double ar;
        double ai;
        double br;
        double bi;
        double cr;
        double ci;
        double dr;
        double di;
        double er;
        double ei;
        double fr;
        double fi;
        double gr;
        double gi;
        double tau4r;
        double tau4i;
        double tau5r;
        double tau5i;
        double tau6r;
        double tau6i;
        double tau7r;
        double tau7i;
        double c1;
        double c2;
        double c3;
        double s1;
        double s2;
        double s3;
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

        ar = op[0].real();
        ai = op[0].imag();

        br = op[tkm1].real();
        bi = op[tkm1].imag();

        cr = op[tkm2].real();
        ci = op[tkm2].imag();

        dr = op[tkm3].real();
        di = op[tkm3].imag();

        er = op[tkm4].real();
        ei = op[tkm4].imag();

        fr = op[tkm5].real();
        fi = op[tkm5].imag();

        gr = op[tkm6].real();
        gi = op[tkm6].imag();

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

        op[0].real(ar + tau0r + tau1r + tau2r);
        op[0].imag(ai + tau0i + tau1i + tau2i);

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

        op[tkm1].real(tau6r - tau7i);
        op[tkm1].imag(tau6i + tau7r);

        op[tkm6].real(tau6r + tau7i);
        op[tkm6].imag(tau6i - tau7r);

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

        op[tkm2].real(tau6r - tau7i);
        op[tkm2].imag(tau6i + tau7r);

        op[tkm5].real(tau6r + tau7i);
        op[tkm5].imag(tau6i - tau7r);

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

        op[tkm3].real(tau6r - tau7i);
        op[tkm3].imag(tau6i + tau7r);

        op[tkm4].real(tau6r + tau7i);
        op[tkm4].imag(tau6i - tau7r);

        for (k = 1; k < m; k++) {
            ind = m - 1 + 6 * k;
            wlr = (obj->data.get() + ind)->real();
            wli = (obj->data.get() + ind)->imag();
            ind++;
            wl2r = (obj->data.get() + ind)->real();
            wl2i = (obj->data.get() + ind)->imag();
            ind++;
            wl3r = (obj->data.get() + ind)->real();
            wl3i = (obj->data.get() + ind)->imag();
            ind++;
            wl4r = (obj->data.get() + ind)->real();
            wl4i = (obj->data.get() + ind)->imag();
            ind++;
            wl5r = (obj->data.get() + ind)->real();
            wl5i = (obj->data.get() + ind)->imag();
            ind++;
            wl6r = (obj->data.get() + ind)->real();
            wl6i = (obj->data.get() + ind)->imag();

            tkm1 = k + m;
            tkm2 = tkm1 + m;
            tkm3 = tkm2 + m;
            tkm4 = tkm3 + m;
            tkm5 = tkm4 + m;
            tkm6 = tkm5 + m;

            ar = op[k].real();
            ai = op[k].imag();

            br = op[tkm1].real() * wlr - op[tkm1].imag() * wli;
            bi = op[tkm1].imag() * wlr + op[tkm1].real() * wli;

            cr = op[tkm2].real() * wl2r - op[tkm2].imag() * wl2i;
            ci = op[tkm2].imag() * wl2r + op[tkm2].real() * wl2i;

            dr = op[tkm3].real() * wl3r - op[tkm3].imag() * wl3i;
            di = op[tkm3].imag() * wl3r + op[tkm3].real() * wl3i;

            er = op[tkm4].real() * wl4r - op[tkm4].imag() * wl4i;
            ei = op[tkm4].imag() * wl4r + op[tkm4].real() * wl4i;

            fr = op[tkm5].real() * wl5r - op[tkm5].imag() * wl5i;
            fi = op[tkm5].imag() * wl5r + op[tkm5].real() * wl5i;

            gr = op[tkm6].real() * wl6r - op[tkm6].imag() * wl6i;
            gi = op[tkm6].imag() * wl6r + op[tkm6].real() * wl6i;

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

            op[k].real(ar + tau0r + tau1r + tau2r);
            op[k].imag(ai + tau0i + tau1i + tau2i);

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

            op[tkm1].real(tau6r - tau7i);
            op[tkm1].imag(tau6i + tau7r);

            op[tkm6].real(tau6r + tau7i);
            op[tkm6].imag(tau6i - tau7r);

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

            op[tkm2].real(tau6r - tau7i);
            op[tkm2].imag(tau6i + tau7r);

            op[tkm5].real(tau6r + tau7i);
            op[tkm5].imag(tau6i - tau7r);

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

            op[tkm3].real(tau6r - tau7i);
            op[tkm3].imag(tau6i + tau7r);

            op[tkm4].real(tau6r + tau7i);
            op[tkm4].imag(tau6i - tau7r);
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
        double wlr;
        double wli;
        double wl2r;
        double wl2i;
        double wl3r;
        double wl3i;
        double wl4r;
        double wl4i;
        double wl5r;
        double wl5i;
        double wl6r;
        double wl6i;
        double wl7r;
        double wl7i;
        double tau0r;
        double tau0i;
        double tau1r;
        double tau1i;
        double tau2r;
        double tau2i;
        double tau3r;
        double tau3i;
        double ar;
        double ai;
        double br;
        double bi;
        double cr;
        double ci;
        double dr;
        double di;
        double er;
        double ei;
        double fr;
        double fi;
        double gr;
        double gi;
        double hr;
        double hi;
        double tau4r;
        double tau4i;
        double tau5r;
        double tau5i;
        double tau6r;
        double tau6i;
        double tau7r;
        double tau7i;
        double tau8r;
        double tau8i;
        double tau9r;
        double tau9i;
        double c1;
        double s1;
        double temp1r;
        double temp1i;
        double temp2r;
        double temp2i;

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
            wlr = (obj->data.get() + ind)->real();
            wli = (obj->data.get() + ind)->imag();
            ind++;
            wl2r = (obj->data.get() + ind)->real();
            wl2i = (obj->data.get() + ind)->imag();
            ind++;
            wl3r = (obj->data.get() + ind)->real();
            wl3i = (obj->data.get() + ind)->imag();
            ind++;
            wl4r = (obj->data.get() + ind)->real();
            wl4i = (obj->data.get() + ind)->imag();
            ind++;
            wl5r = (obj->data.get() + ind)->real();
            wl5i = (obj->data.get() + ind)->imag();
            ind++;
            wl6r = (obj->data.get() + ind)->real();
            wl6i = (obj->data.get() + ind)->imag();
            ind++;
            wl7r = (obj->data.get() + ind)->real();
            wl7i = (obj->data.get() + ind)->imag();

            tkm1 = k + m;
            tkm2 = tkm1 + m;
            tkm3 = tkm2 + m;
            tkm4 = tkm3 + m;
            tkm5 = tkm4 + m;
            tkm6 = tkm5 + m;
            tkm7 = tkm6 + m;

            ar = op[k].real();
            ai = op[k].imag();

            br = op[tkm1].real() * wlr - op[tkm1].imag() * wli;
            bi = op[tkm1].imag() * wlr + op[tkm1].real() * wli;

            cr = op[tkm2].real() * wl2r - op[tkm2].imag() * wl2i;
            ci = op[tkm2].imag() * wl2r + op[tkm2].real() * wl2i;

            dr = op[tkm3].real() * wl3r - op[tkm3].imag() * wl3i;
            di = op[tkm3].imag() * wl3r + op[tkm3].real() * wl3i;

            er = op[tkm4].real() * wl4r - op[tkm4].imag() * wl4i;
            ei = op[tkm4].imag() * wl4r + op[tkm4].real() * wl4i;

            fr = op[tkm5].real() * wl5r - op[tkm5].imag() * wl5i;
            fi = op[tkm5].imag() * wl5r + op[tkm5].real() * wl5i;

            gr = op[tkm6].real() * wl6r - op[tkm6].imag() * wl6i;
            gi = op[tkm6].imag() * wl6r + op[tkm6].real() * wl6i;

            hr = op[tkm7].real() * wl7r - op[tkm7].imag() * wl7i;
            hi = op[tkm7].imag() * wl7r + op[tkm7].real() * wl7i;

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

            op[k].real(tau0r + tau1r + tau2r + tau3r);
            op[k].imag(tau0i + tau1i + tau2i + tau3i);

            op[tkm4].real(tau0r - tau1r - tau2r + tau3r);
            op[tkm4].imag(tau0i - tau1i - tau2i + tau3i);

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

            op[tkm1].real(tau8r - tau9i);
            op[tkm1].imag(tau8i + tau9r);

            op[tkm7].real(tau8r + tau9i);
            op[tkm7].imag(tau8i - tau9r);

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

            op[tkm2].real(tau8r - tau9i);
            op[tkm2].imag(tau8i + tau9r);

            op[tkm6].real(tau8r + tau9i);
            op[tkm6].imag(tau8i - tau9r);

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

            op[tkm3].real(tau8r - tau9i);
            op[tkm3].imag(tau8i + tau9r);

            op[tkm5].real(tau8r + tau9i);
            op[tkm5].imag(tau8i - tau9r);
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
    //     double temp1r;
    //     double temp1i;
    //     double temp2r;
    //     double temp2i;
    //     auto wlr = std::make_unique<double[]>(radix - 1);
    //     auto wli = std::make_unique<double[]>(radix - 1);
    //     auto taur = std::make_unique<double[]>(radix - 1);
    //     auto taui = std::make_unique<double[]>(radix - 1);
    //     auto c1 = std::make_unique<double[]>(radix - 1);
    //     auto s1 = std::make_unique<double[]>(radix - 1);
    //     auto yr = std::make_unique<double[]>(radix);
    //     auto yi = std::make_unique<double[]>(radix);

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
    //         yr[0] = op[k].real();
    //         yi[0] = op[k].imag();
    //         for (i = 0; i < radix - 1; ++i) {
    //             wlr[i] = (obj->data.get() + ind)->real();
    //             wli[i] = (obj->data.get() + ind)->imag();
    //             tkm = k + (i + 1) * m;
    //             yr[i + 1] = op[tkm].real() * wlr[i] - op[tkm].imag() * wli[i];
    //             yi[i + 1] = op[tkm].imag() * wlr[i] + op[tkm].real() * wli[i];
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

    //         op[k].real( temp1r);
    //         op[k].imag( temp1i);

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

    //             op[k + (u + 1) * m].real( temp1r - temp2i);
    //             op[k + (u + 1) * m].imag( temp1i + temp2r);

    //             op[k + (radix - u - 1) * m].real( temp1r + temp2i);
    //             op[k + (radix - u - 1) * m].imag( temp1i - temp2r);
    //         }
    //     }
    // }
}

static auto bluesteinExp(Complex<double>* hl, Complex<double>* hlt, int len, int m) -> void
{
    double pi;
    double theta;
    double angle;
    int l2;
    int len2;
    int i;
    pi = 3.1415926535897932384626433832795;
    theta = pi / len;
    l2 = 0;
    len2 = 2 * len;

    for (i = 0; i < len; ++i) {
        angle = theta * l2;
        hlt[i].real(cos(angle));
        hlt[i].imag(sin(angle));
        hl[i].real(hlt[i].real());
        hl[i].imag(hlt[i].imag());
        l2 += 2 * i + 1;
        while (l2 > len2) {
            l2 -= len2;
        }
    }

    for (i = len; i < m - len + 1; i++) {
        hl[i].real(0.0);
        hl[i].imag(0.0);
    }

    for (i = m - len + 1; i < m; i++) {
        hl[i].real(hlt[m - i].real());
        hl[i].imag(hlt[m - i].imag());
    }
}

static auto bluesteinFft(Complex<double>* data, Complex<double>* oup, FFT* obj, int sgn, int n) -> void
{

    int m;
    int ii;
    int i;
    double scale;
    double temp;

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

    auto yn = std::make_unique<Complex<double>[]>(m);
    auto hk = std::make_unique<Complex<double>[]>(m);
    auto tempop = std::make_unique<Complex<double>[]>(m);
    auto yno = std::make_unique<Complex<double>[]>(m);
    auto hlt = std::make_unique<Complex<double>[]>(n);

    bluesteinExp(tempop.get(), hlt.get(), n, m);
    scale = 1.0 / m;

    for (ii = 0; ii < m; ++ii) {
        tempop[ii].imag(tempop[ii].imag() * scale);
        tempop[ii].real(tempop[ii].real() * scale);
    }

    //fft_set* obj = initialize_fft2(M,1);
    obj->perform(tempop.get(), hk.get());

    if (sgn == 1) {
        for (i = 0; i < n; i++) {
            tempop[i].real(data[i].real() * hlt[i].real() + data[i].imag() * hlt[i].imag());
            tempop[i].imag(-data[i].real() * hlt[i].imag() + data[i].imag() * hlt[i].real());
        }
    } else {
        for (i = 0; i < n; i++) {
            tempop[i].real(data[i].real() * hlt[i].real() - data[i].imag() * hlt[i].imag());
            tempop[i].imag(data[i].real() * hlt[i].imag() + data[i].imag() * hlt[i].real());
        }
    }

    for (i = n; i < m; i++) {
        tempop[i].real(0.0);
        tempop[i].imag(0.0);
    }

    obj->perform(tempop.get(), yn.get());

    if (sgn == 1) {
        for (i = 0; i < m; i++) {
            temp = yn[i].real() * hk[i].real() - yn[i].imag() * hk[i].imag();
            yn[i].imag(yn[i].real() * hk[i].imag() + yn[i].imag() * hk[i].real());
            yn[i].real(temp);
        }
    } else {
        for (i = 0; i < m; i++) {
            temp = yn[i].real() * hk[i].real() + yn[i].imag() * hk[i].imag();
            yn[i].imag(-yn[i].real() * hk[i].imag() + yn[i].imag() * hk[i].real());
            yn[i].real(temp);
        }
    }

    //IFFT

    for (ii = 0; ii < m; ++ii) {
        (obj->data.get() + ii)->imag(-(obj->data.get() + ii)->imag());
    }

    obj->sgn = -1 * sgn;

    obj->perform(yn.get(), yno.get());

    if (sgn == 1) {
        for (i = 0; i < n; i++) {
            oup[i].real(yno[i].real() * hlt[i].real() + yno[i].imag() * hlt[i].imag());
            oup[i].imag(-yno[i].real() * hlt[i].imag() + yno[i].imag() * hlt[i].real());
        }
    } else {
        for (i = 0; i < n; i++) {
            oup[i].real(yno[i].real() * hlt[i].real() - yno[i].imag() * hlt[i].imag());
            oup[i].imag(yno[i].real() * hlt[i].imag() + yno[i].imag() * hlt[i].real());
        }
    }

    obj->sgn = defSgn;
    obj->N = defN;
    obj->lt = defLt;
    for (ii = 0; ii < m; ++ii) {
        (obj->data.get() + ii)->imag(-(obj->data.get() + ii)->imag());
    }
}

auto FFT::perform(Complex<double>* inp, Complex<double>* oup) -> void
{
    if (lt == 0) {
        mixedRadixDitRec(oup, inp, this, sgn, N, 1, 0);
        return;
    }

    assert(lt == 1);
    bluesteinFft(inp, oup, this, sgn, N);
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

auto fftR2cExec(FftRealSet* obj, double const* inp, Complex<double>* oup) -> void
{
    int i;
    int n2;
    int n;
    double temp1;
    double temp2;
    n2 = obj->cobj->N;
    n = n2 * 2;

    auto cinp = std::make_unique<Complex<double>[]>(n2);
    auto coup = std::make_unique<Complex<double>[]>(n2);

    for (i = 0; i < n2; ++i) {
        cinp[i].real(inp[2 * i]);
        cinp[i].imag(inp[2 * i + 1]);
    }

    obj->cobj->perform(cinp.get(), coup.get());

    oup[0].real(coup[0].real() + coup[0].imag());
    oup[0].imag(0.0);

    for (i = 1; i < n2; ++i) {
        temp1 = coup[i].imag() + coup[n2 - i].imag();
        temp2 = coup[n2 - i].real() - coup[i].real();
        oup[i].real((coup[i].real() + coup[n2 - i].real() + (temp1 * obj->data[i].real()) + (temp2 * obj->data[i].imag())) / 2.0);
        oup[i].imag((coup[i].imag() - coup[n2 - i].imag() + (temp2 * obj->data[i].real()) - (temp1 * obj->data[i].imag())) / 2.0);
    }

    oup[n2].real(coup[0].real() - coup[0].imag());
    oup[n2].imag(0.0);

    for (i = 1; i < n2; ++i) {
        oup[n - i].real(oup[i].real());
        oup[n - i].imag(-oup[i].imag());
    }
}

auto fftC2rExec(FftRealSet* obj, Complex<double>* inp, double* oup) -> void
{
    int i;
    int n2;
    double temp1;
    double temp2;
    n2 = obj->cobj->N;

    auto cinp = std::make_unique<Complex<double>[]>(n2);
    auto coup = std::make_unique<Complex<double>[]>(n2);

    for (i = 0; i < n2; ++i) {
        temp1 = -inp[i].imag() - inp[n2 - i].imag();
        temp2 = -inp[n2 - i].real() + inp[i].real();
        cinp[i].real(inp[i].real() + inp[n2 - i].real() + (temp1 * obj->data[i].real()) - (temp2 * obj->data[i].imag()));
        cinp[i].imag(inp[i].imag() - inp[n2 - i].imag() + (temp2 * obj->data[i].real()) + (temp1 * obj->data[i].imag()));
    }

    obj->cobj->perform(cinp.get(), coup.get());

    for (i = 0; i < n2; ++i) {
        oup[2 * i] = coup[i].real();
        oup[2 * i + 1] = coup[i].imag();
    }
}
