/// \copyright Copyright (c) 2018, Rafat Hussain

#include "wtmath.h"

void dwtPerStride(double const* inp, int n, double const* lpd, double const* hpd, int lpdLen, double* cA, int lenCA, double* cD, int istride, int ostride)
{
    int l;
    int l2;
    int isodd;
    int i;
    int t;
    int lenAvg;
    int is;
    int os;

    lenAvg = lpdLen;
    l2 = lenAvg / 2;
    isodd = n % 2;

    for (i = 0; i < lenCA; ++i) {
        t = 2 * i + l2;
        os = i * ostride;
        cA[os] = 0.0;
        cD[os] = 0.0;
        for (l = 0; l < lenAvg; ++l) {
            if ((t - l) >= l2 && (t - l) < n) {
                is = (t - l) * istride;
                cA[os] += lpd[l] * inp[is];
                cD[os] += hpd[l] * inp[is];
            } else if ((t - l) < l2 && (t - l) >= 0) {
                is = (t - l) * istride;
                cA[os] += lpd[l] * inp[is];
                cD[os] += hpd[l] * inp[is];
            } else if ((t - l) < 0 && isodd == 0) {
                is = (t - l + n) * istride;
                cA[os] += lpd[l] * inp[is];
                cD[os] += hpd[l] * inp[is];
            } else if ((t - l) < 0 && isodd == 1) {
                if ((t - l) != -1) {
                    is = (t - l + n + 1) * istride;
                    cA[os] += lpd[l] * inp[is];
                    cD[os] += hpd[l] * inp[is];
                } else {
                    is = (n - 1) * istride;
                    cA[os] += lpd[l] * inp[is];
                    cD[os] += hpd[l] * inp[is];
                }
            } else if ((t - l) >= n && isodd == 0) {
                is = (t - l - n) * istride;
                cA[os] += lpd[l] * inp[is];
                cD[os] += hpd[l] * inp[is];
            } else if ((t - l) >= n && isodd == 1) {
                is = (t - l - (n + 1)) * istride;
                if (t - l != n) {
                    cA[os] += lpd[l] * inp[is];
                    cD[os] += hpd[l] * inp[is];
                } else {
                    is = (n - 1) * istride;
                    cA[os] += lpd[l] * inp[is];
                    cD[os] += hpd[l] * inp[is];
                }
            }
        }
    }
}

void dwtSymStride(double const* inp, int n, double const* lpd, double const* hpd, int lpdLen, double* cA, int lenCA, double* cD, int istride, int ostride)
{
    int i;
    int l;
    int t;
    int lenAvg;
    int is;
    int os;
    lenAvg = lpdLen;

    for (i = 0; i < lenCA; ++i) {
        t = 2 * i + 1;
        os = i * ostride;
        cA[os] = 0.0;
        cD[os] = 0.0;
        for (l = 0; l < lenAvg; ++l) {
            if ((t - l) >= 0 && (t - l) < n) {
                is = (t - l) * istride;
                cA[os] += lpd[l] * inp[is];
                cD[os] += hpd[l] * inp[is];
            } else if ((t - l) < 0) {
                is = (-t + l - 1) * istride;
                cA[os] += lpd[l] * inp[is];
                cD[os] += hpd[l] * inp[is];
            } else if ((t - l) >= n) {
                is = (2 * n - t + l - 1) * istride;
                cA[os] += lpd[l] * inp[is];
                cD[os] += hpd[l] * inp[is];
            }
        }
    }
}

void modwtPerStride(int m, double const* inp, int /*N*/, double const* filt, int lpdLen, double* cA, int lenCA, double* cD, int istride, int ostride)
{
    int l;
    int i;
    int t;
    int lenAvg;
    int is;
    int os;
    lenAvg = lpdLen;

    for (i = 0; i < lenCA; ++i) {
        t = i;
        os = i * ostride;
        is = t * istride;
        cA[os] = filt[0] * inp[is];
        cD[os] = filt[lenAvg] * inp[is];
        for (l = 1; l < lenAvg; l++) {
            t -= m;
            while (t >= lenCA) {
                t -= lenCA;
            }
            while (t < 0) {
                t += lenCA;
            }
            os = i * ostride;
            is = t * istride;
            cA[os] += filt[l] * inp[is];
            cD[os] += filt[lenAvg + l] * inp[is];
        }
    }
}

void swtPerStride(int m, double const* inp, int n, double const* lpd, double const* hpd, int lpdLen, double* cA, int lenCA, double* cD, int istride, int ostride)
{
    int l;
    int l2;
    int isodd;
    int i;
    int t;
    int lenAvg;
    int j;
    int is;
    int os;
    lenAvg = m * lpdLen;
    l2 = lenAvg / 2;
    isodd = n % 2;

    for (i = 0; i < lenCA; ++i) {
        t = i + l2;
        os = i * ostride;
        cA[os] = 0.0;
        cD[os] = 0.0;
        l = -1;
        for (j = 0; j < lenAvg; j += m) {
            l++;
            while (j >= lenCA) {
                j -= lenCA;
            }
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

void idwtPerStride(double const* cA, int lenCA, double const* cD, double const* lpr, double const* hpr, int lprLen, double* x, int istride, int ostride)
{
    int lenAvg;
    int i;
    int l;
    int m;
    int n;
    int t;
    int l2;
    int is;
    int ms;
    int ns;

    lenAvg = lprLen;
    l2 = lenAvg / 2;
    m = -2;
    n = -1;

    for (i = 0; i < lenCA + l2 - 1; ++i) {
        m += 2;
        n += 2;
        ms = m * ostride;
        ns = n * ostride;
        x[ms] = 0.0;
        x[ns] = 0.0;
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

void idwtSymStride(double const* cA, int lenCA, double const* cD, double const* lpr, double const* hpr, int lprLen, double* x, int istride, int ostride)
{
    int lenAvg;
    int i;
    int l;
    int m;
    int n;
    int t;
    int v;
    int ms;
    int ns;
    int is;
    lenAvg = lprLen;
    m = -2;
    n = -1;

    for (v = 0; v < lenCA; ++v) {
        i = v;
        m += 2;
        n += 2;
        ms = m * ostride;
        ns = n * ostride;
        x[ms] = 0.0;
        x[ns] = 0.0;
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

void imodwtPerStride(int m, double const* cA, int lenCA, double const* cD, double const* filt, int lf, double* x, int istride, int ostride)
{
    int lenAvg;
    int i;
    int l;
    int t;
    int is;
    int os;

    lenAvg = lf;

    for (i = 0; i < lenCA; ++i) {
        t = i;
        os = i * ostride;
        is = t * istride;
        x[os] = (filt[0] * cA[is]) + (filt[lenAvg] * cD[is]);
        for (l = 1; l < lenAvg; l++) {
            t += m;
            while (t >= lenCA) {
                t -= lenCA;
            }
            while (t < 0) {
                t += lenCA;
            }
            is = t * istride;
            x[os] += (filt[l] * cA[is]) + (filt[lenAvg + l] * cD[is]);
        }
    }
}

void idwt2Shift(int shift, int rows, int cols, double const* lpr, double const* hpr, int lf, double* a, double* h, double* v, double* d, double* oup)
{
    auto const n = rows > cols ? 2 * rows : 2 * cols;
    auto const dim1 = 2 * rows;
    auto const dim2 = 2 * cols;

    auto xLp = makeZeros<double>(n + 2 * lf - 1);
    auto cL = makeZeros<double>(dim1 * dim2);
    auto cH = makeZeros<double>(dim1 * dim2);

    auto ir = rows;
    auto ic = cols;
    auto istride = ic;
    auto ostride = 1;

    for (auto i = 0; i < ic; ++i) {
        idwtPerStride(a + i, ir, h + i, lpr, hpr, lf, xLp.get(), istride, ostride);

        for (auto k = lf / 2 - 1; k < 2 * ir + lf / 2 - 1; ++k) {
            cL[(k - lf / 2 + 1) * ic + i] = xLp[k];
        }

        idwtPerStride(v + i, ir, d + i, lpr, hpr, lf, xLp.get(), istride, ostride);

        for (auto k = lf / 2 - 1; k < 2 * ir + lf / 2 - 1; ++k) {
            cH[(k - lf / 2 + 1) * ic + i] = xLp[k];
        }
    }

    ir *= 2;
    istride = 1;
    ostride = 1;

    for (auto i = 0; i < ir; ++i) {
        idwtPerStride(cL.get() + i * ic, ic, cH.get() + i * ic, lpr, hpr, lf, xLp.get(), istride, ostride);

        for (auto k = lf / 2 - 1; k < 2 * ic + lf / 2 - 1; ++k) {
            oup[(k - lf / 2 + 1) + i * ic * 2] = xLp[k];
        }
    }

    ic *= 2;

    if (shift == -1) {
        //Save the last column
        for (auto i = 0; i < ir; ++i) {
            cL[i] = oup[(i + 1) * ic - 1];
        }
        // Save the last row
        memcpy(cH.get(), oup + (ir - 1) * ic, sizeof(double) * ic);
        for (auto i = ir - 1; i > 0; --i) {
            memcpy(oup + i * ic + 1, oup + (i - 1) * ic, sizeof(double) * (ic - 1));
        }
        oup[0] = cL[ir - 1];
        for (auto i = 1; i < ir; ++i) {
            oup[i * ic] = cL[i - 1];
        }

        for (auto i = 1; i < ic; ++i) {
            oup[i] = cH[i - 1];
        }
    }
}

auto upsamp(double const* x, int lenx, int m, double* y) -> int
{
    int n;
    int i;
    int j;
    int k;

    if (m < 0) {
        return -1;
    }

    if (m == 0) {
        for (i = 0; i < lenx; ++i) {
            y[i] = x[i];
        }
        return lenx;
    }

    n = m * (lenx - 1) + 1;
    j = 1;
    k = 0;

    for (i = 0; i < n; ++i) {
        j--;
        y[i] = 0.0;
        if (j == 0) {
            y[i] = x[k];
            k++;
            j = m;
        }
    }

    return n;
}

auto upsamp2(double const* x, int lenx, int m, double* y) -> int
{
    int n;
    int i;
    int j;
    int k;
    // upsamp2 returns even numbered output. Last value is set to zero
    if (m < 0) {
        return -1;
    }

    if (m == 0) {
        for (i = 0; i < lenx; ++i) {
            y[i] = x[i];
        }
        return lenx;
    }

    n = m * lenx;
    j = 1;
    k = 0;

    for (i = 0; i < n; ++i) {
        j--;
        y[i] = 0.0;
        if (j == 0) {
            y[i] = x[k];
            k++;
            j = m;
        }
    }

    return n;
}

auto downsamp(double const* x, int lenx, int m, double* y) -> int
{
    int n;
    int i;

    if (m < 0) {
        return -1;
    }
    if (m == 0) {
        for (i = 0; i < lenx; ++i) {
            y[i] = x[i];
        }
        return lenx;
    }

    n = (lenx - 1) / m + 1;

    for (i = 0; i < n; ++i) {
        y[i] = x[i * m];
    }

    return n;
}
/*
int per_ext(double *sig, int len, int a,double *oup) {
	int i,len2;
	// oup is of length len + (len%2) + 2 * a
	for (i = 0; i < len; ++i) {
		oup[a + i] = sig[i];
	}
	len2 = len;
	if ((len % 2) != 0) {
		len2 = len + 1;
		oup[a + len] = sig[len - 1];
	}
	for (i = 0; i < a; ++i) {
		oup[a-1-i] = sig[len - 1 - i];
		oup[len2 + a + i] = sig[i];
	}

	return len2;

}
*/

auto perExt(double const* sig, int len, int a, double* oup) -> int
{
    int i;
    int len2;
    double temp1;
    double temp2;
    for (i = 0; i < len; ++i) {
        oup[a + i] = sig[i];
    }
    len2 = len;
    if ((len % 2) != 0) {
        len2 = len + 1;
        oup[a + len] = sig[len - 1];
    }
    for (i = 0; i < a; ++i) {
        temp1 = oup[a + i];
        temp2 = oup[a + len2 - 1 - i];
        oup[a - 1 - i] = temp2;
        oup[len2 + a + i] = temp1;
    }
    return len2;
}
/*
int symm_ext(double *sig, int len, int a, double *oup) {
	int i, len2;
	// oup is of length len + 2 * a
	for (i = 0; i < len; ++i) {
		oup[a + i] = sig[i];
	}
	len2 = len;
	for (i = 0; i < a; ++i) {
		oup[a - 1 - i] = sig[i];
		oup[len2 + a + i] = sig[len - 1 - i];
	}

	return len2;

}
*/

auto symmExt(double const* sig, int len, int a, double* oup) -> int
{
    int i;
    int len2;
    double temp1;
    double temp2;
    // oup is of length len + 2 * a
    for (i = 0; i < len; ++i) {
        oup[a + i] = sig[i];
    }
    len2 = len;
    for (i = 0; i < a; ++i) {
        temp1 = oup[a + i];
        temp2 = oup[a + len2 - 1 - i];
        oup[a - 1 - i] = temp1;
        oup[len2 + a + i] = temp2;
    }

    return len2;
}

static auto isign(int n) -> int
{
    int m;
    if (n >= 0) {
        m = 1;
    } else {
        m = -1;
    }

    return m;
}

static auto iabs(int n) -> int
{
    if (n >= 0) {
        return n;
    }
    return -n;
}

void circshift(double* array, int n, int l)
{
    if (iabs(l) > n) {
        l = isign(l) * (iabs(l) % n);
    }
    if (l < 0) {
        l = (n + l) % n;
    }

    auto temp = makeZeros<double>(l);
    for (auto i = 0; i < l; ++i) {
        temp[i] = array[i];
    }
    for (auto i = 0; i < n - l; ++i) {
        array[i] = array[i + l];
    }
    for (auto i = 0; i < l; ++i) {
        array[n - l + i] = temp[i];
    }
}

auto testSWTlength(int n, int j) -> int
{
    int ret;
    int div;
    int i;
    ret = 1;

    div = 1;
    for (i = 0; i < j; ++i) {
        div *= 2;
    }

    if ((n % div) != 0) {
        ret = 0;
    }

    return ret;
}

auto wmaxiter(int sigLen, int filtLen) -> int
{
    int lev;
    double temp;

    temp = std::log((double)sigLen / ((double)filtLen - 1.0)) / std::log(2.0);
    lev = (int)temp;

    return lev;
}

static auto entropyS(double const* x, int n) -> double
{
    int i;
    double val;
    double x2;

    val = 0.0;

    for (i = 0; i < n; ++i) {
        if (x[i] != 0) {
            x2 = x[i] * x[i];
            val -= x2 * std::log(x2);
        }
    }
    return val;
}

static auto entropyT(double* x, int n, double t) -> double
{
    int i;
    double val;
    double x2;
    if (t < 0) {
        std::printf("Threshold value must be >= 0");
        std::exit(1);
    }
    val = 0.0;

    for (i = 0; i < n; ++i) {
        x2 = fabs(x[i]);
        if (x2 > t) {
            val += 1;
        }
    }

    return val;
}

static auto entropyN(double* x, int n, double p) -> double
{
    int i;
    double val;
    double x2;
    if (p < 1) {
        std::printf("Norm power value must be >= 1");
        std::exit(1);
    }
    val = 0.0;
    for (i = 0; i < n; ++i) {
        x2 = fabs(x[i]);
        val += std::pow(x2, (double)p);
    }

    return val;
}

static auto entropyL(double const* x, int n) -> double
{
    int i;
    double val;
    double x2;

    val = 0.0;

    for (i = 0; i < n; ++i) {
        if (x[i] != 0) {
            x2 = x[i] * x[i];
            val += std::log(x2);
        }
    }
    return val;
}

auto costfunc(double* x, int n, char const* entropy, double p) -> double
{
    double val;

    if (strcmp(entropy, "shannon") == 0) {
        val = entropyS(x, n);
    } else if (strcmp(entropy, "threshold") == 0) {
        val = entropyT(x, n, p);
    } else if (strcmp(entropy, "norm") == 0) {
        val = entropyN(x, n, p);
    } else if ((strcmp(entropy, "logenergy") == 0) || (strcmp(entropy, "log energy") == 0) || (strcmp(entropy, "energy") == 0)) {
        val = entropyL(x, n);
    } else {
        std::printf("Entropy must be one of shannon, threshold, norm or energy");
        std::exit(-1);
    }

    return val;
}
