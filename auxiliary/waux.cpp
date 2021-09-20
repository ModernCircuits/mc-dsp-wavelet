#include "waux.h"
#include "wauxlib.h"

#include <algorithm>
#include <memory>
#include <numeric>

auto mean(double const* vec, int N) -> double
{
    return std::accumulate(vec, vec + N, 0.0) / static_cast<double>(N);
}

auto var(double const* vec, int N) -> double
{
    double v;
    double temp;
    double m;
    v = 0.0;
    m = mean(vec, N);

    for (auto i = 0; i < N; ++i) {
        temp = vec[i] - m;
        v += temp * temp;
    }

    v = v / N;

    return v;
}

auto median(double* const x, int N) -> double
{
    std::sort(x, x + N, std::less<double> {});

    double sigma;
    if ((N % 2) == 0) {
        sigma = (x[N / 2 - 1] + x[N / 2]) / 2.0;
    } else {
        sigma = x[N / 2];
    }

    return sigma;
}

auto mad(double* x, int N) -> double
{
    double sigma;

    sigma = median(x, N);

    for (auto i = 0; i < N; ++i) {
        x[i] = (x[i] - sigma) > 0 ? (x[i] - sigma) : -(x[i] - sigma);
    }

    sigma = median(x, N);

    return sigma;
}

auto minindex(double const* arr, int N) -> int
{
    double min;
    int index;

    min = DBL_MAX;
    index = 0;
    for (auto i = 0; i < N; ++i) {
        if (arr[i] < min) {
            min = arr[i];
            index = i;
        }
    }

    return index;
}

void getDWTAppx(wt_set* wt, double* appx, int N)
{
    /*
	Wavelet decomposition is stored as
	[A(J) D(J) D(J-1) ..... D(1)] in wt->output vector

	Length of A(J) , N = wt->length[0]
	*/

    for (auto i = 0; i < N; ++i) {
        appx[i] = wt->output[i];
    }
}

void getDWTDetail(wt_set* wt, double* detail, int N, int level)
{
    /*
	returns Detail coefficents at the jth level where j = J,J-1,...,1
	and Wavelet decomposition is stored as
	[A(J) D(J) D(J-1) ..... D(1)] in wt->output vector
	Use getDWTAppx() to get A(J)
	Level 1 : Length of D(J), ie N, is stored in wt->length[1]
	Level 2 :Length of D(J-1), ie N, is stored in wt->length[2]
	....
	Level J : Length of D(1), ie N, is stored in wt->length[J]
	*/
    int iter;
    int J;
    J = wt->J;

    if (level > J || level < 1) {
        printf("The decomposition only has 1,..,%d levels", J);
        exit(-1);
    }

    iter = wt->length[0];

    for (auto i = 1; i < J - level; ++i) {
        iter += wt->length[i];
    }

    for (auto i = 0; i < N; ++i) {
        detail[i] = wt->output[i + iter];
    }
}

void getDWTRecCoeff(double const* coeff, int const* length, char const* ctype, char const* ext, int level, int J, double* lpr,
    double* hpr, int lf, int siglength, double* reccoeff)
{

    int j;
    int k;
    int det_len;
    int N;
    int l;
    int m;
    int n;
    int v;
    int t;
    int l2;

    double* filt;
    auto out = std::make_unique<double[]>(siglength + 1);
    l2 = lf / 2;
    m = -2;
    n = -1;
    if (strcmp(ext, "per") == 0) {
        if (strcmp((ctype), "appx") == 0) {
            det_len = length[0];
        } else {
            det_len = length[J - level + 1];
        }

        N = 2 * length[J];

        auto X_lp = std::make_unique<double[]>(N + 2 * lf - 1);

        for (auto i = 0; i < det_len; ++i) {
            out[i] = coeff[i];
        }

        for (j = level; j > 0; --j) {

            //idwt1(wt, temp, cA_up, out, det_len, wt->output + iter, det_len, X_lp, X_hp, out);

            if ((strcmp((ctype), "det") == 0) && j == level) {
                filt = hpr;
            } else {
                filt = lpr;
            }

            //idwt_per(wt,out, det_len, wt->output + iter, det_len, X_lp);
            m = -2;
            n = -1;

            for (auto i = 0; i < det_len + l2 - 1; ++i) {
                m += 2;
                n += 2;
                X_lp[m] = 0.0;
                X_lp[n] = 0.0;
                for (l = 0; l < l2; ++l) {
                    t = 2 * l;
                    if ((i - l) >= 0 && (i - l) < det_len) {
                        X_lp[m] += filt[t] * out[i - l];
                        X_lp[n] += filt[t + 1] * out[i - l];
                    } else if ((i - l) >= det_len && (i - l) < det_len + lf - 1) {
                        X_lp[m] += filt[t] * out[i - l - det_len];
                        X_lp[n] += filt[t + 1] * out[i - l - det_len];
                    } else if ((i - l) < 0 && (i - l) > -l2) {
                        X_lp[m] += filt[t] * out[det_len + i - l];
                        X_lp[n] += filt[t + 1] * out[det_len + i - l];
                    }
                }
            }

            for (k = lf / 2 - 1; k < 2 * det_len + lf / 2 - 1; ++k) {
                out[k - lf / 2 + 1] = X_lp[k];
            }

            if (j != 1) {
                det_len = length[J - j + 2];
            }
        }

    } else if (strcmp(ext, "sym") == 0) {
        if (strcmp((ctype), "appx") == 0) {
            det_len = length[0];
        } else {
            det_len = length[J - level + 1];
        }

        N = 2 * length[J] - 1;

        auto X_lp = std::make_unique<double[]>(N + 2 * lf - 1);
        for (auto i = 0; i < det_len; ++i) {
            out[i] = coeff[i];
        }

        for (j = level; j > 0; --j) {

            //idwt1(wt, temp, cA_up, out, det_len, wt->output + iter, det_len, X_lp, X_hp, out);

            if ((strcmp((ctype), "det") == 0) && j == level) {
                filt = hpr;
            } else {
                filt = lpr;
            }

            //idwt_sym(wt, out, det_len, wt->output + iter, det_len, X_lp);

            m = -2;
            n = -1;

            for (v = 0; v < det_len; ++v) {
                auto i = v;
                m += 2;
                n += 2;
                X_lp[m] = 0.0;
                X_lp[n] = 0.0;
                for (l = 0; l < lf / 2; ++l) {
                    t = 2 * l;
                    if ((i - l) >= 0 && (i - l) < det_len) {
                        X_lp[m] += filt[t] * out[i - l];
                        X_lp[n] += filt[t + 1] * out[i - l];
                    }
                }
            }

            for (k = lf - 2; k < 2 * det_len; ++k) {
                out[k - lf + 2] = X_lp[k];
            }

            if (j != 1) {
                det_len = length[J - j + 2];
            }
        }

    } else {
        printf("Signal extension can be either per or sym");
        exit(-1);
    }

    for (auto i = 0; i < siglength; ++i) {
        reccoeff[i] = out[i];
    }
}

void autocovar(double const* vec, int N, double* acov, int M)
{
    double m;
    double temp1;
    double temp2;
    int t;
    m = mean(vec, N);

    if (M > N) {
        M = N - 1;
        printf("\n Lag is greater than the length N of the input vector. It is automatically set to length N - 1.\n");
        printf("\n The Output Vector only contains N calculated values.");
    } else if (M < 0) {
        M = 0;
    }

    for (auto i = 0; i < M; i++) {
        acov[i] = 0.0;
        for (t = 0; t < N - i; t++) {
            temp1 = vec[t] - m;
            temp2 = vec[t + i] - m;
            acov[i] += temp1 * temp2;
        }
        acov[i] = acov[i] / N;
    }
}

void autocorr(double const* vec, int N, double* acorr, int M)
{
    double var;
    if (M > N) {
        M = N - 1;
        printf("\n Lag is greater than the length N of the input vector. It is automatically set to length N - 1.\n");
        printf("\n The Output Vector only contains N calculated values.");
    } else if (M < 0) {
        M = 0;
    }
    autocovar(vec, N, acorr, M);
    var = acorr[0];
    acorr[0] = 1.0;

    for (auto i = 1; i < M; i++) {
        acorr[i] = acorr[i] / var;
    }
}
