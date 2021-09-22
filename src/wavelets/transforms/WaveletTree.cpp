#include "WaveletTree.hpp"

#include "wavelets/wtmath.h"

#include <cmath>
#include <string_view>

using namespace std::string_view_literals;

auto wtreeInit(Wavelet* wave, int siglength, int j) -> WaveletTree*
{
    auto const size = wave->size();
    auto const maxIter = wmaxiter(siglength, size);
    if (j > 100) {
        printf("\n The Decomposition Iterations Cannot Exceed 100. Exiting \n");
        exit(-1);
    }
    if (j > maxIter) {
        printf("\n Error - The Signal Can only be iterated %d times using this wavelet. Exiting\n", maxIter);
        exit(-1);
    }

    auto temp = 1;
    auto elength = 0;
    auto nodes = 0;
    for (auto i = 0; i < j; ++i) {
        temp *= 2;
        nodes += temp;
        auto const temp2 = (size - 2) * (temp - 1);
        elength += temp2;
    }

    auto obj = std::make_unique<WaveletTree>();
    obj->params = std::make_unique<double[]>(siglength * (j + 1) + elength + nodes + j + 1);
    obj->outlength = siglength * (j + 1) + elength;
    obj->ext = "sym";

    obj->wave = wave;
    obj->siglength = siglength;
    obj->J = j;
    obj->MaxIter = maxIter;
    obj->method = "dwt";

    if (siglength % 2 == 0) {
        obj->even = 1;
    } else {
        obj->even = 0;
    }

    obj->cobj = nullptr;
    obj->nodes = nodes;

    obj->cfftset = 0;
    obj->lenlength = j + 2;
    obj->output = &obj->params[0];
    obj->nodelength = (int*)&obj->params[siglength * (j + 1) + elength];
    obj->coeflength = (int*)&obj->params[siglength * (j + 1) + elength + nodes];

    for (auto i = 0; i < siglength * (j + 1) + elength + nodes + j + 1; ++i) {
        obj->params[i] = 0.0;
    }

    return obj.release();
}

static auto wtreePer(WaveletTree* wt, double const* inp, int n, double* cA, int lenCA, double* cD) -> void
{
    int l;
    int l2;
    int isodd;
    int t;
    int lenAvg;

    lenAvg = wt->wave->lpdLen();
    l2 = lenAvg / 2;
    isodd = n % 2;

    for (auto i = 0; i < lenCA; ++i) {
        t = 2 * i + l2;
        cA[i] = 0.0;
        cD[i] = 0.0;
        for (l = 0; l < lenAvg; ++l) {
            if ((t - l) >= l2 && (t - l) < n) {
                cA[i] += wt->wave->lpd()[l] * inp[t - l];
                cD[i] += wt->wave->hpd()[l] * inp[t - l];
            } else if ((t - l) < l2 && (t - l) >= 0) {
                cA[i] += wt->wave->lpd()[l] * inp[t - l];
                cD[i] += wt->wave->hpd()[l] * inp[t - l];
            } else if ((t - l) < 0 && isodd == 0) {
                cA[i] += wt->wave->lpd()[l] * inp[t - l + n];
                cD[i] += wt->wave->hpd()[l] * inp[t - l + n];
            } else if ((t - l) < 0 && isodd == 1) {
                if ((t - l) != -1) {
                    cA[i] += wt->wave->lpd()[l] * inp[t - l + n + 1];
                    cD[i] += wt->wave->hpd()[l] * inp[t - l + n + 1];
                } else {
                    cA[i] += wt->wave->lpd()[l] * inp[n - 1];
                    cD[i] += wt->wave->hpd()[l] * inp[n - 1];
                }
            } else if ((t - l) >= n && isodd == 0) {
                cA[i] += wt->wave->lpd()[l] * inp[t - l - n];
                cD[i] += wt->wave->hpd()[l] * inp[t - l - n];
            } else if ((t - l) >= n && isodd == 1) {
                if (t - l != n) {
                    cA[i] += wt->wave->lpd()[l] * inp[t - l - (n + 1)];
                    cD[i] += wt->wave->hpd()[l] * inp[t - l - (n + 1)];
                } else {
                    cA[i] += wt->wave->lpd()[l] * inp[n - 1];
                    cD[i] += wt->wave->hpd()[l] * inp[n - 1];
                }
            }
        }
    }
}

static auto wtreeSym(WaveletTree* wt, double const* inp, int n, double* cA, int lenCA, double* cD) -> void
{
    int l;
    int t;
    int lenAvg;

    lenAvg = wt->wave->lpdLen();

    for (auto i = 0; i < lenCA; ++i) {
        t = 2 * i + 1;
        cA[i] = 0.0;
        cD[i] = 0.0;
        for (l = 0; l < lenAvg; ++l) {
            if ((t - l) >= 0 && (t - l) < n) {
                cA[i] += wt->wave->lpd()[l] * inp[t - l];
                cD[i] += wt->wave->hpd()[l] * inp[t - l];
            } else if ((t - l) < 0) {
                cA[i] += wt->wave->lpd()[l] * inp[-t + l - 1];
                cD[i] += wt->wave->hpd()[l] * inp[-t + l - 1];
            } else if ((t - l) >= n) {
                cA[i] += wt->wave->lpd()[l] * inp[2 * n - t + l - 1];
                cD[i] += wt->wave->hpd()[l] * inp[2 * n - t + l - 1];
            }
        }
    }
}

auto wtree(WaveletTree* wt, double const* inp) -> void
{
    int iter;
    int n;
    int lp;
    int p2;
    int k;
    int n2;
    int np;
    int lenCA;
    int t;
    int t2;
    int it1;

    auto tempLen = wt->siglength;
    auto j = wt->J;
    wt->length[j + 1] = tempLen;
    wt->outlength = 0;
    wt->zpad = 0;

    auto orig = std::make_unique<double[]>(tempLen);

    for (auto i = 0; i < wt->siglength; ++i) {
        orig[i] = inp[i];
    }

    if (wt->zpad == 1) {
        orig[tempLen - 1] = orig[tempLen - 2];
    }

    n = tempLen;
    lp = wt->wave->lpdLen();

    if (wt->ext == "per"sv) {
        auto i = j;
        p2 = 2;
        while (i > 0) {
            n = (int)ceil((double)n / 2.0);
            wt->length[i] = n;
            wt->outlength += p2 * (wt->length[i]);
            i--;
            p2 *= 2;
        }
        wt->length[0] = wt->length[1];

        n2 = wt->outlength;
        p2 = 1;
        for (iter = 0; iter < j; ++iter) {
            lenCA = wt->length[j - iter];
            n2 -= 2 * p2 * lenCA;
            n = n2;
            for (k = 0; k < p2; ++k) {
                if (iter == 0) {
                    wtreePer(wt, orig.get(), tempLen, wt->params.get() + n, lenCA, wt->params.get() + n + lenCA);
                } else {
                    wtreePer(wt, wt->params.get() + np + k * tempLen, tempLen, wt->params.get() + n, lenCA, wt->params.get() + n + lenCA);
                }
                n += 2 * lenCA;
            }

            tempLen = wt->length[j - iter];
            p2 = 2 * p2;
            np = n2;
        }
    } else if (wt->ext == "sym"sv) {
        auto i = j;
        p2 = 2;
        while (i > 0) {
            n = n + lp - 2;
            n = (int)ceil((double)n / 2.0);
            wt->length[i] = n;
            wt->outlength += p2 * (wt->length[i]);
            i--;
            p2 *= 2;
        }
        wt->length[0] = wt->length[1];

        n2 = wt->outlength;
        p2 = 1;

        for (iter = 0; iter < j; ++iter) {
            lenCA = wt->length[j - iter];
            n2 -= 2 * p2 * lenCA;
            n = n2;
            for (k = 0; k < p2; ++k) {
                if (iter == 0) {
                    wtreeSym(wt, orig.get(), tempLen, wt->params.get() + n, lenCA, wt->params.get() + n + lenCA);
                } else {
                    wtreeSym(wt, wt->params.get() + np + k * tempLen, tempLen, wt->params.get() + n, lenCA, wt->params.get() + n + lenCA);
                }
                n += 2 * lenCA;
            }

            tempLen = wt->length[j - iter];
            p2 = 2 * p2;
            np = n2;
        }

    } else {
        printf("Signal extension can be either per or sym");
        exit(-1);
    }

    j = wt->J;
    t2 = wt->outlength - 2 * wt->length[j];
    p2 = 2;
    it1 = 0;
    for (auto i = 0; i < j; ++i) {
        t = t2;
        for (k = 0; k < p2; ++k) {
            wt->nodelength[it1] = t;
            it1++;
            t += wt->length[j - i];
        }
        p2 *= 2;
        t2 = t2 - p2 * wt->length[j - i - 1];
    }

    wt->coeflength[0] = wt->siglength;

    for (auto i = 1; i < j + 1; ++i) {
        wt->coeflength[i] = wt->length[j - i + 1];
    }
}

/// X - Level. All Nodes at any level have the same length
auto getWTREENodelength(WaveletTree* wt, int x) -> int
{
    if (x <= 0 || x > wt->J) {
        printf("X co-ordinate must be >= 1 and <= %d", wt->J);
        exit(-1);
    }

    return wt->length[wt->J - x + 1];
}

auto getWTREECoeffs(WaveletTree* wt, int x, int y, double* coeffs, int n) -> void
{
    int ymax;
    int t;
    int t2;

    if (x <= 0 || x > wt->J) {
        printf("X co-ordinate must be >= 1 and <= %d", wt->J);
        exit(-1);
    }
    ymax = 1;
    for (auto i = 0; i < x; ++i) {
        ymax *= 2;
    }

    ymax -= 1;

    if (y < 0 || y > ymax) {
        printf("Y co-ordinate must be >= 0 and <= %d", ymax);
        exit(-1);
    }

    if (x == 1) {
        t = 0;
    } else {
        t = 0;
        t2 = 1;
        for (auto i = 0; i < x - 1; ++i) {
            t2 *= 2;
            t += t2;
        }
    }

    t += y;
    t2 = wt->nodelength[t];
    for (auto i = 0; i < n; ++i) {
        coeffs[i] = wt->output[t2 + i];
    }
}

auto setWTREEExtension(WaveletTree* wt, char const* extension) -> void
{
    if (extension == "sym"sv) {
        wt->ext = "sym";
    } else if (extension == "per"sv) {
        wt->ext = "per";
    } else {
        printf("Signal extension can be either per or sym");
        exit(-1);
    }
}

auto summary(WaveletTree const& wt) -> void
{
    int k;
    int p2;
    int j;
    int t;
    j = wt.J;
    summary(*wt.wave);
    printf("\n");
    printf("Wavelet Transform : %s \n", wt.method.c_str());
    printf("\n");
    printf("Signal Extension : %s \n", wt.ext.c_str());
    printf("\n");
    printf("Number of Decomposition Levels %d \n", wt.J);
    printf("\n");
    printf("Length of Input Signal %d \n", wt.siglength);
    printf("\n");
    printf("Length of WT Output Vector %d \n", wt.outlength);
    printf("\n");
    printf("Wavelet Coefficients are contained in vector : %s \n", "output");
    printf("\n");
    printf("Coefficients Access \n");
    t = 0;
    p2 = 2;
    for (auto i = 0; i < j; ++i) {
        for (k = 0; k < p2; ++k) {
            printf("Node %d %d Access : output[%d] Length : %d \n", i + 1, k, wt.nodelength[t], wt.length[j - i]);
            t++;
        }
        p2 *= 2;
    }
    printf("\n");
}

auto wtreeFree(WaveletTree* object) -> void
{
    delete object;
}
