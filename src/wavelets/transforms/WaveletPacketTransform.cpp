#include "WaveletPacketTransform.hpp"

#include "wavelets/conv.h"
#include "wavelets/hsfft.h"
#include "wavelets/wtmath.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string_view>

using namespace std::string_view_literals;

auto wptInit(Wavelet* wave, int siglength, int j) -> WaveletPacketTransform*
{
    auto const size = wave->size();

    if (j > 100) {
        printf("\n The Decomposition Iterations Cannot Exceed 100. Exiting \n");
        exit(-1);
    }

    auto const maxIter = wmaxiter(siglength, size);
    if (j > maxIter) {
        printf("\n Error - The Signal Can only be iterated %d times using this wavelet. Exiting\n", maxIter);
        exit(-1);
    }
    auto temp = 1;
    auto nodes = 0;
    for (auto i = 0; i < j; ++i) {
        temp *= 2;
        nodes += temp;
    }

    auto idx = j;
    auto p2 = 2;
    auto n = siglength;
    auto lp = size;
    auto elength = 0;
    while (idx > 0) {
        n = n + lp - 2;
        n = (int)ceil((double)n / 2.0);
        elength = p2 * n;
        idx--;
        p2 *= 2;
    }

    auto obj = std::make_unique<WaveletPacketTransform>();
    obj->params = std::make_unique<double[]>(elength + 4 * nodes + 2 * j + 6);
    obj->outlength = siglength + 2 * (j + 1) * (size + 1);
    obj->ext = "sym";
    obj->entropy = "shannon";
    obj->eparam = 0.0;

    obj->wave = wave;
    obj->siglength = siglength;
    obj->J = j;
    obj->MaxIter = maxIter;

    if (siglength % 2 == 0) {
        obj->even = 1;
    } else {
        obj->even = 0;
    }

    obj->cobj = nullptr;
    obj->nodes = nodes;

    obj->lenlength = j + 2;
    obj->output = &obj->params[0];
    obj->costvalues = &obj->params[elength];
    obj->basisvector = &obj->params[elength + nodes + 1];
    obj->nodeindex = (int*)&obj->params[elength + 2 * nodes + 2];
    obj->numnodeslevel = (int*)&obj->params[elength + 4 * nodes + 4];
    obj->coeflength = (int*)&obj->params[elength + 4 * nodes + j + 5];

    for (auto i = 0; i < elength + 4 * nodes + 2 * j + 6; ++i) {
        obj->params[i] = 0.0;
    }

    return obj.release();
}

static auto dwptPer(WaveletPacketTransform* wt, double const* inp, int n, double* cA, int lenCA, double* cD) -> void
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

static auto dwptSym(WaveletPacketTransform* wt, double const* inp, int n, double* cA, int lenCA, double* cD) -> void
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

static constexpr auto ipow2(int n) -> int
{
    auto p = 1;
    for (auto i = 0; i < n; ++i) {
        p *= 2;
    }
    return p;
}

auto dwpt(WaveletPacketTransform* wt, double const* inp) -> void
{
    int iter;
    int p2;
    int k;
    int n2;
    int np;
    int llb;
    double v1;
    double v2;
    int lenCA;
    int t;

    auto tempLen = wt->siglength;
    auto jj = wt->J;
    wt->length[jj + 1] = tempLen;
    wt->outlength = 0;
    auto temp = 1;
    auto elength = 0;
    auto size = wt->wave->size();
    auto nodes = wt->nodes;
    auto n1 = nodes + 1;
    for (auto i = 0; i < jj; ++i) {
        temp *= 2;
        auto const temp2 = (size - 2) * (temp - 1);
        elength += temp2;
    }

    auto eparam = wt->eparam;
    auto orig = std::make_unique<double[]>(tempLen);
    auto tree = std::make_unique<double[]>((tempLen * (jj + 1) + elength));
    auto nodelength = std::make_unique<int[]>(nodes);

    for (auto i = 0; i < wt->siglength; ++i) {
        orig[i] = inp[i];
    }

    for (auto i = 0; i < tempLen * (jj + 1) + elength; ++i) {
        tree[i] = 0.0;
    }

    for (auto i = 0; i < nodes + 1; ++i) {
        wt->basisvector[i] = 0.0;
        wt->costvalues[i] = 0.0;
    }

    auto n = tempLen;
    auto lp = wt->wave->lpdLen();
    // p2 = 1;

    //set eparam value here
    wt->costvalues[0] = costfunc(orig.get(), wt->siglength, wt->entropy.c_str(), eparam);
    auto it2 = 1;
    if (wt->ext == "per"sv) {
        auto i = jj;
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
        for (iter = 0; iter < jj; ++iter) {
            lenCA = wt->length[jj - iter];
            n2 -= 2 * p2 * lenCA;
            n = n2;
            for (k = 0; k < p2; ++k) {
                if (iter == 0) {
                    dwptPer(wt, orig.get(), tempLen, tree.get() + n, lenCA, tree.get() + n + lenCA);
                } else {
                    dwptPer(wt, tree.get() + np + k * tempLen, tempLen, tree.get() + n, lenCA, tree.get() + n + lenCA);
                }
                wt->costvalues[it2] = costfunc(tree.get() + n, lenCA, wt->entropy.c_str(), eparam);
                it2++;
                wt->costvalues[it2] = costfunc(tree.get() + n + lenCA, lenCA, wt->entropy.c_str(), eparam);
                it2++;
                n += 2 * lenCA;
            }

            tempLen = wt->length[jj - iter];
            p2 = 2 * p2;
            np = n2;
        }
    } else if (wt->ext == "sym"sv) {
        auto i = jj;
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

        for (iter = 0; iter < jj; ++iter) {
            lenCA = wt->length[jj - iter];
            n2 -= 2 * p2 * lenCA;
            n = n2;
            for (k = 0; k < p2; ++k) {
                if (iter == 0) {
                    dwptSym(wt, orig.get(), tempLen, tree.get() + n, lenCA, tree.get() + n + lenCA);
                } else {
                    dwptSym(wt, tree.get() + np + k * tempLen, tempLen, tree.get() + n, lenCA, tree.get() + n + lenCA);
                }
                wt->costvalues[it2] = costfunc(tree.get() + n, lenCA, wt->entropy.c_str(), eparam);
                it2++;
                wt->costvalues[it2] = costfunc(tree.get() + n + lenCA, lenCA, wt->entropy.c_str(), eparam);
                it2++;
                n += 2 * lenCA;
            }

            tempLen = wt->length[jj - iter];
            p2 = 2 * p2;
            np = n2;
        }

    } else {
        printf("Signal extension can be either per or sym");
        exit(-1);
    }

    jj = wt->J;
    auto t2 = wt->outlength - 2 * wt->length[jj];
    p2 = 2;
    auto it1 = 0;
    for (auto i = 0; i < jj; ++i) {
        t = t2;
        for (k = 0; k < p2; ++k) {
            nodelength[it1] = t;
            it1++;
            t += wt->length[jj - i];
        }
        p2 *= 2;
        t2 = t2 - p2 * wt->length[jj - i - 1];
    }

    jj = wt->J;
    llb = 1;
    for (auto i = 0; i < jj; ++i) {
        llb *= 2;
    }

    for (auto i = n1 - llb; i < n1; ++i) {
        wt->basisvector[i] = 1;
    }

    for (auto j = jj - 1; j >= 0; --j) {
        for (k = ipow2(j) - 1; k < ipow2(j + 1) - 1; ++k) {
            v1 = wt->costvalues[k];
            v2 = wt->costvalues[2 * k + 1] + wt->costvalues[2 * k + 2];
            if (v1 <= v2) {
                wt->basisvector[k] = 1;
            } else {
                wt->costvalues[k] = v2;
            }
        }
    }

    for (k = 0; k < nodes / 2; ++k) {
        if (wt->basisvector[k] == 1 || wt->basisvector[k] == 2) {
            wt->basisvector[2 * k + 1] = 2;
            wt->basisvector[2 * k + 2] = 2;
        }
    }

    for (k = 0; k < n1; ++k) {
        if (wt->basisvector[k] == 2) {
            wt->basisvector[k] = 0;
        }
    }

    // N2 = 0;
    it1 = n1;
    it2 = 0;
    wt->nodes = 0;
    wt->numnodeslevel[0] = 0;

    if (wt->basisvector[0] == 1) {
        wt->outlength = wt->siglength;
        for (auto i = 0; i < wt->siglength; ++i) {
            wt->output[i] = inp[i];
        }
        wt->nodes = 1;
        wt->nodeindex[0] = 0;
        wt->nodeindex[1] = 0;
        wt->numnodeslevel[0] = 1;
    } else {
        for (auto i = jj; i > 0; --i) {
            llb = ipow2(i);
            it1 -= llb;
            wt->numnodeslevel[i] = 0;
            for (auto j = 0; j < llb; ++j) {
                if (wt->basisvector[it1 + j] == 1) {
                    wt->nodeindex[2 * wt->nodes] = i;
                    wt->nodeindex[2 * wt->nodes + 1] = j;
                    wt->nodes += 1;
                    wt->numnodeslevel[i] += 1;
                    for (k = 0; k < wt->length[jj - i + 1]; ++k) {
                        wt->output[it2 + k] = tree[nodelength[it1 - 1 + j] + k]; // access tree
                    }
                    it2 += wt->length[jj - i + 1];
                }
            }
        }
        wt->outlength = it2;
    }

    wt->coeflength[0] = wt->siglength;

    for (auto i = 1; i < jj + 1; ++i) {
        wt->coeflength[i] = wt->length[jj - i + 1];
    }
}

/// X - Level. All Nodes at any level have the same length
auto getDWPTNodelength(WaveletPacketTransform* wt, int x) -> int
{
    if (x <= 0 || x > wt->J) {
        printf("X co-ordinate must be >= 1 and <= %d", wt->J);
        exit(-1);
    }

    return wt->length[wt->J - x + 1];
}

static auto idwptPer(WaveletPacketTransform* wt, double const* cA, int lenCA, double const* cD, double* x) -> void
{
    int lenAvg;
    int l;
    int m;
    int n;
    int t;
    int l2;

    lenAvg = (wt->wave->lprLen() + wt->wave->hprLen()) / 2;
    l2 = lenAvg / 2;
    m = -2;
    n = -1;

    for (auto i = 0; i < lenCA + l2 - 1; ++i) {
        m += 2;
        n += 2;
        x[m] = 0.0;
        x[n] = 0.0;
        for (l = 0; l < l2; ++l) {
            t = 2 * l;
            if ((i - l) >= 0 && (i - l) < lenCA) {
                x[m] += wt->wave->lpr()[t] * cA[i - l] + wt->wave->hpr()[t] * cD[i - l];
                x[n] += wt->wave->lpr()[t + 1] * cA[i - l] + wt->wave->hpr()[t + 1] * cD[i - l];
            } else if ((i - l) >= lenCA && (i - l) < lenCA + lenAvg - 1) {
                x[m] += wt->wave->lpr()[t] * cA[i - l - lenCA] + wt->wave->hpr()[t] * cD[i - l - lenCA];
                x[n] += wt->wave->lpr()[t + 1] * cA[i - l - lenCA] + wt->wave->hpr()[t + 1] * cD[i - l - lenCA];
            } else if ((i - l) < 0 && (i - l) > -l2) {
                x[m] += wt->wave->lpr()[t] * cA[lenCA + i - l] + wt->wave->hpr()[t] * cD[lenCA + i - l];
                x[n] += wt->wave->lpr()[t + 1] * cA[lenCA + i - l] + wt->wave->hpr()[t + 1] * cD[lenCA + i - l];
            }
        }
    }
}

static auto idwptSym(WaveletPacketTransform* wt, double const* cA, int lenCA, double const* cD, double* x) -> void
{
    auto lenAvg = (wt->wave->lprLen() + wt->wave->hprLen()) / 2;
    auto m = -2;
    auto n = -1;

    for (auto v = 0; v < lenCA; ++v) {
        auto i = v;
        m += 2;
        n += 2;
        x[m] = 0.0;
        x[n] = 0.0;
        for (auto l = 0; l < lenAvg / 2; ++l) {
            auto const t = 2 * l;
            if ((i - l) >= 0 && (i - l) < lenCA) {
                x[m] += wt->wave->lpr()[t] * cA[i - l] + wt->wave->hpr()[t] * cD[i - l];
                x[n] += wt->wave->lpr()[t + 1] * cA[i - l] + wt->wave->hpr()[t + 1] * cD[i - l];
            }
        }
    }
}

auto idwpt(WaveletPacketTransform* wt, double* dwtop) -> void
{
    int k;
    int l;
    int index;

    auto j = wt->J;
    auto appLen = wt->length[0];
    auto powJ = ipow2(j);
    auto lf = (wt->wave->lprLen() + wt->wave->hprLen()) / 2;
    auto xlen = powJ * (appLen + 2 * lf);

    auto xLp = std::make_unique<double[]>(2 * (wt->length[j] + lf));
    auto x = std::make_unique<double[]>(xlen);
    auto out = std::make_unique<double[]>(wt->length[j]);
    auto out2 = std::make_unique<double[]>(wt->length[j]);
    auto prep = makeZeros<int>(powJ);
    auto ptemp = makeZeros<int>(powJ);
    auto n1 = 1;
    auto llb = 1;
    auto index2 = xlen / powJ;
    auto indexp = 0;
    if (wt->basisvector[0] == 1) {
        for (auto i = 0; i < wt->siglength; ++i) {
            dwtop[i] = wt->output[i];
        }

    } else {
        for (auto i = 0; i < j; ++i) {
            llb *= 2;
            n1 += llb;
        }

        for (auto i = 0; i < xlen; ++i) {
            x[i] = 0.0;
        }

        for (auto i = 0; i < llb; ++i) {
            prep[i] = (int)wt->basisvector[n1 - llb + i];
            ptemp[i] = 0;
        }

        if (wt->ext == "per"sv) {
            index = 0;
            for (auto i = 0; i < j; ++i) {
                auto p = ipow2(j - i - 1);
                auto detLen = wt->length[i + 1];
                index2 *= 2;
                auto index3 = 0;
                auto index4 = 0;
                n1 -= llb;
                for (l = 0; l < llb; ++l) {
                    if (ptemp[l] != 2) {
                        prep[l] = (int)wt->basisvector[n1 + l];
                    } else {
                        prep[l] = ptemp[l];
                    }
                    ptemp[l] = 0;
                }

                for (l = 0; l < p; ++l) {
                    if (prep[2 * l] == 1 && prep[2 * l + 1] == 1) {
                        for (k = 0; k < detLen; ++k) {
                            out[k] = wt->output[index + k];
                            out2[k] = wt->output[index + detLen + k];
                        }
                        idwptPer(wt, out.get(), detLen, out2.get(), xLp.get());
                        for (k = lf / 2 - 1; k < 2 * detLen + lf / 2 - 1; ++k) {
                            x[index3 + k - lf / 2 + 1] = xLp[k];
                        }
                        index += 2 * detLen;
                        index3 += index2;
                        index4 += 2 * indexp;
                        ptemp[l] = 2;
                    } else if (prep[2 * l] == 1 && prep[2 * l + 1] == 2) {
                        index4 += indexp;
                        for (k = 0; k < detLen; ++k) {
                            out[k] = wt->output[index + k];
                            out2[k] = x[index4 + k];
                        }
                        idwptPer(wt, out.get(), detLen, out2.get(), xLp.get());
                        for (k = lf / 2 - 1; k < 2 * detLen + lf / 2 - 1; ++k) {
                            x[index3 + k - lf / 2 + 1] = xLp[k];
                        }
                        index += detLen;
                        index3 += index2;
                        index4 += indexp;
                        ptemp[l] = 2;
                    } else if (prep[2 * l] == 2 && prep[2 * l + 1] == 1) {
                        for (k = 0; k < detLen; ++k) {
                            out[k] = x[index4 + k];
                            out2[k] = wt->output[index + k];
                        }
                        idwptPer(wt, out.get(), detLen, out2.get(), xLp.get());
                        for (k = lf / 2 - 1; k < 2 * detLen + lf / 2 - 1; ++k) {
                            x[index3 + k - lf / 2 + 1] = xLp[k];
                        }
                        index += detLen;
                        index3 += index2;
                        index4 += 2 * indexp;
                        ptemp[l] = 2;
                    } else if (prep[2 * l] == 2 && prep[2 * l + 1] == 2) {
                        for (k = 0; k < detLen; ++k) {
                            out[k] = x[index4 + k];
                            out2[k] = x[index4 + indexp + k];
                        }
                        idwptPer(wt, out.get(), detLen, out2.get(), xLp.get());
                        for (k = lf / 2 - 1; k < 2 * detLen + lf / 2 - 1; ++k) {
                            x[index3 + k - lf / 2 + 1] = xLp[k];
                        }
                        index4 += 2 * indexp;
                        index3 += index2;
                        ptemp[l] = 2;
                    } else {
                        index3 += index2;
                        index4 += 2 * indexp;
                    }
                }

                llb /= 2;
                indexp = index2;
            }

        } else if (wt->ext == "sym"sv) {
            index = 0;

            for (auto i = 0; i < j; ++i) {
                auto p = ipow2(j - i - 1);
                auto detLen = wt->length[i + 1];
                index2 *= 2;
                auto index3 = 0;
                auto index4 = 0;
                n1 -= llb;
                for (l = 0; l < llb; ++l) {
                    if (ptemp[l] != 2) {
                        prep[l] = (int)wt->basisvector[n1 + l];
                    } else {
                        prep[l] = ptemp[l];
                    }
                    ptemp[l] = 0;
                }

                for (l = 0; l < p; ++l) {
                    if (prep[2 * l] == 1 && prep[2 * l + 1] == 1) {
                        for (k = 0; k < detLen; ++k) {
                            out[k] = wt->output[index + k];
                            out2[k] = wt->output[index + detLen + k];
                        }
                        idwptSym(wt, out.get(), detLen, out2.get(), xLp.get());
                        for (k = lf - 2; k < 2 * detLen; ++k) {
                            x[index3 + k - lf + 2] = xLp[k];
                        }
                        index += 2 * detLen;
                        index3 += index2;
                        index4 += 2 * indexp;
                        ptemp[l] = 2;
                    } else if (prep[2 * l] == 1 && prep[2 * l + 1] == 2) {
                        index4 += indexp;
                        for (k = 0; k < detLen; ++k) {
                            out[k] = wt->output[index + k];
                            out2[k] = x[index4 + k];
                        }
                        idwptSym(wt, out.get(), detLen, out2.get(), xLp.get());
                        for (k = lf - 2; k < 2 * detLen; ++k) {
                            x[index3 + k - lf + 2] = xLp[k];
                        }
                        index += detLen;
                        index3 += index2;
                        index4 += indexp;
                        ptemp[l] = 2;
                    } else if (prep[2 * l] == 2 && prep[2 * l + 1] == 1) {
                        for (k = 0; k < detLen; ++k) {
                            out[k] = x[index4 + k];
                            out2[k] = wt->output[index + k];
                        }
                        idwptSym(wt, out.get(), detLen, out2.get(), xLp.get());
                        for (k = lf - 2; k < 2 * detLen; ++k) {
                            x[index3 + k - lf + 2] = xLp[k];
                        }
                        index += detLen;
                        index3 += index2;
                        index4 += 2 * indexp;
                        ptemp[l] = 2;
                    } else if (prep[2 * l] == 2 && prep[2 * l + 1] == 2) {
                        for (k = 0; k < detLen; ++k) {
                            out[k] = x[index4 + k];
                            out2[k] = x[index4 + indexp + k];
                        }
                        idwptSym(wt, out.get(), detLen, out2.get(), xLp.get());
                        for (k = lf - 2; k < 2 * detLen; ++k) {
                            x[index3 + k - lf + 2] = xLp[k];
                        }
                        index4 += 2 * indexp;
                        index3 += index2;
                        ptemp[l] = 2;
                    } else {
                        index3 += index2;
                        index4 += 2 * indexp;
                    }
                }

                //idwt1(wt, temp, cA_up, out, det_len, wt->output().data() + iter, det_len, X_lp.get(), X_hp, out);
                /*
				idwpt_sym(wt, out, det_len, wt->output().data() + iter, det_len, X_lp);
				for (k = lf - 2; k < 2 * det_len; ++k) {
				out[k - lf + 2] = X_lp[k];
				}

				iter += det_len;
				det_len = wt->length[i + 2];
				*/
                llb /= 2;
                indexp = index2;
            }

            //free(X_lp);

        } else {
            printf("Signal extension can be either per or sym");
            exit(-1);
        }

        for (auto i = 0; i < wt->siglength; ++i) {
            dwtop[i] = x[i];
        }
    }
}

auto setDWPTExtension(WaveletPacketTransform* wt, char const* extension) -> void
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

auto setDWPTEntropy(WaveletPacketTransform* wt, char const* entropy, double eparam) -> void
{
    if (strcmp(entropy, "shannon") == 0) {
        wt->entropy = "shannon";
    } else if (strcmp(entropy, "threshold") == 0) {
        wt->entropy = "threshold";
        wt->eparam = eparam;
    } else if (strcmp(entropy, "norm") == 0) {
        wt->entropy = "norm";
        wt->eparam = eparam;
    } else if ((strcmp(entropy, "logenergy") == 0) || (strcmp(entropy, "log energy") == 0) || (strcmp(entropy, "energy") == 0)) {
        wt->entropy = "logenergy";
    } else {
        printf("Entropy should be one of shannon, threshold, norm or logenergy");
        exit(-1);
    }
}

auto summary(WaveletPacketTransform const& wt) -> void
{
    int k;
    int p2;
    int j;
    int it1;
    int it2;
    j = wt.J;
    summary(*wt.wave);
    printf("\n");
    printf("Signal Extension : %s \n", wt.ext.c_str());
    printf("\n");
    printf("Entropy : %s \n", wt.entropy.c_str());
    printf("\n");
    printf("Number of Decomposition Levels %d \n", wt.J);
    printf("\n");
    printf("Number of Active Nodes %d \n", wt.nodes);
    printf("\n");
    printf("Length of Input Signal %d \n", wt.siglength);
    printf("\n");
    printf("Length of WT Output Vector %d \n", wt.outlength);
    printf("\n");
    printf("Wavelet Coefficients are contained in vector : %s \n", "output");
    printf("\n");
    printf("Coefficients Access \n");
    it1 = 1;
    it2 = 0;
    for (auto i = 0; i < j; ++i) {
        it1 += ipow2(i + 1);
    }
    for (auto i = j; i > 0; --i) {
        p2 = ipow2(i);
        it1 -= p2;
        for (k = 0; k < p2; ++k) {
            if (wt.basisvector[it1 + k] == 1) {
                printf("Node %d %d Access : output[%d] Length : %d \n", i, k, it2, wt.length[j - i + 1]);
                it2 += wt.length[j - i + 1];
            }
        }
    }

    printf("\n");
}

auto wptFree(WaveletPacketTransform* object) -> void
{
    delete object;
}
