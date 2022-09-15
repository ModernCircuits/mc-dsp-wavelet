#include "WaveletPacketTransform.hpp"

#include <mc/dsp/convolution/FFTConvolver.hpp>
#include <mc/dsp/fft/FFT.hpp>
#include <mc/dsp/wavelets/common.hpp>

#include <mc/core/cassert.hpp>
#include <mc/core/cmath.hpp>
#include <mc/core/cstdlib.hpp>
#include <mc/core/cstring.hpp>
#include <mc/core/format.hpp>
#include <mc/core/string_view.hpp>
#include <mc/core/utility.hpp>

#include <fmt/printf.h>

namespace mc::dsp {

namespace {
auto entropyS(float const* x, int n) -> float
{
    auto val = 0.0F;
    for (auto i = 0; i < n; ++i) {
        if (x[i] != 0) {
            auto x2 = x[i] * x[i];
            val -= x2 * std::log(x2);
        }
    }
    return val;
}

auto entropyT(float* x, int n, float t) -> float
{
    if (t < 0) { throw std::invalid_argument("Threshold value must be >= 0"); }

    auto val = 0.0F;
    for (auto i = 0; i < n; ++i) {
        auto const x2 = std::fabs(x[i]);
        if (x2 > t) { val += 1; }
    }

    return val;
}

auto entropyN(float* x, int n, float p) -> float
{
    if (p < 1) { throw std::invalid_argument("Norm power value must be >= 1"); }

    auto val = 0.0F;
    for (auto i = 0; i < n; ++i) {
        auto const x2 = std::fabs(x[i]);
        val += std::pow(x2, (float)p);
    }

    return val;
}

auto entropyL(float const* x, int n) -> float
{
    auto val = 0.0F;
    for (auto i = 0; i < n; ++i) {
        if (x[i] != 0) {
            auto const x2 = x[i] * x[i];
            val += std::log(x2);
        }
    }
    return val;
}

auto costfunc(float* x, int n, char const* entropy, float p) -> float
{
    auto e = StringView{entropy};
    if (e == "shannon") { return entropyS(x, n); }
    if (e == "threshold") { return entropyT(x, n, p); }
    if (e == "norm") { return entropyN(x, n, p); }
    if ((e == "logenergy") || (e == "log energy") || (e == "energy")) {
        return entropyL(x, n);
    }

    fmt::printf("Entropy must be one of shannon, threshold, norm or energy");
    return 0.0F;
}
}  // namespace

WaveletPacketTransform::WaveletPacketTransform(
    Wavelet* wave,
    std::size_t siglength,
    std::size_t j
)
{
    auto const size = wave->size();

    if (j > 100) {
        throw std::invalid_argument("Decomposition Iterations Cannot Exceed 100");
    }

    auto const maxIter = maxIterations(siglength, size);
    if (j > maxIter) {
        throw std::invalid_argument(
            "Signal Can only be iterated maxIter times using this wavelet"
        );
    }
    auto temp   = 1;
    auto nodess = 0;
    for (std::size_t i = 0; i < j; ++i) {
        temp *= 2;
        nodess += temp;
    }

    auto idx            = j;
    auto p2             = 2;
    auto n              = siglength;
    auto lp             = size;
    std::size_t elength = 0;
    while (idx > 0) {
        n       = n + lp - 2;
        n       = (int)std::ceil((float)n / 2.0F);
        elength = p2 * n;
        idx--;
        p2 *= 2;
    }

    this->params    = std::make_unique<float[]>(elength + 4 * nodess + 2 * j + 6);
    this->outlength = siglength + 2 * (j + 1) * (size + 1);
    this->ext       = "sym";
    this->entropy   = "shannon";
    this->eparam    = 0.0F;

    this->wave_         = wave;
    this->signalLength_ = siglength;
    this->J             = j;
    this->MaxIter       = maxIter;

    this->cobj  = nullptr;
    this->nodes = nodess;

    this->lenlength     = j + 2;
    this->output        = &this->params[0];
    this->costvalues    = &this->params[elength];
    this->basisvector   = &this->params[elength + nodes + 1];
    this->nodeindex     = (int*)&this->params[elength + 2 * nodes + 2];
    this->numnodeslevel = (int*)&this->params[elength + 4 * nodes + 4];
    this->coeflength    = (int*)&this->params[elength + 4 * nodes + j + 5];

    for (std::size_t i = 0; i < elength + 4 * nodes + 2 * j + 6; ++i) {
        this->params[i] = 0.0F;
    }
}

static auto
dwtPer(WaveletPacketTransform& wt, float const* inp, int n, float* cA, int lenCA, float* cD)
    -> void
{
    int l      = 0;
    int l2     = 0;
    int isodd  = 0;
    int t      = 0;
    int lenAvg = 0;

    lenAvg = wt.wave().lpd().size();
    l2     = lenAvg / 2;
    isodd  = n % 2;

    for (auto i = 0; i < lenCA; ++i) {
        t     = 2 * i + l2;
        cA[i] = 0.0F;
        cD[i] = 0.0F;
        for (l = 0; l < lenAvg; ++l) {
            if ((t - l) >= l2 && (t - l) < n) {
                cA[i] += wt.wave().lpd()[l] * inp[t - l];
                cD[i] += wt.wave().hpd()[l] * inp[t - l];
            } else if ((t - l) < l2 && (t - l) >= 0) {
                cA[i] += wt.wave().lpd()[l] * inp[t - l];
                cD[i] += wt.wave().hpd()[l] * inp[t - l];
            } else if ((t - l) < 0 && isodd == 0) {
                cA[i] += wt.wave().lpd()[l] * inp[t - l + n];
                cD[i] += wt.wave().hpd()[l] * inp[t - l + n];
            } else if ((t - l) < 0 && isodd == 1) {
                if ((t - l) != -1) {
                    cA[i] += wt.wave().lpd()[l] * inp[t - l + n + 1];
                    cD[i] += wt.wave().hpd()[l] * inp[t - l + n + 1];
                } else {
                    cA[i] += wt.wave().lpd()[l] * inp[n - 1];
                    cD[i] += wt.wave().hpd()[l] * inp[n - 1];
                }
            } else if ((t - l) >= n && isodd == 0) {
                cA[i] += wt.wave().lpd()[l] * inp[t - l - n];
                cD[i] += wt.wave().hpd()[l] * inp[t - l - n];
            } else if ((t - l) >= n && isodd == 1) {
                if (t - l != n) {
                    cA[i] += wt.wave().lpd()[l] * inp[t - l - (n + 1)];
                    cD[i] += wt.wave().hpd()[l] * inp[t - l - (n + 1)];
                } else {
                    cA[i] += wt.wave().lpd()[l] * inp[n - 1];
                    cD[i] += wt.wave().hpd()[l] * inp[n - 1];
                }
            }
        }
    }
}

static auto
dwtSym(WaveletPacketTransform& wt, float const* inp, int n, float* cA, int lenCA, float* cD)
    -> void
{
    int l      = 0;
    int t      = 0;
    int lenAvg = 0;

    lenAvg = wt.wave().lpd().size();

    for (auto i = 0; i < lenCA; ++i) {
        t     = 2 * i + 1;
        cA[i] = 0.0F;
        cD[i] = 0.0F;
        for (l = 0; l < lenAvg; ++l) {
            if ((t - l) >= 0 && (t - l) < n) {
                cA[i] += wt.wave().lpd()[l] * inp[t - l];
                cD[i] += wt.wave().hpd()[l] * inp[t - l];
            } else if ((t - l) < 0) {
                cA[i] += wt.wave().lpd()[l] * inp[-t + l - 1];
                cD[i] += wt.wave().hpd()[l] * inp[-t + l - 1];
            } else if ((t - l) >= n) {
                cA[i] += wt.wave().lpd()[l] * inp[2 * n - t + l - 1];
                cD[i] += wt.wave().hpd()[l] * inp[2 * n - t + l - 1];
            }
        }
    }
}

static constexpr auto ipow2(int n) -> int
{
    auto p = 1;
    for (auto i = 0; i < n; ++i) { p *= 2; }
    return p;
}

auto dwpt(WaveletPacketTransform& wt, float const* inp) -> void
{
    int iter  = 0;
    int p2    = 0;
    int k     = 0;
    int n2    = 0;
    int np    = 0;
    int llb   = 0;
    float v1  = NAN;
    float v2  = NAN;
    int lenCA = 0;
    int t     = 0;

    auto tempLen      = wt.signalLength();
    auto jj           = wt.J;
    wt.length[jj + 1] = tempLen;
    wt.outlength      = 0;
    auto temp         = 1;
    auto elength      = 0;
    auto size         = wt.wave().size();
    auto nodes        = wt.nodes;
    auto n1           = nodes + 1;
    for (auto i = 0; i < jj; ++i) {
        temp *= 2;
        auto const temp2 = (size - 2) * (temp - 1);
        elength += temp2;
    }

    auto eparam     = wt.eparam;
    auto orig       = std::make_unique<float[]>(tempLen);
    auto tree       = std::make_unique<float[]>((tempLen * (jj + 1) + elength));
    auto nodelength = std::make_unique<int[]>(nodes);

    for (auto i = 0; i < wt.signalLength(); ++i) { orig[i] = inp[i]; }

    for (auto i = 0; i < tempLen * (jj + 1) + elength; ++i) { tree[i] = 0.0F; }

    for (auto i = 0; i < nodes + 1; ++i) {
        wt.basisvector[i] = 0.0F;
        wt.costvalues[i]  = 0.0F;
    }

    auto n  = tempLen;
    auto lp = wt.wave().lpd().size();
    // p2 = 1;

    // set eparam value here
    wt.costvalues[0] = costfunc(orig.get(), wt.signalLength(), wt.entropy.c_str(), eparam);
    auto it2         = 1;
    if (wt.ext == StringView{"per"}) {
        auto i = jj;
        p2     = 2;
        while (i > 0) {
            n            = (int)std::ceil((float)n / 2.0F);
            wt.length[i] = n;
            wt.outlength += p2 * (wt.length[i]);
            i--;
            p2 *= 2;
        }
        wt.length[0] = wt.length[1];

        n2 = wt.outlength;
        p2 = 1;
        for (iter = 0; iter < jj; ++iter) {
            lenCA = wt.length[jj - iter];
            n2 -= 2 * p2 * lenCA;
            n = n2;
            for (k = 0; k < p2; ++k) {
                if (iter == 0) {
                    dwtPer(
                        wt,
                        orig.get(),
                        tempLen,
                        tree.get() + n,
                        lenCA,
                        tree.get() + n + lenCA
                    );
                } else {
                    dwtPer(
                        wt,
                        tree.get() + np + k * tempLen,
                        tempLen,
                        tree.get() + n,
                        lenCA,
                        tree.get() + n + lenCA
                    );
                }
                wt.costvalues[it2]
                    = costfunc(tree.get() + n, lenCA, wt.entropy.c_str(), eparam);
                it2++;
                wt.costvalues[it2]
                    = costfunc(tree.get() + n + lenCA, lenCA, wt.entropy.c_str(), eparam);
                it2++;
                n += 2 * lenCA;
            }

            tempLen = wt.length[jj - iter];
            p2      = 2 * p2;
            np      = n2;
        }
    } else if (wt.ext == StringView{"sym"}) {
        auto i = jj;
        p2     = 2;
        while (i > 0) {
            n            = n + lp - 2;
            n            = (int)std::ceil((float)n / 2.0F);
            wt.length[i] = n;
            wt.outlength += p2 * (wt.length[i]);
            i--;
            p2 *= 2;
        }
        wt.length[0] = wt.length[1];

        n2 = wt.outlength;
        p2 = 1;

        for (iter = 0; iter < jj; ++iter) {
            lenCA = wt.length[jj - iter];
            n2 -= 2 * p2 * lenCA;
            n = n2;
            for (k = 0; k < p2; ++k) {
                if (iter == 0) {
                    dwtSym(
                        wt,
                        orig.get(),
                        tempLen,
                        tree.get() + n,
                        lenCA,
                        tree.get() + n + lenCA
                    );
                } else {
                    dwtSym(
                        wt,
                        tree.get() + np + k * tempLen,
                        tempLen,
                        tree.get() + n,
                        lenCA,
                        tree.get() + n + lenCA
                    );
                }
                wt.costvalues[it2]
                    = costfunc(tree.get() + n, lenCA, wt.entropy.c_str(), eparam);
                it2++;
                wt.costvalues[it2]
                    = costfunc(tree.get() + n + lenCA, lenCA, wt.entropy.c_str(), eparam);
                it2++;
                n += 2 * lenCA;
            }

            tempLen = wt.length[jj - iter];
            p2      = 2 * p2;
            np      = n2;
        }
    } else {
        throw std::invalid_argument("Signal extension can be either per or sym");
    }

    jj       = wt.J;
    auto t2  = wt.outlength - 2 * wt.length[jj];
    p2       = 2;
    auto it1 = 0;
    for (auto i = 0; i < jj; ++i) {
        t = t2;
        for (k = 0; k < p2; ++k) {
            nodelength[it1] = t;
            it1++;
            t += wt.length[jj - i];
        }
        p2 *= 2;
        t2 = t2 - p2 * wt.length[jj - i - 1];
    }

    jj  = wt.J;
    llb = 1;
    for (auto i = 0; i < jj; ++i) { llb *= 2; }

    for (auto i = n1 - llb; i < n1; ++i) { wt.basisvector[i] = 1; }

    for (auto j = jj - 1; j >= 0; --j) {
        for (k = ipow2(j) - 1; k < ipow2(j + 1) - 1; ++k) {
            v1 = wt.costvalues[k];
            v2 = wt.costvalues[2 * k + 1] + wt.costvalues[2 * k + 2];
            if (v1 <= v2) {
                wt.basisvector[k] = 1;
            } else {
                wt.costvalues[k] = v2;
            }
        }
    }

    for (k = 0; k < nodes / 2; ++k) {
        if (wt.basisvector[k] == 1 || wt.basisvector[k] == 2) {
            wt.basisvector[2 * k + 1] = 2;
            wt.basisvector[2 * k + 2] = 2;
        }
    }

    for (k = 0; k < n1; ++k) {
        if (wt.basisvector[k] == 2) { wt.basisvector[k] = 0; }
    }

    // N2 = 0;
    it1                 = n1;
    it2                 = 0;
    wt.nodes            = 0;
    wt.numnodeslevel[0] = 0;

    if (wt.basisvector[0] == 1) {
        wt.outlength = wt.signalLength();
        for (auto i = 0; i < wt.signalLength(); ++i) { wt.output[i] = inp[i]; }
        wt.nodes            = 1;
        wt.nodeindex[0]     = 0;
        wt.nodeindex[1]     = 0;
        wt.numnodeslevel[0] = 1;
    } else {
        for (auto i = jj; i > 0; --i) {
            llb = ipow2(i);
            it1 -= llb;
            wt.numnodeslevel[i] = 0;
            for (auto j = 0; j < llb; ++j) {
                if (wt.basisvector[it1 + j] == 1) {
                    wt.nodeindex[2 * wt.nodes]     = i;
                    wt.nodeindex[2 * wt.nodes + 1] = j;
                    wt.nodes += 1;
                    wt.numnodeslevel[i] += 1;
                    for (k = 0; k < wt.length[jj - i + 1]; ++k) {
                        wt.output[it2 + k]
                            = tree[nodelength[it1 - 1 + j] + k];  // access tree
                    }
                    it2 += wt.length[jj - i + 1];
                }
            }
        }
        wt.outlength = it2;
    }

    wt.coeflength[0] = wt.signalLength();

    for (auto i = 1; i < jj + 1; ++i) { wt.coeflength[i] = wt.length[jj - i + 1]; }
}

/// X - Level. All Nodes at any level have the same length
auto getDWPTNodelength(WaveletPacketTransform& wt, int x) -> int
{
    if (x <= 0 || x > wt.J) {
        throw std::invalid_argument("X co-ordinate must be >= 1 and <= wt.J");
    }

    return wt.length[wt.J - x + 1];
}

static auto
idwtPer(WaveletPacketTransform& wt, float const* cA, int lenCA, float const* cD, float* x)
    -> void
{
    int lenAvg = 0;
    int l      = 0;
    int m      = 0;
    int n      = 0;
    int t      = 0;
    int l2     = 0;

    lenAvg = (wt.wave().lpr().size() + wt.wave().hpr().size()) / 2;
    l2     = lenAvg / 2;
    m      = -2;
    n      = -1;

    for (auto i = 0; i < lenCA + l2 - 1; ++i) {
        m += 2;
        n += 2;
        x[m] = 0.0F;
        x[n] = 0.0F;
        for (l = 0; l < l2; ++l) {
            t = 2 * l;
            if ((i - l) >= 0 && (i - l) < lenCA) {
                x[m] += wt.wave().lpr()[t] * cA[i - l] + wt.wave().hpr()[t] * cD[i - l];
                x[n] += wt.wave().lpr()[t + 1] * cA[i - l]
                      + wt.wave().hpr()[t + 1] * cD[i - l];
            } else if ((i - l) >= lenCA && (i - l) < lenCA + lenAvg - 1) {
                x[m] += wt.wave().lpr()[t] * cA[i - l - lenCA]
                      + wt.wave().hpr()[t] * cD[i - l - lenCA];
                x[n] += wt.wave().lpr()[t + 1] * cA[i - l - lenCA]
                      + wt.wave().hpr()[t + 1] * cD[i - l - lenCA];
            } else if ((i - l) < 0 && (i - l) > -l2) {
                x[m] += wt.wave().lpr()[t] * cA[lenCA + i - l]
                      + wt.wave().hpr()[t] * cD[lenCA + i - l];
                x[n] += wt.wave().lpr()[t + 1] * cA[lenCA + i - l]
                      + wt.wave().hpr()[t + 1] * cD[lenCA + i - l];
            }
        }
    }
}

static auto
idwtSym(WaveletPacketTransform& wt, float const* cA, int lenCA, float const* cD, float* x)
    -> void
{
    auto lenAvg = (wt.wave().lpr().size() + wt.wave().hpr().size()) / 2;
    auto m      = -2;
    auto n      = -1;

    for (auto v = 0; v < lenCA; ++v) {
        auto i = v;
        m += 2;
        n += 2;
        x[m] = 0.0F;
        x[n] = 0.0F;
        for (auto l = 0; l < static_cast<int>(lenAvg) / 2; ++l) {
            auto const t = 2 * l;
            if ((i - l) >= 0 && i - l < lenCA) {
                x[m] += wt.wave().lpr()[t] * cA[i - l] + wt.wave().hpr()[t] * cD[i - l];
                x[n] += wt.wave().lpr()[t + 1] * cA[i - l]
                      + wt.wave().hpr()[t + 1] * cD[i - l];
            }
        }
    }
}

auto idwpt(WaveletPacketTransform& wt, float* dwtop) -> void
{
    int k     = 0;
    int l     = 0;
    int index = 0;

    auto j      = wt.J;
    auto appLen = wt.length[0];
    auto powJ   = ipow2(j);
    auto lf     = (wt.wave().lpr().size() + wt.wave().hpr().size()) / 2;
    auto xlen   = powJ * (appLen + 2 * lf);

    auto xLp    = std::make_unique<float[]>(2 * (wt.length[j] + lf));
    auto x      = std::make_unique<float[]>(xlen);
    auto out    = std::make_unique<float[]>(wt.length[j]);
    auto out2   = std::make_unique<float[]>(wt.length[j]);
    auto prep   = makeZeros<int>(powJ);
    auto ptemp  = makeZeros<int>(powJ);
    auto n1     = 1;
    auto llb    = 1;
    auto index2 = xlen / powJ;
    auto indexp = 0;
    if (wt.basisvector[0] == 1) {
        for (auto i = 0; i < wt.signalLength(); ++i) { dwtop[i] = wt.output[i]; }
    } else {
        for (auto i = 0; i < j; ++i) {
            llb *= 2;
            n1 += llb;
        }

        for (std::size_t i = 0; i < xlen; ++i) { x[i] = 0.0F; }

        for (auto i = 0; i < llb; ++i) {
            prep[i]  = (int)wt.basisvector[n1 - llb + i];
            ptemp[i] = 0;
        }

        if (wt.ext == StringView{"per"}) {
            index = 0;
            for (auto i = 0; i < j; ++i) {
                auto p      = ipow2(j - i - 1);
                auto detLen = wt.length[i + 1];
                index2 *= 2;
                auto index3 = 0;
                auto index4 = 0;
                n1 -= llb;
                for (l = 0; l < llb; ++l) {
                    if (ptemp[l] != 2) {
                        prep[l] = (int)wt.basisvector[n1 + l];
                    } else {
                        prep[l] = ptemp[l];
                    }
                    ptemp[l] = 0;
                }

                for (l = 0; l < p; ++l) {
                    if (prep[2 * l] == 1 && prep[2 * l + 1] == 1) {
                        for (k = 0; k < detLen; ++k) {
                            out[k]  = wt.output[index + k];
                            out2[k] = wt.output[index + detLen + k];
                        }
                        idwtPer(wt, out.get(), detLen, out2.get(), xLp.get());
                        for (k = lf / 2 - 1; cmp_less(k, 2 * detLen + lf / 2 - 1); ++k) {
                            x[index3 + k - lf / 2 + 1] = xLp[k];
                        }
                        index += 2 * detLen;
                        index3 += index2;
                        index4 += 2 * indexp;
                        ptemp[l] = 2;
                    } else if (prep[2 * l] == 1 && prep[2 * l + 1] == 2) {
                        index4 += indexp;
                        for (k = 0; k < detLen; ++k) {
                            out[k]  = wt.output[index + k];
                            out2[k] = x[index4 + k];
                        }
                        idwtPer(wt, out.get(), detLen, out2.get(), xLp.get());
                        for (k = lf / 2 - 1; cmp_less(k, 2 * detLen + lf / 2 - 1); ++k) {
                            x[index3 + k - lf / 2 + 1] = xLp[k];
                        }
                        index += detLen;
                        index3 += index2;
                        index4 += indexp;
                        ptemp[l] = 2;
                    } else if (prep[2 * l] == 2 && prep[2 * l + 1] == 1) {
                        for (k = 0; k < detLen; ++k) {
                            out[k]  = x[index4 + k];
                            out2[k] = wt.output[index + k];
                        }
                        idwtPer(wt, out.get(), detLen, out2.get(), xLp.get());
                        for (k = lf / 2 - 1; cmp_less(k, 2 * detLen + lf / 2 - 1); ++k) {
                            x[index3 + k - lf / 2 + 1] = xLp[k];
                        }
                        index += detLen;
                        index3 += index2;
                        index4 += 2 * indexp;
                        ptemp[l] = 2;
                    } else if (prep[2 * l] == 2 && prep[2 * l + 1] == 2) {
                        for (k = 0; k < detLen; ++k) {
                            out[k]  = x[index4 + k];
                            out2[k] = x[index4 + indexp + k];
                        }
                        idwtPer(wt, out.get(), detLen, out2.get(), xLp.get());
                        for (k = lf / 2 - 1; cmp_less(k, 2 * detLen + lf / 2 - 1); ++k) {
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
        } else if (wt.ext == StringView{"sym"}) {
            index = 0;

            for (auto i = 0; i < j; ++i) {
                auto p      = ipow2(j - i - 1);
                auto detLen = wt.length[i + 1];
                index2 *= 2;
                auto index3 = 0;
                auto index4 = 0;
                n1 -= llb;
                for (l = 0; l < llb; ++l) {
                    if (ptemp[l] != 2) {
                        prep[l] = (int)wt.basisvector[n1 + l];
                    } else {
                        prep[l] = ptemp[l];
                    }
                    ptemp[l] = 0;
                }

                for (l = 0; l < p; ++l) {
                    if (prep[2 * l] == 1 && prep[2 * l + 1] == 1) {
                        for (k = 0; k < detLen; ++k) {
                            out[k]  = wt.output[index + k];
                            out2[k] = wt.output[index + detLen + k];
                        }
                        idwtSym(wt, out.get(), detLen, out2.get(), xLp.get());
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
                            out[k]  = wt.output[index + k];
                            out2[k] = x[index4 + k];
                        }
                        idwtSym(wt, out.get(), detLen, out2.get(), xLp.get());
                        for (k = lf - 2; k < 2 * detLen; ++k) {
                            x[index3 + k - lf + 2] = xLp[k];
                        }
                        index += detLen;
                        index3 += index2;
                        index4 += indexp;
                        ptemp[l] = 2;
                    } else if (prep[2 * l] == 2 && prep[2 * l + 1] == 1) {
                        for (k = 0; k < detLen; ++k) {
                            out[k]  = x[index4 + k];
                            out2[k] = wt.output[index + k];
                        }
                        idwtSym(wt, out.get(), detLen, out2.get(), xLp.get());
                        for (k = lf - 2; k < 2 * detLen; ++k) {
                            x[index3 + k - lf + 2] = xLp[k];
                        }
                        index += detLen;
                        index3 += index2;
                        index4 += 2 * indexp;
                        ptemp[l] = 2;
                    } else if (prep[2 * l] == 2 && prep[2 * l + 1] == 2) {
                        for (k = 0; k < detLen; ++k) {
                            out[k]  = x[index4 + k];
                            out2[k] = x[index4 + indexp + k];
                        }
                        idwtSym(wt, out.get(), detLen, out2.get(), xLp.get());
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

                // idwt1(wt, temp, cA_up, out, det_len, wt.output().data() + iter, det_len,
                // X_lp.get(), X_hp, out);
                /*
                            idwt_sym(wt, out, det_len, wt.output().data() + iter, det_len,
                   X_lp); for (k = lf - 2; k < 2 * det_len; ++k) { out[k - lf + 2] =
                   X_lp[k];
                            }

                            iter += det_len;
                            det_len = wt.length[i + 2];
                            */
                llb /= 2;
                indexp = index2;
            }

            // free(X_lp);
        } else {
            throw std::invalid_argument("Signal extension can be either per or sym");
        }

        for (auto i = 0; i < wt.signalLength(); ++i) { dwtop[i] = x[i]; }
    }
}

auto setDWPTExtension(WaveletPacketTransform& wt, char const* extension) -> void
{
    if (extension == StringView{"sym"}) {
        wt.ext = "sym";
    } else if (extension == StringView{"per"}) {
        wt.ext = "per";
    } else {
        throw std::invalid_argument("Signal extension can be either per or sym");
    }
}

auto setDWPTEntropy(WaveletPacketTransform& wt, char const* entropy, float eparam) -> void
{
    if (strcmp(entropy, "shannon") == 0) {
        wt.entropy = "shannon";
    } else if (strcmp(entropy, "threshold") == 0) {
        wt.entropy = "threshold";
        wt.eparam  = eparam;
    } else if (strcmp(entropy, "norm") == 0) {
        wt.entropy = "norm";
        wt.eparam  = eparam;
    } else if ((strcmp(entropy, "logenergy") == 0) || (strcmp(entropy, "log energy") == 0) || (strcmp(entropy, "energy") == 0)) {
        wt.entropy = "logenergy";
    } else {
        throw std::invalid_argument(
            "Entropy should be one of shannon, threshold, norm or logenergy"
        );
    }
}

auto summary(WaveletPacketTransform const& wt) -> void
{
    int k   = 0;
    int p2  = 0;
    int j   = 0;
    int it1 = 0;
    int it2 = 0;
    j       = wt.J;
    summary(wt.wave());
    fmt::printf("\n");
    fmt::printf("Signal Extension : %s \n", wt.ext.c_str());
    fmt::printf("\n");
    fmt::printf("Entropy : %s \n", wt.entropy.c_str());
    fmt::printf("\n");
    fmt::printf("Number of Decomposition Levels %d \n", wt.J);
    fmt::printf("\n");
    fmt::printf("Number of Active Nodes %d \n", wt.nodes);
    fmt::printf("\n");
    fmt::printf("Length of Input Signal %d \n", wt.signalLength());
    fmt::printf("\n");
    fmt::printf("Length of WT Output Vector %d \n", wt.outlength);
    fmt::printf("\n");
    fmt::printf("Wavelet Coefficients are contained in vector : %s \n", "output");
    fmt::printf("\n");
    fmt::printf("Coefficients Access \n");
    it1 = 1;
    it2 = 0;
    for (auto i = 0; i < j; ++i) { it1 += ipow2(i + 1); }
    for (auto i = j; i > 0; --i) {
        p2 = ipow2(i);
        it1 -= p2;
        for (k = 0; k < p2; ++k) {
            if (wt.basisvector[it1 + k] == 1) {
                fmt::printf(
                    "Node %d %d Access : output[%d] Length : %d \n",
                    i,
                    k,
                    it2,
                    wt.length[j - i + 1]
                );
                it2 += wt.length[j - i + 1];
            }
        }
    }

    fmt::printf("\n");
}
}  // namespace mc::dsp
