#include "WaveletTree.hpp"

#include "mc/dsp/wavelets/common.hpp"

#include "mc/cassert.hpp"
#include "mc/cmath.hpp"
#include "mc/format.hpp"
#include "mc/string_view.hpp"
#include "mc/utility.hpp"

namespace mc::dsp
{

WaveletTree::WaveletTree(Wavelet* waveIn, std::size_t signalLength, std::size_t j) : wave{waveIn}
{
    auto const size    = wave->size();
    auto const maxIter = maxIterations(signalLength, size);

    if (j > 100)
    {
        fmt::print("\n The Decomposition Iterations Cannot Exceed 100. Exiting \n");
        exit(-1);
    }
    if (j > maxIter)
    {

        fmt::print("\n Error - The Signal Can only be iterated {0} times using this wavelet. Exiting\n", maxIter);
        exit(-1);
    }

    std::size_t temp    = 1;
    std::size_t elength = 0;
    std::size_t nodess  = 0;
    for (std::size_t i = 0; i < j; ++i)
    {
        temp *= 2;
        nodess += temp;
        auto const temp2 = (size - 2) * (temp - 1);
        elength += temp2;
    }

    this->params    = std::make_unique<float[]>(signalLength * (j + 1) + elength + nodess + j + 1);
    this->outlength = signalLength * (j + 1) + elength;
    this->ext_      = "sym";

    this->siglength = signalLength;
    this->J         = j;
    this->MaxIter   = maxIter;
    this->method    = "dwt";

    this->nodes = nodess;

    this->cfftset     = 0;
    this->lenlength   = j + 2;
    this->output      = &this->params[0];
    this->nodeLength_ = (unsigned*)&this->params[signalLength * (j + 1) + elength];
    this->coeflength  = (unsigned*)&this->params[signalLength * (j + 1) + elength + nodess];

    for (auto i = 0; cmp_less(i, signalLength * (j + 1) + elength + nodess + j + 1); ++i) { this->params[i] = 0.0F; }
}

static auto wtreePer(WaveletTree& wt, float const* inp, int n, float* cA, int lenCA, float* cD) -> void
{
    int l      = 0;
    int l2     = 0;
    int isodd  = 0;
    int t      = 0;
    int lenAvg = 0;

    lenAvg = wt.wave->lpd().size();
    l2     = lenAvg / 2;
    isodd  = n % 2;

    for (auto i = 0; i < lenCA; ++i)
    {
        t     = 2 * i + l2;
        cA[i] = 0.0F;
        cD[i] = 0.0F;
        for (l = 0; l < lenAvg; ++l)
        {
            if ((t - l) >= l2 && (t - l) < n)
            {
                cA[i] += wt.wave->lpd()[l] * inp[t - l];
                cD[i] += wt.wave->hpd()[l] * inp[t - l];
            }
            else if ((t - l) < l2 && (t - l) >= 0)
            {
                cA[i] += wt.wave->lpd()[l] * inp[t - l];
                cD[i] += wt.wave->hpd()[l] * inp[t - l];
            }
            else if ((t - l) < 0 && isodd == 0)
            {
                cA[i] += wt.wave->lpd()[l] * inp[t - l + n];
                cD[i] += wt.wave->hpd()[l] * inp[t - l + n];
            }
            else if ((t - l) < 0 && isodd == 1)
            {
                if ((t - l) != -1)
                {
                    cA[i] += wt.wave->lpd()[l] * inp[t - l + n + 1];
                    cD[i] += wt.wave->hpd()[l] * inp[t - l + n + 1];
                }
                else
                {
                    cA[i] += wt.wave->lpd()[l] * inp[n - 1];
                    cD[i] += wt.wave->hpd()[l] * inp[n - 1];
                }
            }
            else if ((t - l) >= n && isodd == 0)
            {
                cA[i] += wt.wave->lpd()[l] * inp[t - l - n];
                cD[i] += wt.wave->hpd()[l] * inp[t - l - n];
            }
            else if ((t - l) >= n && isodd == 1)
            {
                if (t - l != n)
                {
                    cA[i] += wt.wave->lpd()[l] * inp[t - l - (n + 1)];
                    cD[i] += wt.wave->hpd()[l] * inp[t - l - (n + 1)];
                }
                else
                {
                    cA[i] += wt.wave->lpd()[l] * inp[n - 1];
                    cD[i] += wt.wave->hpd()[l] * inp[n - 1];
                }
            }
        }
    }
}

static auto wtreeSym(WaveletTree& wt, float const* inp, int n, float* cA, int lenCA, float* cD) -> void
{
    int l      = 0;
    int t      = 0;
    int lenAvg = 0;

    lenAvg = wt.wave->lpd().size();

    for (auto i = 0; i < lenCA; ++i)
    {
        t     = 2 * i + 1;
        cA[i] = 0.0F;
        cD[i] = 0.0F;
        for (l = 0; l < lenAvg; ++l)
        {
            if ((t - l) >= 0 && (t - l) < n)
            {
                cA[i] += wt.wave->lpd()[l] * inp[t - l];
                cD[i] += wt.wave->hpd()[l] * inp[t - l];
            }
            else if ((t - l) < 0)
            {
                cA[i] += wt.wave->lpd()[l] * inp[-t + l - 1];
                cD[i] += wt.wave->hpd()[l] * inp[-t + l - 1];
            }
            else if ((t - l) >= n)
            {
                cA[i] += wt.wave->lpd()[l] * inp[2 * n - t + l - 1];
                cD[i] += wt.wave->hpd()[l] * inp[2 * n - t + l - 1];
            }
        }
    }
}

auto wtree(WaveletTree& wt, float const* inp) -> void
{
    std::size_t iter = 0;
    int n            = 0;
    int lp           = 0;
    int p2           = 0;
    int k            = 0;
    int n2           = 0;
    int np           = 0;
    int lenCA        = 0;
    int t            = 0;
    int t2           = 0;
    int it1          = 0;

    auto tempLen     = wt.siglength;
    auto j           = wt.J;
    wt.length[j + 1] = tempLen;
    wt.outlength     = 0;
    wt.zpad          = 0;

    auto orig = std::make_unique<float[]>(tempLen);

    for (std::size_t i = 0; i < wt.siglength; ++i) { orig[i] = inp[i]; }

    if (wt.zpad == 1) { orig[tempLen - 1] = orig[tempLen - 2]; }

    n  = tempLen;
    lp = wt.wave->lpd().size();

    if (wt.extension() == string_view{"per"})
    {
        auto i = j;
        p2     = 2;
        while (i > 0)
        {
            n            = (int)std::ceil((float)n / 2.0F);
            wt.length[i] = n;
            wt.outlength += p2 * (wt.length[i]);
            i--;
            p2 *= 2;
        }
        wt.length[0] = wt.length[1];

        n2 = wt.outlength;
        p2 = 1;
        for (iter = 0; iter < j; ++iter)
        {
            lenCA = wt.length[j - iter];
            n2 -= 2 * p2 * lenCA;
            n = n2;
            for (k = 0; k < p2; ++k)
            {
                if (iter == 0)
                {
                    wtreePer(wt, orig.get(), tempLen, wt.params.get() + n, lenCA, wt.params.get() + n + lenCA);
                }
                else
                {
                    wtreePer(wt, wt.params.get() + np + k * tempLen, tempLen, wt.params.get() + n, lenCA,
                             wt.params.get() + n + lenCA);
                }
                n += 2 * lenCA;
            }

            tempLen = wt.length[j - iter];
            p2      = 2 * p2;
            np      = n2;
        }
    }
    else if (wt.extension() == string_view{"sym"})
    {
        auto i = j;
        p2     = 2;
        while (i > 0)
        {
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

        for (iter = 0; iter < j; ++iter)
        {
            lenCA = wt.length[j - iter];
            n2 -= 2 * p2 * lenCA;
            n = n2;
            for (k = 0; k < p2; ++k)
            {
                if (iter == 0)
                {
                    wtreeSym(wt, orig.get(), tempLen, wt.params.get() + n, lenCA, wt.params.get() + n + lenCA);
                }
                else
                {
                    wtreeSym(wt, wt.params.get() + np + k * tempLen, tempLen, wt.params.get() + n, lenCA,
                             wt.params.get() + n + lenCA);
                }
                n += 2 * lenCA;
            }

            tempLen = wt.length[j - iter];
            p2      = 2 * p2;
            np      = n2;
        }
    }
    else
    {
        fmt::printf("Signal extension can be either per or sym");
        exit(-1);
    }

    j   = wt.J;
    t2  = wt.outlength - 2 * wt.length[j];
    p2  = 2;
    it1 = 0;
    for (std::size_t i = 0; i < j; ++i)
    {
        t = t2;
        for (k = 0; k < p2; ++k)
        {
            wt.nodeLength_[it1] = t;
            it1++;
            t += wt.length[j - i];
        }
        p2 *= 2;
        t2 = t2 - p2 * wt.length[j - i - 1];
    }

    wt.coeflength[0] = wt.siglength;

    for (std::size_t i = 1; i < j + 1; ++i) { wt.coeflength[i] = wt.length[j - i + 1]; }
}

/// X - Level. All Nodes at any level have the same length
auto WaveletTree::nodeLength(std::size_t x) -> std::size_t
{
    if (x <= 0 || x > J)
    {
        fmt::printf("X co-ordinate must be >= 1 and <= %d", J);
        exit(-1);
    }

    return length[J - x + 1];
}

auto WaveletTree::coeffs(std::size_t x, std::size_t y, float* coeffs, std::size_t n) const -> void
{
    std::size_t ymax = 0;
    int t            = 0;
    int t2           = 0;

    if (x <= 0 || x > J)
    {
        fmt::printf("X co-ordinate must be >= 1 and <= %d", J);
        exit(-1);
    }
    ymax = 1;
    for (std::size_t i = 0; i < x; ++i) { ymax *= 2; }

    ymax -= 1;

    if (y > ymax)
    {
        fmt::printf("Y co-ordinate must be >= 0 and <= %d", ymax);
        exit(-1);
    }

    if (x == 1) { t = 0; }
    else
    {
        t  = 0;
        t2 = 1;
        for (std::size_t i = 0; i < x - 1; ++i)
        {
            t2 *= 2;
            t += t2;
        }
    }

    t += y;
    t2 = nodeLength_[t];
    for (std::size_t i = 0; i < n; ++i) { coeffs[i] = output[t2 + i]; }
}

auto WaveletTree::extension(char const* newExtension) noexcept -> void
{
    MC_ASSERT((newExtension == string_view{"sym"}) || (newExtension == string_view{"per"}));
    ext_ = newExtension;
}

auto WaveletTree::extension() const noexcept -> std::string const& { return ext_; }

auto summary(WaveletTree const& wt) -> void
{
    int k  = 0;
    int p2 = 0;
    int j  = 0;
    int t  = 0;
    j      = wt.J;
    summary(*wt.wave);
    fmt::printf("\n");
    fmt::printf("Wavelet Transform : %s \n", wt.method.c_str());
    fmt::printf("\n");
    fmt::printf("Signal Extension : %s \n", wt.extension().c_str());
    fmt::printf("\n");
    fmt::printf("Number of Decomposition Levels %d \n", wt.J);
    fmt::printf("\n");
    fmt::printf("Length of Input Signal %d \n", wt.siglength);
    fmt::printf("\n");
    fmt::printf("Length of WT Output Vector %d \n", wt.outlength);
    fmt::printf("\n");
    fmt::printf("Wavelet Coefficients are contained in vector : %s \n", "output");
    fmt::printf("\n");
    fmt::printf("Coefficients Access \n");
    t  = 0;
    p2 = 2;
    for (auto i = 0; i < j; ++i)
    {
        for (k = 0; k < p2; ++k)
        {
            fmt::printf("Node %d %d Access : output[%d] Length : %d \n", i + 1, k, wt.nodeLength_[t], wt.length[j - i]);
            t++;
        }
        p2 *= 2;
    }
    fmt::printf("\n");
}
}  // namespace mc::dsp