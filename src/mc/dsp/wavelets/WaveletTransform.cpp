#include "WaveletTransform.hpp"

#include "mc/dsp/convolution/FFTConvolver.hpp"
#include "mc/dsp/convolution/convolute.hpp"
#include "mc/dsp/fft/FFT.hpp"
#include "mc/dsp/wavelets/common.hpp"

#include "mc/cassert.hpp"
#include "mc/cmath.hpp"
#include "mc/cstdlib.hpp"
#include "mc/cstring.hpp"
#include "mc/format.hpp"
#include "mc/iterator.hpp"
#include "mc/string_view.hpp"

namespace mc::dsp
{

namespace
{

auto upSample(float const* const x, std::size_t const lenx, std::size_t const m, float* const y) -> void
{
    if (m == 0) { std::copy(x, x + static_cast<std::size_t>(lenx), y); }

    auto j       = std::size_t{1};
    auto k       = std::size_t{0};
    auto const n = m * (lenx - 1U) + 1U;
    for (std::size_t i = 0; i < n; ++i)
    {
        j--;
        y[i] = 0.0F;
        if (j == 0)
        {
            y[i] = x[k];
            k++;
            j = m;
        }
    }
}

// Returns even numbered output. Last value is set to zero
template<typename SrcIt, typename DestIt>
auto upSampleEven(SrcIt srcF, SrcIt srcL, DestIt destF, std::size_t m) -> void
{
    using T = typename std::iterator_traits<SrcIt>::value_type;
    if (m == 0) { std::copy(srcF, srcL, destF); }

    auto it   = destF;
    auto last = destF + (m * static_cast<std::size_t>(std::distance(srcF, srcL)));
    auto j    = std::size_t{1};
    for (; it != last; ++it)
    {
        j--;
        *it = T(0);
        if (j == 0)
        {
            *it = *srcF;
            ++srcF;
            j = m;
        }
    }
}

auto isign(int n) -> int { return n >= 0 ? 1 : (-1); }

auto circshift(float* array, int n, int l) -> void
{
    if (std::abs(l) > n) { l = isign(l) * (std::abs(l) % n); }
    if (l < 0) { l = (n + l) % n; }

    auto temp = makeZeros<float>(l);
    std::copy(array, array + static_cast<std::size_t>(l), temp.get());

    for (auto i = 0; i < n - l; ++i) { array[i] = array[i + l]; }
    for (auto i = 0; i < l; ++i) { array[n - l + i] = temp[i]; }
}

auto downsamp(float const* x, std::size_t lenx, std::size_t m, float* y) -> void
{
    if (m == 0) { std::copy(x, x + static_cast<std::size_t>(lenx), y); }

    auto const n = (lenx - 1U) / m + 1U;
    for (std::size_t i = 0; i < n; ++i) { y[i] = x[i * m]; }
}

}  // namespace

WaveletTransform::WaveletTransform(Wavelet& w, char const* method, std::size_t siglength, std::size_t j)
    : wave_{&w}
    , levels_{j}
    , signalLength_{siglength}
    , method_{method}
    , convolver{nullptr}
    , modwtsiglength{siglength}
    , lenlength{levels_ + 2}
    , MaxIter{maxIterations(siglength, w.size())}

{
    auto const size = w.size();

    if (levels_ > 100) { throw std::invalid_argument("decomposition Iterations Cannot Exceed 100."); }

    if (levels_ > MaxIter)
    {
        throw std::invalid_argument("signal Can only be iterated maxIter times using this wavelet");
    }

    if (method == nullptr)
    {
        this->params    = std::make_unique<float[]>(siglength + 2 * levels_ * (size + 1));
        this->outlength = siglength + 2 * levels_ * (size + 1);
        ext_            = SignalExtension::symmetric;
    }
    else if ((method == string_view{"dwt"}) || (method == string_view{"DWT"}))
    {
        this->params    = std::make_unique<float[]>(siglength + 2 * levels_ * (size + 1));
        this->outlength = siglength + 2 * levels_ * (size + 1);
        ext_            = SignalExtension::symmetric;
    }
    else if ((method == string_view{"swt"}) || (method == string_view{"SWT"}))
    {
        if (testSWTlength(siglength, levels_) == 0)
        {
            throw std::invalid_argument("For SWT the signal length must be a multiple of 2^levels");
        }

        this->params    = std::make_unique<float[]>(siglength * (levels_ + 1));
        this->outlength = siglength * (levels_ + 1);
        ext_            = SignalExtension::periodic;
    }
    else if ((method == string_view{"MODWT"}) || (method == string_view{"modwt"}))
    {

        if (strstr(w.name().c_str(), "haar") == nullptr)
        {
            if (strstr(w.name().c_str(), "db") == nullptr)
            {
                if (strstr(w.name().c_str(), "sym") == nullptr)
                {
                    if (strstr(w.name().c_str(), "coif") == nullptr)
                    {
                        throw std::invalid_argument(
                            "MODWT is only implemented for orthogonal wavelet families - db, sym and coif");
                    }
                }
            }
        }

        this->params    = std::make_unique<float[]>(siglength * 2 * (levels_ + 1));
        this->outlength = siglength * (levels_ + 1);
        ext_            = SignalExtension::periodic;
    }

    this->output_ = &this->params[0];
    if ((method == string_view{"dwt"}) || (method == string_view{"DWT"}))
    {
        for (std::size_t i = 0; i < siglength + 2 * levels() * (size + 1); ++i) { this->params[i] = 0.0F; }
    }
    else if ((method == string_view{"swt"}) || (method == string_view{"SWT"}))
    {
        for (std::size_t i = 0; i < siglength * (levels() + 1); ++i) { this->params[i] = 0.0F; }
    }
    else if ((method == string_view{"MODWT"}) || (method == string_view{"modwt"}))
    {
        for (std::size_t i = 0; i < siglength * 2 * (levels() + 1); ++i) { this->params[i] = 0.0F; }
    }
}

auto WaveletTransform::convMethod(ConvolutionMethod method) -> void { cmethod_ = method; }

auto WaveletTransform::extension(SignalExtension ext) -> void
{
    MC_ASSERT((ext == SignalExtension::periodic) || (ext == SignalExtension::symmetric));
    ext_ = ext;
}

auto WaveletTransform::output() const -> span<float>
{
    return span<float>{output_, static_cast<std::size_t>(outlength)};
}

auto WaveletTransform::approx() const -> span<float>
{
    /*
    Wavelet decomposition is stored as
    [A(J) D(J) D(J-1) ..... D(1)] in wt->output vector

    Length of A(J) , N = wt->length[0]
    */

    return {output_, static_cast<size_t>(length[0])};
}

auto WaveletTransform::detail(std::size_t level) const -> span<float>
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

    if (level > static_cast<std::size_t>(levels()) || level < 1U)
    {
        throw std::invalid_argument("The decomposition only has 1,..,%d levels");
    }

    auto iter = length[0];
    for (auto i = 1U; i < levels() - level; ++i) { iter += length[i]; }

    return {&output_[iter], static_cast<size_t>(length[level])};
}

static auto wconv(WaveletTransform& wt, float* sig, std::size_t n, float const* filt, std::size_t l, float* oup) -> void
{
    if (wt.convMethod() == ConvolutionMethod::direct)
    {
        convolute(sig, n, filt, l, oup);
        return;
    }

    MC_ASSERT(wt.convMethod() == ConvolutionMethod::fft);
    if (wt.cfftset == 0)
    {
        wt.convolver = std::make_unique<FFTConvolver>(n, l);
        convolute(*wt.convolver, sig, filt, oup);
    }
    else { convolute(*wt.convolver, sig, filt, oup); }
}

static auto dwtPer(WaveletTransform& wt, float* inp, int n, float* cA, int lenCA, float* cD) -> void
{

    dwtPerStride(inp, n, wt.wave().lpd().data(), wt.wave().hpd().data(), wt.wave().lpd().size(), cA, lenCA, cD, 1, 1);
}

static auto dwtSym(WaveletTransform& wt, float* inp, int n, float* cA, int lenCA, float* cD) -> void
{

    dwtSymStride(inp, n, wt.wave().lpd().data(), wt.wave().hpd().data(), wt.wave().lpd().size(), cA, lenCA, cD, 1, 1);
}

static auto dwt1(WaveletTransform& wt, float* sig, int lenSig, float* cA, float* cD) -> void
{
    if (wt.extension() == SignalExtension::periodic)
    {
        auto lenAvg  = (wt.wave().lpd().size() + wt.wave().hpd().size()) / 2;
        auto signal  = std::make_unique<float[]>(lenSig + lenAvg + (lenSig % 2));
        lenSig       = dsp::periodicExtension(sig, lenSig, lenAvg / 2, signal.get());
        auto cAUndec = std::make_unique<float[]>(lenSig + lenAvg + wt.wave().lpd().size() - 1);

        if (wt.wave().lpd().size() == wt.wave().hpd().size() && (wt.convMethod() == ConvolutionMethod::fft))
        {
            wt.convolver = std::make_unique<FFTConvolver>(lenSig + lenAvg, wt.wave().lpd().size());
            wt.cfftset   = 1;
        }
        else if (!(wt.wave().lpd().size() == wt.wave().hpd().size()))
        {
            throw std::invalid_argument("decomposition filters must have the same length.");
        }

        wconv(wt, signal.get(), lenSig + lenAvg, wt.wave().lpd().data(), wt.wave().lpd().size(), cAUndec.get());
        downsamp(cAUndec.get() + lenAvg, static_cast<std::size_t>(lenSig), 2U, cA);
        wconv(wt, signal.get(), lenSig + lenAvg, wt.wave().hpd().data(), wt.wave().hpd().size(), cAUndec.get());
        downsamp(cAUndec.get() + lenAvg, static_cast<std::size_t>(lenSig), 2U, cD);
    }
    else if (wt.extension() == SignalExtension::symmetric)
    {
        auto lf      = wt.wave().lpd().size();  // lpd and hpd have the same length
        auto signal  = std::make_unique<float[]>(lenSig + 2 * (lf - 1));
        lenSig       = dsp::symmetricExtension(sig, lenSig, lf - 1, signal.get());
        auto cAUndec = std::make_unique<float[]>(lenSig + 3 * (lf - 1));

        if (wt.wave().lpd().size() == wt.wave().hpd().size() && (wt.convMethod() == ConvolutionMethod::fft))
        {
            wt.convolver = std::make_unique<FFTConvolver>(lenSig + 2 * (lf - 1), lf);
            wt.cfftset   = 1;
        }
        else if (!(wt.wave().lpd().size() == wt.wave().hpd().size()))
        {
            throw std::invalid_argument("decomposition filters must have the same length.");
        }

        wconv(wt, signal.get(), lenSig + 2 * (lf - 1), wt.wave().lpd().data(), wt.wave().lpd().size(), cAUndec.get());
        downsamp(cAUndec.get() + lf, static_cast<std::size_t>(lenSig + lf - 2), 2, cA);
        wconv(wt, signal.get(), lenSig + 2 * (lf - 1), wt.wave().hpd().data(), wt.wave().hpd().size(), cAUndec.get());
        downsamp(cAUndec.get() + lf, static_cast<std::size_t>(lenSig + lf - 2), 2, cD);
    }
    else { throw std::invalid_argument("Signal extension can be either per or sym"); }

    if (wt.wave().lpd().size() == wt.wave().hpd().size() && (wt.convMethod() == ConvolutionMethod::fft))
    {

        wt.cfftset = 0;
    }
}

auto dwt(WaveletTransform& wt, float const* inp) -> void
{

    auto tempLen = wt.signalLength();
    auto const j = wt.levels();

    wt.length[j + 1] = tempLen;
    wt.outlength     = 0;
    wt.zpad          = 0;

    auto orig2 = std::make_unique<float[]>(tempLen);
    auto orig  = std::make_unique<float[]>(tempLen);

    std::copy(inp, inp + wt.signalLength(), orig.get());

    if (wt.zpad == 1) { orig[tempLen - 1] = orig[tempLen - 2]; }

    auto n  = tempLen;
    auto lp = wt.wave().lpd().size();

    if (wt.extension() == SignalExtension::periodic)
    {
        auto idx = j;
        while (idx > 0)
        {
            n              = (int)std::ceil((float)n / 2.0F);
            wt.length[idx] = n;
            wt.outlength += wt.length[idx];
            idx--;
        }
        wt.length[0] = wt.length[1];
        wt.outlength += wt.length[0];
        n = wt.outlength;

        for (auto iter = 0; iter < j; ++iter)
        {
            auto const lenCA = wt.length[j - iter];
            n -= lenCA;
            if (wt.convMethod() == ConvolutionMethod::fft)
            {
                dwt1(wt, orig.get(), tempLen, orig2.get(), wt.params.get() + n);
            }
            else { dwtPer(wt, orig.get(), tempLen, orig2.get(), lenCA, wt.params.get() + n); }
            tempLen = wt.length[j - iter];
            if (iter == j - 1)
            {
                for (std::size_t i = 0; i < lenCA; ++i) { wt.params[i] = orig2[i]; }
            }
            else
            {
                for (std::size_t i = 0; i < lenCA; ++i) { orig[i] = orig2[i]; }
            }
        }
    }
    else if (wt.extension() == SignalExtension::symmetric)
    {
        auto idx = j;
        while (idx > 0)
        {
            n              = n + lp - 2;
            n              = (int)std::ceil((float)n / 2.0F);
            wt.length[idx] = n;
            wt.outlength += wt.length[idx];
            idx--;
        }
        wt.length[0] = wt.length[1];
        wt.outlength += wt.length[0];
        n = wt.outlength;

        for (auto iter = 0; iter < j; ++iter)
        {
            auto const lenCA = wt.length[j - iter];
            n -= lenCA;
            if (wt.convMethod() == ConvolutionMethod::fft)
            {
                dwt1(wt, orig.get(), tempLen, orig2.get(), wt.params.get() + n);
            }
            else { dwtSym(wt, orig.get(), tempLen, orig2.get(), lenCA, wt.params.get() + n); }
            tempLen = wt.length[j - iter];

            if (iter == j - 1)
            {
                for (std::size_t i = 0; i < lenCA; ++i) { wt.params[i] = orig2[i]; }
            }
            else
            {
                for (std::size_t i = 0; i < lenCA; ++i) { orig[i] = orig2[i]; }
            }
        }
    }
    else { throw std::invalid_argument("Signal extension can be either per or sym"); }
}

static auto idwt1(WaveletTransform& wt, float* temp, float* cAUp, float* cA, int lenCA, float* cD, int lenCD,
                  float* xLp, float* xHp, float* x) -> void
{
    auto lenAvg = (wt.wave().lpr().size() + wt.wave().hpr().size()) / 2;
    auto n      = 2 * lenCD;
    auto u      = 2;

    upSampleEven(cA, cA + lenCA, cAUp, u);

    dsp::periodicExtension(cAUp, 2 * lenCA, lenAvg / 2, temp);

    auto n2 = 2 * lenCA + lenAvg;

    if (wt.wave().lpr().size() == wt.wave().hpr().size() && (wt.convMethod() == ConvolutionMethod::fft))
    {
        wt.convolver = std::make_unique<FFTConvolver>(n2, lenAvg);
        wt.cfftset   = 1;
    }
    else if (!(wt.wave().lpr().size() == wt.wave().hpr().size()))
    {
        throw std::invalid_argument("Decomposition Filters must have the same length");
    }

    wconv(wt, temp, n2, wt.wave().lpr().data(), lenAvg, xLp);

    upSampleEven(cD, cD + lenCD, cAUp, u);

    dsp::periodicExtension(cAUp, 2 * lenCD, lenAvg / 2, temp);

    n2 = 2 * lenCD + lenAvg;

    wconv(wt, temp, n2, wt.wave().hpr().data(), lenAvg, xHp);

    for (auto i = lenAvg - 1; i < n + lenAvg - 1; ++i) { x[i - lenAvg + 1] = xLp[i] + xHp[i]; }

    if (wt.wave().lpr().size() == wt.wave().hpr().size() && (wt.convMethod() == ConvolutionMethod::fft))
    {

        wt.cfftset = 0;
    }
}

static auto idwtPer(WaveletTransform& wt, float* cA, int lenCA, float* cD, float* x) -> void
{
    idwtPerStride(cA, lenCA, cD, wt.wave().lpr().data(), wt.wave().hpr().data(), wt.wave().lpr().size(), x, 1, 1);
}

static auto idwtSym(WaveletTransform& wt, float* cA, int lenCA, float* cD, float* x) -> void
{
    idwtSymStride(cA, lenCA, cD, wt.wave().lpr().data(), wt.wave().hpr().data(), wt.wave().lpr().size(), x, 1, 1);
}

auto idwt(WaveletTransform& wt, float* dwtop) -> void
{
    std::size_t lf     = 0;
    std::size_t n      = 0;
    std::size_t n2     = 0;
    std::size_t iter   = 0;
    std::size_t k      = 0;
    std::size_t detLen = 0;

    auto j      = wt.levels();
    auto u      = 2;
    auto appLen = wt.length[0];
    auto out    = std::make_unique<float[]>(wt.signalLength() + 1);

    if ((wt.extension() == SignalExtension::periodic) && (wt.convMethod() == ConvolutionMethod::fft))
    {
        appLen = wt.length[0];
        detLen = wt.length[1];
        n      = 2 * wt.length[j];
        lf     = (wt.wave().lpr().size() + wt.wave().hpr().size()) / 2;

        auto cAUp = std::make_unique<float[]>(n);
        auto temp = std::make_unique<float[]>((n + lf));
        auto xLp  = std::make_unique<float[]>((n + 2 * lf - 1));
        auto xHp  = std::make_unique<float[]>((n + 2 * lf - 1));
        iter      = appLen;

        for (std::size_t i = 0; i < appLen; ++i) { out[i] = wt.output()[i]; }

        for (auto i = 0; i < j; ++i)
        {
            idwt1(wt, temp.get(), cAUp.get(), out.get(), detLen, wt.output().data() + iter, detLen, xLp.get(),
                  xHp.get(), out.get());
            iter += detLen;
            detLen = wt.length[i + 2];
        }
    }
    else if ((wt.extension() == SignalExtension::periodic) && (wt.convMethod() == ConvolutionMethod::direct))
    {
        appLen = wt.length[0];
        detLen = wt.length[1];
        n      = 2 * wt.length[j];
        lf     = (wt.wave().lpr().size() + wt.wave().hpr().size()) / 2;

        auto xLp = std::make_unique<float[]>((n + 2 * lf - 1));
        iter     = appLen;

        for (std::size_t i = 0; i < appLen; ++i) { out[i] = wt.output()[i]; }

        for (auto i = 0; i < j; ++i)
        {
            idwtPer(wt, out.get(), detLen, wt.output().data() + iter, xLp.get());
            for (k = lf / 2 - 1; k < 2 * detLen + lf / 2 - 1; ++k) { out[k - lf / 2 + 1] = xLp[k]; }

            iter += detLen;
            detLen = wt.length[i + 2];
        }
    }
    else if ((wt.extension() == SignalExtension::symmetric) && (wt.convMethod() == ConvolutionMethod::direct))
    {
        appLen = wt.length[0];
        detLen = wt.length[1];
        n      = 2 * wt.length[j] - 1;
        lf     = (wt.wave().lpr().size() + wt.wave().hpr().size()) / 2;

        auto xLp = std::make_unique<float[]>((n + 2 * lf - 1));
        iter     = appLen;

        for (std::size_t i = 0; i < appLen; ++i) { out[i] = wt.output()[i]; }

        for (auto i = 0; i < j; ++i)
        {
            idwtSym(wt, out.get(), detLen, wt.output().data() + iter, xLp.get());
            for (k = lf - 2; k < 2 * detLen; ++k) { out[k - lf + 2] = xLp[k]; }

            iter += detLen;
            detLen = wt.length[i + 2];
        }
    }
    else if ((wt.extension() == SignalExtension::symmetric) && (wt.convMethod() == ConvolutionMethod::fft))
    {
        lf = wt.wave().lpd().size();  // lpd and hpd have the same length

        n         = 2 * wt.length[j] - 1;
        auto cAUp = std::make_unique<float[]>(n);
        auto xLp  = std::make_unique<float[]>((n + lf - 1));
        auto xHp  = std::make_unique<float[]>((n + lf - 1));

        for (std::size_t i = 0; i < appLen; ++i) { out[i] = wt.output()[i]; }

        iter = appLen;

        for (auto i = 0; i < j; ++i)
        {
            detLen = wt.length[i + 1];
            upSample(out.get(), detLen, u, cAUp.get());
            n2 = 2 * wt.length[i + 1] - 1;

            if (wt.wave().lpr().size() == wt.wave().hpr().size() && (wt.convMethod() == ConvolutionMethod::fft))
            {
                wt.convolver = std::make_unique<FFTConvolver>(n2, lf);
                wt.cfftset   = 1;
            }
            else if (!(wt.wave().lpr().size() == wt.wave().hpr().size()))
            {
                throw std::invalid_argument("Decomposition Filters must have the same length");
            }

            wconv(wt, cAUp.get(), n2, wt.wave().lpr().data(), lf, xLp.get());
            upSample(wt.output().data() + iter, detLen, u, cAUp.get());
            wconv(wt, cAUp.get(), n2, wt.wave().hpr().data(), lf, xHp.get());

            for (k = lf - 2; k < n2 + 1; ++k) { out[k - lf + 2] = xLp[k] + xHp[k]; }
            iter += detLen;
            if (wt.wave().lpr().size() == wt.wave().hpr().size() && (wt.convMethod() == ConvolutionMethod::fft))
            {

                wt.cfftset = 0;
            }
        }
    }
    else { throw std::invalid_argument("Signal extension can be either per or sym"); }

    std::copy(out.get(), out.get() + wt.signalLength(), dwtop);
}

static auto swtPer(WaveletTransform& wt, int m, float* inp, int n, float* cA, int lenCA, float* cD) -> void
{

    swtPerStride(m, inp, n, wt.wave().lpd().data(), wt.wave().hpd().data(), wt.wave().lpd().size(), cA, lenCA, cD, 1,
                 1);
}

static auto swtFft(WaveletTransform& wt, float const* inp) -> void
{

    auto tempLen = wt.signalLength();
    auto j       = wt.levels();
    wt.length[0] = wt.length[j] = tempLen;
    wt.outlength = wt.length[j + 1] = (j + 1) * tempLen;
    auto m                          = 1;
    for (auto iter = 1; iter < j; ++iter)
    {
        m               = 2 * m;
        wt.length[iter] = tempLen;
    }

    auto const lenFilt = wt.wave().size();

    auto lowPass  = std::make_unique<float[]>(m * lenFilt);
    auto highPass = std::make_unique<float[]>(m * lenFilt);
    auto sig      = std::make_unique<float[]>((m * lenFilt + tempLen + (tempLen % 2)));
    auto cA       = std::make_unique<float[]>((2 * m * lenFilt + tempLen + (tempLen % 2)) - 1);
    auto cD       = std::make_unique<float[]>((2 * m * lenFilt + tempLen + (tempLen % 2)) - 1);

    m = 1;

    for (std::size_t i = 0; i < tempLen; ++i) { wt.params[i] = inp[i]; }

    auto lenacc = wt.outlength;
    std::size_t n{0};

    for (auto iter = 0; iter < j; ++iter)
    {
        lenacc -= tempLen;
        if (iter > 0)
        {
            m = 2 * m;
            n = m * lenFilt;
            upSampleEven(begin(wt.wave().lpd()), end(wt.wave().lpd()), lowPass.get(), m);
            upSampleEven(begin(wt.wave().hpd()), end(wt.wave().hpd()), highPass.get(), m);
        }
        else
        {
            n = lenFilt;
            for (std::size_t i = 0; i < n; ++i)
            {
                lowPass[i]  = wt.wave().lpd()[i];
                highPass[i] = wt.wave().hpd()[i];
            }
        }

        // swt_per(wt,M, wt.params.get(), temp_len, cA, temp_len, cD,temp_len);

        dsp::periodicExtension(wt.params.get(), tempLen, n / 2, sig.get());

        if (wt.wave().lpd().size() == wt.wave().hpd().size() && (wt.convMethod() == ConvolutionMethod::fft))
        {
            wt.convolver = std::make_unique<FFTConvolver>(n + tempLen + (tempLen % 2), n);
            wt.cfftset   = 1;
        }
        else if (!(wt.wave().lpd().size() == wt.wave().hpd().size()))
        {
            throw std::invalid_argument("Decomposition Filters must have the same length");
        }

        wconv(wt, sig.get(), n + tempLen + (tempLen % 2), lowPass.get(), n, cA.get());

        wconv(wt, sig.get(), n + tempLen + (tempLen % 2), highPass.get(), n, cD.get());

        if (wt.wave().lpd().size() == wt.wave().hpd().size() && (wt.convMethod() == ConvolutionMethod::fft))
        {

            wt.cfftset = 0;
        }

        for (std::size_t i = 0; i < tempLen; ++i)
        {
            wt.params[i]          = cA[n + i];
            wt.params[lenacc + i] = cD[n + i];
        }
    }
}

static auto swtDirect(WaveletTransform& wt, float const* inp) -> void
{
    auto tempLen = wt.signalLength();
    auto j       = wt.levels();

    wt.length[0] = wt.length[j] = tempLen;
    wt.outlength = wt.length[j + 1] = (j + 1) * tempLen;

    auto m = 1;
    for (auto iter = 1; iter < j; ++iter)
    {
        m               = 2 * m;
        wt.length[iter] = tempLen;
    }

    auto cA = std::make_unique<float[]>(tempLen);
    auto cD = std::make_unique<float[]>(tempLen);

    m = 1;

    for (std::size_t i = 0; i < tempLen; ++i) { wt.params[i] = inp[i]; }

    auto lenacc = wt.outlength;

    for (auto iter = 0; iter < j; ++iter)
    {
        lenacc -= tempLen;
        if (iter > 0) { m = 2 * m; }

        swtPer(wt, m, wt.params.get(), tempLen, cA.get(), tempLen, cD.get());

        for (std::size_t i = 0; i < tempLen; ++i)
        {
            wt.params[i]          = cA[i];
            wt.params[lenacc + i] = cD[i];
        }
    }
}

auto swt(WaveletTransform& wt, float const* inp) -> void
{
    if ((wt.method() == string_view{"swt"}) && (wt.convMethod() == ConvolutionMethod::direct)) { swtDirect(wt, inp); }
    else if ((wt.method() == string_view{"swt"}) && (wt.convMethod() == ConvolutionMethod::fft)) { swtFft(wt, inp); }
    else { throw std::invalid_argument("SWT Only accepts two methods - direct and fft"); }
}

auto iswt(WaveletTransform& wt, float* swtop) -> void
{
    auto n  = wt.signalLength();
    auto j  = static_cast<std::size_t>(wt.levels());
    auto u  = 2;
    auto lf = wt.wave().lpr().size();

    auto appxSig = std::make_unique<float[]>(n);
    auto detSig  = std::make_unique<float[]>(n);
    auto appx1   = std::make_unique<float[]>(n);
    auto det1    = std::make_unique<float[]>(n);
    auto appx2   = std::make_unique<float[]>(n);
    auto det2    = std::make_unique<float[]>(n);
    auto tempx   = std::make_unique<float[]>(n);
    auto cL0     = std::make_unique<float[]>((n + (n % 2) + lf));
    auto cH0     = std::make_unique<float[]>((n + (n % 2) + lf));
    auto oup00L  = std::make_unique<float[]>((n + 2 * lf));
    auto oup00H  = std::make_unique<float[]>((n + 2 * lf));
    auto oup00   = std::make_unique<float[]>(n);
    auto oup01   = std::make_unique<float[]>(n);

    for (std::size_t iter = 0; iter < j; ++iter)
    {
        for (std::size_t i = 0; i < n; ++i) { swtop[i] = 0.0F; }
        if (iter == 0)
        {
            for (std::size_t i = 0; i < n; ++i)
            {
                appxSig[i] = wt.output()[i];
                detSig[i]  = wt.output()[n + i];
            }
        }
        else
        {
            for (std::size_t i = 0; i < n; ++i) { detSig[i] = wt.output()[(iter + 1) * n + i]; }
        }

        auto const value = (int)std::pow(2.0F, (float)(j - 1 - iter));

        for (auto count = 0; count < value; count++)
        {
            auto len = 0;
            for (std::size_t index = count; index < n; index += value)
            {
                appx1[len] = appxSig[index];
                det1[len]  = detSig[index];
                len++;
            }

            // SHIFT 0
            auto len0 = 0;

            for (auto indexShift = 0; indexShift < len; indexShift += 2)
            {
                appx2[len0] = appx1[indexShift];
                det2[len0]  = det1[indexShift];
                len0++;
            }
            upSampleEven(appx2.get(), appx2.get() + len0, tempx.get(), u);
            dsp::periodicExtension(tempx.get(), 2 * len0, lf / 2, cL0.get());

            upSampleEven(det2.get(), det2.get() + len0, tempx.get(), u);
            dsp::periodicExtension(tempx.get(), 2 * len0, lf / 2, cH0.get());

            auto n1 = 2 * len0 + lf;

            if (wt.wave().lpr().size() == wt.wave().hpr().size() && (wt.convMethod() == ConvolutionMethod::fft))
            {
                wt.convolver = std::make_unique<FFTConvolver>(n1, lf);
                wt.cfftset   = 1;
            }
            else if (!(wt.wave().lpd().size() == wt.wave().hpd().size()))
            {
                throw std::invalid_argument("Decomposition Filters must have the same length");
            }

            wconv(wt, cL0.get(), n1, wt.wave().lpr().data(), lf, oup00L.get());

            wconv(wt, cH0.get(), n1, wt.wave().hpr().data(), lf, oup00H.get());

            for (auto i = lf - 1; i < 2 * len0 + lf - 1; ++i) { oup00[i - lf + 1] = oup00L[i] + oup00H[i]; }

            // SHIFT 1

            len0 = 0;

            for (auto indexShift = 1; indexShift < len; indexShift += 2)
            {
                appx2[len0] = appx1[indexShift];
                det2[len0]  = det1[indexShift];
                len0++;
            }

            upSampleEven(appx2.get(), appx2.get() + len0, tempx.get(), u);
            dsp::periodicExtension(tempx.get(), 2 * len0, lf / 2, cL0.get());

            upSampleEven(det2.get(), det2.get() + len0, tempx.get(), u);
            dsp::periodicExtension(tempx.get(), 2 * len0, lf / 2, cH0.get());

            n1 = 2 * len0 + lf;

            wconv(wt, cL0.get(), n1, wt.wave().lpr().data(), lf, oup00L.get());
            wconv(wt, cH0.get(), n1, wt.wave().hpr().data(), lf, oup00H.get());

            for (auto i = lf - 1; i < 2 * len0 + lf - 1; ++i) { oup01[i - lf + 1] = oup00L[i] + oup00H[i]; }

            circshift(oup01.get(), 2 * len0, -1);

            auto index2 = 0;

            for (auto index = static_cast<std::size_t>(count); index < n; index += value)
            {
                swtop[index] = (oup00[index2] + oup01[index2]) / 2.0F;
                index2++;
            }
        }
        for (std::size_t i = 0; i < n; ++i) { appxSig[i] = swtop[i]; }
    }
}

static auto modwtPer(WaveletTransform& wt, int m, float const* inp, float* cA, int lenCA, float* cD) -> void
{
    auto const lenAvg = wt.wave().lpd().size();
    auto filt         = std::make_unique<float[]>(2 * lenAvg);
    auto s            = std::sqrt(2.0F);

    for (std::size_t i = 0; i < lenAvg; ++i)
    {
        filt[i]          = wt.wave().lpd()[i] / s;
        filt[lenAvg + i] = wt.wave().hpd()[i] / s;
    }

    for (auto i = 0; i < lenCA; ++i)
    {
        auto t = i;
        cA[i]  = filt[0] * inp[t];
        cD[i]  = filt[lenAvg] * inp[t];
        for (std::size_t l = 1; l < lenAvg; l++)
        {
            t -= m;
            while (t >= lenCA) { t -= lenCA; }
            while (t < 0) { t += lenCA; }

            cA[i] += filt[l] * inp[t];
            cD[i] += filt[lenAvg + l] * inp[t];
        }
    }
}

static auto modwtDirect(WaveletTransform& wt, float const* inp) -> void
{
    if (wt.extension() != SignalExtension::periodic)
    {
        throw std::invalid_argument("MODWT direct method only uses periodic extension per.");
    }

    auto tempLen = wt.signalLength();
    auto j       = static_cast<std::size_t>(wt.levels());
    wt.length[0] = wt.length[j] = tempLen;
    wt.outlength = wt.length[j + 1] = (j + 1) * tempLen;
    auto m                          = 1;
    for (std::size_t iter = 1; iter < j; ++iter)
    {
        m               = 2 * m;
        wt.length[iter] = tempLen;
    }

    auto cA = std::make_unique<float[]>(tempLen);
    auto cD = std::make_unique<float[]>(tempLen);

    m = 1;

    for (std::size_t i = 0; i < tempLen; ++i) { wt.params[i] = inp[i]; }

    auto lenacc = wt.outlength;

    for (std::size_t iter = 0; iter < j; ++iter)
    {
        lenacc -= tempLen;
        if (iter > 0) { m = 2 * m; }

        modwtPer(wt, m, wt.params.get(), cA.get(), tempLen, cD.get());

        for (std::size_t i = 0; i < tempLen; ++i)
        {
            wt.params[i]          = cA[i];
            wt.params[lenacc + i] = cD[i];
        }
    }
}

static auto modwtFft(WaveletTransform& wt, float const* inp) -> void
{
    int j      = 0;
    int iter   = 0;
    int m      = 0;
    int lenacc = 0;
    float s    = NAN;
    float tmp1 = NAN;
    float tmp2 = NAN;

    auto tempLen = wt.signalLength();
    auto lenAvg  = static_cast<std::size_t>(wt.wave().lpd().size());
    std::size_t n{0};
    if (wt.extension() == SignalExtension::symmetric) { n = 2 * tempLen; }
    else if (wt.extension() == SignalExtension::periodic) { n = tempLen; }
    j                 = wt.levels();
    wt.modwtsiglength = n;
    wt.length[0] = wt.length[j] = n;
    wt.outlength = wt.length[j + 1] = (j + 1) * n;

    s = std::sqrt(2.0F);
    for (iter = 1; iter < j; ++iter) { wt.length[iter] = n; }

    auto fftFd = std::make_unique<FFT<float, KissFFT>>(n, FFTDirection::forward);
    auto fftBd = std::make_unique<FFT<float, KissFFT>>(n, FFTDirection::backward);

    auto sig      = std::make_unique<Complex<float>[]>(n);
    auto cA       = std::make_unique<Complex<float>[]>(n);
    auto cD       = std::make_unique<Complex<float>[]>(n);
    auto lowPass  = std::make_unique<Complex<float>[]>(n);
    auto highPass = std::make_unique<Complex<float>[]>(n);
    auto index    = std::make_unique<int[]>(n);

    // N-point FFT of low pass and high pass filters

    // Low Pass Filter

    for (std::size_t i = 0; i < lenAvg; ++i)
    {
        sig[i].real((float)wt.wave().lpd()[i] / s);
        sig[i].imag(0.0F);
    }
    for (auto i = lenAvg; i < n; ++i)
    {
        sig[i].real(0.0F);
        sig[i].imag(0.0F);
    }

    fftFd->perform(sig.get(), lowPass.get());

    // High Pass Filter

    for (std::size_t i = 0; i < lenAvg; ++i)
    {
        sig[i].real((float)wt.wave().hpd()[i] / s);
        sig[i].imag(0.0F);
    }
    for (std::size_t i = lenAvg; i < n; ++i)
    {
        sig[i].real(0.0F);
        sig[i].imag(0.0F);
    }

    fftFd->perform(sig.get(), highPass.get());

    // symmetric extension
    for (std::size_t i = 0; i < tempLen; ++i)
    {
        sig[i].real((float)inp[i]);
        sig[i].imag(0.0F);
    }
    for (std::size_t i = tempLen; i < n; ++i)
    {
        sig[i].real((float)inp[n - i - 1]);
        sig[i].imag(0.0F);
    }

    // FFT of data

    fftFd->perform(sig.get(), cA.get());

    lenacc = wt.outlength;

    m = 1;

    for (iter = 0; iter < j; ++iter)
    {
        lenacc -= n;

        for (std::size_t i = 0; i < n; ++i) { index[i] = (m * i) % n; }

        for (std::size_t i = 0; i < n; ++i)
        {
            tmp1 = cA[i].real();
            tmp2 = cA[i].imag();
            cA[i].real(lowPass[index[i]].real() * tmp1 - lowPass[index[i]].imag() * tmp2);
            cA[i].imag(lowPass[index[i]].real() * tmp2 + lowPass[index[i]].imag() * tmp1);

            cD[i].real(highPass[index[i]].real() * tmp1 - highPass[index[i]].imag() * tmp2);
            cD[i].imag(highPass[index[i]].real() * tmp2 + highPass[index[i]].imag() * tmp1);
        }

        fftBd->perform(cD.get(), sig.get());

        for (std::size_t i = 0; i < n; ++i) { wt.params[lenacc + i] = sig[i].real() / static_cast<float>(n); }

        m *= 2;
    }

    fftBd->perform(cA.get(), sig.get());

    for (std::size_t i = 0; i < n; ++i) { wt.params[i] = sig[i].real() / static_cast<float>(n); }
}

auto modwt(WaveletTransform& wt, float const* inp) -> void
{
    if (wt.convMethod() == ConvolutionMethod::direct)
    {
        modwtDirect(wt, inp);
        return;
    }

    modwtFft(wt, inp);
}

static auto conjComplex(Complex<float>* x, int n) -> void
{
    for (auto i = 0; i < n; ++i) { x[i].imag(x[i].imag() * -1.0F); }
}

auto imodwtFft(WaveletTransform& wt, float* oup) -> void
{
    auto n      = wt.modwtsiglength;
    auto lenAvg = static_cast<std::size_t>(wt.wave().lpd().size());
    auto j      = static_cast<std::size_t>(wt.levels());

    auto s     = std::sqrt(2.0F);
    auto fftFd = std::make_unique<FFT<float, KissFFT>>(n, FFTDirection::forward);
    auto fftBd = std::make_unique<FFT<float, KissFFT>>(n, FFTDirection::backward);

    auto sig      = std::make_unique<Complex<float>[]>(n);
    auto cA       = std::make_unique<Complex<float>[]>(n);
    auto cD       = std::make_unique<Complex<float>[]>(n);
    auto lowPass  = std::make_unique<Complex<float>[]>(n);
    auto highPass = std::make_unique<Complex<float>[]>(n);
    auto index    = std::make_unique<std::size_t[]>(n);

    // N-point FFT of low pass and high pass filters

    // Low Pass Filter

    for (std::size_t i = 0; i < lenAvg; ++i)
    {
        sig[i].real((float)wt.wave().lpd()[i] / s);
        sig[i].imag(0.0F);
    }
    for (std::size_t i = lenAvg; i < n; ++i)
    {
        sig[i].real(0.0F);
        sig[i].imag(0.0F);
    }

    fftFd->perform(sig.get(), lowPass.get());

    // High Pass Filter

    for (std::size_t i = 0; i < lenAvg; ++i)
    {
        sig[i].real((float)wt.wave().hpd()[i] / s);
        sig[i].imag(0.0F);
    }
    for (std::size_t i = lenAvg; i < n; ++i)
    {
        sig[i].real(0.0F);
        sig[i].imag(0.0F);
    }

    fftFd->perform(sig.get(), highPass.get());

    // Complex conjugate of the two filters

    conjComplex(lowPass.get(), static_cast<int>(n));
    conjComplex(highPass.get(), static_cast<int>(n));

    auto m      = (int)std::pow(2.0F, (float)j - 1.0F);
    auto lenacc = n;

    //
    for (std::size_t i = 0; i < n; ++i)
    {
        sig[i].real((float)wt.output()[i]);
        sig[i].imag(0.0F);
    }

    for (std::size_t iter = 0; iter < j; ++iter)
    {
        fftFd->perform(sig.get(), cA.get());
        for (std::size_t i = 0; i < n; ++i)
        {
            sig[i].real(wt.output()[lenacc + i]);
            sig[i].imag(0.0F);
        }
        fftFd->perform(sig.get(), cD.get());

        for (std::size_t i = 0; i < n; ++i) { index[i] = (m * i) % n; }

        for (std::size_t i = 0; i < n; ++i)
        {
            auto const tmp1 = cA[i].real();
            auto const tmp2 = cA[i].imag();
            cA[i].real(lowPass[index[i]].real() * tmp1 - lowPass[index[i]].imag() * tmp2
                       + highPass[index[i]].real() * cD[i].real() - highPass[index[i]].imag() * cD[i].imag());
            cA[i].imag(lowPass[index[i]].real() * tmp2 + lowPass[index[i]].imag() * tmp1
                       + highPass[index[i]].real() * cD[i].imag() + highPass[index[i]].imag() * cD[i].real());
        }

        fftBd->perform(cA.get(), sig.get());

        for (std::size_t i = 0; i < n; ++i)
        {
            sig[i].real(sig[i].real() / static_cast<float>(n));
            sig[i].imag(sig[i].imag() / static_cast<float>(n));
        }
        m /= 2;
        lenacc += n;
    }

    std::transform(sig.get(), sig.get() + wt.signalLength(), oup, [](auto c) { return c.real(); });
}

static auto imodwtPer(WaveletTransform& wt, int m, float const* cA, int lenCA, float const* cD, float* x) -> void
{
    auto const lenAvg = wt.wave().lpd().size();
    auto filt         = std::make_unique<float[]>(2 * lenAvg);
    auto s            = std::sqrt(2.0F);

    for (std::size_t i = 0; i < lenAvg; ++i)
    {
        filt[i]          = wt.wave().lpd()[i] / s;
        filt[lenAvg + i] = wt.wave().hpd()[i] / s;
    }

    for (auto i = 0; i < lenCA; ++i)
    {
        auto t = i;
        x[i]   = (filt[0] * cA[t]) + (filt[lenAvg] * cD[t]);
        for (std::size_t l = 1; l < lenAvg; l++)
        {
            t += m;
            while (t >= lenCA) { t -= lenCA; }
            while (t < 0) { t += lenCA; }

            x[i] += (filt[l] * cA[t]) + (filt[lenAvg + l] * cD[t]);
        }
    }
}

static auto imodwtDirect(WaveletTransform& wt, float* dwtop) -> void
{
    auto n      = wt.signalLength();
    auto lenacc = n;

    auto j = static_cast<std::size_t>(wt.levels());

    auto x = std::make_unique<float[]>(n);

    for (std::size_t i = 0; i < n; ++i) { dwtop[i] = wt.output()[i]; }

    auto m = static_cast<int>(std::pow(2.0F, (float)j - 1.0F));
    for (std::size_t iter = 0; iter < j; ++iter)
    {
        if (iter > 0) { m = m / 2; }
        imodwtPer(wt, m, dwtop, static_cast<int>(n), wt.params.get() + lenacc, x.get());
        /*
            for (std::size_t j = lf - 1; j < N; ++j) {
                    dwtop[j - lf + 1] = X[j];
            }
            for (std::size_t j = 0; j < lf - 1; ++j) {
                    dwtop[N - lf + 1 + j] = X[j];
            }
            */
        for (std::size_t jj = 0; jj < n; ++jj) { dwtop[jj] = x[jj]; }

        lenacc += n;
    }
}

auto imodwt(WaveletTransform& wt, float* oup) -> void
{
    if (wt.convMethod() == ConvolutionMethod::direct)
    {
        imodwtDirect(wt, oup);
        return;
    }
    imodwtFft(wt, oup);
}

auto summary(WaveletTransform const& wt) -> void
{
    auto j = wt.levels();
    summary(wt.wave());
    fmt::printf("\n");
    fmt::printf("Wavelet Transform : %s \n", wt.method().c_str());
    fmt::printf("Signal Extension : %s \n", toString(wt.extension()).c_str());
    fmt::printf("Convolutional Method : %s \n", toString(wt.convMethod()).c_str());
    fmt::printf("Number of Decomposition Levels %d \n", wt.levels());
    fmt::printf("Length of Input Signal %d \n", wt.signalLength());
    fmt::printf("Length of WT Output Vector %d \n", wt.outlength);
    fmt::printf("Wavelet Coefficients are contained in vector : %s \n", "output");
    fmt::printf("Approximation Coefficients \n");
    fmt::printf("Level %d Access : output[%d] Length : %d \n", j, 0, wt.length[0]);
    fmt::printf("Detail Coefficients \n");
    auto t = wt.length[0];
    for (auto i = 0; i < j; ++i)
    {
        fmt::printf("Level %d Access : output[%d] Length : %d \n", j - i, t, wt.length[i + 1]);
        t += wt.length[i + 1];
    }
    fmt::printf("\n");
}

}  // namespace mc::dsp