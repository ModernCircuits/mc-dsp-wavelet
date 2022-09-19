// SPDX-License-Identifier: BSL-1.0

#include "wavelet_transform.hpp"

#include <mc/dsp/algorithm/down_sample.hpp>
#include <mc/dsp/algorithm/up_sample.hpp>
#include <mc/dsp/algorithm/up_sample_even.hpp>
#include <mc/dsp/convolution/convolute.hpp>
#include <mc/dsp/convolution/FFTConvolver.hpp>
#include <mc/dsp/fft/FFT.hpp>
#include <mc/dsp/wavelet/transform/common.hpp>

#include <mc/core/cassert.hpp>
#include <mc/core/cmath.hpp>
#include <mc/core/cstdlib.hpp>
#include <mc/core/cstring.hpp>
#include <mc/core/format.hpp>
#include <mc/core/iterator.hpp>
#include <mc/core/raise.hpp>
#include <mc/core/stdexcept.hpp>
#include <mc/core/string_view.hpp>

namespace mc::dsp {

WaveletTransform::WaveletTransform(
    Wavelet& w,
    char const* method,
    std::size_t siglength,
    std::size_t j
)
    : _wave{&w}
    , _levels{j}
    , _signalLength{siglength}
    , _method{method}
    , convolver{nullptr}
    , modwtsiglength{siglength}
    , lenlength{_levels + 2}
    , MaxIter{maxIterations(siglength, w.size())}

{
    auto const size = w.size();

    if (_levels > 100) {
        raise<InvalidArgument>("decomposition Iterations Cannot Exceed 100.");
    }

    if (_levels > MaxIter) {
        raise<InvalidArgument>(
            "signal Can only be iterated maxIter times using this wavelet"
        );
    }

    if (method == nullptr) {
        this->params    = makeUnique<float[]>(siglength + 2 * _levels * (size + 1));
        this->outlength = siglength + 2 * _levels * (size + 1);
        _ext            = SignalExtension::symmetric;
    } else if ((method == StringView{"dwt"}) || (method == StringView{"DWT"})) {
        this->params    = makeUnique<float[]>(siglength + 2 * _levels * (size + 1));
        this->outlength = siglength + 2 * _levels * (size + 1);
        _ext            = SignalExtension::symmetric;
    } else if ((method == StringView{"swt"}) || (method == StringView{"SWT"})) {
        if (testSWTlength(static_cast<int>(siglength), static_cast<int>(_levels)) == 0) {
            raise<InvalidArgument>(
                "For SWT the signal length must be a multiple of 2^levels"
            );
        }

        this->params    = makeUnique<float[]>(siglength * (_levels + 1));
        this->outlength = siglength * (_levels + 1);
        _ext            = SignalExtension::periodic;
    } else if ((method == StringView{"MODWT"}) || (method == StringView{"modwt"})) {

        if (strstr(w.name().c_str(), "haar") == nullptr) {
            if (strstr(w.name().c_str(), "db") == nullptr) {
                if (strstr(w.name().c_str(), "sym") == nullptr) {
                    if (strstr(w.name().c_str(), "coif") == nullptr) {
                        raise<InvalidArgument>(
                            "MODWT is only implemented for orthogonal wavelet families - "
                            "db, sym and coif"
                        );
                    }
                }
            }
        }

        this->params    = makeUnique<float[]>(siglength * 2 * (_levels + 1));
        this->outlength = siglength * (_levels + 1);
        _ext            = SignalExtension::periodic;
    }

    this->_output = &this->params[0];
    if ((method == StringView{"dwt"}) || (method == StringView{"DWT"})) {
        for (std::size_t i = 0; i < siglength + 2 * levels() * (size + 1); ++i) {
            this->params[i] = 0.0F;
        }
    } else if ((method == StringView{"swt"}) || (method == StringView{"SWT"})) {
        for (std::size_t i = 0; i < siglength * (levels() + 1); ++i) {
            this->params[i] = 0.0F;
        }
    } else if ((method == StringView{"MODWT"}) || (method == StringView{"modwt"})) {
        for (std::size_t i = 0; i < siglength * 2 * (levels() + 1); ++i) {
            this->params[i] = 0.0F;
        }
    }
}

auto WaveletTransform::wave() const noexcept -> Wavelet const& { return *_wave; }

auto WaveletTransform::levels() const noexcept -> int { return static_cast<int>(_levels); }

auto WaveletTransform::signalLength() const noexcept -> std::size_t
{
    return _signalLength;
}

auto WaveletTransform::method() const noexcept -> String const& { return _method; }

auto WaveletTransform::extension() const noexcept -> SignalExtension { return _ext; }

auto WaveletTransform::convMethod() const noexcept -> ConvolutionMethod { return _cmethod; }

auto WaveletTransform::convMethod(ConvolutionMethod method) -> void { _cmethod = method; }

auto WaveletTransform::extension(SignalExtension ext) -> void
{
    MC_ASSERT((ext == SignalExtension::periodic) || (ext == SignalExtension::symmetric));
    _ext = ext;
}

auto WaveletTransform::output() const -> Span<float>
{
    return Span<float>{_output, static_cast<std::size_t>(outlength)};
}

auto WaveletTransform::approx() const -> Span<float>
{
    /*
    Wavelet decomposition is stored as
    [A(J) D(J) D(J-1) ..... D(1)] in wt->output vector

    Length of A(J) , N = wt->length[0]
    */

    return {_output, static_cast<size_t>(length[0])};
}

auto WaveletTransform::detail(std::size_t level) const -> Span<float>
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

    if (level > static_cast<std::size_t>(levels()) || level < 1U) {
        raisef<InvalidArgument>("The decomposition only has 1,..,{:d} levels", levels());
    }

    auto iter = length[0];
    for (auto i = 1U; i < levels() - level; ++i) { iter += length[i]; }

    return {&_output[iter], static_cast<size_t>(length[level])};
}

static auto wconv(WaveletTransform& wt, Span<float> sig, Span<float const> filt, float* oup)
    -> void
{
    if (wt.convMethod() == ConvolutionMethod::direct) {
        convolute<float>(sig, filt, oup);
        return;
    }

    MC_ASSERT(wt.convMethod() == ConvolutionMethod::fft);
    if (wt.cfftset == 0) {
        wt.convolver = makeUnique<FFTConvolver>(mc::size(sig), mc::size(filt));
        convolute(*wt.convolver, sig, filt, oup);
    } else {
        convolute(*wt.convolver, sig, filt, oup);
    }
}

static auto dwtPer(WaveletTransform& wt, float* inp, int n, float* cA, int lenCA, float* cD)
    -> void
{

    dwtPerStride(
        inp,
        n,
        wt.wave().lpd().data(),
        wt.wave().hpd().data(),
        wt.wave().lpd().size(),
        cA,
        lenCA,
        cD,
        1,
        1
    );
}

static auto dwtSym(WaveletTransform& wt, float* inp, int n, float* cA, int lenCA, float* cD)
    -> void
{

    dwtSymStride(
        inp,
        n,
        wt.wave().lpd().data(),
        wt.wave().hpd().data(),
        wt.wave().lpd().size(),
        cA,
        lenCA,
        cD,
        1,
        1
    );
}

static auto dwt1(WaveletTransform& wt, float* sig, std::size_t lenSig, float* cA, float* cD)
    -> void
{
    if (wt.extension() == SignalExtension::periodic) {
        auto lenAvg  = (wt.wave().lpd().size() + wt.wave().hpd().size()) / 2;
        auto signal  = makeUnique<float[]>(lenSig + lenAvg + (lenSig % 2));
        lenSig       = dsp::periodicExtension({sig, lenSig}, lenAvg / 2, signal.get());
        auto cAUndec = makeUnique<float[]>(lenSig + lenAvg + wt.wave().lpd().size() - 1);

        if (wt.wave().lpd().size() == wt.wave().hpd().size()
            && (wt.convMethod() == ConvolutionMethod::fft)) {
            wt.convolver
                = makeUnique<FFTConvolver>(lenSig + lenAvg, wt.wave().lpd().size());
            wt.cfftset = 1;
        } else if (!(wt.wave().lpd().size() == wt.wave().hpd().size())) {
            raise<InvalidArgument>("decomposition filters must have the same length.");
        }

        wconv(wt, {signal.get(), (lenSig + lenAvg)}, wt.wave().lpd(), cAUndec.get());
        downSample<float>(cAUndec.get() + lenAvg, static_cast<std::size_t>(lenSig), 2U, cA);

        wconv(wt, {signal.get(), (lenSig + lenAvg)}, wt.wave().hpd(), cAUndec.get());
        downSample<float>(cAUndec.get() + lenAvg, static_cast<std::size_t>(lenSig), 2U, cD);

    } else if (wt.extension() == SignalExtension::symmetric) {
        auto lf      = wt.wave().lpd().size();  // lpd and hpd have the same length
        auto signal  = makeUnique<float[]>(lenSig + 2 * (lf - 1));
        lenSig       = dsp::symmetricExtension({sig, (size_t)lenSig}, lf - 1, signal.get());
        auto cAUndec = makeUnique<float[]>(lenSig + 3 * (lf - 1));

        if (wt.wave().lpd().size() == wt.wave().hpd().size()
            && (wt.convMethod() == ConvolutionMethod::fft)) {
            wt.convolver = makeUnique<FFTConvolver>(lenSig + 2 * (lf - 1), lf);
            wt.cfftset   = 1;
        } else if (!(wt.wave().lpd().size() == wt.wave().hpd().size())) {
            raise<InvalidArgument>("decomposition filters must have the same length.");
        }

        wconv(wt, {signal.get(), lenSig + 2U * (lf - 1U)}, wt.wave().lpd(), cAUndec.get());
        downSample<float>(cAUndec.get() + lf, lenSig + lf - 2U, 2, cA);

        wconv(wt, {signal.get(), lenSig + 2 * (lf - 1)}, wt.wave().hpd(), cAUndec.get());
        downSample<float>(cAUndec.get() + lf, lenSig + lf - 2U, 2, cD);

    } else {
        raise<InvalidArgument>("Signal extension can be either per or sym");
    }

    if (wt.wave().lpd().size() == wt.wave().hpd().size()
        && (wt.convMethod() == ConvolutionMethod::fft)) {

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

    auto orig2 = makeUnique<float[]>(tempLen);
    auto orig  = makeUnique<float[]>(tempLen);

    std::copy(inp, inp + wt.signalLength(), orig.get());

    if (wt.zpad == 1) { orig[tempLen - 1] = orig[tempLen - 2]; }

    auto n  = tempLen;
    auto lp = wt.wave().lpd().size();

    if (wt.extension() == SignalExtension::periodic) {
        auto idx = j;
        while (idx > 0) {
            n              = (int)std::ceil((float)n / 2.0F);
            wt.length[idx] = n;
            wt.outlength += wt.length[idx];
            idx--;
        }
        wt.length[0] = wt.length[1];
        wt.outlength += wt.length[0];
        n = wt.outlength;

        for (auto iter = 0; iter < j; ++iter) {
            auto const lenCA = wt.length[j - iter];
            n -= lenCA;
            if (wt.convMethod() == ConvolutionMethod::fft) {
                dwt1(wt, orig.get(), tempLen, orig2.get(), wt.params.get() + n);
            } else {
                dwtPer(wt, orig.get(), tempLen, orig2.get(), lenCA, wt.params.get() + n);
            }
            tempLen = wt.length[j - iter];
            if (iter == j - 1) {
                for (std::size_t i = 0; i < lenCA; ++i) { wt.params[i] = orig2[i]; }
            } else {
                for (std::size_t i = 0; i < lenCA; ++i) { orig[i] = orig2[i]; }
            }
        }
    } else if (wt.extension() == SignalExtension::symmetric) {
        auto idx = j;
        while (idx > 0) {
            n              = n + lp - 2;
            n              = (int)std::ceil((float)n / 2.0F);
            wt.length[idx] = n;
            wt.outlength += wt.length[idx];
            idx--;
        }
        wt.length[0] = wt.length[1];
        wt.outlength += wt.length[0];
        n = wt.outlength;

        for (auto iter = 0; iter < j; ++iter) {
            auto const lenCA = wt.length[j - iter];
            n -= lenCA;
            if (wt.convMethod() == ConvolutionMethod::fft) {
                dwt1(wt, orig.get(), tempLen, orig2.get(), wt.params.get() + n);
            } else {
                dwtSym(wt, orig.get(), tempLen, orig2.get(), lenCA, wt.params.get() + n);
            }
            tempLen = wt.length[j - iter];

            if (iter == j - 1) {
                for (std::size_t i = 0; i < lenCA; ++i) { wt.params[i] = orig2[i]; }
            } else {
                for (std::size_t i = 0; i < lenCA; ++i) { orig[i] = orig2[i]; }
            }
        }
    } else {
        raise<InvalidArgument>("Signal extension can be either per or sym");
    }
}

static auto idwt1(
    WaveletTransform& wt,
    float* temp,
    float* cAUp,
    float* cA,
    int lenCA,
    float* cD,
    int lenCD,
    float* xLp,
    float* xHp,
    float* x
) -> void
{
    auto lenAvg = (wt.wave().lpr().size() + wt.wave().hpr().size()) / 2UL;
    auto n      = 2 * lenCD;
    auto u      = 2;

    upSampleEven(cA, cA + lenCA, cAUp, u);

    dsp::periodicExtension({cAUp, static_cast<std::size_t>(2 * lenCA)}, lenAvg / 2, temp);

    auto n2 = 2 * lenCA + lenAvg;

    if (wt.wave().lpr().size() == wt.wave().hpr().size()
        && (wt.convMethod() == ConvolutionMethod::fft)) {
        wt.convolver = makeUnique<FFTConvolver>(n2, lenAvg);
        wt.cfftset   = 1;
    } else if (!(wt.wave().lpr().size() == wt.wave().hpr().size())) {
        raise<InvalidArgument>("Decomposition Filters must have the same length");
    }

    wconv(wt, {temp, n2}, {wt.wave().lpr().data(), lenAvg}, xLp);

    upSampleEven(cD, cD + lenCD, cAUp, u);

    dsp::periodicExtension({cAUp, static_cast<std::size_t>(2 * lenCD)}, lenAvg / 2, temp);

    n2 = 2 * lenCD + lenAvg;

    wconv(wt, {temp, n2}, {wt.wave().hpr().data(), lenAvg}, xHp);

    for (auto i = lenAvg - 1; i < n + lenAvg - 1; ++i) {
        x[i - lenAvg + 1] = xLp[i] + xHp[i];
    }

    if (wt.wave().lpr().size() == wt.wave().hpr().size()
        && (wt.convMethod() == ConvolutionMethod::fft)) {

        wt.cfftset = 0;
    }
}

static auto idwtPer(WaveletTransform& wt, float* cA, int lenCA, float* cD, float* x) -> void
{
    idwtPerStride(
        cA,
        lenCA,
        cD,
        wt.wave().lpr().data(),
        wt.wave().hpr().data(),
        wt.wave().lpr().size(),
        x,
        1,
        1
    );
}

static auto idwtSym(WaveletTransform& wt, float* cA, int lenCA, float* cD, float* x) -> void
{
    idwtSymStride(
        cA,
        lenCA,
        cD,
        wt.wave().lpr().data(),
        wt.wave().hpr().data(),
        wt.wave().lpr().size(),
        x,
        1,
        1
    );
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
    auto out    = makeUnique<float[]>(wt.signalLength() + 1);

    if ((wt.extension() == SignalExtension::periodic)
        && (wt.convMethod() == ConvolutionMethod::fft)) {
        appLen = wt.length[0];
        detLen = wt.length[1];
        n      = 2 * wt.length[j];
        lf     = (wt.wave().lpr().size() + wt.wave().hpr().size()) / 2;

        auto cAUp = makeUnique<float[]>(n);
        auto temp = makeUnique<float[]>((n + lf));
        auto xLp  = makeUnique<float[]>((n + 2 * lf - 1));
        auto xHp  = makeUnique<float[]>((n + 2 * lf - 1));
        iter      = appLen;

        for (std::size_t i = 0; i < appLen; ++i) { out[i] = wt.output()[i]; }

        for (auto i = 0; i < j; ++i) {
            idwt1(
                wt,
                temp.get(),
                cAUp.get(),
                out.get(),
                detLen,
                wt.output().data() + iter,
                detLen,
                xLp.get(),
                xHp.get(),
                out.get()
            );
            iter += detLen;
            detLen = wt.length[i + 2];
        }
    } else if ((wt.extension() == SignalExtension::periodic) && (wt.convMethod() == ConvolutionMethod::direct)) {
        appLen = wt.length[0];
        detLen = wt.length[1];
        n      = 2 * wt.length[j];
        lf     = (wt.wave().lpr().size() + wt.wave().hpr().size()) / 2;

        auto xLp = makeUnique<float[]>((n + 2 * lf - 1));
        iter     = appLen;

        for (std::size_t i = 0; i < appLen; ++i) { out[i] = wt.output()[i]; }

        for (auto i = 0; i < j; ++i) {
            idwtPer(wt, out.get(), detLen, wt.output().data() + iter, xLp.get());
            for (k = lf / 2 - 1; k < 2 * detLen + lf / 2 - 1; ++k) {
                out[k - lf / 2 + 1] = xLp[k];
            }

            iter += detLen;
            detLen = wt.length[i + 2];
        }
    } else if ((wt.extension() == SignalExtension::symmetric) && (wt.convMethod() == ConvolutionMethod::direct)) {
        appLen = wt.length[0];
        detLen = wt.length[1];
        n      = 2 * wt.length[j] - 1;
        lf     = (wt.wave().lpr().size() + wt.wave().hpr().size()) / 2;

        auto xLp = makeUnique<float[]>((n + 2 * lf - 1));
        iter     = appLen;

        for (std::size_t i = 0; i < appLen; ++i) { out[i] = wt.output()[i]; }

        for (auto i = 0; i < j; ++i) {
            idwtSym(wt, out.get(), detLen, wt.output().data() + iter, xLp.get());
            for (k = lf - 2; k < 2 * detLen; ++k) { out[k - lf + 2] = xLp[k]; }

            iter += detLen;
            detLen = wt.length[i + 2];
        }
    } else if ((wt.extension() == SignalExtension::symmetric) && (wt.convMethod() == ConvolutionMethod::fft)) {
        lf = wt.wave().lpd().size();  // lpd and hpd have the same length

        n         = 2 * wt.length[j] - 1;
        auto cAUp = makeUnique<float[]>(n);
        auto xLp  = makeUnique<float[]>((n + lf - 1));
        auto xHp  = makeUnique<float[]>((n + lf - 1));

        for (std::size_t i = 0; i < appLen; ++i) { out[i] = wt.output()[i]; }

        iter = appLen;

        for (auto i = 0; i < j; ++i) {
            detLen = wt.length[i + 1];
            upSample<float>(out.get(), detLen, u, cAUp.get());
            n2 = 2 * wt.length[i + 1] - 1;

            if (wt.wave().lpr().size() != wt.wave().hpr().size()) {
                raise<InvalidArgument>("Decomposition Filters must have the same length");
            }

            wt.convolver = makeUnique<FFTConvolver>(n2, lf);
            wt.cfftset   = 1;

            wconv(wt, {cAUp.get(), n2}, {wt.wave().lpr().data(), lf}, xLp.get());
            upSample<float>(wt.output().data() + iter, detLen, u, cAUp.get());
            wconv(wt, {cAUp.get(), n2}, {wt.wave().hpr().data(), lf}, xHp.get());

            for (k = lf - 2; k < n2 + 1; ++k) { out[k - lf + 2] = xLp[k] + xHp[k]; }
            iter += detLen;

            MC_ASSERT(wt.convMethod() == ConvolutionMethod::fft);
            MC_ASSERT(wt.wave().lpr().size() == wt.wave().hpr().size());

            wt.cfftset = 0;
        }
    } else {
        raise<InvalidArgument>("Signal extension can be either per or sym");
    }

    std::copy(out.get(), out.get() + wt.signalLength(), dwtop);
}

static auto
swtPer(WaveletTransform& wt, int m, float* inp, int n, float* cA, int lenCA, float* cD)
    -> void
{

    swtPerStride(
        m,
        inp,
        n,
        wt.wave().lpd().data(),
        wt.wave().hpd().data(),
        wt.wave().lpd().size(),
        cA,
        lenCA,
        cD,
        1,
        1
    );
}

static auto swtFft(WaveletTransform& wt, float const* inp) -> void
{

    auto tempLen = wt.signalLength();
    auto j       = wt.levels();
    wt.length[0] = wt.length[j] = tempLen;
    wt.outlength = wt.length[j + 1] = (j + 1) * tempLen;
    auto m                          = 1;
    for (auto iter = 1; iter < j; ++iter) {
        m               = 2 * m;
        wt.length[iter] = tempLen;
    }

    auto const lenFilt = wt.wave().size();

    auto lowPass  = makeUnique<float[]>(m * lenFilt);
    auto highPass = makeUnique<float[]>(m * lenFilt);
    auto sig      = makeUnique<float[]>((m * lenFilt + tempLen + (tempLen % 2)));
    auto cA       = makeUnique<float[]>((2 * m * lenFilt + tempLen + (tempLen % 2)) - 1);
    auto cD       = makeUnique<float[]>((2 * m * lenFilt + tempLen + (tempLen % 2)) - 1);

    m = 1;

    for (std::size_t i = 0; i < tempLen; ++i) { wt.params[i] = inp[i]; }

    auto lenacc = wt.outlength;
    std::size_t n{0};

    for (auto iter = 0; iter < j; ++iter) {
        lenacc -= tempLen;
        if (iter > 0) {
            m = 2 * m;
            n = m * lenFilt;
            upSampleEven(begin(wt.wave().lpd()), end(wt.wave().lpd()), lowPass.get(), m);
            upSampleEven(begin(wt.wave().hpd()), end(wt.wave().hpd()), highPass.get(), m);
        } else {
            n = lenFilt;
            for (std::size_t i = 0; i < n; ++i) {
                lowPass[i]  = wt.wave().lpd()[i];
                highPass[i] = wt.wave().hpd()[i];
            }
        }

        // swt_per(wt,M, wt.params.get(), temp_len, cA, temp_len, cD,temp_len);

        dsp::periodicExtension({wt.params.get(), tempLen}, n / 2, sig.get());

        if (wt.wave().lpd().size() == wt.wave().hpd().size()
            && (wt.convMethod() == ConvolutionMethod::fft)) {
            wt.convolver = makeUnique<FFTConvolver>(n + tempLen + (tempLen % 2), n);
            wt.cfftset   = 1;
        } else if (!(wt.wave().lpd().size() == wt.wave().hpd().size())) {
            raise<InvalidArgument>("Decomposition Filters must have the same length");
        }

        wconv(wt, {sig.get(), n + tempLen + (tempLen % 2)}, {lowPass.get(), n}, cA.get());

        wconv(wt, {sig.get(), n + tempLen + (tempLen % 2)}, {highPass.get(), n}, cD.get());

        if (wt.wave().lpd().size() == wt.wave().hpd().size()
            && (wt.convMethod() == ConvolutionMethod::fft)) {

            wt.cfftset = 0;
        }

        for (std::size_t i = 0; i < tempLen; ++i) {
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
    for (auto iter = 1; iter < j; ++iter) {
        m               = 2 * m;
        wt.length[iter] = tempLen;
    }

    auto cA = makeUnique<float[]>(tempLen);
    auto cD = makeUnique<float[]>(tempLen);

    m = 1;

    for (std::size_t i = 0; i < tempLen; ++i) { wt.params[i] = inp[i]; }

    auto lenacc = wt.outlength;

    for (auto iter = 0; iter < j; ++iter) {
        lenacc -= tempLen;
        if (iter > 0) { m = 2 * m; }

        swtPer(wt, m, wt.params.get(), tempLen, cA.get(), tempLen, cD.get());

        for (std::size_t i = 0; i < tempLen; ++i) {
            wt.params[i]          = cA[i];
            wt.params[lenacc + i] = cD[i];
        }
    }
}

auto swt(WaveletTransform& wt, float const* inp) -> void
{
    if ((wt.method() == StringView{"swt"})
        && (wt.convMethod() == ConvolutionMethod::direct)) {
        swtDirect(wt, inp);
    } else if ((wt.method() == StringView{"swt"}) && (wt.convMethod() == ConvolutionMethod::fft)) {
        swtFft(wt, inp);
    } else {
        raise<InvalidArgument>("SWT Only accepts two methods - direct and fft");
    }
}

auto iswt(WaveletTransform& wt, float* swtop) -> void
{
    auto n  = wt.signalLength();
    auto j  = static_cast<std::size_t>(wt.levels());
    auto u  = 2;
    auto lf = wt.wave().lpr().size();

    auto appxSig = makeUnique<float[]>(n);
    auto detSig  = makeUnique<float[]>(n);
    auto appx1   = makeUnique<float[]>(n);
    auto det1    = makeUnique<float[]>(n);
    auto appx2   = makeUnique<float[]>(n);
    auto det2    = makeUnique<float[]>(n);
    auto tempx   = makeUnique<float[]>(n);
    auto cL0     = makeUnique<float[]>((n + (n % 2) + lf));
    auto cH0     = makeUnique<float[]>((n + (n % 2) + lf));
    auto oup00L  = makeUnique<float[]>((n + 2 * lf));
    auto oup00H  = makeUnique<float[]>((n + 2 * lf));
    auto oup00   = makeUnique<float[]>(n);
    auto oup01   = makeUnique<float[]>(n);

    for (std::size_t iter = 0; iter < j; ++iter) {
        for (std::size_t i = 0; i < n; ++i) { swtop[i] = 0.0F; }
        if (iter == 0) {
            for (std::size_t i = 0; i < n; ++i) {
                appxSig[i] = wt.output()[i];
                detSig[i]  = wt.output()[n + i];
            }
        } else {
            for (std::size_t i = 0; i < n; ++i) {
                detSig[i] = wt.output()[(iter + 1) * n + i];
            }
        }

        auto const value = (int)std::pow(2.0F, (float)(j - 1 - iter));

        for (auto count = 0; count < value; count++) {
            auto len = 0;
            for (std::size_t index = count; index < n; index += value) {
                appx1[len] = appxSig[index];
                det1[len]  = detSig[index];
                len++;
            }

            // SHIFT 0
            auto len0 = 0;

            for (auto indexShift = 0; indexShift < len; indexShift += 2) {
                appx2[len0] = appx1[indexShift];
                det2[len0]  = det1[indexShift];
                len0++;
            }
            upSampleEven(appx2.get(), appx2.get() + len0, tempx.get(), u);
            dsp::periodicExtension(
                {tempx.get(), static_cast<std::size_t>(2 * len0)},
                lf / 2,
                cL0.get()
            );

            upSampleEven(det2.get(), det2.get() + len0, tempx.get(), u);
            dsp::periodicExtension(
                {tempx.get(), static_cast<std::size_t>(2 * len0)},
                lf / 2,
                cH0.get()
            );

            auto n1 = 2 * len0 + lf;

            if (wt.wave().lpr().size() == wt.wave().hpr().size()
                && (wt.convMethod() == ConvolutionMethod::fft)) {
                wt.convolver = makeUnique<FFTConvolver>(n1, lf);
                wt.cfftset   = 1;
            } else if (!(wt.wave().lpd().size() == wt.wave().hpd().size())) {
                raise<InvalidArgument>("Decomposition Filters must have the same length");
            }

            wconv(wt, {cL0.get(), n1}, {wt.wave().lpr().data(), lf}, oup00L.get());

            wconv(wt, {cH0.get(), n1}, {wt.wave().hpr().data(), lf}, oup00H.get());

            for (auto i = lf - 1; i < 2 * len0 + lf - 1; ++i) {
                oup00[i - lf + 1] = oup00L[i] + oup00H[i];
            }

            // SHIFT 1

            len0 = 0;

            for (auto indexShift = 1; indexShift < len; indexShift += 2) {
                appx2[len0] = appx1[indexShift];
                det2[len0]  = det1[indexShift];
                len0++;
            }

            upSampleEven(appx2.get(), appx2.get() + len0, tempx.get(), u);
            dsp::periodicExtension(
                {tempx.get(), static_cast<std::size_t>(2 * len0)},
                lf / 2,
                cL0.get()
            );

            upSampleEven(det2.get(), det2.get() + len0, tempx.get(), u);
            dsp::periodicExtension(
                {tempx.get(), static_cast<std::size_t>(2 * len0)},
                lf / 2,
                cH0.get()
            );

            n1 = 2 * len0 + lf;

            wconv(wt, {cL0.get(), n1}, {wt.wave().lpr().data(), lf}, oup00L.get());
            wconv(wt, {cH0.get(), n1}, {wt.wave().hpr().data(), lf}, oup00H.get());

            for (auto i = lf - 1; i < 2 * len0 + lf - 1; ++i) {
                oup01[i - lf + 1] = oup00L[i] + oup00H[i];
            }

            // Rotate right by 1
            auto view = Span<float>{oup01.get(), static_cast<size_t>(2 * len0)};
            std::rotate(view.rbegin(), view.rbegin() + 1, view.rend());

            auto index2 = 0;

            for (auto index = static_cast<std::size_t>(count); index < n; index += value) {
                swtop[index] = (oup00[index2] + oup01[index2]) / 2.0F;
                index2++;
            }
        }
        for (std::size_t i = 0; i < n; ++i) { appxSig[i] = swtop[i]; }
    }
}

static auto
modwtPer(WaveletTransform& wt, int m, float const* inp, float* cA, int lenCA, float* cD)
    -> void
{
    auto const lenAvg = wt.wave().lpd().size();
    auto filt         = makeUnique<float[]>(2 * lenAvg);
    auto s            = sqrt(2.0F);

    for (std::size_t i = 0; i < lenAvg; ++i) {
        filt[i]          = wt.wave().lpd()[i] / s;
        filt[lenAvg + i] = wt.wave().hpd()[i] / s;
    }

    for (auto i = 0; i < lenCA; ++i) {
        auto t = i;
        cA[i]  = filt[0] * inp[t];
        cD[i]  = filt[lenAvg] * inp[t];
        for (std::size_t l = 1; l < lenAvg; l++) {
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
    if (wt.extension() != SignalExtension::periodic) {
        raise<InvalidArgument>("MODWT direct method only uses periodic extension per.");
    }

    auto tempLen = wt.signalLength();
    auto j       = static_cast<std::size_t>(wt.levels());
    wt.length[0] = wt.length[j] = tempLen;
    wt.outlength = wt.length[j + 1] = (j + 1) * tempLen;
    auto m                          = 1;
    for (std::size_t iter = 1; iter < j; ++iter) {
        m               = 2 * m;
        wt.length[iter] = tempLen;
    }

    auto cA = makeUnique<float[]>(tempLen);
    auto cD = makeUnique<float[]>(tempLen);

    m = 1;

    for (std::size_t i = 0; i < tempLen; ++i) { wt.params[i] = inp[i]; }

    auto lenacc = wt.outlength;

    for (std::size_t iter = 0; iter < j; ++iter) {
        lenacc -= tempLen;
        if (iter > 0) { m = 2 * m; }

        modwtPer(wt, m, wt.params.get(), cA.get(), tempLen, cD.get());

        for (std::size_t i = 0; i < tempLen; ++i) {
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

    auto tempLen = wt.signalLength();
    auto lenAvg  = static_cast<std::size_t>(wt.wave().lpd().size());
    std::size_t n{0};
    if (wt.extension() == SignalExtension::symmetric) {
        n = 2 * tempLen;
    } else if (wt.extension() == SignalExtension::periodic) {
        n = tempLen;
    }
    j                 = wt.levels();
    wt.modwtsiglength = n;
    wt.length[0] = wt.length[j] = n;
    wt.outlength = wt.length[j + 1] = (j + 1) * n;

    auto const s = sqrt(2.0F);
    for (iter = 1; iter < j; ++iter) { wt.length[iter] = n; }

    auto fftFd = makeUnique<FFT<float, KissFFT>>(n);
    auto fftBd = makeUnique<FFT<float, KissFFT>>(n);

    auto sig      = makeUnique<Complex<float>[]>(n);
    auto cA       = makeUnique<Complex<float>[]>(n);
    auto cD       = makeUnique<Complex<float>[]>(n);
    auto lowPass  = makeUnique<Complex<float>[]>(n);
    auto highPass = makeUnique<Complex<float>[]>(n);
    auto index    = makeUnique<int[]>(n);

    // N-point FFT of low pass and high pass filters

    // Low Pass Filter

    for (std::size_t i = 0; i < lenAvg; ++i) {
        sig[i].real((float)wt.wave().lpd()[i] / s);
        sig[i].imag(0.0F);
    }
    for (auto i = lenAvg; i < n; ++i) {
        sig[i].real(0.0F);
        sig[i].imag(0.0F);
    }

    fft(*fftFd, sig.get(), lowPass.get());

    // High Pass Filter

    for (std::size_t i = 0; i < lenAvg; ++i) {
        sig[i].real((float)wt.wave().hpd()[i] / s);
        sig[i].imag(0.0F);
    }
    for (std::size_t i = lenAvg; i < n; ++i) {
        sig[i].real(0.0F);
        sig[i].imag(0.0F);
    }

    fft(*fftFd, sig.get(), highPass.get());

    // symmetric extension
    for (std::size_t i = 0; i < tempLen; ++i) {
        sig[i].real((float)inp[i]);
        sig[i].imag(0.0F);
    }
    for (std::size_t i = tempLen; i < n; ++i) {
        sig[i].real((float)inp[n - i - 1]);
        sig[i].imag(0.0F);
    }

    // FFT of data

    fft(*fftFd, sig.get(), cA.get());

    lenacc = wt.outlength;

    m = 1;

    for (iter = 0; iter < j; ++iter) {
        lenacc -= n;

        for (std::size_t i = 0; i < n; ++i) { index[i] = (m * i) % n; }

        for (std::size_t i = 0; i < n; ++i) {
            auto const tmp1 = cA[i].real();
            auto const tmp2 = cA[i].imag();
            cA[i].real(lowPass[index[i]].real() * tmp1 - lowPass[index[i]].imag() * tmp2);
            cA[i].imag(lowPass[index[i]].real() * tmp2 + lowPass[index[i]].imag() * tmp1);

            cD[i].real(highPass[index[i]].real() * tmp1 - highPass[index[i]].imag() * tmp2);
            cD[i].imag(highPass[index[i]].real() * tmp2 + highPass[index[i]].imag() * tmp1);
        }

        ifft(*fftBd, cD.get(), sig.get());

        for (std::size_t i = 0; i < n; ++i) {
            wt.params[lenacc + i] = sig[i].real() / static_cast<float>(n);
        }

        m *= 2;
    }

    ifft(*fftBd, cA.get(), sig.get());

    for (std::size_t i = 0; i < n; ++i) {
        wt.params[i] = sig[i].real() / static_cast<float>(n);
    }
}

auto modwt(WaveletTransform& wt, float const* inp) -> void
{
    if (wt.convMethod() == ConvolutionMethod::direct) {
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

    auto s     = sqrt(2.0F);
    auto fftFd = makeUnique<FFT<float, KissFFT>>(n);
    auto fftBd = makeUnique<FFT<float, KissFFT>>(n);

    auto sig      = makeUnique<Complex<float>[]>(n);
    auto cA       = makeUnique<Complex<float>[]>(n);
    auto cD       = makeUnique<Complex<float>[]>(n);
    auto lowPass  = makeUnique<Complex<float>[]>(n);
    auto highPass = makeUnique<Complex<float>[]>(n);
    auto index    = makeUnique<std::size_t[]>(n);

    // N-point FFT of low pass and high pass filters

    // Low Pass Filter

    for (std::size_t i = 0; i < lenAvg; ++i) {
        sig[i].real((float)wt.wave().lpd()[i] / s);
        sig[i].imag(0.0F);
    }
    for (std::size_t i = lenAvg; i < n; ++i) {
        sig[i].real(0.0F);
        sig[i].imag(0.0F);
    }

    fft(*fftFd, sig.get(), lowPass.get());

    // High Pass Filter

    for (std::size_t i = 0; i < lenAvg; ++i) {
        sig[i].real((float)wt.wave().hpd()[i] / s);
        sig[i].imag(0.0F);
    }
    for (std::size_t i = lenAvg; i < n; ++i) {
        sig[i].real(0.0F);
        sig[i].imag(0.0F);
    }

    fft(*fftFd, sig.get(), highPass.get());

    // Complex conjugate of the two filters

    conjComplex(lowPass.get(), static_cast<int>(n));
    conjComplex(highPass.get(), static_cast<int>(n));

    auto m      = (int)std::pow(2.0F, (float)j - 1.0F);
    auto lenacc = n;

    //
    for (std::size_t i = 0; i < n; ++i) {
        sig[i].real((float)wt.output()[i]);
        sig[i].imag(0.0F);
    }

    for (std::size_t iter = 0; iter < j; ++iter) {
        fft(*fftFd, sig.get(), cA.get());
        for (std::size_t i = 0; i < n; ++i) {
            sig[i].real(wt.output()[lenacc + i]);
            sig[i].imag(0.0F);
        }
        fft(*fftFd, sig.get(), cD.get());

        for (std::size_t i = 0; i < n; ++i) { index[i] = (m * i) % n; }

        for (std::size_t i = 0; i < n; ++i) {
            auto const tmp1 = cA[i].real();
            auto const tmp2 = cA[i].imag();
            cA[i].real(
                lowPass[index[i]].real() * tmp1 - lowPass[index[i]].imag() * tmp2
                + highPass[index[i]].real() * cD[i].real()
                - highPass[index[i]].imag() * cD[i].imag()
            );
            cA[i].imag(
                lowPass[index[i]].real() * tmp2 + lowPass[index[i]].imag() * tmp1
                + highPass[index[i]].real() * cD[i].imag()
                + highPass[index[i]].imag() * cD[i].real()
            );
        }

        ifft(*fftBd, cA.get(), sig.get());

        for (std::size_t i = 0; i < n; ++i) {
            sig[i].real(sig[i].real() / static_cast<float>(n));
            sig[i].imag(sig[i].imag() / static_cast<float>(n));
        }
        m /= 2;
        lenacc += n;
    }

    std::transform(sig.get(), sig.get() + wt.signalLength(), oup, [](auto c) {
        return c.real();
    });
}

static auto imodwtPer(
    WaveletTransform& wt,
    int m,
    float const* cA,
    int lenCA,
    float const* cD,
    float* x
) -> void
{
    auto const lenAvg = wt.wave().lpd().size();
    auto filt         = makeUnique<float[]>(2 * lenAvg);
    auto s            = sqrt(2.0F);

    for (std::size_t i = 0; i < lenAvg; ++i) {
        filt[i]          = wt.wave().lpd()[i] / s;
        filt[lenAvg + i] = wt.wave().hpd()[i] / s;
    }

    for (auto i = 0; i < lenCA; ++i) {
        auto t = i;
        x[i]   = (filt[0] * cA[t]) + (filt[lenAvg] * cD[t]);
        for (std::size_t l = 1; l < lenAvg; l++) {
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

    auto x = makeUnique<float[]>(n);

    for (std::size_t i = 0; i < n; ++i) { dwtop[i] = wt.output()[i]; }

    auto m = static_cast<int>(std::pow(2.0F, (float)j - 1.0F));
    for (std::size_t iter = 0; iter < j; ++iter) {
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
    if (wt.convMethod() == ConvolutionMethod::direct) {
        imodwtDirect(wt, oup);
        return;
    }
    imodwtFft(wt, oup);
}

}  // namespace mc::dsp
