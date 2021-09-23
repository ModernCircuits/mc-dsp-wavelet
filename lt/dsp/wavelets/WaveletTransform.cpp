#include "WaveletTransform.hpp"

#include "lt/dsp/convolution/FFTConvolver.hpp"
#include "lt/dsp/convolution/convolute.hpp"
#include "lt/dsp/fft/FFT.hpp"
#include "lt/dsp/wavelets/common.hpp"

#include "lt/cassert.hpp"
#include "lt/cmath.hpp"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string_view>

using namespace std::string_view_literals;

namespace {

auto upsamp(double const* x, int lenx, int m, double* y) -> int
{
    int n = 0;
    int i = 0;
    int j = 0;
    int k = 0;

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
    int n = 0;
    int i = 0;
    int j = 0;
    int k = 0;
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

auto isign(int n) -> int
{
    int m = 0;
    if (n >= 0) {
        m = 1;
    } else {
        m = -1;
    }

    return m;
}

auto circshift(double* array, int n, int l) -> void
{
    if (std::abs(l) > n) {
        l = isign(l) * (std::abs(l) % n);
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

auto perExt(double const* sig, int len, int a, double* oup) -> int
{
    int i = 0;
    int len2 = 0;
    double temp1 = NAN;
    double temp2 = NAN;
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

auto symmExt(double const* sig, int len, int a, double* oup) -> int
{
    int i = 0;
    int len2 = 0;
    double temp1 = NAN;
    double temp2 = NAN;
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

auto downsamp(double const* x, int lenx, int m, double* y) -> int
{
    int n = 0;
    int i = 0;

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

}

WaveletTransform::WaveletTransform(Wavelet& w, char const* method, int siglength, int j)
    : wave_ { &w }
    , levels_ { j }

{
    auto const size = w.size();
    auto const maxIter = wmaxiter(siglength, w.size());

    if (levels_ > 100) {
        printf("\n The Decomposition Iterations Cannot Exceed 100. Exiting \n");
        exit(-1);
    }

    if (levels_ > maxIter) {
        printf("\n Error - The Signal Can only be iterated %d times using this wavelet. Exiting\n", maxIter);
        exit(-1);
    }

    if (method == nullptr) {
        this->params = std::make_unique<double[]>(siglength + 2 * levels_ * (size + 1));
        this->outlength = siglength + 2 * levels_ * (size + 1);
        ext_ = SignalExtension::symmetric;
    } else if ((method == "dwt"sv) || (method == "DWT"sv)) {
        this->params = std::make_unique<double[]>(siglength + 2 * levels_ * (size + 1));
        this->outlength = siglength + 2 * levels_ * (size + 1);
        ext_ = SignalExtension::symmetric;
    } else if ((method == "swt"sv) || (method == "SWT"sv)) {
        if (testSWTlength(siglength, levels_) == 0) {
            printf("\n For SWT the signal length must be a multiple of 2^levels_. \n");
            exit(-1);
        }

        this->params = std::make_unique<double[]>(siglength * (levels_ + 1));
        this->outlength = siglength * (levels_ + 1);
        ext_ = SignalExtension::periodic;
    } else if ((method == "modwt"sv) || (method == "MODWT"sv)) {

        if (strstr(w.name().c_str(), "haar") == nullptr) {
            if (strstr(w.name().c_str(), "db") == nullptr) {
                if (strstr(w.name().c_str(), "sym") == nullptr) {
                    if (strstr(w.name().c_str(), "coif") == nullptr) {
                        printf("\n MODWT is only implemented for orthogonal wavelet families - db, sym and coif \n");
                        exit(-1);
                    }
                }
            }
        }

        this->params = std::make_unique<double[]>(siglength * 2 * (levels_ + 1));
        this->outlength = siglength * (levels_ + 1);
        ext_ = SignalExtension::periodic;
    }

    this->signalLength_ = siglength;
    this->modwtsiglength = siglength;
    this->MaxIter = maxIter;
    method_ = method;

    if (siglength % 2 == 0) {
        this->even = 1;
    } else {
        this->even = 0;
    }

    this->convolver = nullptr;

    this->cmethod_ = ConvolutionMethod::direct;
    this->cfftset = 0;
    this->lenlength = levels_ + 2;
    this->output_ = &this->params[0];
    if ((method == "dwt"sv) || (method == "DWT"sv)) {
        for (auto i = 0; i < siglength + 2 * levels() * (size + 1); ++i) {
            this->params[i] = 0.0;
        }
    } else if ((method == "swt"sv) || (method == "SWT"sv)) {
        for (auto i = 0; i < siglength * (levels() + 1); ++i) {
            this->params[i] = 0.0;
        }
    } else if ((method == "modwt"sv) || (method == "MODWT"sv)) {
        for (auto i = 0; i < siglength * 2 * (levels() + 1); ++i) {
            this->params[i] = 0.0;
        }
    }
}

auto WaveletTransform::convMethod(ConvolutionMethod method) -> void
{
    cmethod_ = method;
}

auto WaveletTransform::extension(SignalExtension ext) -> void
{
    LT_ASSERT((ext == SignalExtension::periodic) || (ext == SignalExtension::symmetric));
    ext_ = ext;
}

auto WaveletTransform::output() const noexcept -> lt::span<double>
{
    return lt::span<double> { output_, static_cast<std::size_t>(outlength) };
}

auto WaveletTransform::approx() const noexcept -> lt::span<double>
{
    /*
	Wavelet decomposition is stored as
	[A(J) D(J) D(J-1) ..... D(1)] in wt->output vector

	Length of A(J) , N = wt->length[0]
	*/

    return lt::span<double>(output_, static_cast<size_t>(length[0]));
}

auto WaveletTransform::detail(std::size_t level) const noexcept -> lt::span<double>
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
        printf("The decomposition only has 1,..,%d levels", levels());
        exit(-1);
    }

    auto iter = length[0];
    for (auto i = 1U; i < levels() - level; ++i) {
        iter += length[i];
    }

    return lt::span<double>(&output_[iter], static_cast<size_t>(length[level]));
}

static auto wconv(WaveletTransform& wt, double* sig, int n, double const* filt, int l, double* oup) -> void
{
    if (wt.convMethod() == ConvolutionMethod::direct) {
        convolute(sig, n, filt, l, oup);
        return;
    }

    LT_ASSERT(wt.convMethod() == ConvolutionMethod::fft);
    if (wt.cfftset == 0) {
        wt.convolver = std::make_unique<FFTConvolver>(n, l);
        convolute(*wt.convolver, sig, filt, oup);
    } else {
        convolute(*wt.convolver, sig, filt, oup);
    }
}

static auto dwtPer(WaveletTransform& wt, double* inp, int n, double* cA, int lenCA, double* cD) -> void
{

    dwtPerStride(inp, n, wt.wave().lpd(), wt.wave().hpd(), wt.wave().lpdLen(), cA, lenCA, cD, 1, 1);
}

static auto dwtSym(WaveletTransform& wt, double* inp, int n, double* cA, int lenCA, double* cD) -> void
{

    dwtSymStride(inp, n, wt.wave().lpd(), wt.wave().hpd(), wt.wave().lpdLen(), cA, lenCA, cD, 1, 1);
}

static auto dwt1(WaveletTransform& wt, double* sig, int lenSig, double* cA, double* cD) -> void
{
    constexpr auto d = 2;

    if (wt.extension() == SignalExtension::periodic) {
        auto lenAvg = (wt.wave().lpdLen() + wt.wave().hpdLen()) / 2;
        auto signal = std::make_unique<double[]>(lenSig + lenAvg + (lenSig % 2));
        lenSig = perExt(sig, lenSig, lenAvg / 2, signal.get());
        auto cAUndec = std::make_unique<double[]>(lenSig + lenAvg + wt.wave().lpdLen() - 1);

        if (wt.wave().lpdLen() == wt.wave().hpdLen() && (wt.convMethod() == ConvolutionMethod::fft)) {
            wt.convolver = std::make_unique<FFTConvolver>(lenSig + lenAvg, wt.wave().lpdLen());
            wt.cfftset = 1;
        } else if (!(wt.wave().lpdLen() == wt.wave().hpdLen())) {
            printf("Decomposition Filters must have the same length.");
            exit(-1);
        }

        wconv(wt, signal.get(), lenSig + lenAvg, wt.wave().lpd(), wt.wave().lpdLen(), cAUndec.get());
        downsamp(cAUndec.get() + lenAvg, lenSig, d, cA);
        wconv(wt, signal.get(), lenSig + lenAvg, wt.wave().hpd(), wt.wave().hpdLen(), cAUndec.get());
        downsamp(cAUndec.get() + lenAvg, lenSig, d, cD);

    } else if (wt.extension() == SignalExtension::symmetric) {
        auto lf = wt.wave().lpdLen(); // lpd and hpd have the same length
        auto signal = std::make_unique<double[]>(lenSig + 2 * (lf - 1));
        lenSig = symmExt(sig, lenSig, lf - 1, signal.get());
        auto cAUndec = std::make_unique<double[]>(lenSig + 3 * (lf - 1));

        if (wt.wave().lpdLen() == wt.wave().hpdLen() && (wt.convMethod() == ConvolutionMethod::fft)) {
            wt.convolver = std::make_unique<FFTConvolver>(lenSig + 2 * (lf - 1), lf);
            wt.cfftset = 1;
        } else if (!(wt.wave().lpdLen() == wt.wave().hpdLen())) {
            printf("Decomposition Filters must have the same length.");
            exit(-1);
        }

        wconv(wt, signal.get(), lenSig + 2 * (lf - 1), wt.wave().lpd(), wt.wave().lpdLen(), cAUndec.get());
        downsamp(cAUndec.get() + lf, lenSig + lf - 2, d, cA);
        wconv(wt, signal.get(), lenSig + 2 * (lf - 1), wt.wave().hpd(), wt.wave().hpdLen(), cAUndec.get());
        downsamp(cAUndec.get() + lf, lenSig + lf - 2, d, cD);
    } else {
        printf("Signal extension can be either per or sym");
        exit(-1);
    }

    if (wt.wave().lpdLen() == wt.wave().hpdLen() && (wt.convMethod() == ConvolutionMethod::fft)) {

        wt.cfftset = 0;
    }
}

auto dwt(WaveletTransform& wt, double const* inp) -> void
{

    auto tempLen = wt.signalLength();
    auto const j = wt.levels();

    wt.length[j + 1] = tempLen;
    wt.outlength = 0;
    wt.zpad = 0;

    auto orig2 = std::make_unique<double[]>(tempLen);
    auto orig = std::make_unique<double[]>(tempLen);

    for (auto i = 0; i < wt.signalLength(); ++i) {
        orig[i] = inp[i];
    }

    if (wt.zpad == 1) {
        orig[tempLen - 1] = orig[tempLen - 2];
    }

    auto n = tempLen;
    auto lp = wt.wave().lpdLen();

    if (wt.extension() == SignalExtension::periodic) {
        auto idx = j;
        while (idx > 0) {
            n = (int)ceil((double)n / 2.0);
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
                for (auto i = 0; i < lenCA; ++i) {
                    wt.params[i] = orig2[i];
                }
            } else {
                for (auto i = 0; i < lenCA; ++i) {
                    orig[i] = orig2[i];
                }
            }
        }
    } else if (wt.extension() == SignalExtension::symmetric) {
        auto idx = j;
        while (idx > 0) {
            n = n + lp - 2;
            n = (int)ceil((double)n / 2.0);
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
                for (auto i = 0; i < lenCA; ++i) {
                    wt.params[i] = orig2[i];
                }
            } else {
                for (auto i = 0; i < lenCA; ++i) {
                    orig[i] = orig2[i];
                }
            }
        }
    } else {
        printf("Signal extension can be either per or sym");
        exit(-1);
    }
}

static auto idwt1(WaveletTransform& wt, double* temp, double* cAUp, double* cA, int lenCA, double* cD, int lenCD, double* xLp, double* xHp, double* x) -> void
{
    auto lenAvg = (wt.wave().lprLen() + wt.wave().hprLen()) / 2;
    auto n = 2 * lenCD;
    auto u = 2;

    upsamp2(cA, lenCA, u, cAUp);

    perExt(cAUp, 2 * lenCA, lenAvg / 2, temp);

    auto n2 = 2 * lenCA + lenAvg;

    if (wt.wave().lprLen() == wt.wave().hprLen() && (wt.convMethod() == ConvolutionMethod::fft)) {
        wt.convolver = std::make_unique<FFTConvolver>(n2, lenAvg);
        wt.cfftset = 1;
    } else if (!(wt.wave().lprLen() == wt.wave().hprLen())) {
        printf("Decomposition Filters must have the same length.");
        exit(-1);
    }

    wconv(wt, temp, n2, wt.wave().lpr(), lenAvg, xLp);

    upsamp2(cD, lenCD, u, cAUp);

    perExt(cAUp, 2 * lenCD, lenAvg / 2, temp);

    n2 = 2 * lenCD + lenAvg;

    wconv(wt, temp, n2, wt.wave().hpr(), lenAvg, xHp);

    for (auto i = lenAvg - 1; i < n + lenAvg - 1; ++i) {
        x[i - lenAvg + 1] = xLp[i] + xHp[i];
    }

    if (wt.wave().lprLen() == wt.wave().hprLen() && (wt.convMethod() == ConvolutionMethod::fft)) {

        wt.cfftset = 0;
    }
}

static auto idwtPer(WaveletTransform& wt, double* cA, int lenCA, double* cD, double* x) -> void
{
    idwtPerStride(cA, lenCA, cD, wt.wave().lpr(), wt.wave().hpr(), wt.wave().lprLen(), x, 1, 1);
}

static auto idwtSym(WaveletTransform& wt, double* cA, int lenCA, double* cD, double* x) -> void
{
    idwtSymStride(cA, lenCA, cD, wt.wave().lpr(), wt.wave().hpr(), wt.wave().lprLen(), x, 1, 1);
}

auto idwt(WaveletTransform& wt, double* dwtop) -> void
{

    int lf = 0;
    int n = 0;
    int n2 = 0;
    int iter = 0;
    int k = 0;
    int detLen = 0;

    auto j = wt.levels();
    auto u = 2;
    auto appLen = wt.length[0];
    auto out = std::make_unique<double[]>(wt.signalLength() + 1);
    if ((wt.extension() == SignalExtension::periodic) && (wt.convMethod() == ConvolutionMethod::fft)) {
        appLen = wt.length[0];
        detLen = wt.length[1];
        n = 2 * wt.length[j];
        lf = (wt.wave().lprLen() + wt.wave().hprLen()) / 2;

        auto cAUp = std::make_unique<double[]>(n);
        auto temp = std::make_unique<double[]>((n + lf));
        auto xLp = std::make_unique<double[]>((n + 2 * lf - 1));
        auto xHp = std::make_unique<double[]>((n + 2 * lf - 1));
        iter = appLen;

        for (auto i = 0; i < appLen; ++i) {
            out[i] = wt.output()[i];
        }

        for (auto i = 0; i < j; ++i) {

            idwt1(wt, temp.get(), cAUp.get(), out.get(), detLen, wt.output().data() + iter, detLen, xLp.get(), xHp.get(), out.get());
            /*
			idwt_per(wt,out.get(), det_len, wt.output().data() + iter, det_len, X_lp);
			for (k = lf/2 - 1; k < 2 * det_len + lf/2 - 1; ++k) {
				out[k - lf/2 + 1] = X_lp[k];
			}
			*/
            iter += detLen;
            detLen = wt.length[i + 2];
        }

    } else if ((wt.extension() == SignalExtension::periodic) && (wt.convMethod() == ConvolutionMethod::direct)) {
        appLen = wt.length[0];
        detLen = wt.length[1];
        n = 2 * wt.length[j];
        lf = (wt.wave().lprLen() + wt.wave().hprLen()) / 2;

        auto xLp = std::make_unique<double[]>((n + 2 * lf - 1));
        iter = appLen;

        for (auto i = 0; i < appLen; ++i) {
            out[i] = wt.output()[i];
        }

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
        n = 2 * wt.length[j] - 1;
        lf = (wt.wave().lprLen() + wt.wave().hprLen()) / 2;

        auto xLp = std::make_unique<double[]>((n + 2 * lf - 1));
        iter = appLen;

        for (auto i = 0; i < appLen; ++i) {
            out[i] = wt.output()[i];
        }

        for (auto i = 0; i < j; ++i) {
            idwtSym(wt, out.get(), detLen, wt.output().data() + iter, xLp.get());
            for (k = lf - 2; k < 2 * detLen; ++k) {
                out[k - lf + 2] = xLp[k];
            }

            iter += detLen;
            detLen = wt.length[i + 2];
        }

    } else if ((wt.extension() == SignalExtension::symmetric) && (wt.convMethod() == ConvolutionMethod::fft)) {
        lf = wt.wave().lpdLen(); // lpd and hpd have the same length

        n = 2 * wt.length[j] - 1;
        auto cAUp = std::make_unique<double[]>(n);
        auto xLp = std::make_unique<double[]>((n + lf - 1));
        auto xHp = std::make_unique<double[]>((n + lf - 1));

        for (auto i = 0; i < appLen; ++i) {
            out[i] = wt.output()[i];
        }

        iter = appLen;

        for (auto i = 0; i < j; ++i) {
            detLen = wt.length[i + 1];
            upsamp(out.get(), detLen, u, cAUp.get());
            n2 = 2 * wt.length[i + 1] - 1;

            if (wt.wave().lprLen() == wt.wave().hprLen() && (wt.convMethod() == ConvolutionMethod::fft)) {
                wt.convolver = std::make_unique<FFTConvolver>(n2, lf);
                wt.cfftset = 1;
            } else if (!(wt.wave().lprLen() == wt.wave().hprLen())) {
                printf("Decomposition Filters must have the same length.");
                exit(-1);
            }

            wconv(wt, cAUp.get(), n2, wt.wave().lpr(), lf, xLp.get());
            upsamp(wt.output().data() + iter, detLen, u, cAUp.get());
            wconv(wt, cAUp.get(), n2, wt.wave().hpr(), lf, xHp.get());

            for (k = lf - 2; k < n2 + 1; ++k) {
                out[k - lf + 2] = xLp[k] + xHp[k];
            }
            iter += detLen;
            if (wt.wave().lprLen() == wt.wave().hprLen() && (wt.convMethod() == ConvolutionMethod::fft)) {

                wt.cfftset = 0;
            }
        }

    } else {
        printf("Signal extension can be either per or sym");
        exit(-1);
    }

    for (auto i = 0; i < wt.signalLength(); ++i) {
        dwtop[i] = out[i];
    }
}

static auto swtPer(WaveletTransform& wt, int m, double* inp, int n, double* cA, int lenCA, double* cD) -> void
{

    swtPerStride(m, inp, n, wt.wave().lpd(), wt.wave().hpd(), wt.wave().lpdLen(), cA, lenCA, cD, 1, 1);
}

static auto swtFft(WaveletTransform& wt, double const* inp) -> void
{
    int n { 0 };

    auto tempLen = wt.signalLength();
    auto j = wt.levels();
    wt.length[0] = wt.length[j] = tempLen;
    wt.outlength = wt.length[j + 1] = (j + 1) * tempLen;
    auto m = 1;
    for (auto iter = 1; iter < j; ++iter) {
        m = 2 * m;
        wt.length[iter] = tempLen;
    }

    auto const lenFilt = wt.wave().size();

    auto lowPass = std::make_unique<double[]>(m * lenFilt);
    auto highPass = std::make_unique<double[]>(m * lenFilt);
    auto sig = std::make_unique<double[]>((m * lenFilt + tempLen + (tempLen % 2)));
    auto cA = std::make_unique<double[]>((2 * m * lenFilt + tempLen + (tempLen % 2)) - 1);
    auto cD = std::make_unique<double[]>((2 * m * lenFilt + tempLen + (tempLen % 2)) - 1);

    m = 1;

    for (auto i = 0; i < tempLen; ++i) {
        wt.params[i] = inp[i];
    }

    auto lenacc = wt.outlength;

    for (auto iter = 0; iter < j; ++iter) {
        lenacc -= tempLen;
        if (iter > 0) {
            m = 2 * m;
            n = m * lenFilt;
            upsamp2(wt.wave().lpd(), wt.wave().lpdLen(), m, lowPass.get());
            upsamp2(wt.wave().hpd(), wt.wave().hpdLen(), m, highPass.get());
        } else {
            n = lenFilt;
            for (auto i = 0; i < n; ++i) {
                lowPass[i] = wt.wave().lpd()[i];
                highPass[i] = wt.wave().hpd()[i];
            }
        }

        //swt_per(wt,M, wt.params.get(), temp_len, cA, temp_len, cD,temp_len);

        perExt(wt.params.get(), tempLen, n / 2, sig.get());

        if (wt.wave().lpdLen() == wt.wave().hpdLen() && (wt.convMethod() == ConvolutionMethod::fft)) {
            wt.convolver = std::make_unique<FFTConvolver>(n + tempLen + (tempLen % 2), n);
            wt.cfftset = 1;
        } else if (!(wt.wave().lpdLen() == wt.wave().hpdLen())) {
            printf("Decomposition Filters must have the same length.");
            exit(-1);
        }

        wconv(wt, sig.get(), n + tempLen + (tempLen % 2), lowPass.get(), n, cA.get());

        wconv(wt, sig.get(), n + tempLen + (tempLen % 2), highPass.get(), n, cD.get());

        if (wt.wave().lpdLen() == wt.wave().hpdLen() && (wt.convMethod() == ConvolutionMethod::fft)) {

            wt.cfftset = 0;
        }

        for (auto i = 0; i < tempLen; ++i) {
            wt.params[i] = cA[n + i];
            wt.params[lenacc + i] = cD[n + i];
        }
    }
}

static auto swtDirect(WaveletTransform& wt, double const* inp) -> void
{
    int j = 0;
    int tempLen = 0;
    int iter = 0;
    int m = 0;
    int lenacc = 0;

    tempLen = wt.signalLength();
    j = wt.levels();
    wt.length[0] = wt.length[j] = tempLen;
    wt.outlength = wt.length[j + 1] = (j + 1) * tempLen;
    m = 1;
    for (iter = 1; iter < j; ++iter) {
        m = 2 * m;
        wt.length[iter] = tempLen;
    }

    auto cA = std::make_unique<double[]>(tempLen);
    auto cD = std::make_unique<double[]>(tempLen);

    m = 1;

    for (auto i = 0; i < tempLen; ++i) {
        wt.params[i] = inp[i];
    }

    lenacc = wt.outlength;

    for (iter = 0; iter < j; ++iter) {
        lenacc -= tempLen;
        if (iter > 0) {
            m = 2 * m;
        }

        swtPer(wt, m, wt.params.get(), tempLen, cA.get(), tempLen, cD.get());

        for (auto i = 0; i < tempLen; ++i) {
            wt.params[i] = cA[i];
            wt.params[lenacc + i] = cD[i];
        }
    }
}

auto swt(WaveletTransform& wt, double const* inp) -> void
{
    if ((wt.method() == "swt"sv) && (wt.convMethod() == ConvolutionMethod::direct)) {
        swtDirect(wt, inp);
    } else if ((wt.method() == "swt"sv) && (wt.convMethod() == ConvolutionMethod::fft)) {
        swtFft(wt, inp);
    } else {
        printf("SWT Only accepts two methods - direct and fft");
        exit(-1);
    }
}

auto iswt(WaveletTransform& wt, double* swtop) -> void
{
    auto n = wt.signalLength();
    auto j = wt.levels();
    auto u = 2;
    auto lf = wt.wave().lprLen();

    auto appxSig = std::make_unique<double[]>(n);
    auto detSig = std::make_unique<double[]>(n);
    auto appx1 = std::make_unique<double[]>(n);
    auto det1 = std::make_unique<double[]>(n);
    auto appx2 = std::make_unique<double[]>(n);
    auto det2 = std::make_unique<double[]>(n);
    auto tempx = std::make_unique<double[]>(n);
    auto cL0 = std::make_unique<double[]>((n + (n % 2) + lf));
    auto cH0 = std::make_unique<double[]>((n + (n % 2) + lf));
    auto oup00L = std::make_unique<double[]>((n + 2 * lf));
    auto oup00H = std::make_unique<double[]>((n + 2 * lf));
    auto oup00 = std::make_unique<double[]>(n);
    auto oup01 = std::make_unique<double[]>(n);

    for (auto iter = 0; iter < j; ++iter) {
        for (auto i = 0; i < n; ++i) {
            swtop[i] = 0.0;
        }
        if (iter == 0) {
            for (auto i = 0; i < n; ++i) {
                appxSig[i] = wt.output()[i];
                detSig[i] = wt.output()[n + i];
            }
        } else {
            for (auto i = 0; i < n; ++i) {
                detSig[i] = wt.output()[(iter + 1) * n + i];
            }
        }

        auto const value = (int)std::pow(2.0, (double)(j - 1 - iter));

        for (auto count = 0; count < value; count++) {
            auto len = 0;
            for (auto index = count; index < n; index += value) {
                appx1[len] = appxSig[index];
                det1[len] = detSig[index];
                len++;
            }

            //SHIFT 0
            auto len0 = 0;

            for (auto indexShift = 0; indexShift < len; indexShift += 2) {
                appx2[len0] = appx1[indexShift];
                det2[len0] = det1[indexShift];
                len0++;
            }
            upsamp2(appx2.get(), len0, u, tempx.get());
            perExt(tempx.get(), 2 * len0, lf / 2, cL0.get());

            upsamp2(det2.get(), len0, u, tempx.get());
            perExt(tempx.get(), 2 * len0, lf / 2, cH0.get());

            auto n1 = 2 * len0 + lf;

            if (wt.wave().lprLen() == wt.wave().hprLen() && (wt.convMethod() == ConvolutionMethod::fft)) {
                wt.convolver = std::make_unique<FFTConvolver>(n1, lf);
                wt.cfftset = 1;
            } else if (!(wt.wave().lpdLen() == wt.wave().hpdLen())) {
                printf("Decomposition Filters must have the same length.");
                exit(-1);
            }

            wconv(wt, cL0.get(), n1, wt.wave().lpr(), lf, oup00L.get());

            wconv(wt, cH0.get(), n1, wt.wave().hpr(), lf, oup00H.get());

            for (auto i = lf - 1; i < 2 * len0 + lf - 1; ++i) {
                oup00[i - lf + 1] = oup00L[i] + oup00H[i];
            }

            //SHIFT 1

            len0 = 0;

            for (auto indexShift = 1; indexShift < len; indexShift += 2) {
                appx2[len0] = appx1[indexShift];
                det2[len0] = det1[indexShift];
                len0++;
            }

            upsamp2(appx2.get(), len0, u, tempx.get());
            perExt(tempx.get(), 2 * len0, lf / 2, cL0.get());

            upsamp2(det2.get(), len0, u, tempx.get());
            perExt(tempx.get(), 2 * len0, lf / 2, cH0.get());

            n1 = 2 * len0 + lf;

            wconv(wt, cL0.get(), n1, wt.wave().lpr(), lf, oup00L.get());
            wconv(wt, cH0.get(), n1, wt.wave().hpr(), lf, oup00H.get());

            for (auto i = lf - 1; i < 2 * len0 + lf - 1; ++i) {
                oup01[i - lf + 1] = oup00L[i] + oup00H[i];
            }

            circshift(oup01.get(), 2 * len0, -1);

            auto index2 = 0;

            for (auto index = count; index < n; index += value) {
                swtop[index] = (oup00[index2] + oup01[index2]) / 2.0;
                index2++;
            }
        }
        for (auto i = 0; i < n; ++i) {
            appxSig[i] = swtop[i];
        }
    }
}

static auto modwtPer(WaveletTransform& wt, int m, double const* inp, double* cA, int lenCA, double* cD) -> void
{
    auto const lenAvg = wt.wave().lpdLen();
    auto filt = std::make_unique<double[]>(2 * lenAvg);
    auto s = std::sqrt(2.0);

    for (auto i = 0; i < lenAvg; ++i) {
        filt[i] = wt.wave().lpd()[i] / s;
        filt[lenAvg + i] = wt.wave().hpd()[i] / s;
    }

    for (auto i = 0; i < lenCA; ++i) {
        auto t = i;
        cA[i] = filt[0] * inp[t];
        cD[i] = filt[lenAvg] * inp[t];
        for (auto l = 1; l < lenAvg; l++) {
            t -= m;
            while (t >= lenCA) {
                t -= lenCA;
            }
            while (t < 0) {
                t += lenCA;
            }

            cA[i] += filt[l] * inp[t];
            cD[i] += filt[lenAvg + l] * inp[t];
        }
    }
}

static auto modwtDirect(WaveletTransform& wt, double const* inp) -> void
{
    if (wt.extension() != SignalExtension::periodic) {
        printf("MODWT direct method only uses periodic extension per. \n");
        printf(" Use MODWT fft method for symmetric extension sym \n");
        exit(-1);
    }

    auto tempLen = wt.signalLength();
    auto j = wt.levels();
    wt.length[0] = wt.length[j] = tempLen;
    wt.outlength = wt.length[j + 1] = (j + 1) * tempLen;
    auto m = 1;
    for (auto iter = 1; iter < j; ++iter) {
        m = 2 * m;
        wt.length[iter] = tempLen;
    }

    auto cA = std::make_unique<double[]>(tempLen);
    auto cD = std::make_unique<double[]>(tempLen);

    m = 1;

    for (auto i = 0; i < tempLen; ++i) {
        wt.params[i] = inp[i];
    }

    auto lenacc = wt.outlength;

    for (auto iter = 0; iter < j; ++iter) {
        lenacc -= tempLen;
        if (iter > 0) {
            m = 2 * m;
        }

        modwtPer(wt, m, wt.params.get(), cA.get(), tempLen, cD.get());

        for (auto i = 0; i < tempLen; ++i) {
            wt.params[i] = cA[i];
            wt.params[lenacc + i] = cD[i];
        }
    }
}

static auto modwtFft(WaveletTransform& wt, double const* inp) -> void
{
    int j = 0;
    int iter = 0;
    int m = 0;
    int lenacc = 0;
    double s = NAN;
    double tmp1 = NAN;
    double tmp2 = NAN;

    auto tempLen = wt.signalLength();
    auto lenAvg = wt.wave().lpdLen();
    int n { 0 };
    if (wt.extension() == SignalExtension::symmetric) {
        n = 2 * tempLen;
    } else if (wt.extension() == SignalExtension::periodic) {
        n = tempLen;
    }
    j = wt.levels();
    wt.modwtsiglength = n;
    wt.length[0] = wt.length[j] = n;
    wt.outlength = wt.length[j + 1] = (j + 1) * n;

    s = std::sqrt(2.0);
    for (iter = 1; iter < j; ++iter) {
        wt.length[iter] = n;
    }

    auto fftFd = std::make_unique<FFT>(n, FFT::forward);
    auto fftBd = std::make_unique<FFT>(n, FFT::backward);

    auto sig = std::make_unique<Complex<double>[]>(n);
    auto cA = std::make_unique<Complex<double>[]>(n);
    auto cD = std::make_unique<Complex<double>[]>(n);
    auto lowPass = std::make_unique<Complex<double>[]>(n);
    auto highPass = std::make_unique<Complex<double>[]>(n);
    auto index = std::make_unique<int[]>(n);

    // N-point FFT of low pass and high pass filters

    // Low Pass Filter

    for (auto i = 0; i < lenAvg; ++i) {
        sig[i].real((double)wt.wave().lpd()[i] / s);
        sig[i].imag(0.0);
    }
    for (auto i = lenAvg; i < n; ++i) {
        sig[i].real(0.0);
        sig[i].imag(0.0);
    }

    fftFd->perform(sig.get(), lowPass.get());

    // High Pass Filter

    for (auto i = 0; i < lenAvg; ++i) {
        sig[i].real((double)wt.wave().hpd()[i] / s);
        sig[i].imag(0.0);
    }
    for (auto i = lenAvg; i < n; ++i) {
        sig[i].real(0.0);
        sig[i].imag(0.0);
    }

    fftFd->perform(sig.get(), highPass.get());

    // symmetric extension
    for (auto i = 0; i < tempLen; ++i) {
        sig[i].real((double)inp[i]);
        sig[i].imag(0.0);
    }
    for (auto i = tempLen; i < n; ++i) {
        sig[i].real((double)inp[n - i - 1]);
        sig[i].imag(0.0);
    }

    // FFT of data

    fftFd->perform(sig.get(), cA.get());

    lenacc = wt.outlength;

    m = 1;

    for (iter = 0; iter < j; ++iter) {
        lenacc -= n;

        for (auto i = 0; i < n; ++i) {
            index[i] = (m * i) % n;
        }

        for (auto i = 0; i < n; ++i) {
            tmp1 = cA[i].real();
            tmp2 = cA[i].imag();
            cA[i].real(lowPass[index[i]].real() * tmp1 - lowPass[index[i]].imag() * tmp2);
            cA[i].imag(lowPass[index[i]].real() * tmp2 + lowPass[index[i]].imag() * tmp1);

            cD[i].real(highPass[index[i]].real() * tmp1 - highPass[index[i]].imag() * tmp2);
            cD[i].imag(highPass[index[i]].real() * tmp2 + highPass[index[i]].imag() * tmp1);
        }

        fftBd->perform(cD.get(), sig.get());

        for (auto i = 0; i < n; ++i) {
            wt.params[lenacc + i] = sig[i].real() / n;
        }

        m *= 2;
    }

    fftBd->perform(cA.get(), sig.get());

    for (auto i = 0; i < n; ++i) {
        wt.params[i] = sig[i].real() / n;
    }
}

auto modwt(WaveletTransform& wt, double const* inp) -> void
{
    if (wt.convMethod() == ConvolutionMethod::direct) {
        modwtDirect(wt, inp);
        return;
    }

    modwtFft(wt, inp);
}

static auto conjComplex(Complex<double>* x, int n) -> void
{
    for (auto i = 0; i < n; ++i) {
        x[i].imag(x[i].imag() * -1.0);
    }
}

auto imodwtFft(WaveletTransform& wt, double* oup) -> void
{
    auto n = wt.modwtsiglength;
    auto lenAvg = wt.wave().lpdLen();
    auto j = wt.levels();

    auto s = std::sqrt(2.0);
    auto fftFd = std::make_unique<FFT>(n, FFT::forward);
    auto fftBd = std::make_unique<FFT>(n, FFT::backward);

    auto sig = std::make_unique<Complex<double>[]>(n);
    auto cA = std::make_unique<Complex<double>[]>(n);
    auto cD = std::make_unique<Complex<double>[]>(n);
    auto lowPass = std::make_unique<Complex<double>[]>(n);
    auto highPass = std::make_unique<Complex<double>[]>(n);
    auto index = std::make_unique<int[]>(n);

    // N-point FFT of low pass and high pass filters

    // Low Pass Filter

    for (auto i = 0; i < lenAvg; ++i) {
        sig[i].real((double)wt.wave().lpd()[i] / s);
        sig[i].imag(0.0);
    }
    for (auto i = lenAvg; i < n; ++i) {
        sig[i].real(0.0);
        sig[i].imag(0.0);
    }

    fftFd->perform(sig.get(), lowPass.get());

    // High Pass Filter

    for (auto i = 0; i < lenAvg; ++i) {
        sig[i].real((double)wt.wave().hpd()[i] / s);
        sig[i].imag(0.0);
    }
    for (auto i = lenAvg; i < n; ++i) {
        sig[i].real(0.0);
        sig[i].imag(0.0);
    }

    fftFd->perform(sig.get(), highPass.get());

    // Complex conjugate of the two filters

    conjComplex(lowPass.get(), n);
    conjComplex(highPass.get(), n);

    auto m = (int)std::pow(2.0, (double)j - 1.0);
    auto lenacc = n;

    //
    for (auto i = 0; i < n; ++i) {
        sig[i].real((double)wt.output()[i]);
        sig[i].imag(0.0);
    }

    for (auto iter = 0; iter < j; ++iter) {
        fftFd->perform(sig.get(), cA.get());
        for (auto i = 0; i < n; ++i) {
            sig[i].real(wt.output()[lenacc + i]);
            sig[i].imag(0.0);
        }
        fftFd->perform(sig.get(), cD.get());

        for (auto i = 0; i < n; ++i) {
            index[i] = (m * i) % n;
        }

        for (auto i = 0; i < n; ++i) {
            auto const tmp1 = cA[i].real();
            auto const tmp2 = cA[i].imag();
            cA[i].real(lowPass[index[i]].real() * tmp1 - lowPass[index[i]].imag() * tmp2 + highPass[index[i]].real() * cD[i].real() - highPass[index[i]].imag() * cD[i].imag());
            cA[i].imag(lowPass[index[i]].real() * tmp2 + lowPass[index[i]].imag() * tmp1 + highPass[index[i]].real() * cD[i].imag() + highPass[index[i]].imag() * cD[i].real());
        }

        fftBd->perform(cA.get(), sig.get());

        for (auto i = 0; i < n; ++i) {
            sig[i].real(sig[i].real() / n);
            sig[i].imag(sig[i].imag() / n);
        }
        m /= 2;
        lenacc += n;
    }

    for (auto i = 0; i < wt.signalLength(); ++i) {
        oup[i] = sig[i].real();
    }
}

static auto imodwtPer(WaveletTransform& wt, int m, double const* cA, int lenCA, double const* cD, double* x) -> void
{
    auto const lenAvg = wt.wave().lpdLen();
    auto filt = std::make_unique<double[]>(2 * lenAvg);
    auto s = std::sqrt(2.0);

    for (auto i = 0; i < lenAvg; ++i) {
        filt[i] = wt.wave().lpd()[i] / s;
        filt[lenAvg + i] = wt.wave().hpd()[i] / s;
    }

    for (auto i = 0; i < lenCA; ++i) {
        auto t = i;
        x[i] = (filt[0] * cA[t]) + (filt[lenAvg] * cD[t]);
        for (auto l = 1; l < lenAvg; l++) {
            t += m;
            while (t >= lenCA) {
                t -= lenCA;
            }
            while (t < 0) {
                t += lenCA;
            }

            x[i] += (filt[l] * cA[t]) + (filt[lenAvg + l] * cD[t]);
        }
    }
}

static auto imodwtDirect(WaveletTransform& wt, double* dwtop) -> void
{
    auto n = wt.signalLength();
    auto lenacc = n;

    auto j = wt.levels();
    auto m = (int)std::pow(2.0, (double)j - 1.0);

    auto x = std::make_unique<double[]>(n);

    for (auto i = 0; i < n; ++i) {
        dwtop[i] = wt.output()[i];
    }

    for (auto iter = 0; iter < j; ++iter) {
        if (iter > 0) {
            m = m / 2;
        }
        imodwtPer(wt, m, dwtop, n, wt.params.get() + lenacc, x.get());
        /*
		for (auto j = lf - 1; j < N; ++j) {
			dwtop[j - lf + 1] = X[j];
		}
		for (auto j = 0; j < lf - 1; ++j) {
			dwtop[N - lf + 1 + j] = X[j];
		}
		*/
        for (auto jj = 0; jj < n; ++jj) {
            dwtop[jj] = x[jj];
        }

        lenacc += n;
    }
}

auto imodwt(WaveletTransform& wt, double* oup) -> void
{
    if (wt.convMethod() == ConvolutionMethod::direct) {
        imodwtDirect(wt, oup);
        return;
    }
    imodwtFft(wt, oup);
}

auto summary(WaveletTransform const& wt) -> void
{
    int j = 0;
    int t = 0;
    j = wt.levels();
    summary(wt.wave());
    printf("\n");
    printf("Wavelet Transform : %s \n", wt.method().c_str());
    printf("Signal Extension : %s \n", toString(wt.extension()).c_str());
    printf("Convolutional Method : %s \n", toString(wt.convMethod()).c_str());
    printf("Number of Decomposition Levels %d \n", wt.levels());
    printf("Length of Input Signal %d \n", wt.signalLength());
    printf("Length of WT Output Vector %d \n", wt.outlength);
    printf("Wavelet Coefficients are contained in vector : %s \n", "output");
    printf("Approximation Coefficients \n");
    printf("Level %d Access : output[%d] Length : %d \n", j, 0, wt.length[0]);
    printf("Detail Coefficients \n");
    t = wt.length[0];
    for (auto i = 0; i < j; ++i) {
        printf("Level %d Access : output[%d] Length : %d \n", j - i, t, wt.length[i + 1]);
        t += wt.length[i + 1];
    }
    printf("\n");
}