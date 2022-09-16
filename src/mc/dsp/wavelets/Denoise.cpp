#include "Denoise.hpp"

#include <mc/core/algorithm.hpp>
#include <mc/core/cmath.hpp>
#include <mc/core/cstdlib.hpp>
#include <mc/core/cstring.hpp>
#include <mc/core/format.hpp>
#include <mc/core/memory.hpp>
#include <mc/core/numeric.hpp>
#include <mc/core/string_view.hpp>
#include <mc/core/utility.hpp>

#include <fmt/printf.h>

namespace mc::dsp {

DenoiseSet::DenoiseSet(int length, int j, char const* name)
{

    N = length;
    J = j;

    wname = name;

    // Set Default Values
    dmethod = "sureshrink";
    ext     = "sym";
    level   = "all";
    thresh  = "soft";
    wmethod = "dwt";
    cmethod = "direct";
}

auto visushrink(
    float* signal,
    std::size_t n,
    std::size_t j,
    char const* wname,
    char const* method,
    char const* ext,
    char const* thresh,
    char const* level,
    float* denoised
) -> void
{
    int dwtLen  = 0;
    int sgn     = 0;
    int it      = 0;
    float sigma = NAN;
    float td    = NAN;
    float tmp   = NAN;

    auto wave    = Wavelet{wname};
    auto filtLen = wave.size();
    auto maxIter = (int)(std::log((float)n / ((float)filtLen - 1.0F)) / std::log(2.0F));

    if (cmp_greater(j, maxIter)) {
        fmt::printf(
            "\n Error - The Signal Can only be iterated %d times using this Wavelet. "
            "Exiting\n",
            maxIter
        );
        std::exit(-1);
    }

    auto wt = WaveletTransform{wave, method, n, j};
    if (method == StringView{"dwt"}) {
        wt.extension(
            ext == StringView{"per"} ? SignalExtension::periodic
                                     : SignalExtension::symmetric
        );
        dwt(wt, signal);
    } else if (method == StringView{"swt"}) {
        swt(wt, signal);
    } else {
        fmt::printf("Acceptable WT methods are - dwt,swt and modwt\n");
        std::exit(-1);
    }

    auto lnoise = makeUnique<float[]>(j);

    // Set sigma

    auto iter = wt.length[0];
    auto dlen = wt.length[j];

    auto dout = makeUnique<float[]>(dlen);

    if (level == StringView{"first"}) {
        for (std::size_t i = 1; i < j; ++i) { iter += wt.length[i]; }

        for (std::size_t i = 0; i < dlen; ++i) {
            dout[i] = std::abs(wt.output()[iter + i]);
        }

        sigma = median(dout.get(), dlen) / 0.6745F;
        for (it = 0; cmp_less(it, j); ++it) { lnoise[it] = sigma; }
    } else if (level == StringView{"all"}) {
        for (it = 0; cmp_less(it, j); ++it) {
            dlen = wt.length[it + 1];
            for (std::size_t i = 0; i < dlen; ++i) {
                dout[i] = std::abs(wt.output()[iter + i]);
            }
            sigma      = median(dout.get(), dlen) / 0.6745F;
            lnoise[it] = sigma;
            iter += dlen;
        }
    } else {
        fmt::printf("Acceptable Noise estimation level values are - first and all \n");
        std::exit(-1);
    }

    dwtLen = wt.outlength;
    iter   = wt.length[0];
    for (it = 0; cmp_less(it, j); ++it) {
        sigma = lnoise[it];
        dlen  = wt.length[it + 1];
        td    = sqrt(2.0F * std::log(dwtLen)) * sigma;

        if (thresh == StringView{"hard"}) {
            for (std::size_t i = 0; i < dlen; ++i) {
                if (std::abs(wt.output()[iter + i]) < td) { wt.output()[iter + i] = 0; }
            }
        } else if (thresh == StringView{"soft"}) {
            for (std::size_t i = 0; i < dlen; ++i) {
                if (std::abs(wt.output()[iter + i]) < td) {
                    wt.output()[iter + i] = 0;
                } else {
                    sgn                   = wt.output()[iter + i] >= 0 ? 1 : -1;
                    tmp                   = sgn * (std::abs(wt.output()[iter + i]) - td);
                    wt.output()[iter + i] = tmp;
                }
            }
        }

        iter += wt.length[it + 1];
    }

    if (method == StringView{"dwt"}) {
        idwt(wt, denoised);
    } else if (method == StringView{"swt"}) {
        iswt(wt, denoised);
    }
}

auto sureshrink(
    float* signal,
    std::size_t n,
    std::size_t j,
    char const* wname,
    char const* method,
    char const* ext,
    char const* thresh,
    char const* level,
    float* denoised
) -> void
{
    int filtLen  = 0;
    int it       = 0;
    int len      = 0;
    int dlen     = 0;
    int dwtLen   = 0;
    int minIndex = 0;
    int sgn      = 0;
    int maxIter  = 0;
    int iter     = 0;
    float sigma  = NAN;
    float norm   = NAN;
    float td     = NAN;
    float tv     = NAN;
    float te     = NAN;
    float ct     = NAN;
    float thr    = NAN;
    float temp   = NAN;
    float xSum   = NAN;

    auto wave = Wavelet{wname};
    filtLen   = wave.size();

    maxIter = (int)(std::log((float)n / ((float)filtLen - 1.0F)) / std::log(2.0F));
    // Depends on J
    if (cmp_greater(j, maxIter)) {
        fmt::printf(
            "\n Error - The Signal Can only be iterated %d times using this Wavelet. "
            "Exiting\n",
            maxIter
        );
        std::exit(-1);
    }

    auto wt = WaveletTransform(wave, method, n, j);

    if (method == StringView{"dwt"}) {
        wt.extension(
            ext == StringView{"per"} ? SignalExtension::periodic
                                     : SignalExtension::symmetric
        );
        dwt(wt, signal);
    } else if (method == StringView{"swt"}) {
        swt(wt, signal);
    } else {
        fmt::printf("Acceptable WT methods are - dwt and swt\n");
        std::exit(-1);
    }

    len  = wt.length[0];
    dlen = wt.length[j];

    auto dout   = makeUnique<float[]>(dlen);
    auto risk   = makeUnique<float[]>(dlen);
    auto dsum   = makeUnique<float[]>(dlen);
    auto lnoise = makeUnique<float[]>(j);

    iter = wt.length[0];

    if (level == StringView{"first"}) {
        for (std::size_t i = 1; i < j; ++i) { iter += wt.length[i]; }

        for (auto i = 0; i < dlen; ++i) { dout[i] = std::abs(wt.output()[iter + i]); }

        sigma = median(dout.get(), dlen) / 0.6745F;
        for (it = 0; cmp_less(it, j); ++it) { lnoise[it] = sigma; }
    } else if (level == StringView{"all"}) {
        for (it = 0; cmp_less(it, j); ++it) {
            dlen = wt.length[it + 1];
            for (auto i = 0; i < dlen; ++i) { dout[i] = std::abs(wt.output()[iter + i]); }
            sigma      = median(dout.get(), dlen) / 0.6745F;
            lnoise[it] = sigma;
            iter += dlen;
        }
    } else {
        fmt::printf("Acceptable Noise estimation level values are - first and all \n");
        std::exit(-1);
    }

    for (it = 0; cmp_less(it, j); ++it) {
        dwtLen = wt.length[it + 1];
        sigma  = lnoise[it];

        if (sigma < 0.00000001) {
            td = 0;
        } else {
            tv   = mc::sqrt(2.0F * mc::log((float)dwtLen));
            norm = 0.0F;
            for (auto i = 0; i < dwtLen; ++i) {
                norm += (wt.output()[len + i] * wt.output()[len + i] / (sigma * sigma));
            }
            te = (norm - (float)dwtLen) / (float)dwtLen;
            ct = mc::pow(mc::log((float)dwtLen) / mc::log(2.0F), 1.5F)
               / mc::sqrt((float)dwtLen);

            if (te < ct) {
                td = tv;
            } else {
                xSum = 0.0F;

                for (auto i = 0; i < dwtLen; ++i) {
                    dout[i] = std::abs(wt.output()[len + i] / sigma);
                }

                std::sort(dout.get(), dout.get() + dwtLen, std::less<float>{});
                for (auto i = 0; i < dwtLen; ++i) {
                    dout[i] = (dout[i] * dout[i]);
                    xSum += dout[i];
                    dsum[i] = xSum;
                }

                for (auto i = 0; i < dwtLen; ++i) {
                    risk[i] = ((float)dwtLen - 2 * ((float)i + 1) + dsum[i]
                               + dout[i] * ((float)dwtLen - 1 - (float)i))
                            / (float)dwtLen;
                }
                minIndex = minindex(risk.get(), dwtLen);
                thr      = sqrt(dout[minIndex]);
                td       = thr < tv ? thr : tv;
            }
        }

        td = td * sigma;

        if (thresh == StringView{"hard"}) {
            for (auto i = 0; i < dwtLen; ++i) {
                if (std::abs(wt.output()[len + i]) < td) { wt.output()[len + i] = 0; }
            }
        } else if (thresh == StringView{"soft"}) {
            for (auto i = 0; i < dwtLen; ++i) {
                if (std::abs(wt.output()[len + i]) < td) {
                    wt.output()[len + i] = 0;
                } else {
                    sgn                  = wt.output()[len + i] >= 0 ? 1 : -1;
                    temp                 = sgn * (std::abs(wt.output()[len + i]) - td);
                    wt.output()[len + i] = temp;
                }
            }
        }

        len += wt.length[it + 1];
    }

    if (method == StringView{"dwt"}) {
        idwt(wt, denoised);
    } else if (method == StringView{"swt"}) {
        iswt(wt, denoised);
    }
}

auto modwtshrink(
    float* signal,
    std::size_t n,
    std::size_t j,
    char const* wname,
    char const* cmethod,
    char const* ext,
    char const* thresh,
    float* denoised
) -> void
{
    int sgn     = 0;
    int it      = 0;
    float sigma = NAN;
    float td    = NAN;
    float tmp   = NAN;
    float m     = NAN;
    float llen  = NAN;

    auto wave    = Wavelet{wname};
    auto filtLen = wave.size();

    auto maxIter = (int)(std::log((float)n / ((float)filtLen - 1.0F)) / std::log(2.0F));

    if (cmp_greater(j, maxIter)) {
        fmt::printf(
            "\n Error - The Signal Can only be iterated %d times using this Wavelet. "
            "Exiting\n",
            maxIter
        );
        std::exit(-1);
    }

    auto wt = WaveletTransform(wave, "modwt", n, j);

    if ((ext == StringView{"sym"}) && (cmethod == StringView{"fft"})) {
        wt.convMethod(ConvolutionMethod::fft);
        wt.extension(SignalExtension::symmetric);
    } else if ((ext == StringView{"sym"}) && (cmethod == StringView{"direct"})) {
        fmt::printf("Symmetric Extension is not available for direct method");
        std::exit(-1);
    } else if ((ext == StringView{"per"}) && (cmethod == StringView{"direct"})) {
        wt.convMethod(ConvolutionMethod::direct);
        wt.extension(SignalExtension::periodic);
    } else if ((ext == StringView{"per"}) && (cmethod == StringView{"fft"})) {
        wt.convMethod(ConvolutionMethod::fft);
        wt.extension(SignalExtension::periodic);
    } else {
        fmt::printf("Signal extension can be either per or sym");
        std::exit(-1);
    }

    modwt(wt, signal);

    auto lnoise = makeUnique<float[]>(j);

    // Set sigma

    auto iter = wt.length[0];
    auto dlen = wt.length[j];
    auto dout = makeUnique<float[]>(dlen);

    for (it = 0; cmp_less(it, j); ++it) {
        dlen = wt.length[it + 1];
        for (std::size_t i = 0; i < dlen; ++i) {
            dout[i] = std::abs(wt.output()[iter + i]);
        }

        sigma      = sqrt(2.0F) * median(dout.get(), dlen) / 0.6745F;
        lnoise[it] = sigma;
        iter += dlen;
    }

    m    = pow(2.0F, j);
    llen = std::log((float)wt.modwtsiglength);
    // Thresholding

    iter = wt.length[0];
    for (it = 0; cmp_less(it, j); ++it) {
        sigma = lnoise[it];
        dlen  = wt.length[it + 1];
        td    = sqrt(2.0F * llen / m) * sigma;

        if (thresh == StringView{"hard"}) {
            for (std::size_t i = 0; i < dlen; ++i) {
                if (std::abs(wt.output()[iter + i]) < td) { wt.output()[iter + i] = 0; }
            }
        } else if (thresh == StringView{"soft"}) {
            for (std::size_t i = 0; i < dlen; ++i) {
                if (std::abs(wt.output()[iter + i]) < td) {
                    wt.output()[iter + i] = 0;
                } else {
                    sgn                   = wt.output()[iter + i] >= 0 ? 1 : -1;
                    tmp                   = sgn * (std::abs(wt.output()[iter + i]) - td);
                    wt.output()[iter + i] = tmp;
                }
            }
        }

        iter += wt.length[it + 1];
        m /= 2.0F;
    }

    imodwt(wt, denoised);
}

auto denoise(DenoiseSet& obj, float* signal, float* denoised) -> void
{
    if (obj.dmethod == "sureshrink") {
        if (obj.wmethod == StringView{"MODWT"}) {
            fmt::printf("sureshrink method only works with swt and dwt. Please use "
                        "setDenoiseWTMethod to set the "
                        "correct method\n");
            std::exit(-1);
        }
        sureshrink(
            signal,
            obj.N,
            obj.J,
            obj.wname.c_str(),
            obj.wmethod.c_str(),
            obj.ext.c_str(),
            obj.thresh.c_str(),
            obj.level.c_str(),
            denoised
        );
    } else if (obj.dmethod == "visushrink") {
        if (obj.wmethod == StringView{"MODWT"}) {
            fmt::printf("visushrink method only works with swt and dwt. Please use "
                        "setDenoiseWTMethod to set the "
                        "correct method\n");
            std::exit(-1);
        }
        visushrink(
            signal,
            obj.N,
            obj.J,
            obj.wname.c_str(),
            obj.wmethod.c_str(),
            obj.ext.c_str(),
            obj.thresh.c_str(),
            obj.level.c_str(),
            denoised
        );
    } else if (obj.dmethod == "modwtshrink") {
        if (obj.wmethod != StringView{"MODWT"}) {
            fmt::printf("modwtshrink method only works with modwt. Please use "
                        "setDenoiseWTMethod to set the correct method\n");
            std::exit(-1);
        }
        modwtshrink(
            signal,
            obj.N,
            obj.J,
            obj.wname.c_str(),
            obj.cmethod.c_str(),
            obj.ext.c_str(),
            obj.thresh.c_str(),
            denoised
        );
    } else {
        fmt::printf("Acceptable Denoising methods are - sureshrink and visushrink\n");
        std::exit(-1);
    }
}

auto setDenoiseMethod(DenoiseSet& obj, char const* dmethod) -> void
{
    if (strcmp(dmethod, "sureshrink") == 0) {
        obj.dmethod = "sureshrink";
    } else if (strcmp(dmethod, "visushrink") == 0) {
        obj.dmethod = "visushrink";
    } else if (strcmp(dmethod, "modwtshrink") == 0) {
        obj.dmethod = "modwtshrink";
    } else {
        fmt::printf(
            "Acceptable Denoising methods are - sureshrink, visushrink and modwtshrink\n"
        );
        std::exit(-1);
    }
}

auto setDenoiseWTMethod(DenoiseSet& obj, char const* wmethod) -> void
{
    obj.wmethod = wmethod;
    if (!((wmethod == StringView{"dwt"}) || (wmethod == StringView{"swt"})
          || (wmethod == StringView{"MODWT"}))) {
        fmt::printf("Wavelet decomposition method can be one of dwt, modwt or swt.\n");
        std::exit(-1);
    }
}

auto setDenoiseWTExtension(DenoiseSet& obj, char const* extension) -> void
{
    if (strcmp(extension, "sym") == 0) {
        obj.ext = "sym";
    } else if (strcmp(extension, "per") == 0) {
        obj.ext = "per";
    } else {
        fmt::printf("Signal extension can be either per or sym");
        std::exit(-1);
    }
}

auto setDenoiseParameters(DenoiseSet& obj, char const* thresh, char const* level) -> void
{

    // Set thresholding
    if (thresh == StringView{"soft"}) {
        obj.thresh = "soft";
    } else if (thresh == StringView{"hard"}) {
        obj.thresh = "hard";
    } else {
        fmt::printf("Thresholding Method - soft or hard");
        std::exit(-1);
    }

    // Set Noise estimation at the first level or at all levels

    if (level == StringView{"first"}) {
        obj.level = "first";
    } else if (level == StringView{"all"}) {
        obj.level = "all";
    } else {
        fmt::printf("Noise Estimation at level - first or all");
        std::exit(-1);
    }
}

auto median(float* const x, int n) -> float
{
    std::sort(x, x + n, std::less<float>{});

    float sigma = NAN;
    if ((n % 2) == 0) {
        sigma = (x[n / 2 - 1] + x[n / 2]) / 2.0F;
    } else {
        sigma = x[n / 2];
    }

    return sigma;
}

auto minindex(float const* arr, int n) -> int
{
    auto min   = std::numeric_limits<float>::max();
    auto index = 0;
    for (auto i = 0; i < n; ++i) {
        if (arr[i] < min) {
            min   = arr[i];
            index = i;
        }
    }

    return index;
}

}  // namespace mc::dsp
