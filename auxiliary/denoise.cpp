#include "waux.h"
#include "wauxlib.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string_view>

using namespace std::string_view_literals;

auto denoise_init(int length, int J, char const* wname) -> denoise_set*
{

    auto obj = std::make_unique<denoise_set>();

    obj->N = length;
    obj->J = J;

    obj->wname = wname;

    //Set Default Values
    obj->dmethod = "sureshrink";
    obj->ext = "sym";
    obj->level = "all";
    obj->thresh = "soft";
    obj->wmethod = "dwt";
    obj->cmethod = "direct";

    return obj.release();
}

void visushrink(double* signal, int N, int J, char const* wname, char const* method, char const* ext, char const* thresh, char const* level, double* denoised)
{
    int dwt_len;
    int sgn;
    int it;
    double sigma;
    double td;
    double tmp;

    auto wave = wavelet { wname };
    auto filt_len = wave.filtlength;
    auto MaxIter = (int)(std::log((double)N / ((double)filt_len - 1.0)) / std::log(2.0));

    if (J > MaxIter) {
        std::printf("\n Error - The Signal Can only be iterated %d times using this wavelet. Exiting\n", MaxIter);
        std::exit(-1);
    }

    wt_set* wt = wt_init(wave, method, N, J);
    if (method == "dwt"sv) {
        setDWTExtension(wt, ext);
        dwt(wt, signal);
    } else if (method == "swt"sv) {
        swt(wt, signal);
    } else {
        std::printf("Acceptable WT methods are - dwt,swt and modwt\n");
        std::exit(-1);
    }

    auto lnoise = std::make_unique<double[]>(J);

    //Set sigma

    auto iter = wt->length[0];
    auto dlen = wt->length[J];

    auto dout = std::make_unique<double[]>(dlen);

    if (level == "first"sv) {
        for (auto i = 1; i < J; ++i) {
            iter += wt->length[i];
        }

        for (auto i = 0; i < dlen; ++i) {
            dout[i] = fabs(wt->output[iter + i]);
        }

        sigma = median(dout.get(), dlen) / 0.6745;
        for (it = 0; it < J; ++it) {
            lnoise[it] = sigma;
        }
    } else if (level == "all"sv) {
        for (it = 0; it < J; ++it) {
            dlen = wt->length[it + 1];
            for (auto i = 0; i < dlen; ++i) {
                dout[i] = fabs(wt->output[iter + i]);
            }
            sigma = median(dout.get(), dlen) / 0.6745;
            lnoise[it] = sigma;
            iter += dlen;
        }

    } else {
        std::printf("Acceptable Noise estimation level values are - first and all \n");
        std::exit(-1);
    }

    dwt_len = wt->outlength;
    iter = wt->length[0];
    for (it = 0; it < J; ++it) {
        sigma = lnoise[it];
        dlen = wt->length[it + 1];
        td = std::sqrt(2.0 * std::log(dwt_len)) * sigma;

        if (thresh == "hard"sv) {
            for (auto i = 0; i < dlen; ++i) {
                if (fabs(wt->output[iter + i]) < td) {
                    wt->output[iter + i] = 0;
                }
            }
        } else if (thresh == "soft"sv) {
            for (auto i = 0; i < dlen; ++i) {
                if (fabs(wt->output[iter + i]) < td) {
                    wt->output[iter + i] = 0;
                } else {
                    sgn = wt->output[iter + i] >= 0 ? 1 : -1;
                    tmp = sgn * (fabs(wt->output[iter + i]) - td);
                    wt->output[iter + i] = tmp;
                }
            }
        }

        iter += wt->length[it + 1];
    }

    if (method == "dwt"sv) {
        idwt(wt, denoised);
    } else if (method == "swt"sv) {
        iswt(wt, denoised);
    }

    wt_free(wt);
}

void sureshrink(double* signal, int N, int J, char const* wname, char const* method, char const* ext, char const* thresh, char const* level, double* denoised)
{
    int filt_len;
    int it;
    int len;
    int dlen;
    int dwt_len;
    int min_index;
    int sgn;
    int MaxIter;
    int iter;
    double sigma;
    double norm;
    double td;
    double tv;
    double te;
    double ct;
    double thr;
    double temp;
    double x_sum;
    wt_set* wt;

    auto wave = wavelet { wname };
    filt_len = wave.filtlength;

    MaxIter = (int)(std::log((double)N / ((double)filt_len - 1.0)) / std::log(2.0));
    // Depends on J
    if (J > MaxIter) {
        std::printf("\n Error - The Signal Can only be iterated %d times using this wavelet. Exiting\n", MaxIter);
        std::exit(-1);
    }

    wt = wt_init(wave, method, N, J);

    if (method == "dwt"sv) {
        setDWTExtension(wt, ext);
        dwt(wt, signal);
    } else if (method == "swt"sv) {
        swt(wt, signal);
    } else {
        std::printf("Acceptable WT methods are - dwt and swt\n");
        std::exit(-1);
    }

    len = wt->length[0];
    dlen = wt->length[J];

    auto dout = std::make_unique<double[]>(dlen);
    auto risk = std::make_unique<double[]>(dlen);
    auto dsum = std::make_unique<double[]>(dlen);
    auto lnoise = std::make_unique<double[]>(J);

    iter = wt->length[0];

    if (level == "first"sv) {
        for (auto i = 1; i < J; ++i) {
            iter += wt->length[i];
        }

        for (auto i = 0; i < dlen; ++i) {
            dout[i] = fabs(wt->output[iter + i]);
        }

        sigma = median(dout.get(), dlen) / 0.6745;
        for (it = 0; it < J; ++it) {
            lnoise[it] = sigma;
        }
    } else if (level == "all"sv) {
        for (it = 0; it < J; ++it) {
            dlen = wt->length[it + 1];
            for (auto i = 0; i < dlen; ++i) {
                dout[i] = fabs(wt->output[iter + i]);
            }
            sigma = median(dout.get(), dlen) / 0.6745;
            lnoise[it] = sigma;
            iter += dlen;
        }

    } else {
        std::printf("Acceptable Noise estimation level values are - first and all \n");
        std::exit(-1);
    }

    for (it = 0; it < J; ++it) {
        dwt_len = wt->length[it + 1];
        sigma = lnoise[it];

        if (sigma < 0.00000001) {
            td = 0;
        } else {
            tv = std::sqrt(2.0 * std::log(dwt_len));
            norm = 0.0;
            for (auto i = 0; i < dwt_len; ++i) {
                norm += (wt->output[len + i] * wt->output[len + i] / (sigma * sigma));
            }
            te = (norm - (double)dwt_len) / (double)dwt_len;
            ct = pow(std::log((double)dwt_len) / std::log(2.0), 1.5) / std::sqrt((double)dwt_len);

            if (te < ct) {
                td = tv;
            } else {
                x_sum = 0.0;

                for (auto i = 0; i < dwt_len; ++i) {
                    dout[i] = fabs(wt->output[len + i] / sigma);
                }

                std::sort(dout.get(), dout.get() + dwt_len, std::less<double> {});
                for (auto i = 0; i < dwt_len; ++i) {
                    dout[i] = (dout[i] * dout[i]);
                    x_sum += dout[i];
                    dsum[i] = x_sum;
                }

                for (auto i = 0; i < dwt_len; ++i) {
                    risk[i] = ((double)dwt_len - 2 * ((double)i + 1) + dsum[i] + dout[i] * ((double)dwt_len - 1 - (double)i)) / (double)dwt_len;
                }
                min_index = minindex(risk.get(), dwt_len);
                thr = std::sqrt(dout[min_index]);
                td = thr < tv ? thr : tv;
            }
        }

        td = td * sigma;

        if (thresh == "hard"sv) {
            for (auto i = 0; i < dwt_len; ++i) {
                if (fabs(wt->output[len + i]) < td) {
                    wt->output[len + i] = 0;
                }
            }
        } else if (thresh == "soft"sv) {
            for (auto i = 0; i < dwt_len; ++i) {
                if (fabs(wt->output[len + i]) < td) {
                    wt->output[len + i] = 0;
                } else {
                    sgn = wt->output[len + i] >= 0 ? 1 : -1;
                    temp = sgn * (fabs(wt->output[len + i]) - td);
                    wt->output[len + i] = temp;
                }
            }
        }

        len += wt->length[it + 1];
    }

    if (method == "dwt"sv) {
        idwt(wt, denoised);
    } else if (method == "swt"sv) {
        iswt(wt, denoised);
    }

    wt_free(wt);
}

void modwtshrink(double* signal, int N, int J, char const* wname, char const* cmethod, char const* ext, char const* thresh, double* denoised)
{
    int sgn;
    int it;
    double sigma;
    double td;
    double tmp;
    double M;
    double llen;
    wt_set* wt;

    auto wave = wavelet { wname };
    auto filt_len = wave.filtlength;

    auto MaxIter = (int)(std::log((double)N / ((double)filt_len - 1.0)) / std::log(2.0));

    if (J > MaxIter) {
        std::printf("\n Error - The Signal Can only be iterated %d times using this wavelet. Exiting\n", MaxIter);
        std::exit(-1);
    }

    wt = wt_init(wave, "modwt", N, J);

    if ((ext == "sym"sv) && (cmethod == "fft"sv)) {
        setWTConv(wt, "fft");
        setDWTExtension(wt, "sym");
    } else if ((ext == "sym"sv) && (cmethod == "direct"sv)) {
        std::printf("Symmetric Extension is not available for direct method");
        std::exit(-1);
    } else if ((ext == "per"sv) && (cmethod == "direct"sv)) {
        setWTConv(wt, "direct");
        setDWTExtension(wt, "per");
    } else if ((ext == "per"sv) && (cmethod == "fft"sv)) {
        setWTConv(wt, "fft");
        setDWTExtension(wt, "per");
    } else {
        std::printf("Signal extension can be either per or sym");
        std::exit(-1);
    }

    modwt(wt, signal);

    auto lnoise = std::make_unique<double[]>(J);

    //Set sigma

    auto iter = wt->length[0];
    auto dlen = wt->length[J];
    auto dout = std::make_unique<double[]>(dlen);

    for (it = 0; it < J; ++it) {
        dlen = wt->length[it + 1];
        for (auto i = 0; i < dlen; ++i) {
            dout[i] = fabs(wt->output[iter + i]);
        }

        sigma = std::sqrt(2.0) * median(dout.get(), dlen) / 0.6745;
        lnoise[it] = sigma;
        iter += dlen;
    }

    M = pow(2.0, J);
    llen = std::log((double)wt->modwtsiglength);
    // Thresholding

    iter = wt->length[0];
    for (it = 0; it < J; ++it) {
        sigma = lnoise[it];
        dlen = wt->length[it + 1];
        td = std::sqrt(2.0 * llen / M) * sigma;

        if (thresh == "hard"sv) {
            for (auto i = 0; i < dlen; ++i) {
                if (fabs(wt->output[iter + i]) < td) {
                    wt->output[iter + i] = 0;
                }
            }
        } else if (thresh == "soft"sv) {
            for (auto i = 0; i < dlen; ++i) {
                if (fabs(wt->output[iter + i]) < td) {
                    wt->output[iter + i] = 0;
                } else {
                    sgn = wt->output[iter + i] >= 0 ? 1 : -1;
                    tmp = sgn * (fabs(wt->output[iter + i]) - td);
                    wt->output[iter + i] = tmp;
                }
            }
        }

        iter += wt->length[it + 1];
        M /= 2.0;
    }

    imodwt(wt, denoised);

    wt_free(wt);
}

void denoise(denoise_set* obj, double* signal, double* denoised)
{
    if (obj->dmethod == "sureshrink"sv) {
        if (obj->wmethod == "modwt"sv) {
            std::printf("sureshrink method only works with swt and dwt. Please use setDenoiseWTMethod to set the correct method\n");
            std::exit(-1);
        }
        sureshrink(signal, obj->N, obj->J, obj->wname.c_str(), obj->wmethod.c_str(), obj->ext.c_str(), obj->thresh.c_str(), obj->level.c_str(), denoised);
    } else if (obj->dmethod == "visushrink"sv) {
        if (obj->wmethod == "modwt"sv) {
            std::printf("visushrink method only works with swt and dwt. Please use setDenoiseWTMethod to set the correct method\n");
            std::exit(-1);
        }
        visushrink(signal, obj->N, obj->J, obj->wname.c_str(), obj->wmethod.c_str(), obj->ext.c_str(), obj->thresh.c_str(), obj->level.c_str(), denoised);
    } else if (obj->dmethod == "modwtshrink"sv) {
        if (obj->wmethod != "modwt"sv) {
            std::printf("modwtshrink method only works with modwt. Please use setDenoiseWTMethod to set the correct method\n");
            std::exit(-1);
        }
        modwtshrink(signal, obj->N, obj->J, obj->wname.c_str(), obj->cmethod.c_str(), obj->ext.c_str(), obj->thresh.c_str(), denoised);
    } else {
        std::printf("Acceptable Denoising methods are - sureshrink and visushrink\n");
        std::exit(-1);
    }
}

void setDenoiseMethod(denoise_set* obj, char const* dmethod)
{
    if (strcmp(dmethod, "sureshrink") == 0) {
        obj->dmethod = "sureshrink";
    } else if (strcmp(dmethod, "visushrink") == 0) {
        obj->dmethod = "visushrink";
    } else if (strcmp(dmethod, "modwtshrink") == 0) {
        obj->dmethod = "modwtshrink";
    } else {
        std::printf("Acceptable Denoising methods are - sureshrink, visushrink and modwtshrink\n");
        std::exit(-1);
    }
}

void setDenoiseWTMethod(denoise_set* obj, char const* wmethod)
{
    obj->wmethod = wmethod;
    if (!((wmethod == "dwt"sv) || (wmethod == "swt"sv) || (wmethod == "modwt"sv))) {
        std::printf("Wavelet decomposition method can be one of dwt, modwt or swt.\n");
        std::exit(-1);
    }
}

void setDenoiseWTExtension(denoise_set* obj, char const* extension)
{
    if (strcmp(extension, "sym") == 0) {
        obj->ext = "sym";
    } else if (strcmp(extension, "per") == 0) {
        obj->ext = "per";
    } else {
        std::printf("Signal extension can be either per or sym");
        std::exit(-1);
    }
}

void setDenoiseParameters(denoise_set* obj, char const* thresh, char const* level)
{

    //Set thresholding
    if (thresh == "soft"sv) {
        obj->thresh = "soft";
    } else if (thresh == "hard"sv) {
        obj->thresh = "hard";
    } else {
        std::printf("Thresholding Method - soft or hard");
        std::exit(-1);
    }

    // Set Noise estimation at the first level or at all levels

    if (level == "first"sv) {
        obj->level = "first";
    } else if (level == "all"sv) {
        obj->level = "all";
    } else {
        std::printf("Noise Estimation at level - first or all");
        std::exit(-1);
    }
}

void denoise_free(denoise_set* object)
{
    delete object;
}
