#include "ContinuousWaveletTransform.hpp"

#include "lt/dsp/convolution/FFTConvolver.hpp"
#include "lt/dsp/fft/FFT.hpp"
#include "lt/dsp/wavelets/common.hpp"

#include "lt/cassert.hpp"
#include "lt/cmath.hpp"
#include "lt/cstdlib.hpp"
#include "lt/cstring.hpp"
#include "lt/format.hpp"
#include "lt/memory.hpp"
#include "lt/string_view.hpp"

namespace lt {
namespace dsp {
    namespace {
        auto factorial(int n) -> float
        {
            static const float fact[] = {
                1.0F,
                1.0F,
                2.0F,
                6.0F,
                24.0F,
                120.0F,
                720.0F,
                5040.0F,
                40320.0F,
                362880.0F,
                3628800.0F,
                39916800.0F,
                479001600.0F,
                6227020800.0F,
                87178291200.0F,
                1307674368000.0F,
                20922789888000.0F,
                355687428096000.0F,
                6402373705728000.0F,
                121645100408832000.0F,
                2432902008176640000.0F,
                51090942171709440000.0F,
                1124000727777607680000.0F,
                25852016738884976640000.0F,
                620448401733239439360000.0F,
                15511210043330985984000000.0F,
                403291461126605635584000000.0F,
                10888869450418352160768000000.0F,
                304888344611713860501504000000.0F,
                8841761993739701954543616000000.0F,
                265252859812191058636308480000000.0F,
                8222838654177922817725562880000000.0F,
                263130836933693530167218012160000000.0F,
                8683317618811886495518194401280000000.0F,
                295232799039604140847618609643520000000.0F,
                // 10333147966386144929666651337523200000000.0f,
                // 371993326789901217467999448150835200000000.0f,
                // 13763753091226345046315979581580902400000000.0f,
                // 523022617466601111760007224100074291200000000.0f,
                // 20397882081197443358640281739902897356800000000.0f,
                // 815915283247897734345611269596115894272000000000.0f,
            };

            if (n > 40 - 6 || n < 0) {
                fmt::print("This program is only valid for 0 <= N <= 40-6 \n");
                return -1.0F;
            }

            return fact[n];
        }

        auto nsfftFd(FFT* obj, Complex<float>* inp, Complex<float>* oup, float lb, float ub, float* w) -> void
        {
            auto const n = obj->size();
            auto const l = n / 2;

            auto temp1 = std::make_unique<float[]>(l);
            auto temp2 = std::make_unique<float[]>(l);

            auto const delta = (ub - lb) / n;
            auto const den = 2 * (ub - lb);
            auto j = -n;

            for (auto i = 0; i < n; ++i) {
                w[i] = (float)j / den;
                j += 2;
            }

            obj->perform(inp, oup);

            for (auto i = 0; i < l; ++i) {
                temp1[i] = oup[i].real();
                temp2[i] = oup[i].imag();
            }

            for (auto i = 0; i < n - l; ++i) {
                oup[i].real(oup[i + l].real());
                oup[i].imag(oup[i + l].imag());
            }

            for (auto i = 0; i < l; ++i) {
                oup[n - l + i].real(temp1[i]);
                oup[n - l + i].imag(temp2[i]);
            }

            auto const plb = pi2 * lb;

            for (auto i = 0; i < n; ++i) {
                auto const tempr = oup[i].real();
                auto const tempi = oup[i].imag();
                auto const theta = w[i] * plb;

                oup[i].real(delta * (tempr * cos(theta) + tempi * sin(theta)));
                oup[i].imag(delta * (tempi * cos(theta) - tempr * sin(theta)));
            }
        }

        auto nsfftBk(FFT* obj, Complex<float>* inp, Complex<float>* oup, float lb, float ub, float* t) -> void
        {

            auto const n = obj->size();
            auto const l = n / 2;

            auto const m = divideby(n, 2);

            if (m == 0) {
                fmt::printf("The Non-Standard FFT Length must be a power of 2");
                exit(1);
            }

            auto temp1 = std::make_unique<float[]>(l);
            auto temp2 = std::make_unique<float[]>(l);
            auto w = std::make_unique<float[]>(n);
            auto inpt = std::make_unique<Complex<float>[]>(n);

            auto const delta = (ub - lb) / n;
            auto const den = 2 * (ub - lb);
            auto j = -n;

            for (auto i = 0; i < n; ++i) {
                w[i] = (float)j / den;
                j += 2;
            }

            auto const plb = pi2 * lb;

            for (auto i = 0; i < n; ++i) {
                auto const theta = w[i] * plb;
                inpt[i].real((inp[i].real() * cos(theta) - inp[i].imag() * sin(theta)) / delta);
                inpt[i].imag((inp[i].imag() * cos(theta) + inp[i].real() * sin(theta)) / delta);
            }

            for (auto i = 0; i < l; ++i) {
                temp1[i] = inpt[i].real();
                temp2[i] = inpt[i].imag();
            }

            for (auto i = 0; i < n - l; ++i) {
                inpt[i].real(inpt[i + l].real());
                inpt[i].imag(inpt[i + l].imag());
            }

            for (auto i = 0; i < l; ++i) {
                inpt[n - l + i].real(temp1[i]);
                inpt[n - l + i].imag(temp2[i]);
            }

            obj->perform(inpt.get(), oup);

            for (auto i = 0; i < n; ++i) {
                t[i] = lb + i * delta;
            }
        }

        auto nsfftExec(FFT* fft, Complex<float>* inp, Complex<float>* oup, float lb, float ub, float* w) -> void
        {
            if (fft->direction() == FFT::forward) {
                nsfftFd(fft, inp, oup, lb, ub, w);
            } else {
                nsfftBk(fft, inp, oup, lb, ub, w);
            }
        }

        auto roundTowardsZero(float x) -> float
        {
            return x >= 0.0F ? std::floor(x) : std::ceil(x);
        }

        auto cwtGamma(float x) -> float
        {
            /*
	 * This C program code is based on  W J Cody's fortran code.
	 * http://www.netlib.org/specfun/gamma
	 *
	 * References:
   "An Overview of Software Development for Special Functions",
	W. J. Cody, Lecture Notes in Mathematics, 506,
	Numerical Analysis Dundee, 1975, G. A. Watson (ed.),
	Springer Verlag, Berlin, 1976.

   Computer Approximations, Hart, Et. Al., Wiley and sons, New York, 1968.
   */

            // numerator and denominator coefficients for 1 <= x <= 2

            float y = NAN;
            float oup = NAN;
            float fact = NAN;
            float sum = NAN;
            float y2 = NAN;
            float yi = NAN;
            float z = NAN;
            float nsum = NAN;
            float dsum = NAN;
            int swi = 0;
            int n = 0;

            constexpr float spi = 0.9189385332046727417803297F;
            constexpr float pi = 3.1415926535897932384626434F;
            constexpr float xmax = 171.624e+0F;
            constexpr float xinf = 1.79e38F;
            constexpr float eps = 2.22e-16F;
            constexpr float xninf = 1.41e-45F;

            float num[8] = {
                -1.71618513886549492533811e+0F,
                2.47656508055759199108314e+1F,
                -3.79804256470945635097577e+2F,
                6.29331155312818442661052e+2F,
                8.66966202790413211295064e+2F,
                -3.14512729688483675254357e+4F,
                -3.61444134186911729807069e+4F,
                6.64561438202405440627855e+4F,
            };

            float den[8] = {
                -3.08402300119738975254353e+1F,
                3.15350626979604161529144e+2F,
                -1.01515636749021914166146e+3F,
                -3.10777167157231109440444e+3F,
                2.25381184209801510330112e+4F,
                4.75584627752788110767815e+3F,
                -1.34659959864969306392456e+5F,
                -1.15132259675553483497211e+5F,
            };

            // Coefficients for Hart's Minimax approximation x >= 12

            float c[7] = {
                -1.910444077728e-03F,
                8.4171387781295e-04F,
                -5.952379913043012e-04F,
                7.93650793500350248e-04F,
                -2.777777777777681622553e-03F,
                8.333333333333333331554247e-02F,
                5.7083835261e-03F,
            };

            y = x;
            swi = 0;
            fact = 1.0F;
            n = 0;

            if (y < 0.0F) {
                // Negative x
                y = -x;
                yi = roundTowardsZero(y);
                oup = y - yi;

                if (oup != 0.0F) {
                    if (yi != roundTowardsZero(yi * 0.5F) * 2.0F) {
                        swi = 1;
                    }
                    fact = -pi / std::sin(pi * oup);
                    y += 1.0F;
                } else {
                    return xinf;
                }
            }

            if (y < eps) {
                if (y >= xninf) {
                    oup = 1.0F / y;
                } else {
                    return xinf;
                }

            } else if (y < 12.0F) {
                yi = y;
                if (y < 1.0F) {
                    z = y;
                    y += 1.0F;
                } else {
                    n = (int)y - 1;
                    y -= (float)n;
                    z = y - 1.0F;
                }
                nsum = 0.0F;
                dsum = 1.0F;
                for (auto i = 0; i < 8; ++i) {
                    nsum = (nsum + num[i]) * z;
                    dsum = dsum * z + den[i];
                }
                oup = nsum / dsum + 1.0F;

                if (yi < y) {

                    oup /= yi;
                } else if (yi > y) {

                    for (auto i = 0; i < n; ++i) {
                        oup *= y;
                        y += 1.0F;
                    }
                }

            } else {
                if (y <= xmax) {
                    y2 = y * y;
                    sum = c[6];
                    for (auto i = 0; i < 6; ++i) {
                        sum = sum / y2 + c[i];
                    }
                    sum = sum / y - y + spi;
                    sum += (y - .5) * std::log(y);
                    oup = std::exp(sum);
                } else {
                    return (xinf);
                }
            }

            if (swi != 0) {
                oup = -oup;
            }
            if (fact != 1.0F) {
                oup = fact / oup;
            }

            return oup;
        }

        auto waveFunction(int nk, float dt, int mother, float param, float scale1, float const* kwave, float pi, float* period1,
            float* coi1, Complex<float>* daughter) -> void
        {

            float norm = NAN;
            float expnt = NAN;
            float fourierFactor = NAN;
            int k = 0;
            int m = 0;
            float temp = NAN;
            int sign = 0;
            int re = 0;

            if (mother == 0) {
                //MORLET
                if (param < 0.0F) {
                    param = 6.0F;
                }
                norm = std::sqrt(2.0F * pi * scale1 / dt) * std::pow(pi, -0.25);

                for (k = 1; k <= nk / 2 + 1; ++k) {
                    temp = (scale1 * kwave[k - 1] - param);
                    expnt = -0.5 * temp * temp;
                    daughter[k - 1].real(norm * std::exp(expnt));
                    daughter[k - 1].imag(0.0F);
                }
                for (k = nk / 2 + 2; k <= nk; ++k) {
                    daughter[k - 1].real(0.0F);
                    daughter[k - 1].imag(0.0F);
                }
                fourierFactor = (4.0F * pi) / (param + std::sqrt(2.0F + param * param));
                *period1 = scale1 * fourierFactor;
                *coi1 = fourierFactor / std::sqrt(2.0F);
            } else if (mother == 1) {
                // PAUL
                if (param < 0.0F) {
                    param = 4.0F;
                }
                m = (int)param;
                norm = std::sqrt(2.0F * pi * scale1 / dt) * (std::pow(2.0F, (float)m) / std::sqrt((float)(m * factorial(2 * m - 1))));
                for (k = 1; k <= nk / 2 + 1; ++k) {
                    temp = scale1 * kwave[k - 1];
                    expnt = -temp;
                    daughter[k - 1].real(norm * std::pow(temp, (float)m) * std::exp(expnt));
                    daughter[k - 1].imag(0.0F);
                }
                for (k = nk / 2 + 2; k <= nk; ++k) {
                    daughter[k - 1].real(0.0F);
                    daughter[k - 1].imag(0.0F);
                }
                fourierFactor = (4.0F * pi) / (2.0F * m + 1.0F);
                *period1 = scale1 * fourierFactor;
                *coi1 = fourierFactor * std::sqrt(2.0F);
            } else if (mother == 2) {
                if (param < 0.0F) {
                    param = 2.0F;
                }
                m = (int)param;

                if (m % 2 == 0) {
                    re = 1;
                } else {
                    re = 0;
                }

                if (m % 4 == 0 || m % 4 == 1) {
                    sign = -1;
                } else {
                    sign = 1;
                }

                norm = std::sqrt(2.0F * pi * scale1 / dt) * std::sqrt(1.0F / cwtGamma(m + 0.50));
                norm *= sign;

                if (re == 1) {
                    for (k = 1; k <= nk; ++k) {
                        temp = scale1 * kwave[k - 1];
                        daughter[k - 1].real(norm * std::pow(temp, (float)m) * std::exp(-0.50 * std::pow(temp, 2.0F)));
                        daughter[k - 1].imag(0.0F);
                    }
                } else if (re == 0) {
                    for (k = 1; k <= nk; ++k) {
                        temp = scale1 * kwave[k - 1];
                        daughter[k - 1].real(0.0F);
                        daughter[k - 1].imag(norm * std::pow(temp, (float)m) * std::exp(-0.50 * std::pow(temp, 2.0F)));
                    }
                }
                fourierFactor = (2.0F * pi) * std::sqrt(2.0F / (2.0F * m + 1.0F));
                *period1 = scale1 * fourierFactor;
                *coi1 = fourierFactor / std::sqrt(2.0F);
            }
        }

    }

    auto cwavelet(float const* y, int n, float dt, int mother, float param, float s0, float dj, int jtot, int npad,
        float* wave, float const* scale, float* period, float* coi) -> void
    {

        int j = 0;
        int k = 0;
        int iter = 0;
        float ymean = NAN;
        float freq1 = NAN;
        float pi = NAN;
        float period1 = NAN;
        float coi1 = NAN;
        float tmp1 = NAN;
        float tmp2 = NAN;
        float scale1 = NAN;

        (void)s0;
        (void)dj; /* yes, we need these parameters unused */

        pi = 4.0F * std::atan(1.0F);

        if (npad < n) {
            fmt::printf("npad must be >= N \n");
            exit(-1);
        }

        auto obj = std::make_unique<FFT>(npad, FFT::forward);
        auto iobj = std::make_unique<FFT>(npad, FFT::backward);

        auto ypad = std::make_unique<Complex<float>[]>(npad);
        auto yfft = std::make_unique<Complex<float>[]>(npad);
        auto daughter = std::make_unique<Complex<float>[]>(npad);
        auto kwave = std::make_unique<float[]>(npad);
        ymean = 0.0F;

        for (auto i = 0; i < n; ++i) {
            ymean += y[i];
        }

        ymean /= n;

        for (auto i = 0; i < n; ++i) {
            ypad[i].real(y[i] - ymean);
            ypad[i].imag(0.0F);
        }

        for (auto i = n; i < npad; ++i) {
            ypad[i].real(0.0F);
            ypad[i].imag(0.0F);
        }

        // Find FFT of the input y (ypad)

        obj->perform(ypad.get(), yfft.get());

        for (auto i = 0; i < npad; ++i) {
            yfft[i].real(yfft[i].real() / (float)npad);
            yfft[i].imag(yfft[i].imag() / (float)npad);
        }

        //Construct the wavenumber array

        freq1 = 2.0F * pi / ((float)npad * dt);
        kwave[0] = 0.0F;

        for (auto i = 1; i < npad / 2 + 1; ++i) {
            kwave[i] = i * freq1;
        }

        for (auto i = npad / 2 + 1; i < npad; ++i) {
            kwave[i] = -kwave[npad - i];
        }

        // Main loop

        for (j = 1; j <= jtot; ++j) {
            scale1 = scale[j - 1]; // = s0*std::pow(2.0F, (float)(j - 1)*dj);
            waveFunction(npad, dt, mother, param, scale1, kwave.get(), pi, &period1, &coi1, daughter.get());
            period[j - 1] = period1;
            for (k = 0; k < npad; ++k) {
                tmp1 = daughter[k].real() * yfft[k].real() - daughter[k].imag() * yfft[k].imag();
                tmp2 = daughter[k].real() * yfft[k].imag() + daughter[k].imag() * yfft[k].real();
                daughter[k].real(tmp1);
                daughter[k].imag(tmp2);
            }
            iobj->perform(daughter.get(), ypad.get());
            iter = 2 * (j - 1) * n;
            for (auto i = 0; i < n; ++i) {
                wave[iter + 2 * i] = ypad[i].real();
                wave[iter + 2 * i + 1] = ypad[i].imag();
            }
        }

        for (auto i = 1; i <= (n + 1) / 2; ++i) {
            coi[i - 1] = coi1 * dt * ((float)i - 1.0F);
            coi[n - i] = coi[i - 1];
        }
    }

    auto psi0(int mother, float param, float* val, int* real) -> void
    {
        float pi = NAN;
        float coeff = NAN;
        int m = 0;
        int sign = 0;

        m = (int)param;
        pi = 4.0F * std::atan(1.0F);

        if (mother == 0) {
            // Morlet
            *val = 1.0F / std::sqrt(std::sqrt(pi));
            *real = 1;
        } else if (mother == 1) {
            //Paul
            if (m % 2 == 0) {
                *real = 1;
            } else {
                *real = 0;
            }

            if (m % 4 == 0 || m % 4 == 1) {
                sign = 1;
            } else {
                sign = -1;
            }
            *val = sign * std::pow(2.0F, (float)m) * factorial(m) / (std::sqrt(pi * factorial(2 * m)));

        } else if (mother == 2) {
            // D.O.G
            *real = 1;

            if (m % 2 == 0) {
                if (m % 4 == 0) {
                    sign = -1;
                } else {
                    sign = 1;
                }
                coeff = sign * std::pow(2.0F, (float)m / 2) / cwtGamma(0.5);
                *val = coeff * cwtGamma(((float)m + 1.0F) / 2.0F) / std::sqrt(cwtGamma(m + 0.50));
            } else {
                *val = 0;
            }
        }
    }

    static auto maxabs(float* array, int n) -> int
    {
        float maxval = NAN;
        float temp = NAN;

        int index = 0;
        maxval = 0.0F;
        index = -1;

        for (auto i = 0; i < n; ++i) {
            temp = std::fabs(array[i]);
            if (temp >= maxval) {
                maxval = temp;
                index = i;
            }
        }

        return index;
    }

    auto cdelta(int mother, float param, float psi0) -> float
    {
        auto s0 { 0.0F };
        auto n { 1 };
        auto subscale = 8.0F;
        auto dt = 0.25;

        if (mother == 0) {
            n = 16;
            s0 = dt / 4;
        } else if (mother == 1) {
            n = 16;
            s0 = dt / 4.0F;
        } else if (mother == 2) {
            s0 = dt / 8.0F;
            n = 256;
            if (param == 2.0F) {
                subscale = 16.0F;
                s0 = dt / 16.0F;
                n = 2048;
            }
        }

        auto const dj = 1.0F / subscale;
        auto const jtot = 16 * (int)subscale;

        auto delta = std::make_unique<float[]>(n);
        auto wave = std::make_unique<float[]>(2 * n * jtot);
        auto coi = std::make_unique<float[]>(n);
        auto scale = std::make_unique<float[]>(jtot);
        auto period = std::make_unique<float[]>(jtot);
        auto mval = std::make_unique<float[]>(n);

        delta[0] = 1;

        for (auto i = 1; i < n; ++i) {
            delta[i] = 0;
        }

        for (auto i = 0; i < jtot; ++i) {
            scale[i] = s0 * std::pow(2.0F, (float)(i)*dj);
        }

        cwavelet(delta.get(), n, dt, mother, param, s0, dj, jtot, n, wave.get(), scale.get(), period.get(), coi.get());

        for (auto i = 0; i < n; ++i) {
            mval[i] = 0;
        }

        for (auto j = 0; j < jtot; ++j) {
            auto const iter = 2 * j * n;
            auto const den = std::sqrt(scale[j]);
            for (auto i = 0; i < n; ++i) {
                mval[i] += wave[iter + 2 * i] / den;
            }
        }

        auto const maxarr = maxabs(mval.get(), n);
        auto const cdel = std::sqrt(dt) * dj * mval[maxarr] / psi0;
        return cdel;
    }

    auto icwavelet(float const* wave, int n, float* scale, int jtot, float dt, float dj, float cdelta, float psi0, float* oup) -> void
    {

        int j = 0;
        int iter = 0;
        float den = NAN;
        float coeff = NAN;

        coeff = std::sqrt(dt) * dj / (cdelta * psi0);

        for (auto i = 0; i < n; ++i) {
            oup[i] = 0.0F;
        }

        for (j = 0; j < jtot; ++j) {
            iter = 2 * j * n;
            den = std::sqrt(scale[j]);
            for (auto i = 0; i < n; ++i) {
                oup[i] += wave[iter + 2 * i] / den;
            }
        }

        for (auto i = 0; i < n; ++i) {
            oup[i] *= coeff;
        }
    }

    ContinuousWaveletTransform::ContinuousWaveletTransform(char const* wave, float param, int siglength, float dtIn, int j)
    {

        float s0Tmp {};
        float djTmp {};
        char const* pdefault = "pow";

        auto const mm = (int)param;
        auto odd = 1;
        if (2 * (mm / 2) == mm) {
            odd = 0;
        }

        auto n = siglength;
        auto nj2 = 2 * n * j;
        params = std::make_unique<float[]>(nj2 + 2 * j + n);

        int motherTmp { 0 };
        if ((wave == string_view { "morlet" }) || (wave == string_view { "morl" })) {
            s0Tmp = 2 * dtIn;
            djTmp = 0.4875;
            motherTmp = 0;
            if (param < 0.0F) {
                fmt::printf("\n Morlet Wavelet Parameter should be >= 0 \n");
                exit(-1);
            }
            if (param == 0) {
                param = 6.0F;
            }
            wave_ = "morlet";

        } else if (wave == string_view { "paul" }) {
            s0Tmp = 2 * dtIn;
            djTmp = 0.4875;
            motherTmp = 1;
            if (param < 0 || param > 20) {
                fmt::printf("\n Paul Wavelet Parameter should be > 0 and <= 20 \n");
                exit(-1);
            }
            if (param == 0) {
                param = 4.0F;
            }
            wave_ = "paul";

        } else if ((wave == string_view { "dgauss" }) || (wave == string_view { "dog" })) {
            s0Tmp = 2 * dtIn;
            djTmp = 0.4875;
            motherTmp = 2;
            if (param < 0 || odd == 1) {
                fmt::printf("\n DOG Wavelet Parameter should be > 0 and even \n");
                exit(-1);
            }
            if (param == 0) {
                param = 2.0F;
            }
            wave_ = "dog";
        }

        pow = 2;
        type = pdefault;

        s0 = s0Tmp;
        dj = djTmp;
        dt = dtIn;
        J = j;
        signalLength = siglength;
        sflag = 0;
        pflag = 1;
        mother = motherTmp;
        m = param;

        auto const t1 = 0.499999 + std::log((float)n) / std::log(2.0F);
        auto const ibase2 = 1 + (int)t1;

        npad = (int)std::pow(2.0F, (float)ibase2);

        output = (Complex<float>*)&params[0];
        scale = &params[nj2];
        period = &params[nj2 + j];
        coi = &params[nj2 + 2 * j];

        for (auto i = 0; i < nj2 + 2 * j + n; ++i) {
            params[i] = 0.0F;
        }
    }

    auto ContinuousWaveletTransform::scales(float newS0, float newDj, char const* newType, int power) -> void
    {
        type = newType;
        //newS0*std::pow(2.0F, (float)(j - 1)*newDj);
        if ((type == "pow") || (type == "power")) {
            for (auto i = 0; i < J; ++i) {
                scale[i] = newS0 * std::pow((float)power, (float)(i)*newDj);
            }
            sflag = 1;
            pow = power;

        } else if ((type == "lin") || (type == "linear")) {
            for (auto i = 0; i < J; ++i) {
                scale[i] = newS0 + (float)i * newDj;
            }
            sflag = 1;
        } else {
            fmt::printf("\n Type accepts only two values : pow and lin\n");
            exit(-1);
        }
        s0 = newS0;
        dj = newDj;
    }

    auto cwt(ContinuousWaveletTransform& wt, float const* inp) -> void
    {
        int n = 0;
        int npad = 0;
        int nj2 = 0;
        int j = 0;
        int j2 = 0;
        n = wt.signalLength;
        if (wt.sflag == 0) {
            for (auto i = 0; i < wt.J; ++i) {
                wt.scale[i] = wt.s0 * std::pow(2.0F, (float)(i)*wt.dj);
            }
            wt.sflag = 1;
        }

        if (wt.pflag == 0) {
            npad = n;
        } else {
            npad = wt.npad;
        }

        nj2 = 2 * n * wt.J;
        j = wt.J;
        j2 = 2 * j;

        wt.smean = 0.0F;

        for (auto i = 0; i < n; ++i) {
            wt.smean += inp[i];
        }
        wt.smean /= n;

        cwavelet(inp, n, wt.dt, wt.mother, wt.m, wt.s0, wt.dj, wt.J, npad, wt.params.get(), wt.params.get() + nj2, wt.params.get() + nj2 + j, wt.params.get() + nj2 + j2);
    }

    auto icwt(ContinuousWaveletTransform& wt, float* cwtop) -> void
    {
        float psi = NAN;
        float cdel = NAN;
        int real = 0;
        int n = 0;
        int nj2 = 0;

        n = wt.signalLength;
        nj2 = n * 2 * wt.J;

        psi0(wt.mother, wt.m, &psi, &real);
        cdel = cdelta(wt.mother, wt.m, psi);

        if (((wt.type == "pow") || (wt.type == "power")) && wt.pow == 2) {
            icwavelet(wt.params.get(), n, wt.params.get() + nj2, wt.J, wt.dt, wt.dj, cdel, psi, cwtop);
        } else {
            fmt::printf("Inverse CWT is only available for power of 2.0F scales \n");
            exit(-1);
        }
        for (auto i = 0; i < n; ++i) {
            cwtop[i] += wt.smean;
        }
    }

    auto summary(ContinuousWaveletTransform const& wt) -> void
    {
        fmt::printf("\nWavelet : %s Parameter %lf \n", wt.wave().c_str(), wt.m);
        fmt::printf("\nLength of Input Signal : %d \n", wt.signalLength);
        fmt::printf("\nSampling Rate : %g \n", wt.dt);
        fmt::printf("\nTotal Number of Scales : %d \n", wt.J);
        fmt::printf("\nSmallest Scale (s0) : %lf \n", wt.s0);
        fmt::printf("\nSeparation Between Scales (dj) %lf \n", wt.dj);
        fmt::printf("\nScale Type %s \n", wt.type.c_str());
        fmt::printf("\nComplex CWT Output Vector is of size %d  const& %d stored in Row Major format \n", wt.J, wt.signalLength);
        fmt::printf("\nThe ith real value can be accessed using wt.output()[i].real() and imaginary value by wt.output()[i].imag() \n");
        fmt::printf("\n");
    }

    auto meyer(int n, float lb, float ub, float* phi, float* psi, float* tgrid) -> void
    {
        auto const m = divideby(n, 2);
        if (m == 0) {
            fmt::printf("Size of Wavelet must be a power of 2");
            exit(1);
        }
        if (lb >= ub) {
            fmt::printf("upper bound must be greater than lower bound");
            exit(1);
        }

        auto obj = std::make_unique<FFT>(n, FFT::backward);
        auto w = std::make_unique<float[]>(n);
        auto phiw = std::make_unique<Complex<float>[]>(n);
        auto psiw = std::make_unique<Complex<float>[]>(n);
        auto oup = std::make_unique<Complex<float>[]>(n);

        auto const delta = 2 * (ub - lb) / pi2;

        auto j = (float)n;
        j *= -1.0F;

        for (auto i = 0; i < n; ++i) {
            w[i] = j / delta;
            j += 2.0F;
            psiw[i].real(0.0F);
            psiw[i].imag(0.0F);
            phiw[i].real(0.0F);
            phiw[i].imag(0.0F);
        }

        for (auto i = 0; i < n; ++i) {
            auto const wf = fabs(w[i]);
            if (wf <= pi2 / 3.0F) {
                phiw[i].real(1.0F);
            }
            if (wf > pi2 / 3.0F && wf <= 2 * pi2 / 3.0F) {
                auto const x = (3 * wf / pi2) - 1.0F;
                auto const x2 = x * x;
                auto const x3 = x2 * x;
                auto const x4 = x3 * x;
                auto const v = x4 * (35 - 84 * x + 70 * x2 - 20 * x3);
                auto const theta = v * pi2 / 4.0F;
                auto const cs = std::cos(theta);
                auto const sn = std::sin(theta);

                phiw[i].real(cs);
                psiw[i].real(std::cos(w[i] / 2.0F) * sn);
                psiw[i].imag(std::sin(w[i] / 2.0F) * sn);
            }
            if (wf > 2.0F * pi2 / 3.0F && wf <= 4 * pi2 / 3.0F) {
                auto const x = (1.5 * wf / pi2) - 1.0F;
                auto const x2 = x * x;
                auto const x3 = x2 * x;
                auto const x4 = x3 * x;
                auto const v = x4 * (35 - 84 * x + 70 * x2 - 20 * x3);
                auto const theta = v * pi2 / 4.0F;
                auto const cs = std::cos(theta);

                psiw[i].real(std::cos(w[i] / 2.0F) * cs);
                psiw[i].imag(std::sin(w[i] / 2.0F) * cs);
            }
        }

        nsfftExec(obj.get(), phiw.get(), oup.get(), lb, ub, tgrid);

        for (auto i = 0; i < n; ++i) {
            phi[i] = oup[i].real() / n;
        }

        nsfftExec(obj.get(), psiw.get(), oup.get(), lb, ub, tgrid);

        for (auto i = 0; i < n; ++i) {
            psi[i] = oup[i].real() / n;
        }
    }

    auto gauss(int n, int p, float lb, float ub, float* psi, float* t) -> void
    {
        float delta = NAN;
        float num = NAN;
        float den = NAN;
        float t2 = NAN;
        float t4 = NAN;
        int i = 0;

        if (lb >= ub) {
            fmt::printf("upper bound must be greater than lower bound");
            exit(1);
        }

        t[0] = lb;
        t[n - 1] = ub;
        delta = (ub - lb) / (n - 1);
        for (i = 1; i < n - 1; ++i) {
            t[i] = lb + delta * i;
        }

        den = std::sqrt(cwtGamma(p + 0.5));

        if ((p + 1) % 2 == 0) {
            num = 1.0F;
        } else {
            num = -1.0F;
        }

        num /= den;

        //fmt::printf("\n%g\n",num);

        if (p == 1) {
            for (i = 0; i < n; ++i) {
                psi[i] = -t[i] * std::exp(-t[i] * t[i] / 2.0F) * num;
            }
        } else if (p == 2) {
            for (i = 0; i < n; ++i) {
                t2 = t[i] * t[i];
                psi[i] = (-1.0F + t2) * std::exp(-t2 / 2.0F) * num;
            }
        } else if (p == 3) {
            for (i = 0; i < n; ++i) {
                t2 = t[i] * t[i];
                psi[i] = t[i] * (3.0F - t2) * std::exp(-t2 / 2.0F) * num;
            }
        } else if (p == 4) {
            for (i = 0; i < n; ++i) {
                t2 = t[i] * t[i];
                psi[i] = (t2 * t2 - 6.0F * t2 + 3.0F) * std::exp(-t2 / 2.0F) * num;
            }
        } else if (p == 5) {
            for (i = 0; i < n; ++i) {
                t2 = t[i] * t[i];
                psi[i] = t[i] * (-t2 * t2 + 10.0F * t2 - 15.0F) * std::exp(-t2 / 2.0F) * num;
            }
        } else if (p == 6) {
            for (i = 0; i < n; ++i) {
                t2 = t[i] * t[i];
                psi[i] = (t2 * t2 * t2 - 15.0F * t2 * t2 + 45.0F * t2 - 15.0F) * std::exp(-t2 / 2.0F) * num;
            }
        } else if (p == 7) {
            for (i = 0; i < n; ++i) {
                t2 = t[i] * t[i];
                psi[i] = t[i] * (-t2 * t2 * t2 + 21.0F * t2 * t2 - 105.0F * t2 + 105.0F) * std::exp(-t2 / 2.0F) * num;
            }
        } else if (p == 8) {
            for (i = 0; i < n; ++i) {
                t2 = t[i] * t[i];
                t4 = t2 * t2;
                psi[i] = (t4 * t4 - 28.0F * t4 * t2 + 210.0F * t4 - 420.0F * t2 + 105.0F) * std::exp(-t2 / 2.0F) * num;
            }
        } else if (p == 9) {
            for (i = 0; i < n; ++i) {
                t2 = t[i] * t[i];
                t4 = t2 * t2;
                psi[i] = t[i] * (-t4 * t4 + 36.0F * t4 * t2 - 378.0F * t4 + 1260.0F * t2 - 945.0F) * std::exp(-t2 / 2.0F) * num;
            }
        } else if (p == 10) {
            for (i = 0; i < n; ++i) {
                t2 = t[i] * t[i];
                t4 = t2 * t2;
                psi[i] = (t4 * t4 * t2 - 45.0F * t4 * t4 + 630.0F * t4 * t2 - 3150.0F * t4 + 4725.0F * t2 - 945.0F) * std::exp(-t2 / 2.0F) * num;
            }
        } else {
            fmt::printf("\n The Gaussian Derivative Wavelet is only available for Derivatives 1 to 10");
            exit(1);
        }
    }

    auto mexhat(int n, float lb, float ub, float* psi, float* t) -> void
    {
        gauss(n, 2, lb, ub, psi, t);
    }

    auto morlet(int n, float lb, float ub, float* psi, float* t) -> void
    {
        int i = 0;
        float delta = NAN;

        if (lb >= ub) {
            fmt::printf("upper bound must be greater than lower bound");
            exit(1);
        }

        t[0] = lb;
        t[n - 1] = ub;
        delta = (ub - lb) / (n - 1);
        for (i = 1; i < n - 1; ++i) {
            t[i] = lb + delta * i;
        }

        for (i = 0; i < n; ++i) {
            psi[i] = std::exp(-t[i] * t[i] / 2.0F) * std::cos(5 * t[i]);
        }
    }

}
}