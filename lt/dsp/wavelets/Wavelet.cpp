#include "Wavelet.hpp"

#include "lt/dsp/convolution/Convolution.hpp"
#include "lt/dsp/wavelets/filters/coif.hpp"
#include "lt/dsp/wavelets/filters/daubechies.hpp"
#include "lt/dsp/wavelets/filters/h.hpp"
#include "lt/dsp/wavelets/filters/meyer.hpp"
#include "lt/dsp/wavelets/filters/sym.hpp"

#define USE_MATH_DEFINES
#include <cmath>
#include <cstdio>
#include <cstring>
#include <memory>
#include <string_view>

using namespace std::string_view_literals;

namespace {
auto copy(double const* in, int n, double* out)
{
    int count = 0;
    for (count = 0; count < n; count++) {
        out[count] = in[count];
    }
}

auto copyReverse(double const* in, int n, double* out)
{
    int count = 0;
    for (count = 0; count < n; count++) {
        out[count] = in[n - count - 1];
    }
}

auto qmfEven(double const* in, int n, double* out)
{
    int count = 0;
    for (count = 0; count < n; count++) {
        out[count] = in[n - count - 1];
        if (count % 2 != 0) {
            out[count] = -1 * out[count];
        }
    }
}

auto qmfWrev(double const* in, int n, double* out)
{
    auto sigOutTemp = std::make_unique<double[]>(n);

    qmfEven(in, n, sigOutTemp.get());
    copyReverse(sigOutTemp.get(), n, out);
}

auto waveletFilterLength(char const* name) -> int
{
    int len = strlen(name);
    int i = 0;
    int n = 0;
    if (name == "haar"sv || name == "db1"sv) {
        return 2;
    }
    if (len > 2 && strstr(name, "db") != nullptr) {
        auto newStr = std::make_unique<char[]>((len - 2 + 1));
        for (i = 2; i < len + 1; i++) {
            newStr[i - 2] = name[i];
        }

        n = atoi(newStr.get());
        if (n > 38) {
            printf("\n Filter Not in Database \n");
            return -1;
        }

        return n * 2;
    }
    if (name == "bior1.1"sv) {
        return 2;
    }

    if (name == "bior1.3"sv) {
        return 6;
    }

    if (name == "bior1.5"sv) {
        return 10;
    }

    if (name == "bior2.2"sv) {
        return 6;
    }

    if (name == "bior2.4"sv) {
        return 10;
    }

    if (name == "bior2.6"sv) {
        return 14;
    }
    if (name == "bior2.8"sv) {
        return 18;
    }

    if (name == "bior3.1"sv) {
        return 4;
    }
    if (name == "bior3.3"sv) {
        return 8;
    }
    if (name == "bior3.5"sv) {
        return 12;
    }

    if (name == "bior3.7"sv) {
        return 16;
    }
    if (name == "bior3.9"sv) {
        return 20;
    }
    if (name == "bior4.4"sv) {
        return 10;
    }
    if (name == "bior5.5"sv) {
        return 12;
    }
    if (name == "bior6.8"sv) {
        return 18;
    }
    if (name == "rbior1.1"sv) {
        return 2;
    }

    if (name == "rbior1.3"sv) {
        return 6;
    }

    if (name == "rbior1.5"sv) {
        return 10;
    }

    if (name == "rbior2.2"sv) {
        return 6;
    }

    if (name == "rbior2.4"sv) {
        return 10;
    }

    if (name == "rbior2.6"sv) {
        return 14;
    }
    if (name == "rbior2.8"sv) {
        return 18;
    }

    if (name == "rbior3.1"sv) {
        return 4;
    }
    if (name == "rbior3.3"sv) {
        return 8;
    }
    if (name == "rbior3.5"sv) {
        return 12;
    }

    if (name == "rbior3.7"sv) {
        return 16;
    }
    if (name == "rbior3.9"sv) {
        return 20;
    }
    if (name == "rbior4.4"sv) {
        return 10;
    }
    if (name == "rbior5.5"sv) {
        return 12;
    }
    if (name == "rbior6.8"sv) {
        return 18;
    }
    if (len > 4 && strstr(name, "coif") != nullptr) {
        auto newStr = std::make_unique<char[]>((len - 4 + 1));
        for (i = 4; i < len + 1; i++) {
            newStr[i - 4] = name[i];
        }

        n = atoi(newStr.get());
        if (n > 17) {
            printf("\n Filter Not in Database \n");
            return -1;
        }

        return n * 6;
    }
    if (len > 3 && strstr(name, "sym") != nullptr) {
        auto newStr = std::make_unique<char[]>((len - 3 + 1));
        for (i = 3; i < len + 1; i++) {
            newStr[i - 3] = name[i];
        }

        n = atoi(newStr.get());
        if (n > 20 || n < 2) {
            printf("\n Filter Not in Database \n");
            return -1;
        }

        return n * 2;
    }
    if (strcmp(name, "meyer") == 0) {
        return 102;
    }
    printf("\n Filter Not in Database \n");
    return -1;
}

auto waveletFilterCoefficients(char const* name, double* lp1, double* hp1, double* lp2, double* hp2) -> int
{
    int i = 0;
    int n = waveletFilterLength(name);
    if ((strcmp(name, "haar") == 0) || (strcmp(name, "db1") == 0)) {
        copyReverse(daubechies1, n, lp1);
        qmfWrev(daubechies1, n, hp1);
        copy(daubechies1, n, lp2);
        qmfEven(daubechies1, n, hp2);

        return n;
    }
    if (name == "db2"sv) {
        copyReverse(daubechies2, n, lp1);
        qmfWrev(daubechies2, n, hp1);
        copy(daubechies2, n, lp2);
        qmfEven(daubechies2, n, hp2);

        return n;
    }
    if (name == "db3"sv) {
        copyReverse(daubechies3, n, lp1);
        qmfWrev(daubechies3, n, hp1);
        copy(daubechies3, n, lp2);
        qmfEven(daubechies3, n, hp2);

        return n;
    }
    if (name == "db4"sv) {
        copyReverse(daubechies4, n, lp1);
        qmfWrev(daubechies4, n, hp1);
        copy(daubechies4, n, lp2);
        qmfEven(daubechies4, n, hp2);

        return n;
    }
    if (name == "db5"sv) {
        copyReverse(daubechies5, n, lp1);
        qmfWrev(daubechies5, n, hp1);
        copy(daubechies5, n, lp2);
        qmfEven(daubechies5, n, hp2);

        return n;
    }
    if (name == "db6"sv) {
        copyReverse(daubechies6, n, lp1);
        qmfWrev(daubechies6, n, hp1);
        copy(daubechies6, n, lp2);
        qmfEven(daubechies6, n, hp2);

        return n;
    }
    if (name == "db7"sv) {
        copyReverse(daubechies7, n, lp1);
        qmfWrev(daubechies7, n, hp1);
        copy(daubechies7, n, lp2);
        qmfEven(daubechies7, n, hp2);

        return n;
    }
    if (name == "db8"sv) {
        copyReverse(daubechies8, n, lp1);
        qmfWrev(daubechies8, n, hp1);
        copy(daubechies8, n, lp2);
        qmfEven(daubechies8, n, hp2);

        return n;
    }
    if (name == "db9"sv) {
        copyReverse(daubechies9, n, lp1);
        qmfWrev(daubechies9, n, hp1);
        copy(daubechies9, n, lp2);
        qmfEven(daubechies9, n, hp2);

        return n;
    }

    if (name == "db10"sv) {
        copyReverse(daubechies10, n, lp1);
        qmfWrev(daubechies10, n, hp1);
        copy(daubechies10, n, lp2);
        qmfEven(daubechies10, n, hp2);

        return n;
    }

    if (name == "db11"sv) {
        copyReverse(daubechies11, n, lp1);
        qmfWrev(daubechies11, n, hp1);
        copy(daubechies11, n, lp2);
        qmfEven(daubechies11, n, hp2);

        return n;
    }
    if (name == "db12"sv) {
        copyReverse(daubechies12, n, lp1);
        qmfWrev(daubechies12, n, hp1);
        copy(daubechies12, n, lp2);
        qmfEven(daubechies12, n, hp2);

        return n;
    }

    if (name == "db13"sv) {
        copyReverse(daubechies13, n, lp1);
        qmfWrev(daubechies13, n, hp1);
        copy(daubechies13, n, lp2);
        qmfEven(daubechies13, n, hp2);

        return n;
    }

    if (name == "db14"sv) {
        copyReverse(daubechies14, n, lp1);
        qmfWrev(daubechies14, n, hp1);
        copy(daubechies14, n, lp2);
        qmfEven(daubechies14, n, hp2);

        return n;
    }

    if (name == "db15"sv) {
        copyReverse(daubechies15, n, lp1);
        qmfWrev(daubechies15, n, hp1);
        copy(daubechies15, n, lp2);
        qmfEven(daubechies15, n, hp2);

        return n;
    }

    if (name == "db16"sv) {
        copyReverse(daubechies16, n, lp1);
        qmfWrev(daubechies16, n, hp1);
        copy(daubechies16, n, lp2);
        qmfEven(daubechies16, n, hp2);

        return n;
    }

    if (name == "db17"sv) {
        copyReverse(daubechies17, n, lp1);
        qmfWrev(daubechies17, n, hp1);
        copy(daubechies17, n, lp2);
        qmfEven(daubechies17, n, hp2);

        return n;
    }

    if (name == "db18"sv) {
        copyReverse(daubechies18, n, lp1);
        qmfWrev(daubechies18, n, hp1);
        copy(daubechies18, n, lp2);
        qmfEven(daubechies18, n, hp2);

        return n;
    }

    if (name == "db19"sv) {
        copyReverse(daubechies19, n, lp1);
        qmfWrev(daubechies19, n, hp1);
        copy(daubechies19, n, lp2);
        qmfEven(daubechies19, n, hp2);

        return n;
    }

    if (name == "db20"sv) {
        copyReverse(daubechies20, n, lp1);
        qmfWrev(daubechies20, n, hp1);
        copy(daubechies20, n, lp2);
        qmfEven(daubechies20, n, hp2);

        return n;
    }

    if (name == "db21"sv) {
        copyReverse(daubechies21, n, lp1);
        qmfWrev(daubechies21, n, hp1);
        copy(daubechies21, n, lp2);
        qmfEven(daubechies21, n, hp2);

        return n;
    }

    if (name == "db22"sv) {
        copyReverse(daubechies22, n, lp1);
        qmfWrev(daubechies22, n, hp1);
        copy(daubechies22, n, lp2);
        qmfEven(daubechies22, n, hp2);

        return n;
    }

    if (name == "db23"sv) {
        copyReverse(daubechies23, n, lp1);
        qmfWrev(daubechies23, n, hp1);
        copy(daubechies23, n, lp2);
        qmfEven(daubechies23, n, hp2);

        return n;
    }

    if (name == "db24"sv) {
        copyReverse(daubechies24, n, lp1);
        qmfWrev(daubechies24, n, hp1);
        copy(daubechies24, n, lp2);
        qmfEven(daubechies24, n, hp2);

        return n;
    }

    if (name == "db25"sv) {
        copyReverse(daubechies25, n, lp1);
        qmfWrev(daubechies25, n, hp1);
        copy(daubechies25, n, lp2);
        qmfEven(daubechies25, n, hp2);

        return n;
    }

    if (name == "db26"sv) {
        copyReverse(daubechies26, n, lp1);
        qmfWrev(daubechies26, n, hp1);
        copy(daubechies26, n, lp2);
        qmfEven(daubechies26, n, hp2);
        return n;
    }

    if (name == "db27"sv) {
        copyReverse(daubechies27, n, lp1);
        qmfWrev(daubechies27, n, hp1);
        copy(daubechies27, n, lp2);
        qmfEven(daubechies27, n, hp2);

        return n;
    }

    if (name == "db28"sv) {
        copyReverse(daubechies28, n, lp1);
        qmfWrev(daubechies28, n, hp1);
        copy(daubechies28, n, lp2);
        qmfEven(daubechies28, n, hp2);

        return n;
    }

    if (name == "db29"sv) {
        copyReverse(daubechies29, n, lp1);
        qmfWrev(daubechies29, n, hp1);
        copy(daubechies29, n, lp2);
        qmfEven(daubechies29, n, hp2);

        return n;
    }

    if (name == "db30"sv) {
        copyReverse(daubechies30, n, lp1);
        qmfWrev(daubechies30, n, hp1);
        copy(daubechies30, n, lp2);
        qmfEven(daubechies30, n, hp2);

        return n;
    }

    if (name == "db31"sv) {
        copyReverse(daubechies31, n, lp1);
        qmfWrev(daubechies31, n, hp1);
        copy(daubechies31, n, lp2);
        qmfEven(daubechies31, n, hp2);

        return n;
    }

    if (name == "db32"sv) {
        copyReverse(daubechies32, n, lp1);
        qmfWrev(daubechies32, n, hp1);
        copy(daubechies32, n, lp2);
        qmfEven(daubechies32, n, hp2);

        return n;
    }

    if (name == "db33"sv) {
        copyReverse(daubechies33, n, lp1);
        qmfWrev(daubechies33, n, hp1);
        copy(daubechies33, n, lp2);
        qmfEven(daubechies33, n, hp2);

        return n;
    }

    if (name == "db34"sv) {
        copyReverse(daubechies34, n, lp1);
        qmfWrev(daubechies34, n, hp1);
        copy(daubechies34, n, lp2);
        qmfEven(daubechies34, n, hp2);

        return n;
    }

    if (name == "db35"sv) {
        copyReverse(daubechies35, n, lp1);
        qmfWrev(daubechies35, n, hp1);
        copy(daubechies35, n, lp2);
        qmfEven(daubechies35, n, hp2);

        return n;
    }

    if (name == "db36"sv) {
        copyReverse(daubechies36, n, lp1);
        qmfWrev(daubechies36, n, hp1);
        copy(daubechies36, n, lp2);
        qmfEven(daubechies36, n, hp2);

        return n;
    }

    if (name == "db37"sv) {
        copyReverse(daubechies37, n, lp1);
        qmfWrev(daubechies37, n, hp1);
        copy(daubechies37, n, lp2);
        qmfEven(daubechies37, n, hp2);

        return n;
    }

    if (name == "db38"sv) {
        copyReverse(daubechies38, n, lp1);
        qmfWrev(daubechies38, n, hp1);
        copy(daubechies38, n, lp2);
        qmfEven(daubechies38, n, hp2);

        return n;
    }

    if (name == "bior1.1"sv) {
        copyReverse(hm111, n, lp1);
        qmfWrev(h1 + 4, n, hp1);
        copy(h1 + 4, n, lp2);
        qmfEven(hm111, n, hp2);
        return n;
    }

    if (name == "bior1.3"sv) {
        copyReverse(hm113, n, lp1);
        qmfWrev(h1 + 2, n, hp1);
        copy(h1 + 2, n, lp2);
        qmfEven(hm113, n, hp2);
        return n;
    }

    if (name == "bior1.5"sv) {
        copyReverse(hm115, n, lp1);
        qmfWrev(h1, n, hp1);
        copy(h1, n, lp2);
        qmfEven(hm115, n, hp2);
        return n;
    }

    if (name == "bior2.2"sv) {
        copyReverse(hm222, n, lp1);
        qmfWrev(h2 + 6, n, hp1);
        copy(h2 + 6, n, lp2);
        qmfEven(hm222, n, hp2);
        return n;
    }

    if (name == "bior2.4"sv) {
        copyReverse(hm224, n, lp1);
        qmfWrev(h2 + 4, n, hp1);
        copy(h2 + 4, n, lp2);
        qmfEven(hm224, n, hp2);
        return n;
    }

    if (name == "bior2.6"sv) {
        copyReverse(hm226, n, lp1);
        qmfWrev(h2 + 2, n, hp1);
        copy(h2 + 2, n, lp2);
        qmfEven(hm226, n, hp2);
        return n;
    }

    if (name == "bior2.8"sv) {
        copyReverse(hm228, n, lp1);
        qmfWrev(h2, n, hp1);
        copy(h2, n, lp2);
        qmfEven(hm228, n, hp2);
        return n;
    }

    if (name == "bior3.1"sv) {
        copyReverse(hm331, n, lp1);
        qmfWrev(h3 + 8, n, hp1);
        copy(h3 + 8, n, lp2);
        qmfEven(hm331, n, hp2);
        return n;
    }

    if (name == "bior3.3"sv) {
        copyReverse(hm333, n, lp1);
        qmfWrev(h3 + 6, n, hp1);
        copy(h3 + 6, n, lp2);
        qmfEven(hm333, n, hp2);
        return n;
    }

    if (name == "bior3.5"sv) {
        copyReverse(hm335, n, lp1);
        qmfWrev(h3 + 4, n, hp1);
        copy(h3 + 4, n, lp2);
        qmfEven(hm335, n, hp2);
        return n;
    }

    if (name == "bior3.7"sv) {
        copyReverse(hm337, n, lp1);
        qmfWrev(h3 + 2, n, hp1);
        copy(h3 + 2, n, lp2);
        qmfEven(hm337, n, hp2);
        return n;
    }

    if (name == "bior3.9"sv) {
        copyReverse(hm339, n, lp1);
        qmfWrev(h3, n, hp1);
        copy(h3, n, lp2);
        qmfEven(hm339, n, hp2);
        return n;
    }

    if (name == "bior4.4"sv) {
        copyReverse(hm444, n, lp1);
        qmfWrev(h4, n, hp1);
        copy(h4, n, lp2);
        qmfEven(hm444, n, hp2);
        return n;
    }

    if (name == "bior5.5"sv) {
        copyReverse(hm555, n, lp1);
        qmfWrev(h5, n, hp1);
        copy(h5, n, lp2);
        qmfEven(hm555, n, hp2);
        return n;
    }

    if (name == "bior6.8"sv) {
        copyReverse(hm668, n, lp1);
        qmfWrev(h6, n, hp1);
        copy(h6, n, lp2);
        qmfEven(hm668, n, hp2);
        return n;
    }

    if (name == "rbior1.1"sv) {
        copyReverse(h1 + 4, n, lp1);
        qmfWrev(hm111, n, hp1);
        copy(hm111, n, lp2);
        qmfEven(h1 + 4, n, hp2);
        return n;
    }

    if (name == "rbior1.3"sv) {
        copyReverse(h1 + 2, n, lp1);
        qmfWrev(hm113, n, hp1);
        copy(hm113, n, lp2);
        qmfEven(h1 + 2, n, hp2);
        return n;
    }

    if (name == "rbior1.5"sv) {
        copyReverse(h1, n, lp1);
        qmfWrev(hm115, n, hp1);
        copy(hm115, n, lp2);
        qmfEven(h1, n, hp2);
        return n;
    }

    if (name == "rbior2.2"sv) {
        copyReverse(h2 + 6, n, lp1);
        qmfWrev(hm222, n, hp1);
        copy(hm222, n, lp2);
        qmfEven(h2 + 6, n, hp2);
        return n;
    }

    if (name == "rbior2.4"sv) {
        copyReverse(h2 + 4, n, lp1);
        qmfWrev(hm224, n, hp1);
        copy(hm224, n, lp2);
        qmfEven(h2 + 4, n, hp2);
        return n;
    }

    if (name == "rbior2.6"sv) {
        copyReverse(h2 + 2, n, lp1);
        qmfWrev(hm226, n, hp1);
        copy(hm226, n, lp2);
        qmfEven(h2 + 2, n, hp2);
        return n;
    }

    if (name == "rbior2.8"sv) {
        copyReverse(h2, n, lp1);
        qmfWrev(hm228, n, hp1);
        copy(hm228, n, lp2);
        qmfEven(h2, n, hp2);
        return n;
    }

    if (name == "rbior3.1"sv) {
        copyReverse(h3 + 8, n, lp1);
        qmfWrev(hm331, n, hp1);
        copy(hm331, n, lp2);
        qmfEven(h3 + 8, n, hp2);
        return n;
    }

    if (name == "rbior3.3"sv) {
        copyReverse(h3 + 6, n, lp1);
        qmfWrev(hm333, n, hp1);
        copy(hm333, n, lp2);
        qmfEven(h3 + 6, n, hp2);
        return n;
    }

    if (name == "rbior3.5"sv) {
        copyReverse(h3 + 4, n, lp1);
        qmfWrev(hm335, n, hp1);
        copy(hm335, n, lp2);
        qmfEven(h3 + 4, n, hp2);
        return n;
    }

    if (name == "rbior3.7"sv) {
        copyReverse(h3 + 2, n, lp1);
        qmfWrev(hm337, n, hp1);
        copy(hm337, n, lp2);
        qmfEven(h3 + 2, n, hp2);
        return n;
    }

    if (name == "rbior3.9"sv) {
        copyReverse(h3, n, lp1);
        qmfWrev(hm339, n, hp1);
        copy(hm339, n, lp2);
        qmfEven(h3, n, hp2);
        return n;
    }

    if (name == "rbior4.4"sv) {
        copyReverse(h4, n, lp1);
        qmfWrev(hm444, n, hp1);
        copy(hm444, n, lp2);
        qmfEven(h4, n, hp2);
        return n;
    }

    if (name == "rbior5.5"sv) {
        copyReverse(h5, n, lp1);
        qmfWrev(hm555, n, hp1);
        copy(hm555, n, lp2);
        qmfEven(h5, n, hp2);
        return n;
    }

    if (name == "rbior6.8"sv) {
        copyReverse(h6, n, lp1);
        qmfWrev(hm668, n, hp1);
        copy(hm668, n, lp2);
        qmfEven(h6, n, hp2);
        return n;
    }

    if (name == "coif1"sv) {
        auto coeffTemp = std::make_unique<double[]>(n);

        copy(coif1, n, coeffTemp.get());
        for (i = 0; i < n; ++i) {
            coeffTemp[i] *= M_SQRT2;
        }

        copyReverse(coeffTemp.get(), n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        copy(coeffTemp.get(), n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }

    if (name == "coif2"sv) {
        auto coeffTemp = std::make_unique<double[]>(n);

        copy(coif2, n, coeffTemp.get());
        for (i = 0; i < n; ++i) {
            coeffTemp[i] *= M_SQRT2;
        }

        copyReverse(coeffTemp.get(), n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        copy(coeffTemp.get(), n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }

    if (name == "coif3"sv) {
        auto coeffTemp = std::make_unique<double[]>(n);

        copy(coif3, n, coeffTemp.get());
        for (i = 0; i < n; ++i) {
            coeffTemp[i] *= M_SQRT2;
        }

        copyReverse(coeffTemp.get(), n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        copy(coeffTemp.get(), n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }

    if (name == "coif4"sv) {
        auto coeffTemp = std::make_unique<double[]>(n);

        copy(coif4, n, coeffTemp.get());
        for (i = 0; i < n; ++i) {
            coeffTemp[i] *= M_SQRT2;
        }

        copyReverse(coeffTemp.get(), n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        copy(coeffTemp.get(), n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }

    if (name == "coif5"sv) {
        auto coeffTemp = std::make_unique<double[]>(n);

        copy(coif5, n, coeffTemp.get());
        for (i = 0; i < n; ++i) {
            coeffTemp[i] *= M_SQRT2;
        }

        copyReverse(coeffTemp.get(), n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        copy(coeffTemp.get(), n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }

    if (name == "coif6"sv) {
        auto coeffTemp = std::make_unique<double[]>(n);

        copy(coif6, n, coeffTemp.get());
        for (i = 0; i < n; ++i) {
            coeffTemp[i] *= M_SQRT2;
        }

        copyReverse(coeffTemp.get(), n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        copy(coeffTemp.get(), n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }

    if (name == "coif7"sv) {
        auto coeffTemp = std::make_unique<double[]>(n);

        copy(coif7, n, coeffTemp.get());
        for (i = 0; i < n; ++i) {
            coeffTemp[i] *= M_SQRT2;
        }

        copyReverse(coeffTemp.get(), n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        copy(coeffTemp.get(), n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }

    if (name == "coif8"sv) {
        auto coeffTemp = std::make_unique<double[]>(n);

        copy(coif8, n, coeffTemp.get());
        for (i = 0; i < n; ++i) {
            coeffTemp[i] *= M_SQRT2;
        }

        copyReverse(coeffTemp.get(), n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        copy(coeffTemp.get(), n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }

    if (name == "coif9"sv) {
        auto coeffTemp = std::make_unique<double[]>(n);

        copy(coif9, n, coeffTemp.get());
        for (i = 0; i < n; ++i) {
            coeffTemp[i] *= M_SQRT2;
        }

        copyReverse(coeffTemp.get(), n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        copy(coeffTemp.get(), n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }

    if (name == "coif10"sv) {
        auto coeffTemp = std::make_unique<double[]>(n);

        copy(coif10, n, coeffTemp.get());
        for (i = 0; i < n; ++i) {
            coeffTemp[i] *= M_SQRT2;
        }

        copyReverse(coeffTemp.get(), n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        copy(coeffTemp.get(), n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }
    if (name == "coif11"sv) {
        auto coeffTemp = std::make_unique<double[]>(n);

        copy(coif11, n, coeffTemp.get());
        for (i = 0; i < n; ++i) {
            coeffTemp[i] *= M_SQRT2;
        }

        copyReverse(coeffTemp.get(), n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        copy(coeffTemp.get(), n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }
    if (name == "coif12"sv) {
        auto coeffTemp = std::make_unique<double[]>(n);

        copy(coif12, n, coeffTemp.get());
        for (i = 0; i < n; ++i) {
            coeffTemp[i] *= M_SQRT2;
        }

        copyReverse(coeffTemp.get(), n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        copy(coeffTemp.get(), n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }
    if (name == "coif13"sv) {
        auto coeffTemp = std::make_unique<double[]>(n);

        copy(coif13, n, coeffTemp.get());
        for (i = 0; i < n; ++i) {
            coeffTemp[i] *= M_SQRT2;
        }

        copyReverse(coeffTemp.get(), n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        copy(coeffTemp.get(), n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }
    if (name == "coif14"sv) {
        auto coeffTemp = std::make_unique<double[]>(n);

        copy(coif14, n, coeffTemp.get());
        for (i = 0; i < n; ++i) {
            coeffTemp[i] *= M_SQRT2;
        }

        copyReverse(coeffTemp.get(), n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        copy(coeffTemp.get(), n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }
    if (name == "coif15"sv) {
        auto coeffTemp = std::make_unique<double[]>(n);

        copy(coif15, n, coeffTemp.get());
        for (i = 0; i < n; ++i) {
            coeffTemp[i] *= M_SQRT2;
        }

        copyReverse(coeffTemp.get(), n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        copy(coeffTemp.get(), n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }
    if (name == "coif16"sv) {
        auto coeffTemp = std::make_unique<double[]>(n);

        copy(coif16, n, coeffTemp.get());
        for (i = 0; i < n; ++i) {
            coeffTemp[i] *= M_SQRT2;
        }

        copyReverse(coeffTemp.get(), n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        copy(coeffTemp.get(), n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }
    if (name == "coif17"sv) {
        auto coeffTemp = std::make_unique<double[]>(n);

        copy(coif17, n, coeffTemp.get());
        for (i = 0; i < n; ++i) {
            coeffTemp[i] *= M_SQRT2;
        }

        copyReverse(coeffTemp.get(), n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        copy(coeffTemp.get(), n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }
    if (name == "sym2"sv) {
        copyReverse(sym2, n, lp1);
        qmfWrev(sym2, n, hp1);
        copy(sym2, n, lp2);
        qmfEven(sym2, n, hp2);
        return n;
    }

    if (name == "sym3"sv) {
        copyReverse(sym3, n, lp1);
        qmfWrev(sym3, n, hp1);
        copy(sym3, n, lp2);
        qmfEven(sym3, n, hp2);
        return n;
    }

    if (name == "sym4"sv) {
        copyReverse(sym4, n, lp1);
        qmfWrev(sym4, n, hp1);
        copy(sym4, n, lp2);
        qmfEven(sym4, n, hp2);
        return n;
    }

    if (name == "sym5"sv) {
        copyReverse(sym5, n, lp1);
        qmfWrev(sym5, n, hp1);
        copy(sym5, n, lp2);
        qmfEven(sym5, n, hp2);
        return n;
    }

    if (name == "sym6"sv) {
        copyReverse(sym6, n, lp1);
        qmfWrev(sym6, n, hp1);
        copy(sym6, n, lp2);
        qmfEven(sym6, n, hp2);
        return n;
    }

    if (name == "sym7"sv) {
        copyReverse(sym7, n, lp1);
        qmfWrev(sym7, n, hp1);
        copy(sym7, n, lp2);
        qmfEven(sym7, n, hp2);
        return n;
    }

    if (name == "sym8"sv) {
        copyReverse(sym8, n, lp1);
        qmfWrev(sym8, n, hp1);
        copy(sym8, n, lp2);
        qmfEven(sym8, n, hp2);
        return n;
    }

    if (name == "sym9"sv) {
        copyReverse(sym9, n, lp1);
        qmfWrev(sym9, n, hp1);
        copy(sym9, n, lp2);
        qmfEven(sym9, n, hp2);
        return n;
    }

    if (name == "sym10"sv) {
        copyReverse(sym10, n, lp1);
        qmfWrev(sym10, n, hp1);
        copy(sym10, n, lp2);
        qmfEven(sym10, n, hp2);
        return n;
    }
    if (name == "sym11"sv) {
        copyReverse(sym11, n, lp1);
        qmfWrev(sym11, n, hp1);
        copy(sym11, n, lp2);
        qmfEven(sym11, n, hp2);
        return n;
    }
    if (name == "sym12"sv) {
        copyReverse(sym12, n, lp1);
        qmfWrev(sym12, n, hp1);
        copy(sym12, n, lp2);
        qmfEven(sym12, n, hp2);
        return n;
    }
    if (name == "sym13"sv) {
        copyReverse(sym13, n, lp1);
        qmfWrev(sym13, n, hp1);
        copy(sym13, n, lp2);
        qmfEven(sym13, n, hp2);
        return n;
    }
    if (name == "sym14"sv) {
        copyReverse(sym14, n, lp1);
        qmfWrev(sym14, n, hp1);
        copy(sym14, n, lp2);
        qmfEven(sym14, n, hp2);
        return n;
    }
    if (name == "sym15"sv) {
        copyReverse(sym15, n, lp1);
        qmfWrev(sym15, n, hp1);
        copy(sym15, n, lp2);
        qmfEven(sym15, n, hp2);
        return n;
    }
    if (name == "sym16"sv) {
        copyReverse(sym16, n, lp1);
        qmfWrev(sym16, n, hp1);
        copy(sym16, n, lp2);
        qmfEven(sym16, n, hp2);
        return n;
    }
    if (name == "sym17"sv) {
        copyReverse(sym17, n, lp1);
        qmfWrev(sym17, n, hp1);
        copy(sym17, n, lp2);
        qmfEven(sym17, n, hp2);
        return n;
    }
    if (name == "sym18"sv) {
        copyReverse(sym18, n, lp1);
        qmfWrev(sym18, n, hp1);
        copy(sym18, n, lp2);
        qmfEven(sym18, n, hp2);
        return n;
    }
    if (name == "sym19"sv) {
        copyReverse(sym19, n, lp1);
        qmfWrev(sym19, n, hp1);
        copy(sym19, n, lp2);
        qmfEven(sym19, n, hp2);
        return n;
    }
    if (name == "sym20"sv) {
        copyReverse(sym20, n, lp1);
        qmfWrev(sym20, n, hp1);
        copy(sym20, n, lp2);
        qmfEven(sym20, n, hp2);
        return n;
    }
    if (name == "meyer"sv) {
        copyReverse(meyer, n, lp1);
        qmfWrev(meyer, n, hp1);
        copy(meyer, n, lp2);
        qmfEven(meyer, n, hp2);
        return n;
    }

    printf("\n Filter Not in Database \n");
    return -1;
}

}

Wavelet::Wavelet(char const* name)
    : name_ { name }
    , size_ { static_cast<std::size_t>(::waveletFilterLength(name)) }
    , params_ { std::make_unique<double[]>(4 * size_) }
    , lpd_ { &params_[0], size_ }
    , hpd_ { &params_[size_], size_ }
    , lpr_ { &params_[2 * size_], size_ }
    , hpr_ { &params_[3 * size_], size_ }
{
    auto* p = params_.get();
    if (name != nullptr) {
        waveletFilterCoefficients(name, p, p + size_, p + 2 * size_, p + 3 * size_);
    }
}

auto summary(Wavelet const& obj) -> void
{
    auto const n = obj.size();
    printf("\n");
    printf("Wavelet Name: %s \n", obj.name().c_str());
    printf("\n");
    printf("Wavelet Filters \n");
    printf("lpd: [");
    for (auto i = 0; i < n - 1; ++i) {
        printf("%g,", obj.lpd()[i]);
    }
    printf("%g", obj.lpd()[n - 1]);
    printf("] \n");
    printf("hpd: [");
    for (auto i = 0; i < n - 1; ++i) {
        printf("%g,", obj.hpd()[i]);
    }
    printf("%g", obj.hpd()[n - 1]);
    printf("] \n");
    printf("lpr: [");
    for (auto i = 0; i < n - 1; ++i) {
        printf("%g,", obj.lpr()[i]);
    }
    printf("%g", obj.lpr()[n - 1]);
    printf("] \n");
    printf("hpr: [");
    for (auto i = 0; i < n - 1; ++i) {
        printf("%g,", obj.hpr()[i]);
    }
    printf("%g", obj.hpr()[n - 1]);
    printf("] \n");
}