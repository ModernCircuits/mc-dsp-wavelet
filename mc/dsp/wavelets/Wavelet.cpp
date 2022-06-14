#include "Wavelet.hpp"

#include "mc/dsp/convolution/FFTConvolver.hpp"
#include "mc/dsp/wavelets/filters/coif.hpp"
#include "mc/dsp/wavelets/filters/daubechies.hpp"
#include "mc/dsp/wavelets/filters/h.hpp"
#include "mc/dsp/wavelets/filters/meyer.hpp"
#include "mc/dsp/wavelets/filters/sym.hpp"

#include "mc/cmath.hpp"
#include "mc/cstring.hpp"
#include "mc/format.hpp"
#include "mc/memory.hpp"
#include "mc/numbers.hpp"
#include "mc/string_view.hpp"
#include "mc/utility.hpp"

namespace mc::dsp
{

namespace
{
auto copy(float const* in, int n, float* out)
{
    int count = 0;
    for (count = 0; count < n; count++) { out[count] = in[count]; }
}

auto copyReverse(float const* in, int n, float* out)
{
    int count = 0;
    for (count = 0; count < n; count++) { out[count] = in[n - count - 1]; }
}

auto qmfEven(float const* in, int n, float* out)
{
    int count = 0;
    for (count = 0; count < n; count++)
    {
        out[count] = in[n - count - 1];
        if (count % 2 != 0) { out[count] = -1 * out[count]; }
    }
}

auto qmfWrev(float const* in, int n, float* out)
{
    auto sigOutTemp = std::make_unique<float[]>(n);

    qmfEven(in, n, sigOutTemp.get());
    copyReverse(sigOutTemp.get(), n, out);
}

auto waveletFilterLength(char const* name) -> int
{
    int len = strlen(name);
    int i   = 0;
    int n   = 0;
    if (name == string_view{"haar"} || name == string_view{"db1"}) { return 2; }
    if (len > 2 && strstr(name, "db") != nullptr)
    {
        auto newStr = std::make_unique<char[]>((len - 2 + 1));
        for (i = 2; i < len + 1; i++) { newStr[i - 2] = name[i]; }

        n = atoi(newStr.get());
        if (n > 38)
        {
            fmt::printf("\n Filter Not in Database \n");
            return -1;
        }

        return n * 2;
    }
    if (name == string_view{"bior1.1"}) { return 2; }

    if (name == string_view{"bior1.3"}) { return 6; }

    if (name == string_view{"bior1.5"}) { return 10; }

    if (name == string_view{"bior2.2"}) { return 6; }

    if (name == string_view{"bior2.4"}) { return 10; }

    if (name == string_view{"bior2.6"}) { return 14; }
    if (name == string_view{"bior2.8"}) { return 18; }

    if (name == string_view{"bior3.1"}) { return 4; }
    if (name == string_view{"bior3.3"}) { return 8; }
    if (name == string_view{"bior3.5"}) { return 12; }

    if (name == string_view{"bior3.7"}) { return 16; }
    if (name == string_view{"bior3.9"}) { return 20; }
    if (name == string_view{"bior4.4"}) { return 10; }
    if (name == string_view{"bior5.5"}) { return 12; }
    if (name == string_view{"bior6.8"}) { return 18; }
    if (name == string_view{"rbior1.1"}) { return 2; }

    if (name == string_view{"rbior1.3"}) { return 6; }

    if (name == string_view{"rbior1.5"}) { return 10; }

    if (name == string_view{"rbior2.2"}) { return 6; }

    if (name == string_view{"rbior2.4"}) { return 10; }

    if (name == string_view{"rbior2.6"}) { return 14; }
    if (name == string_view{"rbior2.8"}) { return 18; }

    if (name == string_view{"rbior3.1"}) { return 4; }
    if (name == string_view{"rbior3.3"}) { return 8; }
    if (name == string_view{"rbior3.5"}) { return 12; }

    if (name == string_view{"rbior3.7"}) { return 16; }
    if (name == string_view{"rbior3.9"}) { return 20; }
    if (name == string_view{"rbior4.4"}) { return 10; }
    if (name == string_view{"rbior5.5"}) { return 12; }
    if (name == string_view{"rbior6.8"}) { return 18; }
    if (len > 4 && strstr(name, "coif") != nullptr)
    {
        auto newStr = std::make_unique<char[]>((len - 4 + 1));
        for (i = 4; i < len + 1; i++) { newStr[i - 4] = name[i]; }

        n = atoi(newStr.get());
        if (n > 17)
        {
            fmt::printf("\n Filter Not in Database \n");
            return -1;
        }

        return n * 6;
    }
    if (len > 3 && strstr(name, "sym") != nullptr)
    {
        auto newStr = std::make_unique<char[]>((len - 3 + 1));
        for (i = 3; i < len + 1; i++) { newStr[i - 3] = name[i]; }

        n = atoi(newStr.get());
        if (n > 20 || n < 2)
        {
            fmt::printf("\n Filter Not in Database \n");
            return -1;
        }

        return n * 2;
    }
    if (strcmp(name, "meyer") == 0) { return 102; }
    fmt::printf("\n Filter Not in Database \n");
    return -1;
}

auto waveletFilterCoefficients(char const* name, float* lp1, float* hp1, float* lp2, float* hp2) -> int
{
    int i = 0;
    int n = waveletFilterLength(name);
    if ((strcmp(name, "haar") == 0) || (strcmp(name, "db1") == 0))
    {
        copyReverse(daubechies1, n, lp1);
        qmfWrev(daubechies1, n, hp1);
        copy(daubechies1, n, lp2);
        qmfEven(daubechies1, n, hp2);

        return n;
    }
    if (name == string_view{"db2"})
    {
        copyReverse(daubechies2, n, lp1);
        qmfWrev(daubechies2, n, hp1);
        copy(daubechies2, n, lp2);
        qmfEven(daubechies2, n, hp2);

        return n;
    }
    if (name == string_view{"db3"})
    {
        copyReverse(daubechies3, n, lp1);
        qmfWrev(daubechies3, n, hp1);
        copy(daubechies3, n, lp2);
        qmfEven(daubechies3, n, hp2);

        return n;
    }
    if (name == string_view{"db4"})
    {
        copyReverse(daubechies4, n, lp1);
        qmfWrev(daubechies4, n, hp1);
        copy(daubechies4, n, lp2);
        qmfEven(daubechies4, n, hp2);

        return n;
    }
    if (name == string_view{"db5"})
    {
        copyReverse(daubechies5, n, lp1);
        qmfWrev(daubechies5, n, hp1);
        copy(daubechies5, n, lp2);
        qmfEven(daubechies5, n, hp2);

        return n;
    }
    if (name == string_view{"db6"})
    {
        copyReverse(daubechies6, n, lp1);
        qmfWrev(daubechies6, n, hp1);
        copy(daubechies6, n, lp2);
        qmfEven(daubechies6, n, hp2);

        return n;
    }
    if (name == string_view{"db7"})
    {
        copyReverse(daubechies7, n, lp1);
        qmfWrev(daubechies7, n, hp1);
        copy(daubechies7, n, lp2);
        qmfEven(daubechies7, n, hp2);

        return n;
    }
    if (name == string_view{"db8"})
    {
        copyReverse(daubechies8, n, lp1);
        qmfWrev(daubechies8, n, hp1);
        copy(daubechies8, n, lp2);
        qmfEven(daubechies8, n, hp2);

        return n;
    }
    if (name == string_view{"db9"})
    {
        copyReverse(daubechies9, n, lp1);
        qmfWrev(daubechies9, n, hp1);
        copy(daubechies9, n, lp2);
        qmfEven(daubechies9, n, hp2);

        return n;
    }

    if (name == string_view{"db10"})
    {
        copyReverse(daubechies10, n, lp1);
        qmfWrev(daubechies10, n, hp1);
        copy(daubechies10, n, lp2);
        qmfEven(daubechies10, n, hp2);

        return n;
    }

    if (name == string_view{"db11"})
    {
        copyReverse(daubechies11, n, lp1);
        qmfWrev(daubechies11, n, hp1);
        copy(daubechies11, n, lp2);
        qmfEven(daubechies11, n, hp2);

        return n;
    }
    if (name == string_view{"db12"})
    {
        copyReverse(daubechies12, n, lp1);
        qmfWrev(daubechies12, n, hp1);
        copy(daubechies12, n, lp2);
        qmfEven(daubechies12, n, hp2);

        return n;
    }

    if (name == string_view{"db13"})
    {
        copyReverse(daubechies13, n, lp1);
        qmfWrev(daubechies13, n, hp1);
        copy(daubechies13, n, lp2);
        qmfEven(daubechies13, n, hp2);

        return n;
    }

    if (name == string_view{"db14"})
    {
        copyReverse(daubechies14, n, lp1);
        qmfWrev(daubechies14, n, hp1);
        copy(daubechies14, n, lp2);
        qmfEven(daubechies14, n, hp2);

        return n;
    }

    if (name == string_view{"db15"})
    {
        copyReverse(daubechies15, n, lp1);
        qmfWrev(daubechies15, n, hp1);
        copy(daubechies15, n, lp2);
        qmfEven(daubechies15, n, hp2);

        return n;
    }

    if (name == string_view{"db16"})
    {
        copyReverse(daubechies16, n, lp1);
        qmfWrev(daubechies16, n, hp1);
        copy(daubechies16, n, lp2);
        qmfEven(daubechies16, n, hp2);

        return n;
    }

    if (name == string_view{"db17"})
    {
        copyReverse(daubechies17, n, lp1);
        qmfWrev(daubechies17, n, hp1);
        copy(daubechies17, n, lp2);
        qmfEven(daubechies17, n, hp2);

        return n;
    }

    if (name == string_view{"db18"})
    {
        copyReverse(daubechies18, n, lp1);
        qmfWrev(daubechies18, n, hp1);
        copy(daubechies18, n, lp2);
        qmfEven(daubechies18, n, hp2);

        return n;
    }

    if (name == string_view{"db19"})
    {
        copyReverse(daubechies19, n, lp1);
        qmfWrev(daubechies19, n, hp1);
        copy(daubechies19, n, lp2);
        qmfEven(daubechies19, n, hp2);

        return n;
    }

    if (name == string_view{"db20"})
    {
        copyReverse(daubechies20, n, lp1);
        qmfWrev(daubechies20, n, hp1);
        copy(daubechies20, n, lp2);
        qmfEven(daubechies20, n, hp2);

        return n;
    }

    if (name == string_view{"db21"})
    {
        copyReverse(daubechies21, n, lp1);
        qmfWrev(daubechies21, n, hp1);
        copy(daubechies21, n, lp2);
        qmfEven(daubechies21, n, hp2);

        return n;
    }

    if (name == string_view{"db22"})
    {
        copyReverse(daubechies22, n, lp1);
        qmfWrev(daubechies22, n, hp1);
        copy(daubechies22, n, lp2);
        qmfEven(daubechies22, n, hp2);

        return n;
    }

    if (name == string_view{"db23"})
    {
        copyReverse(daubechies23, n, lp1);
        qmfWrev(daubechies23, n, hp1);
        copy(daubechies23, n, lp2);
        qmfEven(daubechies23, n, hp2);

        return n;
    }

    if (name == string_view{"db24"})
    {
        copyReverse(daubechies24, n, lp1);
        qmfWrev(daubechies24, n, hp1);
        copy(daubechies24, n, lp2);
        qmfEven(daubechies24, n, hp2);

        return n;
    }

    if (name == string_view{"db25"})
    {
        copyReverse(daubechies25, n, lp1);
        qmfWrev(daubechies25, n, hp1);
        copy(daubechies25, n, lp2);
        qmfEven(daubechies25, n, hp2);

        return n;
    }

    if (name == string_view{"db26"})
    {
        copyReverse(daubechies26, n, lp1);
        qmfWrev(daubechies26, n, hp1);
        copy(daubechies26, n, lp2);
        qmfEven(daubechies26, n, hp2);
        return n;
    }

    if (name == string_view{"db27"})
    {
        copyReverse(daubechies27, n, lp1);
        qmfWrev(daubechies27, n, hp1);
        copy(daubechies27, n, lp2);
        qmfEven(daubechies27, n, hp2);

        return n;
    }

    if (name == string_view{"db28"})
    {
        copyReverse(daubechies28, n, lp1);
        qmfWrev(daubechies28, n, hp1);
        copy(daubechies28, n, lp2);
        qmfEven(daubechies28, n, hp2);

        return n;
    }

    if (name == string_view{"db29"})
    {
        copyReverse(daubechies29, n, lp1);
        qmfWrev(daubechies29, n, hp1);
        copy(daubechies29, n, lp2);
        qmfEven(daubechies29, n, hp2);

        return n;
    }

    if (name == string_view{"db30"})
    {
        copyReverse(daubechies30, n, lp1);
        qmfWrev(daubechies30, n, hp1);
        copy(daubechies30, n, lp2);
        qmfEven(daubechies30, n, hp2);

        return n;
    }

    if (name == string_view{"db31"})
    {
        copyReverse(daubechies31, n, lp1);
        qmfWrev(daubechies31, n, hp1);
        copy(daubechies31, n, lp2);
        qmfEven(daubechies31, n, hp2);

        return n;
    }

    if (name == string_view{"db32"})
    {
        copyReverse(daubechies32, n, lp1);
        qmfWrev(daubechies32, n, hp1);
        copy(daubechies32, n, lp2);
        qmfEven(daubechies32, n, hp2);

        return n;
    }

    if (name == string_view{"db33"})
    {
        copyReverse(daubechies33, n, lp1);
        qmfWrev(daubechies33, n, hp1);
        copy(daubechies33, n, lp2);
        qmfEven(daubechies33, n, hp2);

        return n;
    }

    if (name == string_view{"db34"})
    {
        copyReverse(daubechies34, n, lp1);
        qmfWrev(daubechies34, n, hp1);
        copy(daubechies34, n, lp2);
        qmfEven(daubechies34, n, hp2);

        return n;
    }

    if (name == string_view{"db35"})
    {
        copyReverse(daubechies35, n, lp1);
        qmfWrev(daubechies35, n, hp1);
        copy(daubechies35, n, lp2);
        qmfEven(daubechies35, n, hp2);

        return n;
    }

    if (name == string_view{"db36"})
    {
        copyReverse(daubechies36, n, lp1);
        qmfWrev(daubechies36, n, hp1);
        copy(daubechies36, n, lp2);
        qmfEven(daubechies36, n, hp2);

        return n;
    }

    if (name == string_view{"db37"})
    {
        copyReverse(daubechies37, n, lp1);
        qmfWrev(daubechies37, n, hp1);
        copy(daubechies37, n, lp2);
        qmfEven(daubechies37, n, hp2);

        return n;
    }

    if (name == string_view{"db38"})
    {
        copyReverse(daubechies38, n, lp1);
        qmfWrev(daubechies38, n, hp1);
        copy(daubechies38, n, lp2);
        qmfEven(daubechies38, n, hp2);

        return n;
    }

    if (name == string_view{"bior1.1"})
    {
        copyReverse(hm111, n, lp1);
        qmfWrev(h1 + 4, n, hp1);
        copy(h1 + 4, n, lp2);
        qmfEven(hm111, n, hp2);
        return n;
    }

    if (name == string_view{"bior1.3"})
    {
        copyReverse(hm113, n, lp1);
        qmfWrev(h1 + 2, n, hp1);
        copy(h1 + 2, n, lp2);
        qmfEven(hm113, n, hp2);
        return n;
    }

    if (name == string_view{"bior1.5"})
    {
        copyReverse(hm115, n, lp1);
        qmfWrev(h1, n, hp1);
        copy(h1, n, lp2);
        qmfEven(hm115, n, hp2);
        return n;
    }

    if (name == string_view{"bior2.2"})
    {
        copyReverse(hm222, n, lp1);
        qmfWrev(h2 + 6, n, hp1);
        copy(h2 + 6, n, lp2);
        qmfEven(hm222, n, hp2);
        return n;
    }

    if (name == string_view{"bior2.4"})
    {
        copyReverse(hm224, n, lp1);
        qmfWrev(h2 + 4, n, hp1);
        copy(h2 + 4, n, lp2);
        qmfEven(hm224, n, hp2);
        return n;
    }

    if (name == string_view{"bior2.6"})
    {
        copyReverse(hm226, n, lp1);
        qmfWrev(h2 + 2, n, hp1);
        copy(h2 + 2, n, lp2);
        qmfEven(hm226, n, hp2);
        return n;
    }

    if (name == string_view{"bior2.8"})
    {
        copyReverse(hm228, n, lp1);
        qmfWrev(h2, n, hp1);
        copy(h2, n, lp2);
        qmfEven(hm228, n, hp2);
        return n;
    }

    if (name == string_view{"bior3.1"})
    {
        copyReverse(hm331, n, lp1);
        qmfWrev(h3 + 8, n, hp1);
        copy(h3 + 8, n, lp2);
        qmfEven(hm331, n, hp2);
        return n;
    }

    if (name == string_view{"bior3.3"})
    {
        copyReverse(hm333, n, lp1);
        qmfWrev(h3 + 6, n, hp1);
        copy(h3 + 6, n, lp2);
        qmfEven(hm333, n, hp2);
        return n;
    }

    if (name == string_view{"bior3.5"})
    {
        copyReverse(hm335, n, lp1);
        qmfWrev(h3 + 4, n, hp1);
        copy(h3 + 4, n, lp2);
        qmfEven(hm335, n, hp2);
        return n;
    }

    if (name == string_view{"bior3.7"})
    {
        copyReverse(hm337, n, lp1);
        qmfWrev(h3 + 2, n, hp1);
        copy(h3 + 2, n, lp2);
        qmfEven(hm337, n, hp2);
        return n;
    }

    if (name == string_view{"bior3.9"})
    {
        copyReverse(hm339, n, lp1);
        qmfWrev(h3, n, hp1);
        copy(h3, n, lp2);
        qmfEven(hm339, n, hp2);
        return n;
    }

    if (name == string_view{"bior4.4"})
    {
        copyReverse(hm444, n, lp1);
        qmfWrev(h4, n, hp1);
        copy(h4, n, lp2);
        qmfEven(hm444, n, hp2);
        return n;
    }

    if (name == string_view{"bior5.5"})
    {
        copyReverse(hm555, n, lp1);
        qmfWrev(h5, n, hp1);
        copy(h5, n, lp2);
        qmfEven(hm555, n, hp2);
        return n;
    }

    if (name == string_view{"bior6.8"})
    {
        copyReverse(hm668, n, lp1);
        qmfWrev(h6, n, hp1);
        copy(h6, n, lp2);
        qmfEven(hm668, n, hp2);
        return n;
    }

    if (name == string_view{"rbior1.1"})
    {
        copyReverse(h1 + 4, n, lp1);
        qmfWrev(hm111, n, hp1);
        copy(hm111, n, lp2);
        qmfEven(h1 + 4, n, hp2);
        return n;
    }

    if (name == string_view{"rbior1.3"})
    {
        copyReverse(h1 + 2, n, lp1);
        qmfWrev(hm113, n, hp1);
        copy(hm113, n, lp2);
        qmfEven(h1 + 2, n, hp2);
        return n;
    }

    if (name == string_view{"rbior1.5"})
    {
        copyReverse(h1, n, lp1);
        qmfWrev(hm115, n, hp1);
        copy(hm115, n, lp2);
        qmfEven(h1, n, hp2);
        return n;
    }

    if (name == string_view{"rbior2.2"})
    {
        copyReverse(h2 + 6, n, lp1);
        qmfWrev(hm222, n, hp1);
        copy(hm222, n, lp2);
        qmfEven(h2 + 6, n, hp2);
        return n;
    }

    if (name == string_view{"rbior2.4"})
    {
        copyReverse(h2 + 4, n, lp1);
        qmfWrev(hm224, n, hp1);
        copy(hm224, n, lp2);
        qmfEven(h2 + 4, n, hp2);
        return n;
    }

    if (name == string_view{"rbior2.6"})
    {
        copyReverse(h2 + 2, n, lp1);
        qmfWrev(hm226, n, hp1);
        copy(hm226, n, lp2);
        qmfEven(h2 + 2, n, hp2);
        return n;
    }

    if (name == string_view{"rbior2.8"})
    {
        copyReverse(h2, n, lp1);
        qmfWrev(hm228, n, hp1);
        copy(hm228, n, lp2);
        qmfEven(h2, n, hp2);
        return n;
    }

    if (name == string_view{"rbior3.1"})
    {
        copyReverse(h3 + 8, n, lp1);
        qmfWrev(hm331, n, hp1);
        copy(hm331, n, lp2);
        qmfEven(h3 + 8, n, hp2);
        return n;
    }

    if (name == string_view{"rbior3.3"})
    {
        copyReverse(h3 + 6, n, lp1);
        qmfWrev(hm333, n, hp1);
        copy(hm333, n, lp2);
        qmfEven(h3 + 6, n, hp2);
        return n;
    }

    if (name == string_view{"rbior3.5"})
    {
        copyReverse(h3 + 4, n, lp1);
        qmfWrev(hm335, n, hp1);
        copy(hm335, n, lp2);
        qmfEven(h3 + 4, n, hp2);
        return n;
    }

    if (name == string_view{"rbior3.7"})
    {
        copyReverse(h3 + 2, n, lp1);
        qmfWrev(hm337, n, hp1);
        copy(hm337, n, lp2);
        qmfEven(h3 + 2, n, hp2);
        return n;
    }

    if (name == string_view{"rbior3.9"})
    {
        copyReverse(h3, n, lp1);
        qmfWrev(hm339, n, hp1);
        copy(hm339, n, lp2);
        qmfEven(h3, n, hp2);
        return n;
    }

    if (name == string_view{"rbior4.4"})
    {
        copyReverse(h4, n, lp1);
        qmfWrev(hm444, n, hp1);
        copy(hm444, n, lp2);
        qmfEven(h4, n, hp2);
        return n;
    }

    if (name == string_view{"rbior5.5"})
    {
        copyReverse(h5, n, lp1);
        qmfWrev(hm555, n, hp1);
        copy(hm555, n, lp2);
        qmfEven(h5, n, hp2);
        return n;
    }

    if (name == string_view{"rbior6.8"})
    {
        copyReverse(h6, n, lp1);
        qmfWrev(hm668, n, hp1);
        copy(hm668, n, lp2);
        qmfEven(h6, n, hp2);
        return n;
    }

    if (name == string_view{"coif1"})
    {
        auto coeffTemp = std::make_unique<float[]>(n);

        copy(coif1, n, coeffTemp.get());
        for (i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        copyReverse(coeffTemp.get(), n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        copy(coeffTemp.get(), n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }

    if (name == string_view{"coif2"})
    {
        auto coeffTemp = std::make_unique<float[]>(n);

        copy(coif2, n, coeffTemp.get());
        for (i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        copyReverse(coeffTemp.get(), n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        copy(coeffTemp.get(), n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }

    if (name == string_view{"coif3"})
    {
        auto coeffTemp = std::make_unique<float[]>(n);

        copy(coif3, n, coeffTemp.get());
        for (i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        copyReverse(coeffTemp.get(), n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        copy(coeffTemp.get(), n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }

    if (name == string_view{"coif4"})
    {
        auto coeffTemp = std::make_unique<float[]>(n);

        copy(coif4, n, coeffTemp.get());
        for (i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        copyReverse(coeffTemp.get(), n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        copy(coeffTemp.get(), n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }

    if (name == string_view{"coif5"})
    {
        auto coeffTemp = std::make_unique<float[]>(n);

        copy(coif5, n, coeffTemp.get());
        for (i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        copyReverse(coeffTemp.get(), n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        copy(coeffTemp.get(), n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }

    if (name == string_view{"coif6"})
    {
        auto coeffTemp = std::make_unique<float[]>(n);

        copy(coif6, n, coeffTemp.get());
        for (i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        copyReverse(coeffTemp.get(), n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        copy(coeffTemp.get(), n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }

    if (name == string_view{"coif7"})
    {
        auto coeffTemp = std::make_unique<float[]>(n);

        copy(coif7, n, coeffTemp.get());
        for (i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        copyReverse(coeffTemp.get(), n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        copy(coeffTemp.get(), n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }

    if (name == string_view{"coif8"})
    {
        auto coeffTemp = std::make_unique<float[]>(n);

        copy(coif8, n, coeffTemp.get());
        for (i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        copyReverse(coeffTemp.get(), n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        copy(coeffTemp.get(), n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }

    if (name == string_view{"coif9"})
    {
        auto coeffTemp = std::make_unique<float[]>(n);

        copy(coif9, n, coeffTemp.get());
        for (i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        copyReverse(coeffTemp.get(), n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        copy(coeffTemp.get(), n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }

    if (name == string_view{"coif10"})
    {
        auto coeffTemp = std::make_unique<float[]>(n);

        copy(coif10, n, coeffTemp.get());
        for (i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        copyReverse(coeffTemp.get(), n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        copy(coeffTemp.get(), n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }
    if (name == string_view{"coif11"})
    {
        auto coeffTemp = std::make_unique<float[]>(n);

        copy(coif11, n, coeffTemp.get());
        for (i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        copyReverse(coeffTemp.get(), n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        copy(coeffTemp.get(), n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }
    if (name == string_view{"coif12"})
    {
        auto coeffTemp = std::make_unique<float[]>(n);

        copy(coif12, n, coeffTemp.get());
        for (i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        copyReverse(coeffTemp.get(), n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        copy(coeffTemp.get(), n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }
    if (name == string_view{"coif13"})
    {
        auto coeffTemp = std::make_unique<float[]>(n);

        copy(coif13, n, coeffTemp.get());
        for (i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        copyReverse(coeffTemp.get(), n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        copy(coeffTemp.get(), n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }
    if (name == string_view{"coif14"})
    {
        auto coeffTemp = std::make_unique<float[]>(n);

        copy(coif14, n, coeffTemp.get());
        for (i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        copyReverse(coeffTemp.get(), n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        copy(coeffTemp.get(), n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }
    if (name == string_view{"coif15"})
    {
        auto coeffTemp = std::make_unique<float[]>(n);

        copy(coif15, n, coeffTemp.get());
        for (i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        copyReverse(coeffTemp.get(), n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        copy(coeffTemp.get(), n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }
    if (name == string_view{"coif16"})
    {
        auto coeffTemp = std::make_unique<float[]>(n);

        copy(coif16, n, coeffTemp.get());
        for (i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        copyReverse(coeffTemp.get(), n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        copy(coeffTemp.get(), n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }
    if (name == string_view{"coif17"})
    {
        auto coeffTemp = std::make_unique<float[]>(n);

        copy(coif17, n, coeffTemp.get());
        for (i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        copyReverse(coeffTemp.get(), n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        copy(coeffTemp.get(), n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }
    if (name == string_view{"sym2"})
    {
        copyReverse(sym2, n, lp1);
        qmfWrev(sym2, n, hp1);
        copy(sym2, n, lp2);
        qmfEven(sym2, n, hp2);
        return n;
    }

    if (name == string_view{"sym3"})
    {
        copyReverse(sym3, n, lp1);
        qmfWrev(sym3, n, hp1);
        copy(sym3, n, lp2);
        qmfEven(sym3, n, hp2);
        return n;
    }

    if (name == string_view{"sym4"})
    {
        copyReverse(sym4, n, lp1);
        qmfWrev(sym4, n, hp1);
        copy(sym4, n, lp2);
        qmfEven(sym4, n, hp2);
        return n;
    }

    if (name == string_view{"sym5"})
    {
        copyReverse(sym5, n, lp1);
        qmfWrev(sym5, n, hp1);
        copy(sym5, n, lp2);
        qmfEven(sym5, n, hp2);
        return n;
    }

    if (name == string_view{"sym6"})
    {
        copyReverse(sym6, n, lp1);
        qmfWrev(sym6, n, hp1);
        copy(sym6, n, lp2);
        qmfEven(sym6, n, hp2);
        return n;
    }

    if (name == string_view{"sym7"})
    {
        copyReverse(sym7, n, lp1);
        qmfWrev(sym7, n, hp1);
        copy(sym7, n, lp2);
        qmfEven(sym7, n, hp2);
        return n;
    }

    if (name == string_view{"sym8"})
    {
        copyReverse(sym8, n, lp1);
        qmfWrev(sym8, n, hp1);
        copy(sym8, n, lp2);
        qmfEven(sym8, n, hp2);
        return n;
    }

    if (name == string_view{"sym9"})
    {
        copyReverse(sym9, n, lp1);
        qmfWrev(sym9, n, hp1);
        copy(sym9, n, lp2);
        qmfEven(sym9, n, hp2);
        return n;
    }

    if (name == string_view{"sym10"})
    {
        copyReverse(sym10, n, lp1);
        qmfWrev(sym10, n, hp1);
        copy(sym10, n, lp2);
        qmfEven(sym10, n, hp2);
        return n;
    }
    if (name == string_view{"sym11"})
    {
        copyReverse(sym11, n, lp1);
        qmfWrev(sym11, n, hp1);
        copy(sym11, n, lp2);
        qmfEven(sym11, n, hp2);
        return n;
    }
    if (name == string_view{"sym12"})
    {
        copyReverse(sym12, n, lp1);
        qmfWrev(sym12, n, hp1);
        copy(sym12, n, lp2);
        qmfEven(sym12, n, hp2);
        return n;
    }
    if (name == string_view{"sym13"})
    {
        copyReverse(sym13, n, lp1);
        qmfWrev(sym13, n, hp1);
        copy(sym13, n, lp2);
        qmfEven(sym13, n, hp2);
        return n;
    }
    if (name == string_view{"sym14"})
    {
        copyReverse(sym14, n, lp1);
        qmfWrev(sym14, n, hp1);
        copy(sym14, n, lp2);
        qmfEven(sym14, n, hp2);
        return n;
    }
    if (name == string_view{"sym15"})
    {
        copyReverse(sym15, n, lp1);
        qmfWrev(sym15, n, hp1);
        copy(sym15, n, lp2);
        qmfEven(sym15, n, hp2);
        return n;
    }
    if (name == string_view{"sym16"})
    {
        copyReverse(sym16, n, lp1);
        qmfWrev(sym16, n, hp1);
        copy(sym16, n, lp2);
        qmfEven(sym16, n, hp2);
        return n;
    }
    if (name == string_view{"sym17"})
    {
        copyReverse(sym17, n, lp1);
        qmfWrev(sym17, n, hp1);
        copy(sym17, n, lp2);
        qmfEven(sym17, n, hp2);
        return n;
    }
    if (name == string_view{"sym18"})
    {
        copyReverse(sym18, n, lp1);
        qmfWrev(sym18, n, hp1);
        copy(sym18, n, lp2);
        qmfEven(sym18, n, hp2);
        return n;
    }
    if (name == string_view{"sym19"})
    {
        copyReverse(sym19, n, lp1);
        qmfWrev(sym19, n, hp1);
        copy(sym19, n, lp2);
        qmfEven(sym19, n, hp2);
        return n;
    }
    if (name == string_view{"sym20"})
    {
        copyReverse(sym20, n, lp1);
        qmfWrev(sym20, n, hp1);
        copy(sym20, n, lp2);
        qmfEven(sym20, n, hp2);
        return n;
    }
    if (name == string_view{"meyer"})
    {
        copyReverse(meyer, n, lp1);
        qmfWrev(meyer, n, hp1);
        copy(meyer, n, lp2);
        qmfEven(meyer, n, hp2);
        return n;
    }

    fmt::printf("\n Filter Not in Database \n");
    return -1;
}

}  // namespace

Wavelet::Wavelet(char const* name)
    : name_{name}
    , size_{static_cast<std::size_t>(waveletFilterLength(name))}
    , params_{std::make_unique<float[]>(4 * size_)}
    , lpd_{&params_[0], size_}
    , hpd_{&params_[size_], size_}
    , lpr_{&params_[2 * size_], size_}
    , hpr_{&params_[3 * size_], size_}
{
    auto* p = params_.get();
    if (name != nullptr) { waveletFilterCoefficients(name, p, p + size_, p + 2 * size_, p + 3 * size_); }
}

auto summary(Wavelet const& obj) -> void
{
    auto const n = obj.size();
    fmt::printf("\n");
    fmt::printf("Wavelet Name: %s \n", obj.name().c_str());
    fmt::printf("\n");
    fmt::printf("Wavelet Filters \n");
    fmt::printf("lpd: [");
    for (std::size_t i = 0; i < n - 1; ++i) { fmt::printf("%g,", obj.lpd()[i]); }
    fmt::printf("%g", obj.lpd()[n - 1]);
    fmt::printf("] \n");
    fmt::printf("hpd: [");
    for (std::size_t i = 0; i < n - 1; ++i) { fmt::printf("%g,", obj.hpd()[i]); }
    fmt::printf("%g", obj.hpd()[n - 1]);
    fmt::printf("] \n");
    fmt::printf("lpr: [");
    for (std::size_t i = 0; i < n - 1; ++i) { fmt::printf("%g,", obj.lpr()[i]); }
    fmt::printf("%g", obj.lpr()[n - 1]);
    fmt::printf("] \n");
    fmt::printf("hpr: [");
    for (std::size_t i = 0; i < n - 1; ++i) { fmt::printf("%g,", obj.hpr()[i]); }
    fmt::printf("%g", obj.hpr()[n - 1]);
    fmt::printf("] \n");
}

}  // namespace mc::dsp