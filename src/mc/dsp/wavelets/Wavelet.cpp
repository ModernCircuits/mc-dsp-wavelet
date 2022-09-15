#include "Wavelet.hpp"

#include <mc/dsp/convolution/FFTConvolver.hpp>
#include <mc/dsp/wavelets/filters/coif.hpp>
#include <mc/dsp/wavelets/filters/daubechies.hpp>
#include <mc/dsp/wavelets/filters/h.hpp>
#include <mc/dsp/wavelets/filters/meyer.hpp>
#include <mc/dsp/wavelets/filters/sym.hpp>

#include <mc/core/cmath.hpp>
#include <mc/core/cstring.hpp>
#include <mc/core/format.hpp>
#include <mc/core/iterator.hpp>
#include <mc/core/memory.hpp>
#include <mc/core/numbers.hpp>
#include <mc/core/string_view.hpp>
#include <mc/core/utility.hpp>

#include <fmt/printf.h>

namespace mc::dsp {

namespace {

template<typename T>
auto qmfEven(T const* in, int n, T* out)
{
    for (auto count = 0; count < n; ++count) {
        auto evenIndex = count % 2 == 0;
        out[count]     = in[n - count - 1] * (evenIndex ? T{1} : T{-1});
    }
}

template<typename T>
auto qmfWrev(T const* in, int n, T* out)
{
    auto tmp = makeUnique<T[]>(n);
    qmfEven(in, n, tmp.get());
    std::reverse_copy(tmp.get(), tmp.get() + n, out);
}

auto waveletFilterLength(char const* name) -> int
{
    int len = strlen(name);
    int n   = 0;
    if (name == StringView{"haar"} || name == StringView{"db1"}) { return 2; }
    if (len > 2 && strstr(name, "db") != nullptr) {
        auto newStr = makeUnique<char[]>((len - 2 + 1));
        for (auto i = 2; i < len + 1; i++) { newStr[i - 2] = name[i]; }

        n = atoi(newStr.get());
        if (n > 38) { throw std::invalid_argument("wavelet filter not in database"); }

        return n * 2;
    }
    if (name == StringView{"bior1.1"}) { return 2; }
    if (name == StringView{"bior1.3"}) { return 6; }
    if (name == StringView{"bior1.5"}) { return 10; }
    if (name == StringView{"bior2.2"}) { return 6; }
    if (name == StringView{"bior2.4"}) { return 10; }
    if (name == StringView{"bior2.6"}) { return 14; }
    if (name == StringView{"bior2.8"}) { return 18; }
    if (name == StringView{"bior3.1"}) { return 4; }
    if (name == StringView{"bior3.3"}) { return 8; }
    if (name == StringView{"bior3.5"}) { return 12; }
    if (name == StringView{"bior3.7"}) { return 16; }
    if (name == StringView{"bior3.9"}) { return 20; }
    if (name == StringView{"bior4.4"}) { return 10; }
    if (name == StringView{"bior5.5"}) { return 12; }
    if (name == StringView{"bior6.8"}) { return 18; }
    if (name == StringView{"rbior1.1"}) { return 2; }
    if (name == StringView{"rbior1.3"}) { return 6; }
    if (name == StringView{"rbior1.5"}) { return 10; }
    if (name == StringView{"rbior2.2"}) { return 6; }
    if (name == StringView{"rbior2.4"}) { return 10; }
    if (name == StringView{"rbior2.6"}) { return 14; }
    if (name == StringView{"rbior2.8"}) { return 18; }
    if (name == StringView{"rbior3.1"}) { return 4; }
    if (name == StringView{"rbior3.3"}) { return 8; }
    if (name == StringView{"rbior3.5"}) { return 12; }
    if (name == StringView{"rbior3.7"}) { return 16; }
    if (name == StringView{"rbior3.9"}) { return 20; }
    if (name == StringView{"rbior4.4"}) { return 10; }
    if (name == StringView{"rbior5.5"}) { return 12; }
    if (name == StringView{"rbior6.8"}) { return 18; }
    if (len > 4 && strstr(name, "coif") != nullptr) {
        auto newStr = makeUnique<char[]>((len - 4 + 1));
        for (auto i = 4; i < len + 1; i++) { newStr[i - 4] = name[i]; }

        n = atoi(newStr.get());
        if (n > 17) { throw std::invalid_argument("wavelet filter not in database"); }

        return n * 6;
    }
    if (len > 3 && strstr(name, "sym") != nullptr) {
        auto newStr = makeUnique<char[]>((len - 3 + 1));
        for (auto i = 3; i < len + 1; i++) { newStr[i - 3] = name[i]; }

        n = atoi(newStr.get());
        if (n > 20 || n < 2) {
            throw std::invalid_argument("wavelet filter not in database");
        }

        return n * 2;
    }
    if (strcmp(name, "meyer") == 0) { return 102; }

    throw std::invalid_argument("wavelet filter not in database");
}

auto fillDauberchies(char const* name, float* lp1, float* hp1, float* lp2, float* hp2)
{
    using namespace std::string_view_literals;

    auto const* coeffs = static_cast<float*>(nullptr);
    if (name == "haar"sv || name == "db1"sv) { coeffs = std::cbegin(daubechies1); }
    if (name == "db2"sv) { coeffs = std::cbegin(daubechies2); }
    if (name == "db3"sv) { coeffs = std::cbegin(daubechies3); }
    if (name == "db4"sv) { coeffs = std::cbegin(daubechies4); }
    if (name == "db5"sv) { coeffs = std::cbegin(daubechies5); }
    if (name == "db6"sv) { coeffs = std::cbegin(daubechies6); }
    if (name == "db7"sv) { coeffs = std::cbegin(daubechies7); }
    if (name == "db8"sv) { coeffs = std::cbegin(daubechies8); }
    if (name == "db9"sv) { coeffs = std::cbegin(daubechies9); }
    if (name == "db10"sv) { coeffs = std::cbegin(daubechies10); }
    if (name == "db11"sv) { coeffs = std::cbegin(daubechies11); }
    if (name == "db12"sv) { coeffs = std::cbegin(daubechies12); }
    if (name == "db13"sv) { coeffs = std::cbegin(daubechies13); }
    if (name == "db14"sv) { coeffs = std::cbegin(daubechies14); }
    if (name == "db15"sv) { coeffs = std::cbegin(daubechies15); }
    if (name == "db16"sv) { coeffs = std::cbegin(daubechies16); }
    if (name == "db17"sv) { coeffs = std::cbegin(daubechies17); }
    if (name == "db18"sv) { coeffs = std::cbegin(daubechies18); }
    if (name == "db19"sv) { coeffs = std::cbegin(daubechies19); }
    if (name == "db20"sv) { coeffs = std::cbegin(daubechies20); }
    if (name == "db21"sv) { coeffs = std::cbegin(daubechies21); }
    if (name == "db22"sv) { coeffs = std::cbegin(daubechies22); }
    if (name == "db23"sv) { coeffs = std::cbegin(daubechies23); }
    if (name == "db24"sv) { coeffs = std::cbegin(daubechies24); }
    if (name == "db25"sv) { coeffs = std::cbegin(daubechies25); }
    if (name == "db26"sv) { coeffs = std::cbegin(daubechies26); }
    if (name == "db27"sv) { coeffs = std::cbegin(daubechies27); }
    if (name == "db28"sv) { coeffs = std::cbegin(daubechies28); }
    if (name == "db29"sv) { coeffs = std::cbegin(daubechies29); }
    if (name == "db30"sv) { coeffs = std::cbegin(daubechies30); }
    if (name == "db31"sv) { coeffs = std::cbegin(daubechies31); }
    if (name == "db32"sv) { coeffs = std::cbegin(daubechies32); }
    if (name == "db33"sv) { coeffs = std::cbegin(daubechies33); }
    if (name == "db34"sv) { coeffs = std::cbegin(daubechies34); }
    if (name == "db35"sv) { coeffs = std::cbegin(daubechies35); }
    if (name == "db36"sv) { coeffs = std::cbegin(daubechies36); }
    if (name == "db37"sv) { coeffs = std::cbegin(daubechies37); }
    if (name == "db38"sv) { coeffs = std::cbegin(daubechies38); }
    if (name == "sym2"sv) { coeffs = std::cbegin(sym2); }
    if (name == "sym3"sv) { coeffs = std::cbegin(sym3); }
    if (name == "sym4"sv) { coeffs = std::cbegin(sym4); }
    if (name == "sym5"sv) { coeffs = std::cbegin(sym5); }
    if (name == "sym6"sv) { coeffs = std::cbegin(sym6); }
    if (name == "sym7"sv) { coeffs = std::cbegin(sym7); }
    if (name == "sym8"sv) { coeffs = std::cbegin(sym8); }
    if (name == "sym9"sv) { coeffs = std::cbegin(sym9); }
    if (name == "sym10"sv) { coeffs = std::cbegin(sym10); }
    if (name == "sym11"sv) { coeffs = std::cbegin(sym11); }
    if (name == "sym12"sv) { coeffs = std::cbegin(sym12); }
    if (name == "sym13"sv) { coeffs = std::cbegin(sym13); }
    if (name == "sym14"sv) { coeffs = std::cbegin(sym14); }
    if (name == "sym15"sv) { coeffs = std::cbegin(sym15); }
    if (name == "sym16"sv) { coeffs = std::cbegin(sym16); }
    if (name == "sym17"sv) { coeffs = std::cbegin(sym17); }
    if (name == "sym18"sv) { coeffs = std::cbegin(sym18); }
    if (name == "sym19"sv) { coeffs = std::cbegin(sym19); }
    if (name == "sym20"sv) { coeffs = std::cbegin(sym20); }
    if (name == "meyer"sv) { coeffs = std::cbegin(meyer); }

    auto const n = waveletFilterLength(name);
    std::reverse_copy(coeffs, coeffs + n, lp1);
    qmfWrev(coeffs, n, hp1);
    std::copy(coeffs, std::next(coeffs, n), lp2);
    qmfEven(coeffs, n, hp2);
    return n;
}

auto waveletFilterCoefficients(
    char const* name,
    float* lp1,
    float* hp1,
    float* lp2,
    float* hp2
) -> int
{
    auto const n        = waveletFilterLength(name);
    auto const nameView = StringView{name};
    if (nameView.find("haar") != StringView::npos) {
        return fillDauberchies(name, lp1, hp1, lp2, hp2);
    }
    if (nameView.find("db") != StringView::npos) {
        return fillDauberchies(name, lp1, hp1, lp2, hp2);
    }
    if (nameView.find("sym") != StringView::npos) {
        return fillDauberchies(name, lp1, hp1, lp2, hp2);
    }
    if (nameView.find("meyer") != StringView::npos) {
        return fillDauberchies(name, lp1, hp1, lp2, hp2);
    }

    if (name == StringView{"bior1.1"}) {
        std::reverse_copy(hm111, hm111 + n, lp1);
        qmfWrev(h1 + 4, n, hp1);
        std::copy(h1 + 4, h1 + 4 + n, lp2);
        qmfEven(hm111, n, hp2);
        return n;
    }

    if (name == StringView{"bior1.3"}) {
        std::reverse_copy(hm113, hm113 + n, lp1);
        qmfWrev(h1 + 2, n, hp1);
        std::copy(h1 + 2, h1 + 2 + n, lp2);
        qmfEven(hm113, n, hp2);
        return n;
    }

    if (name == StringView{"bior1.5"}) {
        std::reverse_copy(hm115, hm115 + n, lp1);
        qmfWrev(h1, n, hp1);
        std::copy(h1, h1 + n, lp2);
        qmfEven(hm115, n, hp2);
        return n;
    }

    if (name == StringView{"bior2.2"}) {
        std::reverse_copy(hm222, hm222 + n, lp1);
        qmfWrev(h2 + 6, n, hp1);
        std::copy(h2 + 6, h2 + 6 + n, lp2);
        qmfEven(hm222, n, hp2);
        return n;
    }

    if (name == StringView{"bior2.4"}) {
        std::reverse_copy(hm224, hm224 + n, lp1);
        qmfWrev(h2 + 4, n, hp1);
        std::copy(h2 + 4, h2 + 4 + n, lp2);
        qmfEven(hm224, n, hp2);
        return n;
    }

    if (name == StringView{"bior2.6"}) {
        std::reverse_copy(hm226, hm226 + n, lp1);
        qmfWrev(h2 + 2, n, hp1);
        std::copy(h2 + 2, h2 + 2 + n, lp2);
        qmfEven(hm226, n, hp2);
        return n;
    }

    if (name == StringView{"bior2.8"}) {
        std::reverse_copy(hm228, hm228 + n, lp1);
        qmfWrev(h2, n, hp1);
        std::copy(h2, h2 + n, lp2);
        qmfEven(hm228, n, hp2);
        return n;
    }

    if (name == StringView{"bior3.1"}) {
        std::reverse_copy(hm331, hm331 + n, lp1);
        qmfWrev(h3 + 8, n, hp1);
        std::copy(h3 + 8, h3 + 8 + n, lp2);
        qmfEven(hm331, n, hp2);
        return n;
    }

    if (name == StringView{"bior3.3"}) {
        std::reverse_copy(hm333, hm333 + n, lp1);
        qmfWrev(h3 + 6, n, hp1);
        std::copy(h3 + 6, h3 + 6 + n, lp2);
        qmfEven(hm333, n, hp2);
        return n;
    }

    if (name == StringView{"bior3.5"}) {
        std::reverse_copy(hm335, hm335 + n, lp1);
        qmfWrev(h3 + 4, n, hp1);
        std::copy(h3 + 4, h3 + 4 + n, lp2);
        qmfEven(hm335, n, hp2);
        return n;
    }

    if (name == StringView{"bior3.7"}) {
        std::reverse_copy(hm337, hm337 + n, lp1);
        qmfWrev(h3 + 2, n, hp1);
        std::copy(h3 + 2, h3 + 2 + n, lp2);
        qmfEven(hm337, n, hp2);
        return n;
    }

    if (name == StringView{"bior3.9"}) {
        std::reverse_copy(hm339, hm339 + n, lp1);
        qmfWrev(h3, n, hp1);
        std::copy(h3, h3 + n, lp2);
        qmfEven(hm339, n, hp2);
        return n;
    }

    if (name == StringView{"bior4.4"}) {
        std::reverse_copy(hm444, hm444 + n, lp1);
        qmfWrev(h4, n, hp1);
        std::copy(h4, h4 + n, lp2);
        qmfEven(hm444, n, hp2);
        return n;
    }

    if (name == StringView{"bior5.5"}) {
        std::reverse_copy(hm555, hm555 + n, lp1);
        qmfWrev(h5, n, hp1);
        std::copy(h5, h5 + n, lp2);
        qmfEven(hm555, n, hp2);
        return n;
    }

    if (name == StringView{"bior6.8"}) {
        std::reverse_copy(hm668, hm668 + n, lp1);
        qmfWrev(h6, n, hp1);
        std::copy(h6, h6 + n, lp2);
        qmfEven(hm668, n, hp2);
        return n;
    }

    if (name == StringView{"rbior1.1"}) {
        std::reverse_copy(h1 + 4, h1 + 4 + n, lp1);
        qmfWrev(hm111, n, hp1);
        std::copy(hm111, hm111 + n, lp2);
        qmfEven(h1 + 4, n, hp2);
        return n;
    }

    if (name == StringView{"rbior1.3"}) {
        std::reverse_copy(h1 + 2, h1 + 2 + n, lp1);
        qmfWrev(hm113, n, hp1);
        std::copy(hm113, hm113 + n, lp2);
        qmfEven(h1 + 2, n, hp2);
        return n;
    }

    if (name == StringView{"rbior1.5"}) {
        std::reverse_copy(h1, h1 + n, lp1);
        qmfWrev(hm115, n, hp1);
        std::copy(hm115, hm115 + n, lp2);
        qmfEven(h1, n, hp2);
        return n;
    }

    if (name == StringView{"rbior2.2"}) {
        std::reverse_copy(h2 + 6, h2 + 6 + n, lp1);
        qmfWrev(hm222, n, hp1);
        std::copy(hm222, hm222 + n, lp2);
        qmfEven(h2 + 6, n, hp2);
        return n;
    }

    if (name == StringView{"rbior2.4"}) {
        std::reverse_copy(h2 + 4, h2 + 4 + n, lp1);
        qmfWrev(hm224, n, hp1);
        std::copy(hm224, hm224 + n, lp2);
        qmfEven(h2 + 4, n, hp2);
        return n;
    }

    if (name == StringView{"rbior2.6"}) {
        std::reverse_copy(h2 + 2, h2 + 2 + n, lp1);
        qmfWrev(hm226, n, hp1);
        std::copy(hm226, hm226 + n, lp2);
        qmfEven(h2 + 2, n, hp2);
        return n;
    }

    if (name == StringView{"rbior2.8"}) {
        std::reverse_copy(h2, h2 + n, lp1);
        qmfWrev(hm228, n, hp1);
        std::copy(hm228, hm228 + n, lp2);
        qmfEven(h2, n, hp2);
        return n;
    }

    if (name == StringView{"rbior3.1"}) {
        std::reverse_copy(h3 + 8, h3 + 8 + n, lp1);
        qmfWrev(hm331, n, hp1);
        std::copy(hm331, hm331 + n, lp2);
        qmfEven(h3 + 8, n, hp2);
        return n;
    }

    if (name == StringView{"rbior3.3"}) {
        std::reverse_copy(h3 + 6, h3 + 6 + n, lp1);
        qmfWrev(hm333, n, hp1);
        std::copy(hm333, hm333 + n, lp2);
        qmfEven(h3 + 6, n, hp2);
        return n;
    }

    if (name == StringView{"rbior3.5"}) {
        std::reverse_copy(h3 + 4, h3 + 4 + n, lp1);
        qmfWrev(hm335, n, hp1);
        std::copy(hm335, hm335 + n, lp2);
        qmfEven(h3 + 4, n, hp2);
        return n;
    }

    if (name == StringView{"rbior3.7"}) {
        std::reverse_copy(h3 + 2, h3 + 2 + n, lp1);
        qmfWrev(hm337, n, hp1);
        std::copy(hm337, hm337 + n, lp2);
        qmfEven(h3 + 2, n, hp2);
        return n;
    }

    if (name == StringView{"rbior3.9"}) {
        std::reverse_copy(h3, h3 + n, lp1);
        qmfWrev(hm339, n, hp1);
        std::copy(hm339, hm339 + n, lp2);
        qmfEven(h3, n, hp2);
        return n;
    }

    if (name == StringView{"rbior4.4"}) {
        std::reverse_copy(h4, h4 + n, lp1);
        qmfWrev(hm444, n, hp1);
        std::copy(hm444, hm444 + n, lp2);
        qmfEven(h4, n, hp2);
        return n;
    }

    if (name == StringView{"rbior5.5"}) {
        std::reverse_copy(h5, h5 + n, lp1);
        qmfWrev(hm555, n, hp1);
        std::copy(hm555, hm555 + n, lp2);
        qmfEven(h5, n, hp2);
        return n;
    }

    if (name == StringView{"rbior6.8"}) {
        std::reverse_copy(h6, h6 + n, lp1);
        qmfWrev(hm668, n, hp1);
        std::copy(hm668, hm668 + n, lp2);
        qmfEven(h6, n, hp2);
        return n;
    }

    if (name == StringView{"coif1"}) {
        auto coeffTemp = makeUnique<float[]>(n);

        std::copy(coif1, coif1 + n, coeffTemp.get());
        for (auto i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        std::reverse_copy(coeffTemp.get(), coeffTemp.get() + n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        std::copy(coeffTemp.get(), coeffTemp.get() + n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }

    if (name == StringView{"coif2"}) {
        auto coeffTemp = makeUnique<float[]>(n);

        std::copy(coif2, coif2 + n, coeffTemp.get());
        for (auto i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        std::reverse_copy(coeffTemp.get(), coeffTemp.get() + n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        std::copy(coeffTemp.get(), coeffTemp.get() + n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }

    if (name == StringView{"coif3"}) {
        auto coeffTemp = makeUnique<float[]>(n);

        std::copy(coif3, coif3 + n, coeffTemp.get());
        for (auto i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        std::reverse_copy(coeffTemp.get(), coeffTemp.get() + n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        std::copy(coeffTemp.get(), coeffTemp.get() + n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }

    if (name == StringView{"coif4"}) {
        auto coeffTemp = makeUnique<float[]>(n);

        std::copy(coif4, coif4 + n, coeffTemp.get());
        for (auto i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        std::reverse_copy(coeffTemp.get(), coeffTemp.get() + n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        std::copy(coeffTemp.get(), coeffTemp.get() + n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }

    if (name == StringView{"coif5"}) {
        auto coeffTemp = makeUnique<float[]>(n);

        std::copy(coif5, coif5 + n, coeffTemp.get());
        for (auto i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        std::reverse_copy(coeffTemp.get(), coeffTemp.get() + n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        std::copy(coeffTemp.get(), coeffTemp.get() + n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }

    if (name == StringView{"coif6"}) {
        auto coeffTemp = makeUnique<float[]>(n);

        std::copy(coif6, coif6 + n, coeffTemp.get());
        for (auto i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        std::reverse_copy(coeffTemp.get(), coeffTemp.get() + n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        std::copy(coeffTemp.get(), coeffTemp.get() + n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }

    if (name == StringView{"coif7"}) {
        auto coeffTemp = makeUnique<float[]>(n);

        std::copy(coif7, coif7 + n, coeffTemp.get());
        for (auto i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        std::reverse_copy(coeffTemp.get(), coeffTemp.get() + n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        std::copy(coeffTemp.get(), coeffTemp.get() + n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }

    if (name == StringView{"coif8"}) {
        auto coeffTemp = makeUnique<float[]>(n);

        std::copy(coif8, coif8 + n, coeffTemp.get());
        for (auto i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        std::reverse_copy(coeffTemp.get(), coeffTemp.get() + n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        std::copy(coeffTemp.get(), coeffTemp.get() + n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }

    if (name == StringView{"coif9"}) {
        auto coeffTemp = makeUnique<float[]>(n);

        std::copy(coif9, coif9 + n, coeffTemp.get());
        for (auto i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        std::reverse_copy(coeffTemp.get(), coeffTemp.get() + n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        std::copy(coeffTemp.get(), coeffTemp.get() + n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }

    if (name == StringView{"coif10"}) {
        auto coeffTemp = makeUnique<float[]>(n);

        std::copy(coif10, coif10 + n, coeffTemp.get());
        for (auto i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        std::reverse_copy(coeffTemp.get(), coeffTemp.get() + n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        std::copy(coeffTemp.get(), coeffTemp.get() + n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }
    if (name == StringView{"coif11"}) {
        auto coeffTemp = makeUnique<float[]>(n);

        std::copy(coif11, coif11 + n, coeffTemp.get());
        for (auto i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        std::reverse_copy(coeffTemp.get(), coeffTemp.get() + n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        std::copy(coeffTemp.get(), coeffTemp.get() + n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }
    if (name == StringView{"coif12"}) {
        auto coeffTemp = makeUnique<float[]>(n);

        std::copy(coif12, coif12 + n, coeffTemp.get());
        for (auto i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        std::reverse_copy(coeffTemp.get(), coeffTemp.get() + n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        std::copy(coeffTemp.get(), coeffTemp.get() + n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }
    if (name == StringView{"coif13"}) {
        auto coeffTemp = makeUnique<float[]>(n);

        std::copy(coif13, coif13 + n, coeffTemp.get());
        for (auto i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        std::reverse_copy(coeffTemp.get(), coeffTemp.get() + n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        std::copy(coeffTemp.get(), coeffTemp.get() + n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }
    if (name == StringView{"coif14"}) {
        auto coeffTemp = makeUnique<float[]>(n);

        std::copy(coif14, coif14 + n, coeffTemp.get());
        for (auto i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        std::reverse_copy(coeffTemp.get(), coeffTemp.get() + n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        std::copy(coeffTemp.get(), coeffTemp.get() + n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }
    if (name == StringView{"coif15"}) {
        auto coeffTemp = makeUnique<float[]>(n);

        std::copy(coif15, coif15 + n, coeffTemp.get());
        for (auto i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        std::reverse_copy(coeffTemp.get(), coeffTemp.get() + n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        std::copy(coeffTemp.get(), coeffTemp.get() + n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }
    if (name == StringView{"coif16"}) {
        auto coeffTemp = makeUnique<float[]>(n);

        std::copy(coif16, coif16 + n, coeffTemp.get());
        for (auto i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        std::reverse_copy(coeffTemp.get(), coeffTemp.get() + n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        std::copy(coeffTemp.get(), coeffTemp.get() + n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }
    if (name == StringView{"coif17"}) {
        auto coeffTemp = makeUnique<float[]>(n);

        std::copy(coif17, coif17 + n, coeffTemp.get());
        for (auto i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        std::reverse_copy(coeffTemp.get(), coeffTemp.get() + n, lp1);
        qmfWrev(coeffTemp.get(), n, hp1);
        std::copy(coeffTemp.get(), coeffTemp.get() + n, lp2);
        qmfEven(coeffTemp.get(), n, hp2);

        return n;
    }

    fmt::printf("\n Filter Not in Database \n");
    return -1;
}

}  // namespace

Wavelet::Wavelet(char const* name)
    : name_{name}
    , size_{static_cast<std::size_t>(waveletFilterLength(name))}
    , params_{makeUnique<float[]>(4 * size_)}
{
    waveletFilterCoefficients(name, data(lpd()), data(hpd()), data(lpr()), data(hpr()));
}

auto Wavelet::size() const noexcept -> std::size_t { return size_; }

auto Wavelet::name() const noexcept -> String const& { return name_; }

auto Wavelet::lpd() const noexcept -> Span<float> { return {&params_[0], size()}; }

auto Wavelet::hpd() const noexcept -> Span<float> { return {&params_[size()], size()}; }

auto Wavelet::lpr() const noexcept -> Span<float>
{
    return {&params_[size() * 2U], size()};
}

auto Wavelet::hpr() const noexcept -> Span<float>
{
    return {&params_[size() * 3U], size()};
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
