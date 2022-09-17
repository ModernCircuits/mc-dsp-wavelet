#include "Wavelet.hpp"

#include <mc/dsp/convolution/FFTConvolver.hpp>
#include <mc/dsp/wavelets/family.hpp>
#include <mc/dsp/wavelets/filters/coif.hpp>
#include <mc/dsp/wavelets/filters/daubechies.hpp>
#include <mc/dsp/wavelets/filters/sym.hpp>

#include <mc/core/algorithm.hpp>
#include <mc/core/cmath.hpp>
#include <mc/core/cstring.hpp>
#include <mc/core/format.hpp>
#include <mc/core/iterator.hpp>
#include <mc/core/memory.hpp>
#include <mc/core/numbers.hpp>
#include <mc/core/raise.hpp>
#include <mc/core/stdexcept.hpp>
#include <mc/core/string_view.hpp>
#include <mc/core/utility.hpp>

namespace mc::dsp {

namespace {

template<typename T>
auto qmfEven(Span<T const> in, T* out)
{
    auto const n = in.size();
    for (std::size_t count = 0; count < n; ++count) {
        auto evenIndex = count % 2 == 0;
        out[count]     = in[n - count - 1] * (evenIndex ? T{1} : T{-1});
    }
}

template<typename T>
auto qmfWrev(Span<T const> in, T* out)
{
    qmfEven(in, out);
    ranges::reverse(Span<T>{out, in.size()});
}

auto filterLength(char const* name) -> std::size_t
{
    auto nameMatches  = [name](auto const& w) { return w.name == name; };
    auto const filter = ranges::find_if(allWavelets<float>, nameMatches);
    if (filter == ranges::end(allWavelets<float>)) {
        raise<InvalidArgument>("wavelet filter not in database");
    }
    return filter->length;
}

auto fillDauberchies(char const* name, float* lp1, float* hp1, float* lp2, float* hp2)
    -> size_t
{
    using namespace std::string_view_literals;

    if (name == "haar"sv) { return fillDauberchies("db1", lp1, hp1, lp2, hp2); }

    auto nameMatches  = [name](auto const& w) { return w.name == name; };
    auto const filter = ranges::find_if(allWavelets<float>, nameMatches);
    if (filter == ranges::end(allWavelets<float>)) {
        raise<InvalidArgument>("wavelet filter not in database");
    }

    ranges::reverse_copy(filter->coefficients, lp1);
    qmfWrev(filter->coefficients, hp1);

    ranges::copy(filter->coefficients, lp2);
    qmfEven(filter->coefficients, hp2);

    return filter->length;
}

auto waveletFilterCoefficients(
    char const* name,
    float* lp1,
    float* hp1,
    float* lp2,
    float* hp2
) -> std::size_t
{
    auto const n        = filterLength(name);
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

    if (name == StringView{"coif1"}) {
        auto coeffTemp = makeUnique<float[]>(n);

        std::copy(cbegin(coif1<float>), cend(coif1<float>), coeffTemp.get());
        for (auto i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        std::reverse_copy(coeffTemp.get(), coeffTemp.get() + n, lp1);
        qmfWrev({coeffTemp.get(), n}, hp1);
        std::copy(coeffTemp.get(), coeffTemp.get() + n, lp2);
        qmfEven({coeffTemp.get(), n}, hp2);

        return n;
    }

    if (name == StringView{"coif2"}) {
        auto coeffTemp = makeUnique<float[]>(n);

        std::copy(cbegin(coif2<float>), cend(coif2<float>), coeffTemp.get());
        for (auto i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        std::reverse_copy(coeffTemp.get(), coeffTemp.get() + n, lp1);
        qmfWrev({coeffTemp.get(), n}, hp1);
        std::copy(coeffTemp.get(), coeffTemp.get() + n, lp2);
        qmfEven({coeffTemp.get(), n}, hp2);

        return n;
    }

    if (name == StringView{"coif3"}) {
        auto coeffTemp = makeUnique<float[]>(n);

        std::copy(cbegin(coif3<float>), cend(coif3<float>), coeffTemp.get());
        for (auto i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        std::reverse_copy(coeffTemp.get(), coeffTemp.get() + n, lp1);
        qmfWrev({coeffTemp.get(), n}, hp1);
        std::copy(coeffTemp.get(), coeffTemp.get() + n, lp2);
        qmfEven({coeffTemp.get(), n}, hp2);

        return n;
    }

    if (name == StringView{"coif4"}) {
        auto coeffTemp = makeUnique<float[]>(n);

        std::copy(cbegin(coif4<float>), cend(coif4<float>), coeffTemp.get());
        for (auto i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        std::reverse_copy(coeffTemp.get(), coeffTemp.get() + n, lp1);
        qmfWrev({coeffTemp.get(), n}, hp1);
        std::copy(coeffTemp.get(), coeffTemp.get() + n, lp2);
        qmfEven({coeffTemp.get(), n}, hp2);

        return n;
    }

    if (name == StringView{"coif5"}) {
        auto coeffTemp = makeUnique<float[]>(n);

        std::copy(cbegin(coif5<float>), cend(coif5<float>), coeffTemp.get());
        for (auto i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        std::reverse_copy(coeffTemp.get(), coeffTemp.get() + n, lp1);
        qmfWrev({coeffTemp.get(), n}, hp1);
        std::copy(coeffTemp.get(), coeffTemp.get() + n, lp2);
        qmfEven({coeffTemp.get(), n}, hp2);

        return n;
    }

    if (name == StringView{"coif6"}) {
        auto coeffTemp = makeUnique<float[]>(n);

        std::copy(cbegin(coif6<float>), cend(coif6<float>), coeffTemp.get());
        for (auto i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        std::reverse_copy(coeffTemp.get(), coeffTemp.get() + n, lp1);
        qmfWrev({coeffTemp.get(), n}, hp1);
        std::copy(coeffTemp.get(), coeffTemp.get() + n, lp2);
        qmfEven({coeffTemp.get(), n}, hp2);

        return n;
    }

    if (name == StringView{"coif7"}) {
        auto coeffTemp = makeUnique<float[]>(n);

        std::copy(cbegin(coif7<float>), cend(coif7<float>), coeffTemp.get());
        for (auto i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        std::reverse_copy(coeffTemp.get(), coeffTemp.get() + n, lp1);
        qmfWrev({coeffTemp.get(), n}, hp1);
        std::copy(coeffTemp.get(), coeffTemp.get() + n, lp2);
        qmfEven({coeffTemp.get(), n}, hp2);

        return n;
    }

    if (name == StringView{"coif8"}) {
        auto coeffTemp = makeUnique<float[]>(n);

        std::copy(cbegin(coif8<float>), cend(coif8<float>), coeffTemp.get());
        for (auto i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        std::reverse_copy(coeffTemp.get(), coeffTemp.get() + n, lp1);
        qmfWrev({coeffTemp.get(), n}, hp1);
        std::copy(coeffTemp.get(), coeffTemp.get() + n, lp2);
        qmfEven({coeffTemp.get(), n}, hp2);

        return n;
    }

    if (name == StringView{"coif9"}) {
        auto coeffTemp = makeUnique<float[]>(n);

        std::copy(cbegin(coif9<float>), cend(coif9<float>), coeffTemp.get());
        for (auto i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        std::reverse_copy(coeffTemp.get(), coeffTemp.get() + n, lp1);
        qmfWrev({coeffTemp.get(), n}, hp1);
        std::copy(coeffTemp.get(), coeffTemp.get() + n, lp2);
        qmfEven({coeffTemp.get(), n}, hp2);

        return n;
    }

    if (name == StringView{"coif10"}) {
        auto coeffTemp = makeUnique<float[]>(n);

        std::copy(cbegin(coif10<float>), cend(coif10<float>), coeffTemp.get());
        for (auto i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        std::reverse_copy(coeffTemp.get(), coeffTemp.get() + n, lp1);
        qmfWrev({coeffTemp.get(), n}, hp1);
        std::copy(coeffTemp.get(), coeffTemp.get() + n, lp2);
        qmfEven({coeffTemp.get(), n}, hp2);

        return n;
    }
    if (name == StringView{"coif11"}) {
        auto coeffTemp = makeUnique<float[]>(n);

        std::copy(cbegin(coif11<float>), cend(coif11<float>), coeffTemp.get());
        for (auto i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        std::reverse_copy(coeffTemp.get(), coeffTemp.get() + n, lp1);
        qmfWrev({coeffTemp.get(), n}, hp1);
        std::copy(coeffTemp.get(), coeffTemp.get() + n, lp2);
        qmfEven({coeffTemp.get(), n}, hp2);

        return n;
    }
    if (name == StringView{"coif12"}) {
        auto coeffTemp = makeUnique<float[]>(n);

        std::copy(cbegin(coif12<float>), cend(coif12<float>), coeffTemp.get());
        for (auto i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        std::reverse_copy(coeffTemp.get(), coeffTemp.get() + n, lp1);
        qmfWrev({coeffTemp.get(), n}, hp1);
        std::copy(coeffTemp.get(), coeffTemp.get() + n, lp2);
        qmfEven({coeffTemp.get(), n}, hp2);

        return n;
    }
    if (name == StringView{"coif13"}) {
        auto coeffTemp = makeUnique<float[]>(n);

        std::copy(cbegin(coif13<float>), cend(coif13<float>), coeffTemp.get());
        for (auto i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        std::reverse_copy(coeffTemp.get(), coeffTemp.get() + n, lp1);
        qmfWrev({coeffTemp.get(), n}, hp1);
        std::copy(coeffTemp.get(), coeffTemp.get() + n, lp2);
        qmfEven({coeffTemp.get(), n}, hp2);

        return n;
    }
    if (name == StringView{"coif14"}) {
        auto coeffTemp = makeUnique<float[]>(n);

        std::copy(cbegin(coif14<float>), cend(coif14<float>), coeffTemp.get());
        for (auto i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        std::reverse_copy(coeffTemp.get(), coeffTemp.get() + n, lp1);
        qmfWrev({coeffTemp.get(), n}, hp1);
        std::copy(coeffTemp.get(), coeffTemp.get() + n, lp2);
        qmfEven({coeffTemp.get(), n}, hp2);

        return n;
    }
    if (name == StringView{"coif15"}) {
        auto coeffTemp = makeUnique<float[]>(n);

        std::copy(cbegin(coif15<float>), cend(coif15<float>), coeffTemp.get());
        for (auto i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        std::reverse_copy(coeffTemp.get(), coeffTemp.get() + n, lp1);
        qmfWrev({coeffTemp.get(), n}, hp1);
        std::copy(coeffTemp.get(), coeffTemp.get() + n, lp2);
        qmfEven({coeffTemp.get(), n}, hp2);

        return n;
    }
    if (name == StringView{"coif16"}) {
        auto coeffTemp = makeUnique<float[]>(n);

        std::copy(cbegin(coif16<float>), cend(coif16<float>), coeffTemp.get());
        for (auto i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        std::reverse_copy(coeffTemp.get(), coeffTemp.get() + n, lp1);
        qmfWrev({coeffTemp.get(), n}, hp1);
        std::copy(coeffTemp.get(), coeffTemp.get() + n, lp2);
        qmfEven({coeffTemp.get(), n}, hp2);

        return n;
    }
    if (name == StringView{"coif17"}) {
        auto coeffTemp = makeUnique<float[]>(n);

        std::copy(cbegin(coif17<float>), cend(coif17<float>), coeffTemp.get());
        for (auto i = 0; i < n; ++i) { coeffTemp[i] *= static_cast<float>(numbers::sqrt2); }

        std::reverse_copy(coeffTemp.get(), coeffTemp.get() + n, lp1);
        qmfWrev({coeffTemp.get(), n}, hp1);
        std::copy(coeffTemp.get(), coeffTemp.get() + n, lp2);
        qmfEven({coeffTemp.get(), n}, hp2);

        return n;
    }

    print("\n Filter Not in Database \n");
    return -1;
}

}  // namespace

Wavelet::Wavelet(char const* name)
    : _name{name}
    , _size{static_cast<std::size_t>(filterLength(name))}
    , _params{makeUnique<float[]>(4 * _size)}
{
    waveletFilterCoefficients(name, data(lpd()), data(hpd()), data(lpr()), data(hpr()));
}

auto Wavelet::size() const noexcept -> std::size_t { return _size; }

auto Wavelet::name() const noexcept -> String const& { return _name; }

auto Wavelet::lpd() const noexcept -> Span<float> { return {&_params[0], size()}; }

auto Wavelet::hpd() const noexcept -> Span<float> { return {&_params[size()], size()}; }

auto Wavelet::lpr() const noexcept -> Span<float>
{
    return {&_params[size() * 2U], size()};
}

auto Wavelet::hpr() const noexcept -> Span<float>
{
    return {&_params[size() * 3U], size()};
}

auto summary(Wavelet const& wavelet) -> String
{
    auto s = String{};
    fmt::format_to(std::back_inserter(s), "Wavelet: {0}\n", wavelet.name());
    fmt::format_to(std::back_inserter(s), "  Filters length: {0}\n", wavelet.size());

    fmt::format_to(std::back_inserter(s), "lpd: [{0}]\n", fmt::join(wavelet.lpd(), ", "));
    fmt::format_to(std::back_inserter(s), "hpd: [{0}]\n", fmt::join(wavelet.hpd(), ", "));
    fmt::format_to(std::back_inserter(s), "lpr: [{0}]\n", fmt::join(wavelet.lpr(), ", "));
    fmt::format_to(std::back_inserter(s), "hpr: [{0}]\n", fmt::join(wavelet.hpr(), ", "));
    return s;
}

}  // namespace mc::dsp
