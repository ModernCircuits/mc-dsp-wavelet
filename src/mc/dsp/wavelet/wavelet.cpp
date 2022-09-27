// SPDX-License-Identifier: BSL-1.0

#include "wavelet.hpp"

#include <mc/dsp/wavelet/family.hpp>

#include <mc/core/algorithm.hpp>
#include <mc/core/cmath.hpp>
#include <mc/core/cstring.hpp>
#include <mc/core/exception.hpp>
#include <mc/core/iterator.hpp>
#include <mc/core/memory.hpp>
#include <mc/core/numbers.hpp>
#include <mc/core/print.hpp>
#include <mc/core/ranges.hpp>
#include <mc/core/stdexcept.hpp>
#include <mc/core/string_view.hpp>
#include <mc/core/utility.hpp>

namespace mc {

template<typename T>
static auto qmfEven(Span<T const> in, T* out)
{
    auto const n = in.size();
    for (std::size_t count = 0; count < n; ++count) {
        auto evenIndex = count % 2 == 0;
        out[count]     = in[n - count - 1] * (evenIndex ? T{1} : T{-1});
    }
}

template<typename T>
static auto qmfWrev(Span<T const> in, T* out)
{
    qmfEven(in, out);
    ranges::reverse(Span<T>{out, in.size()});
}

static auto filterLength(StringView name) -> std::size_t
{
    auto const& filters = allWavelets<float>;
    auto const filter   = ranges::find(filters, name, &WaveletCoefficients<float>::name);
    if (filter == ranges::end(filters)) {
        raise<InvalidArgument>("wavelet filter not in database");
    }
    return filter->length;
}

static auto fillDaubechiesWaveletCoefficients(
    StringView name,
    Span<float> lp1,
    Span<float> hp1,
    Span<float> lp2,
    Span<float> hp2
) -> size_t
{
    using namespace std::string_view_literals;

    if (name == "haar"sv) {
        return fillDaubechiesWaveletCoefficients("db1"sv, lp1, hp1, lp2, hp2);
    }

    auto const& filters = allWavelets<float>;
    auto const filter   = ranges::find(filters, name, &WaveletCoefficients<float>::name);
    if (filter == ranges::end(filters)) {
        raise<InvalidArgument>("wavelet filter not in database");
    }

    ranges::reverse_copy(filter->coefficients, ranges::begin(lp1));
    qmfWrev(filter->coefficients, mc::data(hp1));

    ranges::copy(filter->coefficients, ranges::begin(lp2));
    qmfEven(filter->coefficients, mc::data(hp2));

    return filter->length;
}

static auto fillCoifWaveletCoefficients(
    StringView name,
    Span<float> lp1,
    Span<float> hp1,
    Span<float> lp2,
    Span<float> hp2
) -> size_t
{
    auto const& filters = coifWavelets<float>;
    auto const filter   = ranges::find(filters, name, &WaveletCoefficients<float>::name);
    if (filter == ranges::end(filters)) {
        raise<InvalidArgument>("wavelet filter not in database");
    }

    auto scale  = [](auto c) { return c * static_cast<float>(numbers::sqrt2); };
    auto scaled = filter->coefficients | ranges::views::transform(scale);

    ranges::copy(scaled, ranges::begin(lp2));
    ranges::reverse_copy(scaled, ranges::begin(lp1));
    qmfWrev<float>(lp2, mc::data(hp1));
    qmfEven<float>(lp2, mc::data(hp2));

    return filter->length;
}

static auto fillWaveletFilterCoefficients(
    StringView name,
    Span<float> lp1,
    Span<float> hp1,
    Span<float> lp2,
    Span<float> hp2
) -> std::size_t
{
    if (StringView{name}.find("haar") != StringView::npos) {
        return fillDaubechiesWaveletCoefficients(name, lp1, hp1, lp2, hp2);
    }
    if (StringView{name}.find("db") != StringView::npos) {
        return fillDaubechiesWaveletCoefficients(name, lp1, hp1, lp2, hp2);
    }
    if (StringView{name}.find("sym") != StringView::npos) {
        return fillDaubechiesWaveletCoefficients(name, lp1, hp1, lp2, hp2);
    }
    if (StringView{name}.find("coif") != StringView::npos) {
        return fillCoifWaveletCoefficients(name, lp1, hp1, lp2, hp2);
    }

    raisef<InvalidArgument>("filter not in database: {}", name);
    return -1;
}

Wavelet::Wavelet(StringView name)
    : _name{name}
    , _size{static_cast<std::size_t>(filterLength(name))}
    , _params(_size * 4U)
{
    // We can't use the member functions to access the coefficients,
    // because they return a const Span.
    auto lp1 = Span<float>{&_params[0], size()};
    auto hp1 = Span<float>{&_params[size() * 1U], size()};
    auto lp2 = Span<float>{&_params[size() * 2U], size()};
    auto hp2 = Span<float>{&_params[size() * 3U], size()};
    fillWaveletFilterCoefficients(name, lp1, hp1, lp2, hp2);
}

auto Wavelet::size() const noexcept -> std::size_t { return _size; }

auto Wavelet::name() const noexcept -> String const& { return _name; }

auto Wavelet::lpd() const noexcept -> Span<float const> { return {&_params[0], size()}; }

auto Wavelet::hpd() const noexcept -> Span<float const>
{
    return {&_params[size()], size()};
}

auto Wavelet::lpr() const noexcept -> Span<float const>
{
    return {&_params[size() * 2U], size()};
}

auto Wavelet::hpr() const noexcept -> Span<float const>
{
    return {&_params[size() * 3U], size()};
}

}  // namespace mc
