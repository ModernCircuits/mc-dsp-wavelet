// SPDX-License-Identifier: BSL-1.0

#include <mc/dsp/fft.hpp>

#include <mc/core/format.hpp>

#include <catch2/catch_template_test_macros.hpp>

using namespace mc;

TEMPLATE_TEST_CASE("dsp/fft: Window", "[dsp][fft]", float, double)
{
    using T = TestType;

    auto allEqual = [](auto const& range, auto val) {
        return ranges::all_of(range, [val](auto v) { return v == val; });
    };

    auto allPositive = [](auto const& range) {
        return ranges::all_of(range, [](auto v) { return v >= T(0); });
    };

    auto allLessEqualOne = [](auto const& range) {
        return ranges::all_of(range, [](auto v) { return v <= T(1); });
    };

    auto size   = std::size_t(2048);
    auto window = Vector<T>(size, T(0));

    dsp::fillWindow(begin(window), end(window), dsp::WindowFunction::triangular, true);
    REQUIRE(allPositive(window));
    dsp::fillWindow(begin(window), end(window), dsp::WindowFunction::triangular, false);
    REQUIRE(allPositive(window));
    REQUIRE(allLessEqualOne(window));

    dsp::fillWindow(begin(window), end(window), dsp::WindowFunction::hann, true);
    REQUIRE(allPositive(window));
    dsp::fillWindow(begin(window), end(window), dsp::WindowFunction::hann, false);
    REQUIRE(allPositive(window));
    REQUIRE(allLessEqualOne(window));

    dsp::fillWindow(begin(window), end(window), dsp::WindowFunction::hamming, true);
    REQUIRE(allPositive(window));
    dsp::fillWindow(begin(window), end(window), dsp::WindowFunction::hamming, false);
    REQUIRE(allPositive(window));
    REQUIRE(allLessEqualOne(window));

    dsp::fillWindow(begin(window), end(window), dsp::WindowFunction::blackman, true);
    REQUIRE(allPositive(window));
    dsp::fillWindow(begin(window), end(window), dsp::WindowFunction::blackman, false);
    REQUIRE(allPositive(window));
    REQUIRE(allLessEqualOne(window));

    dsp::fillWindow(begin(window), end(window), dsp::WindowFunction::blackmanHarris, true);
    REQUIRE(allPositive(window));
    dsp::fillWindow(begin(window), end(window), dsp::WindowFunction::blackmanHarris, false);
    REQUIRE(allPositive(window));
    REQUIRE(allLessEqualOne(window));

    dsp::fillWindow(begin(window), end(window), dsp::WindowFunction::rectangular, true);
    REQUIRE(allPositive(window));
    REQUIRE(allLessEqualOne(window));
    dsp::fillWindow(begin(window), end(window), dsp::WindowFunction::rectangular, false);
    REQUIRE(allPositive(window));
    REQUIRE(allEqual(window, T(1.0)));

    auto signal = Vector<T>(size, T(0.5));
    dsp::multiplyWithWindow(begin(signal), end(signal), cbegin(window), cend(window));
    REQUIRE(allPositive(signal));
    REQUIRE(allEqual(signal, T(0.5)));

    REQUIRE(toString(dsp::WindowFunction::rectangular) == "Rectangular");
    REQUIRE(toString(dsp::WindowFunction::triangular) == "Triangular");
    REQUIRE(toString(dsp::WindowFunction::hann) == "Hann");
    REQUIRE(toString(dsp::WindowFunction::hamming) == "Hamming");
    REQUIRE(toString(dsp::WindowFunction::blackman) == "Blackman");
    REQUIRE(toString(dsp::WindowFunction::blackmanHarris) == "Blackman-Harris");
}
