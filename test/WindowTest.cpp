#include "mc/dsp/fft.hpp"

#include "mc/format.hpp"
#include "mc/testing/test.hpp"

namespace dsp = mc::dsp;

template<typename T>
auto testWindow() -> bool
{
    auto allEqual = [](auto const& range, auto val)
    {
        auto equals = [val](auto v) { return v == val; };
        return std::all_of(cbegin(range), cend(range), equals);
    };

    auto allPositive = [](auto const& range)
    {
        auto positive = [](auto v) { return v >= T(0); };
        return std::all_of(cbegin(range), cend(range), positive);
    };

    auto allLessEqualOne = [](auto const& range)
    {
        auto lessEqual = [](auto v) { return v <= T(1); };
        return std::all_of(cbegin(range), cend(range), lessEqual);
    };

    auto size   = std::size_t(2048);
    auto window = std::vector<T>(size, T(0));

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

    auto signal = std::vector<T>(size, T(0.5));
    dsp::multiplyWithWindow(begin(signal), end(signal), cbegin(window), cend(window));
    REQUIRE(allPositive(signal));
    REQUIRE(allEqual(signal, T(0.5)));

    REQUIRE(toString(dsp::WindowFunction::rectangular) == "Rectangular");
    REQUIRE(toString(dsp::WindowFunction::triangular) == "Triangular");
    REQUIRE(toString(dsp::WindowFunction::hann) == "Hann");
    REQUIRE(toString(dsp::WindowFunction::hamming) == "Hamming");
    REQUIRE(toString(dsp::WindowFunction::blackman) == "Blackman");
    REQUIRE(toString(dsp::WindowFunction::blackmanHarris) == "Blackman-Harris");

    return true;
}

auto main() -> int
{
    testWindow<float>();
    testWindow<double>();
    testWindow<long double>();
    return EXIT_SUCCESS;
}