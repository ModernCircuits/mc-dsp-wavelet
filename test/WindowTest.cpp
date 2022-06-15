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
    MC_REQUIRE(allPositive(window));
    dsp::fillWindow(begin(window), end(window), dsp::WindowFunction::triangular, false);
    MC_REQUIRE(allPositive(window));
    MC_REQUIRE(allLessEqualOne(window));

    dsp::fillWindow(begin(window), end(window), dsp::WindowFunction::hann, true);
    MC_REQUIRE(allPositive(window));
    dsp::fillWindow(begin(window), end(window), dsp::WindowFunction::hann, false);
    MC_REQUIRE(allPositive(window));
    MC_REQUIRE(allLessEqualOne(window));

    dsp::fillWindow(begin(window), end(window), dsp::WindowFunction::hamming, true);
    MC_REQUIRE(allPositive(window));
    dsp::fillWindow(begin(window), end(window), dsp::WindowFunction::hamming, false);
    MC_REQUIRE(allPositive(window));
    MC_REQUIRE(allLessEqualOne(window));

    dsp::fillWindow(begin(window), end(window), dsp::WindowFunction::blackman, true);
    MC_REQUIRE(allPositive(window));
    dsp::fillWindow(begin(window), end(window), dsp::WindowFunction::blackman, false);
    MC_REQUIRE(allPositive(window));
    MC_REQUIRE(allLessEqualOne(window));

    dsp::fillWindow(begin(window), end(window), dsp::WindowFunction::blackmanHarris, true);
    MC_REQUIRE(allPositive(window));
    dsp::fillWindow(begin(window), end(window), dsp::WindowFunction::blackmanHarris, false);
    MC_REQUIRE(allPositive(window));
    MC_REQUIRE(allLessEqualOne(window));

    dsp::fillWindow(begin(window), end(window), dsp::WindowFunction::rectangular, true);
    MC_REQUIRE(allPositive(window));
    MC_REQUIRE(allLessEqualOne(window));
    dsp::fillWindow(begin(window), end(window), dsp::WindowFunction::rectangular, false);
    MC_REQUIRE(allPositive(window));
    MC_REQUIRE(allEqual(window, T(1.0)));

    auto signal = std::vector<T>(size, T(0.5));
    dsp::multiplyWithWindow(begin(signal), end(signal), cbegin(window), cend(window));
    MC_REQUIRE(allPositive(signal));
    MC_REQUIRE(allEqual(signal, T(0.5)));

    MC_REQUIRE(toString(dsp::WindowFunction::rectangular) == "Rectangular");
    MC_REQUIRE(toString(dsp::WindowFunction::triangular) == "Triangular");
    MC_REQUIRE(toString(dsp::WindowFunction::hann) == "Hann");
    MC_REQUIRE(toString(dsp::WindowFunction::hamming) == "Hamming");
    MC_REQUIRE(toString(dsp::WindowFunction::blackman) == "Blackman");
    MC_REQUIRE(toString(dsp::WindowFunction::blackmanHarris) == "Blackman-Harris");

    return true;
}

auto main() -> int
{
    testWindow<float>();
    testWindow<double>();
    testWindow<long double>();
    return EXIT_SUCCESS;
}