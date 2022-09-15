#include "mc/dsp/fft.hpp"

#include <mc/core/format.hpp>
#include <mc/core/iterator.hpp>

#include "mc/testing/test.hpp"

#include <catch2/catch_test_macros.hpp>

using namespace mc;

struct FFTRealToComplexTestCase
{
    std::vector<float> input{};
    std::vector<Complex<float>> expected{};
};

auto readFFTRealToComplexTestData(char const* path) -> FFTRealToComplexTestCase
{
    auto rawData = loadTestData(path);
    auto result  = FFTRealToComplexTestCase{};
    result.input = rawData[0];
    for (auto i = std::size_t{0}; i < rawData[1].size(); ++i)
    {
        result.expected.emplace_back(rawData[1][i], rawData[2][i]);
    }
    return result;
}

TEST_CASE("dsp/fft: realToComplex", "[dsp][fft]")
{

    auto const testFiles = std::vector<char const*>{
        "test_data/raw/fft_real_to_complex_01.txt",
    };

    auto closeEnough = [](auto l, auto r) -> bool
    {
        auto const epsilon = std::numeric_limits<float>::epsilon();
        auto const re      = std::fabs(l.real() - r.real()) < epsilon * 4;
        auto const im      = std::fabs(l.imag() - r.imag()) < epsilon * 4;
        return re && im;
    };

    for (auto const* testFile : testFiles)
    {
        auto testCase = readFFTRealToComplexTestData(testFile);

        auto fft    = dsp::RFFT{static_cast<int>(testCase.input.size()), dsp::FFTDirection::forward};
        auto output = std::vector<Complex<float>>(testCase.expected.size());
        fft.performRealToComplex(mc::data(testCase.input), mc::data(output));
        REQUIRE(std::equal(begin(output), end(output), begin(testCase.expected), closeEnough));
    }
}