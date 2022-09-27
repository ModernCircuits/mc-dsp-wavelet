// SPDX-License-Identifier: BSL-1.0

#include <mc/dsp/fft.hpp>

#include <mc/core/iterator.hpp>
#include <mc/testing/test.hpp>

#include <catch2/catch_test_macros.hpp>

using namespace mc;

struct FFTRealToComplexTestCase
{
    Vector<float> input{};
    Vector<Complex<float>> expected{};
};

auto readFFTRealToComplexTestData(char const* path) -> FFTRealToComplexTestCase
{
    auto rawData = loadTestData(path);
    auto result  = FFTRealToComplexTestCase{};
    result.input = rawData[0];
    for (auto i = size_t{0}; i < rawData[1].size(); ++i) {
        result.expected.emplace_back(rawData[1][i], rawData[2][i]);
    }
    return result;
}

TEST_CASE("dsp/fft: rfft", "[dsp][fft]")
{

    auto const testFiles = Vector<char const*>{
        "test_data/raw/fft_real_to_complex_01.txt",
    };

    auto closeEnough = [](auto l, auto r) -> bool {
        auto const epsilon = std::numeric_limits<float>::epsilon();
        auto const re      = std::abs(l.real() - r.real()) < epsilon * 4;
        auto const im      = std::abs(l.imag() - r.imag()) < epsilon * 4;
        return re && im;
    };

    for (auto const* testFile : testFiles) {
        auto testCase = readFFTRealToComplexTestData(testFile);
        auto fft      = makeRFFT(testCase.input.size());
        auto output   = Vector<Complex<float>>(testCase.expected.size());

        rfft(fft, testCase.input, output);
        REQUIRE(ranges::equal(output, testCase.expected, closeEnough));
    }
}
