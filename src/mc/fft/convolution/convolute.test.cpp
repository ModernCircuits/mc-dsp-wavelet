// SPDX-License-Identifier: BSL-1.0

#include <mc/fft/convolution.hpp>
#include <mc/testing/test.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

using namespace mc;

template<typename T>
auto testConvolute(Vector<Vector<T>> const& testData) -> bool
{
    auto const& signal   = testData[0];
    auto const& patch    = testData[1];
    auto const& expected = testData[2];

    auto output = Vector<T>(expected.size());
    convolute<float>(signal, patch, data(output));
    CHECK(approxEqual<T>(output, expected));
    return true;
}

TEST_CASE("fft/convolution: convolute", "[fft][convolution]")
{
    auto const* const testFile = GENERATE(
        "test_data/raw/conv_xcorr_01.txt",
        "test_data/raw/conv_xcorr_02.txt",
        "test_data/raw/conv_xcorr_03.txt",
        "test_data/raw/conv_xcorr_04.txt"
    );

    auto const testData = loadTestData(testFile);
    CHECK(testData.size() == 4U);
    // CHECK(testConvolute<double>(testData));
    CHECK(testConvolute<float>(toFloat(testData)));
}
