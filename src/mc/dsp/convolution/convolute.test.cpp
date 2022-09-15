#include <mc/dsp/convolution.hpp>

#include <mc/core/format.hpp>
#include <mc/testing/test.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

using namespace mc;

template<typename T>
auto testConvolute(std::vector<std::vector<T>> const& testData) -> bool
{
    auto const& signal   = testData[0];
    auto const& patch    = testData[1];
    auto const& expected = testData[2];

    auto output = std::vector<T>(expected.size());
    dsp::convolute(
        mc::data(signal),
        mc::size(signal),
        mc::data(patch),
        mc::size(patch),
        mc::data(output)
    );
    CHECK(approxEqual(output, expected));
    return true;
}

TEST_CASE("dsp/convolution: convolute", "[dsp][convolution]")
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
