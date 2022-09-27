// SPDX-License-Identifier: BSL-1.0

#include <mc/fft/convolution.hpp>

#include <mc/testing/test.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

using namespace mc;

TEST_CASE("fft/convolution: OverlapSaveConvolver", "[fft][convolution]")
{

    auto const* const testFile = GENERATE(
        "test_data/raw/conv_xcorr_01.txt",
        "test_data/raw/conv_xcorr_02.txt",
        "test_data/raw/conv_xcorr_03.txt",
        "test_data/raw/conv_xcorr_04.txt"
    );

    auto testData = loadTestData(testFile);
    CHECK(testData.size() == 4U);

    auto s = FloatSignal{testData[0].data(), testData[0].size()};
    auto p = FloatSignal{testData[1].data(), testData[1].size()};
    auto x = OverlapSaveConvolver{s, p};

    x.convolute();
    CHECK(approxEqual<float>(x.extractResult(), testData[2]));

    x.crossCorrelate();
    CHECK(approxEqual<float>(x.extractResult(), testData[3]));
}
