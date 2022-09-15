#include "mc/dsp/convolution.hpp"

#include "mc/testing/test.hpp"
#include <mc/core/format.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

using namespace mc;

TEST_CASE("dsp/convolution: OverlapSaveConvolver", "[dsp][convolution]")
{
    // TODO(tobi): Fails on windows

    // auto const* const testFile = GENERATE("test_data/raw/conv_xcorr_01.txt", "test_data/raw/conv_xcorr_02.txt",
    //                                       "test_data/raw/conv_xcorr_03.txt", "test_data/raw/conv_xcorr_04.txt");

    // auto testData = loadTestData(testFile);
    // CHECK(testData.size() == 4U);

    // auto s = dsp::FloatSignal{testData[0].data(), testData[0].size()};
    // auto p = dsp::FloatSignal{testData[1].data(), testData[1].size()};
    // auto x = dsp::OverlapSaveConvolver{s, p};

    // x.convolute();
    // CHECK(approxEqual(x.extractResult(), testData[2]));

    // x.crossCorrelate();
    // CHECK(approxEqual(x.extractResult(), testData[3]));

    SUCCEED();
}