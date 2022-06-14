#include "mc/dsp/convolution.hpp"

#include "mc/format.hpp"

#include "mc/testing/test.hpp"

namespace dsp = mc::dsp;

auto main() -> int
{
    auto const testFiles = std::vector<char const*>{
        "testData/conv_xcorr_01.txt",
        "testData/conv_xcorr_02.txt",
        "testData/conv_xcorr_03.txt",
        "testData/conv_xcorr_04.txt",
    };

    for (auto const* testFile : testFiles)
    {
        fmt::print("Testing: {0} ...\n", testFile);
        auto testData = loadTestData(testFile);
        REQUIRE(testData.size() == 4U);

        auto s = dsp::FloatSignal{testData[0].data(), testData[0].size()};
        auto p = dsp::FloatSignal{testData[1].data(), testData[1].size()};
        auto x = dsp::OverlapSaveConvolver{s, p};

        x.convolute();
        REQUIRE(approxEqual(x.extractResult(), testData[2]));

        x.crossCorrelate();
        REQUIRE(approxEqual(x.extractResult(), testData[3]));

        fmt::print("Testing: {0} done\n", testFile);
    }

    return EXIT_SUCCESS;
}