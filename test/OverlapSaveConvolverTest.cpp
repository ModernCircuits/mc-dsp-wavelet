#include "lt/dsp/convolution.hpp"

#include "lt/testing/test.hpp"

auto main() -> int
{
    auto const testFiles = std::vector<char const*> {
        "testData/conv_xcorr_01.txt",
        "testData/conv_xcorr_02.txt",
        "testData/conv_xcorr_03.txt",
        "testData/conv_xcorr_04.txt",
    };

    for (auto const* testFile : testFiles) {
        std::cout << "Testing: " << testFile << "...\n";
        auto testData = loadTestData(testFile);
        REQUIRE(testData.size() == 4U);

        auto s = DoubleSignal { testData[0].data(), testData[0].size() };
        auto p = DoubleSignal { testData[1].data(), testData[1].size() };
        auto x = OverlapSaveConvolver { s, p };

        x.convolute();
        REQUIRE(approxEqual(x.extractResult(), testData[2]));

        x.crossCorrelate();
        REQUIRE(approxEqual(x.extractResult(), testData[3]));

        std::cout << "Testing: " << testFile << " done\n";
    }

    return EXIT_SUCCESS;
}