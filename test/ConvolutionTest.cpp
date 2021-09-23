#include "lt/dsp/convolution.hpp"

#include "testing.hpp"

auto testDirectConvolve(std::vector<std::vector<double>> const& testData) -> bool
{
    auto const& signal = testData[0];
    auto const& patch = testData[1];
    auto const& expected = testData[2];

    auto output = std::vector<double>(expected.size());
    Convolution::direct(data(signal), size(signal), data(patch), size(patch), data(output));
    REQUIRE(approxEqual(output, expected));
    return true;
}

auto main() -> int
{
    auto const testFiles = std::vector<char const*> {
        "testData/conv_xcorr_01.txt",
        "testData/conv_xcorr_02.txt",
        "testData/conv_xcorr_03.txt",
        "testData/conv_xcorr_04.txt",
    };

    for (auto const* testFile : testFiles) {
        std::cout << "Testing: " << testFile << " ...\n";
        auto const testData = loadConvolutionTestData(testFile);
        REQUIRE(testData.size() == 4U);
        REQUIRE(testDirectConvolve(testData));
        std::cout << "Testing: " << testFile << " done\n";
    }

    return EXIT_SUCCESS;
}