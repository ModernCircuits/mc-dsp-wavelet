#include "lt/dsp/convolution.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#define REQUIRE(e)                                                                         \
    do {                                                                                   \
        if (!(e)) {                                                                        \
            std::cerr << "TEST ASSERTION FAILED: " << __FILE__ << ":" << __LINE__ << '\n'; \
            std::exit(EXIT_FAILURE);                                                       \
        }                                                                                  \
    } while (false)

template <typename It1, typename It2>
auto approxEqual(It1 f1, It1 l1, It2 f2, It2 l2, int epsilonFactor = 4) -> bool
{
    // This function only makes sense when comparing two ranges
    // of the same floating-point type.
    using v1_t = typename std::iterator_traits<It1>::value_type;
    using v2_t = typename std::iterator_traits<It2>::value_type;
    static_assert(std::is_same_v<v1_t, v2_t>);

    auto epsilonEqual = [epsilonFactor](auto l, auto r) {
        auto const epsilon = std::numeric_limits<v1_t>::epsilon();
        return std::fabs(l - r) < epsilon * epsilonFactor;
    };

    return std::equal(f1, l1, f2, l2, epsilonEqual);
}

template <typename Container1, typename Container2>
auto approxEqual(Container1 c1, Container2 c2, int epsilonFactor = 4) -> bool
{
    return approxEqual(std::begin(c1), std::end(c1), std::begin(c2), std::end(c2), epsilonFactor);
}

using ConvolutionTestData = std::vector<std::vector<double>>;

auto split(std::string const& s, char delim) -> std::vector<std::string>
{
    auto result = std::vector<std::string> {};
    auto ss = std::stringstream(s);
    auto item = std::string {};

    while (std::getline(ss, item, delim)) {
        result.push_back(item);
    }
    return result;
}

auto loadConvolutionTestData(char const* filePath) -> ConvolutionTestData
{
    auto parseLine = [](auto const& line) {
        auto splits = split(line, ' ');
        auto values = std::vector<double> {};
        for (auto const& s : splits) {
            values.push_back(std::stod(s));
        }
        return values;
    };

    auto file = std::fstream { filePath, std::ios::in };
    auto tmp = std::string {};
    auto result = ConvolutionTestData {};

    if (file.is_open()) {
        while (std::getline(file, tmp)) {
            result.push_back(parseLine(tmp));
        }
        file.close();
    }

    return result;
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
        std::cout << "Testing: " << testFile << "...\n";
        auto testData = loadConvolutionTestData(testFile);
        REQUIRE(testData.size() == 4U);

        auto s = DoubleSignal { testData[0].data(), testData[0].size() };
        auto p = DoubleSignal { testData[1].data(), testData[1].size() };
        auto x = OverlapSave { s, p };

        x.convolute();
        REQUIRE(approxEqual(x.extractResult(), testData[2]));

        x.crossCorrelate();
        REQUIRE(approxEqual(x.extractResult(), testData[3]));

        std::cout << "Testing: " << testFile << " done\n";
    }

    return EXIT_SUCCESS;
}