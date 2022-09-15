#include <mc/dsp/wavelets.hpp>

#include <mc/core/algorithm.hpp>
#include <mc/core/cmath.hpp>
#include <mc/core/cstdlib.hpp>
#include <mc/core/cstring.hpp>
#include <mc/core/format.hpp>
#include <mc/core/memory.hpp>
#include <mc/core/random.hpp>
#include <mc/core/sstream.hpp>
#include <mc/core/vector.hpp>
#include <mc/testing/test.hpp>

#include <catch2/catch_test_macros.hpp>

using namespace mc;

TEST_CASE("dsp/wavelet: dbCoefTests", "[dsp][wavelet]")
{
    constexpr auto epsilon = 1e-6;
    auto waveletNames      = std::vector<std::string>(38);
    std::generate(begin(waveletNames), end(waveletNames), [i = 1]() mutable {
        return std::string("db") + std::to_string(i);
        ++i;
    });

    for (auto const& name : waveletNames) {
        auto obj = dsp::Wavelet{name.c_str()};
        auto t1  = sum1(obj.lpr().data(), obj.lpr().size()) - std::sqrt(2.0F);
        auto t2  = sum2(obj.lpr().data(), obj.lpr().size()) - 1.0F / std::sqrt(2.0F);
        auto t3  = sum3(obj.lpr().data(), obj.lpr().size()) - 1.0F / std::sqrt(2.0F);
        auto t4  = sum4(obj.lpr().data(), obj.lpr().size()) - 1.0F;

        REQUIRE(fabs(t1) <= epsilon);
        REQUIRE(fabs(t2) <= epsilon);
        REQUIRE(fabs(t3) <= epsilon);
        REQUIRE(fabs(t4) <= epsilon);

        for (std::size_t m = 1; m < (obj.lpr().size() / 2) - 1; m++) {
            auto t5 = sum5(obj.lpr().data(), obj.lpr().size(), m);
            REQUIRE(fabs(t5) <= epsilon);
        }
    }
}

TEST_CASE("dsp/wavelet: coifCoefTests", "[dsp][wavelet]")
{
    float epsilon = 1e-6;
    float t1      = NAN;
    float t2      = NAN;
    float t3      = NAN;
    float t4      = NAN;
    float t5      = NAN;
    std::vector<std::string> waveletNames;
    waveletNames.resize(17);
    for (std::size_t i = 0; i < waveletNames.size(); i++) {
        waveletNames[i] = std::string("coif") + std::to_string(i + 1);
    }

    for (auto const& name : waveletNames) {
        auto obj = dsp::Wavelet{name.c_str()};
        t1       = sum1(obj.lpr().data(), obj.lpr().size()) - std::sqrt(2.0F);
        t2       = sum2(obj.lpr().data(), obj.lpr().size()) - 1.0F / std::sqrt(2.0F);
        t3       = sum3(obj.lpr().data(), obj.lpr().size()) - 1.0F / std::sqrt(2.0F);
        t4       = sum4(obj.lpr().data(), obj.lpr().size()) - 1.0F;

        REQUIRE(std::fabs(t1) <= epsilon);
        REQUIRE(std::fabs(t2) <= epsilon);
        REQUIRE(std::fabs(t3) <= epsilon);
        REQUIRE(std::fabs(t4) <= epsilon);
        for (std::size_t m = 1; m < (obj.lpr().size() / 2) - 1; m++) {
            t5 = sum5(obj.lpr().data(), obj.lpr().size(), m);
            REQUIRE(std::fabs(t5) <= epsilon);
        }
    }
}

TEST_CASE("dsp/wavelet: symCoefTests", "[dsp][wavelet]")
{
    float epsilon = 1e-6;
    float t1      = NAN;
    float t2      = NAN;
    float t3      = NAN;
    float t4      = NAN;
    float t5      = NAN;
    std::vector<std::string> waveletNames;
    for (std::size_t i = 1; i < 20; i++) {
        waveletNames.push_back(std::string("sym") + std::to_string(i + 1));
    }

    for (auto const& name : waveletNames) {
        auto obj = dsp::Wavelet{name.c_str()};
        t1       = sum1(obj.lpr().data(), obj.lpr().size()) - std::sqrt(2.0F);
        t2       = sum2(obj.lpr().data(), obj.lpr().size()) - 1.0F / std::sqrt(2.0F);
        t3       = sum3(obj.lpr().data(), obj.lpr().size()) - 1.0F / std::sqrt(2.0F);
        t4       = sum4(obj.lpr().data(), obj.lpr().size()) - 1.0F;

        REQUIRE(std::fabs(t1) <= epsilon);
        REQUIRE(std::fabs(t2) <= epsilon);
        REQUIRE(std::fabs(t3) <= epsilon);
        REQUIRE(std::fabs(t4) <= epsilon);

        for (std::size_t m = 1; m < (obj.lpr().size() / 2) - 1; m++) {
            t5 = sum5(obj.lpr().data(), obj.lpr().size(), m);
            REQUIRE(std::fabs(t5) <= epsilon);
        }
    }
}

TEST_CASE("dsp/wavelet: biorCoefTests", "[dsp][wavelet]")
{
    constexpr float epsilon = 1e-6;
    std::vector<std::string> waveletNames;
    waveletNames.emplace_back("bior1.1");
    waveletNames.emplace_back("bior1.3");
    waveletNames.emplace_back("bior1.5");
    waveletNames.emplace_back("bior2.2");
    waveletNames.emplace_back("bior2.4");
    waveletNames.emplace_back("bior2.6");
    waveletNames.emplace_back("bior2.8");
    waveletNames.emplace_back("bior3.1");
    waveletNames.emplace_back("bior3.3");
    waveletNames.emplace_back("bior3.5");
    waveletNames.emplace_back("bior3.7");
    waveletNames.emplace_back("bior3.9");
    waveletNames.emplace_back("bior4.4");
    waveletNames.emplace_back("bior5.5");
    waveletNames.emplace_back("bior6.8");

    for (auto const& name : waveletNames) {
        auto obj = dsp::Wavelet{name.c_str()};

        auto const t1 = sum1(obj.lpr().data(), obj.lpr().size()) - std::sqrt(2.0F);
        auto const t2 = sum1(obj.lpd().data(), obj.lpd().size()) - std::sqrt(2.0F);

        auto const t3 = sum2(obj.lpr().data(), obj.lpr().size()) - 1.0F / std::sqrt(2.0F);
        auto const t4 = sum2(obj.lpd().data(), obj.lpd().size()) - 1.0F / std::sqrt(2.0F);

        auto const t5 = sum3(obj.lpr().data(), obj.lpr().size()) - 1.0F / std::sqrt(2.0F);
        auto const t6 = sum3(obj.lpd().data(), obj.lpd().size()) - 1.0F / std::sqrt(2.0F);

        REQUIRE(fabs(t1) <= epsilon);
        REQUIRE(fabs(t2) <= epsilon);
        REQUIRE(fabs(t3) <= epsilon);
        REQUIRE(fabs(t4) <= epsilon);
        REQUIRE(fabs(t5) <= epsilon);
        REQUIRE(fabs(t6) <= epsilon);
    }
}

TEST_CASE("dsp/wavelet: rbiorCoefTests", "[dsp][wavelet]")
{
    auto const epsilon = 1e-6F;

    std::vector<std::string> waveletNames;
    waveletNames.emplace_back("rbior1.1");
    waveletNames.emplace_back("rbior1.3");
    waveletNames.emplace_back("rbior1.5");
    waveletNames.emplace_back("rbior2.2");
    waveletNames.emplace_back("rbior2.4");
    waveletNames.emplace_back("rbior2.6");
    waveletNames.emplace_back("rbior2.8");
    waveletNames.emplace_back("rbior3.1");
    waveletNames.emplace_back("rbior3.3");
    waveletNames.emplace_back("rbior3.5");
    waveletNames.emplace_back("rbior3.7");
    waveletNames.emplace_back("rbior3.9");
    waveletNames.emplace_back("rbior4.4");
    waveletNames.emplace_back("rbior5.5");
    waveletNames.emplace_back("rbior6.8");

    for (auto const& name : waveletNames) {
        auto obj = dsp::Wavelet{name.c_str()};

        auto const t1 = sum1(obj.lpr().data(), obj.lpr().size()) - std::sqrt(2.0F);
        auto const t2 = sum1(obj.lpd().data(), obj.lpd().size()) - std::sqrt(2.0F);

        auto const t3 = sum2(obj.lpr().data(), obj.lpr().size()) - 1.0F / std::sqrt(2.0F);
        auto const t4 = sum2(obj.lpd().data(), obj.lpd().size()) - 1.0F / std::sqrt(2.0F);

        auto const t5 = sum3(obj.lpr().data(), obj.lpr().size()) - 1.0F / std::sqrt(2.0F);
        auto const t6 = sum3(obj.lpd().data(), obj.lpd().size()) - 1.0F / std::sqrt(2.0F);

        REQUIRE(std::fabs(t1) <= epsilon);
        REQUIRE(std::fabs(t2) <= epsilon);
        REQUIRE(std::fabs(t3) <= epsilon);
        REQUIRE(std::fabs(t4) <= epsilon);
        REQUIRE(std::fabs(t5) <= epsilon);
        REQUIRE(std::fabs(t6) <= epsilon);
    }
}
