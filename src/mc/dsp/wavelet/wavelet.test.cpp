// SPDX-License-Identifier: BSL-1.0

#include <mc/dsp/wavelet.hpp>

#include <mc/core/algorithm.hpp>
#include <mc/core/cmath.hpp>
#include <mc/core/cstdlib.hpp>
#include <mc/core/cstring.hpp>
#include <mc/core/memory.hpp>
#include <mc/core/random.hpp>
#include <mc/core/sstream.hpp>
#include <mc/core/vector.hpp>
#include <mc/testing/test.hpp>

#include <catch2/catch_test_macros.hpp>

using namespace mc;

static auto sum1(Span<float const> array) -> float
{
    auto sum = 0.0F;
    for (auto i = size_t{0}; i < size(array); ++i) { sum += array[i]; }
    return sum;
}

static auto sum2(Span<float const> array) -> float
{
    auto sum = 0.0F;
    for (size_t i = 0; i < size(array); i += 2) { sum += array[i]; }
    return sum;
}

static auto sum3(Span<float const> array) -> float
{
    auto sum = 0.0F;
    for (size_t i = 1; i < size(array); i += 2) { sum += array[i]; }
    return sum;
}

// np.sum(w[2*m:(2*N+2*m)]*w[0:2*N])
static auto sum4(Span<float const> array) -> float
{
    auto sum = 0.0F;
    for (size_t i = 0; i < size(array); i += 1) { sum += array[i] * array[i]; }
    return sum;
}

// np.sum(w[2 * m:(2 * N)] * w[0:2 * N - 2 * m])
static auto sum5(float const* array, size_t n, size_t m) -> float
{
    auto sum = 0.0F;
    for (size_t i = 2 * m; i < n; i += 1) { sum += array[i] * array[i - 2 * m]; }
    return sum;
}

static constexpr auto epsilon = 1e-6F;

TEST_CASE("dsp/wavelet: dbCoefTests", "[dsp][wavelet]")
{
    auto waveletNames = Vector<String>(38);
    std::generate(begin(waveletNames), end(waveletNames), [i = 1]() mutable {
        return String("db") + std::to_string(i);
        ++i;
    });

    for (auto const& name : waveletNames) {
        auto obj = Wavelet{name.c_str()};
        auto t1  = sum1(obj.lpr()) - mc::sqrt(2.0F);
        auto t2  = sum2(obj.lpr()) - 1.0F / mc::sqrt(2.0F);
        auto t3  = sum3(obj.lpr()) - 1.0F / mc::sqrt(2.0F);
        auto t4  = sum4(obj.lpr()) - 1.0F;

        REQUIRE(std::abs(t1) <= epsilon);
        REQUIRE(std::abs(t2) <= epsilon);
        REQUIRE(std::abs(t3) <= epsilon);
        REQUIRE(std::abs(t4) <= epsilon);

        for (size_t m = 1; m < (obj.lpr().size() / 2) - 1; m++) {
            auto t5 = sum5(obj.lpr().data(), obj.lpr().size(), m);
            REQUIRE(std::abs(t5) <= epsilon);
        }
    }
}

TEST_CASE("dsp/wavelet: coifCoefTests", "[dsp][wavelet]")
{

    Vector<String> waveletNames;
    waveletNames.resize(17);
    for (size_t i = 0; i < waveletNames.size(); i++) {
        waveletNames[i] = String("coif") + std::to_string(i + 1);
    }

    for (auto const& name : waveletNames) {
        auto obj      = Wavelet{name.c_str()};
        auto const t1 = sum1(obj.lpr()) - mc::sqrt(2.0F);
        auto const t2 = sum2(obj.lpr()) - 1.0F / mc::sqrt(2.0F);
        auto const t3 = sum3(obj.lpr()) - 1.0F / mc::sqrt(2.0F);
        auto const t4 = sum4(obj.lpr()) - 1.0F;

        REQUIRE(std::abs(t1) <= epsilon);
        REQUIRE(std::abs(t2) <= epsilon);
        REQUIRE(std::abs(t3) <= epsilon);
        REQUIRE(std::abs(t4) <= epsilon);
        for (size_t m = 1; m < (obj.lpr().size() / 2) - 1; m++) {
            auto const t5 = sum5(obj.lpr().data(), obj.lpr().size(), m);
            REQUIRE(std::abs(t5) <= epsilon);
        }
    }
}

TEST_CASE("dsp/wavelet: symCoefTests", "[dsp][wavelet]")
{
    Vector<String> waveletNames;
    for (size_t i = 1; i < 20; i++) {
        waveletNames.push_back(String("sym") + std::to_string(i + 1));
    }

    for (auto const& name : waveletNames) {
        auto obj      = Wavelet{name.c_str()};
        auto const t1 = sum1(obj.lpr()) - mc::sqrt(2.0F);
        auto const t2 = sum2(obj.lpr()) - 1.0F / mc::sqrt(2.0F);
        auto const t3 = sum3(obj.lpr()) - 1.0F / mc::sqrt(2.0F);
        auto const t4 = sum4(obj.lpr()) - 1.0F;

        REQUIRE(std::abs(t1) <= epsilon);
        REQUIRE(std::abs(t2) <= epsilon);
        REQUIRE(std::abs(t3) <= epsilon);
        REQUIRE(std::abs(t4) <= epsilon);

        for (size_t m = 1; m < (obj.lpr().size() / 2) - 1; m++) {
            auto const t5 = sum5(obj.lpr().data(), obj.lpr().size(), m);
            REQUIRE(std::abs(t5) <= epsilon);
        }
    }
}
