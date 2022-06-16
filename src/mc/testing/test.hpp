#pragma once

#include "mc/algorithm.hpp"
#include "mc/cmath.hpp"
#include "mc/cstdlib.hpp"
#include "mc/format.hpp"
#include "mc/fstream.hpp"
#include "mc/iterator.hpp"
#include "mc/limits.hpp"
#include "mc/sstream.hpp"
#include "mc/string.hpp"
#include "mc/type_traits.hpp"
#include "mc/vector.hpp"

template<typename It1, typename It2>
auto approxEqual(It1 f1, It1 l1, It2 f2, It2 l2, int epsilonFactor = 4) -> bool
{
    // This function only makes sense when comparing two ranges
    // of the same floating-point type.
    using v1_t = typename std::iterator_traits<It1>::value_type;
    using v2_t = typename std::iterator_traits<It2>::value_type;
    static_assert(std::is_same<v1_t, v2_t>::value);

    auto epsilonEqual = [epsilonFactor](auto l, auto r)
    {
        auto const epsilon = std::numeric_limits<v1_t>::epsilon();
        return std::fabs(l - r) < epsilon * epsilonFactor;
    };

    return std::equal(f1, l1, f2, l2, epsilonEqual);
}

template<typename Container1, typename Container2>
auto approxEqual(Container1 c1, Container2 c2, int epsilonFactor = 4) -> bool
{
    return approxEqual(std::begin(c1), std::end(c1), std::begin(c2), std::end(c2), epsilonFactor);
}

auto absmax(float* array, std::size_t n) -> float;
auto sum1(float const* array, std::size_t n) -> float;
auto sum2(float const* array, std::size_t n) -> float;
auto sum3(float const* array, std::size_t n) -> float;

// np.sum(w[2*m:(2*N+2*m)]*w[0:2*N])
auto sum4(float const* array, std::size_t n) -> float;
// np.sum(w[2 * m:(2 * N)] * w[0:2 * N - 2 * m])
auto sum5(float const* array, std::size_t n, std::size_t m) -> float;

auto rmsError(float const* data, float const* rec, std::size_t n) -> float;
auto relError(float const* data, float const* rec, std::size_t n) -> float;

auto generateRnd() -> float;

template<typename T>
using TestData = std::vector<std::vector<T>>;

auto split(std::string const& s, char delim) -> std::vector<std::string>;
auto loadTestData(char const* filePath) -> TestData<float>;
auto toFloat(TestData<float> const& d) -> TestData<float>;

template<typename T>
[[nodiscard]] auto readFileToVector(char const* filePath) -> std::vector<T>
{
    auto in = std::ifstream{filePath};
    if (!in) { throw std::invalid_argument{fmt::format("Cannot Open File: %s\n", filePath)}; }

    auto result = std::vector<T>{};
    result.reserve(8096 * 4);

    auto line = std::string{};
    while (std::getline(in, line))
    {
        if (!line.empty())
        {
            auto i = T{};
            std::istringstream{line} >> i;
            result.push_back(i);
        }
    }

    in.close();
    return result;
}

auto generateRandomTestData(std::size_t n) -> std::vector<float>;

auto corrcoef(int n, float const* x, float const* y) -> float;