#pragma once

#include "lt/cmath.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

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

auto absmax(double* array, int N) -> double;
auto sum1(double const* array, int N) -> double;
auto sum2(double const* array, int N) -> double;
auto sum3(double const* array, int N) -> double;

// np.sum(w[2*m:(2*N+2*m)]*w[0:2*N])
auto sum4(double const* array, int N) -> double;
// np.sum(w[2 * m:(2 * N)] * w[0:2 * N - 2 * m])
auto sum5(double const* array, int N, int m) -> double;

auto rmsError(double const* data, double const* rec, int N) -> double;
auto relError(double const* data, double const* rec, int N) -> double;

auto generateRnd() -> double;

template <typename T>
using TestData = std::vector<std::vector<T>>;

auto split(std::string const& s, char delim) -> std::vector<std::string>;
auto loadTestData(char const* filePath) -> TestData<double>;
auto toFloat(TestData<double> const& d) -> TestData<float>;