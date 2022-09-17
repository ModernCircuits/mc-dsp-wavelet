#pragma once

#include <mc/core/algorithm.hpp>
#include <mc/core/cmath.hpp>
#include <mc/core/cstdlib.hpp>
#include <mc/core/format.hpp>
#include <mc/core/fstream.hpp>
#include <mc/core/iterator.hpp>
#include <mc/core/limits.hpp>
#include <mc/core/raise.hpp>
#include <mc/core/span.hpp>
#include <mc/core/sstream.hpp>
#include <mc/core/stdexcept.hpp>
#include <mc/core/string.hpp>
#include <mc/core/type_traits.hpp>
#include <mc/core/vector.hpp>

namespace mc {
template<typename It1, typename It2>
auto approxEqual(It1 f1, It1 l1, It2 f2, It2 l2, int epsilonFactor = 4) -> bool
{
    // This function only makes sense when comparing two ranges
    // of the same floating-point type.
    using v1_t = typename std::iterator_traits<It1>::value_type;
    using v2_t = typename std::iterator_traits<It2>::value_type;
    static_assert(std::is_same<v1_t, v2_t>::value);

    auto epsilonEqual = [epsilonFactor](auto l, auto r) {
        auto const epsilon = std::numeric_limits<v1_t>::epsilon();
        return std::abs(l - r) < epsilon * epsilonFactor;
    };

    return std::equal(f1, l1, f2, l2, epsilonEqual);
}

template<typename T>
auto approxEqual(Span<T const> lhs, Span<T const> rhs, int epsilonFactor = 4) -> bool
{
    return approxEqual(begin(lhs), end(lhs), begin(rhs), end(rhs), epsilonFactor);
}

auto absmax(float const* array, std::size_t n) -> float;
auto rmsError(float const* data, float const* rec, std::size_t n) -> float;
auto relError(float const* data, float const* rec, std::size_t n) -> float;

auto generateRnd() -> float;

template<typename T>
using TestData = Vector<Vector<T>>;

auto split(String const& s, char delim) -> Vector<String>;
auto loadTestData(char const* filePath) -> TestData<float>;
auto toFloat(TestData<float> const& d) -> TestData<float>;

template<typename T>
[[nodiscard]] auto readFileToVector(char const* filePath) -> Vector<T>
{
    auto in = std::ifstream{filePath};
    if (!in) { raisef<InvalidArgument>("cannot open file: {:s}", filePath); }

    auto result = Vector<T>{};
    result.reserve(8096 * 4);

    auto line = String{};
    while (std::getline(in, line)) {
        if (!line.empty()) {
            auto i = T{};
            std::istringstream{line} >> i;
            result.push_back(i);
        }
    }

    in.close();
    return result;
}

auto generateRandomTestData(std::size_t n) -> Vector<float>;

auto corrcoef(int n, float const* x, float const* y) -> float;

}  // namespace mc
