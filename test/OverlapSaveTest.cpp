#include "lt/dsp/convolution.hpp"

#include <iostream>
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

auto main() -> int
{
    // auto si = std::array {1.0, 2.0, 3.0};
    // auto pi = std::array {0.0, 1.0, 0.5};

    auto si = std::array { 1.0, 1.0, 1.0, 1.0, 1.0 };
    auto pi = std::array { 0.0, 1.0, 0.5, 1.0, 1.0 };

    auto s = DoubleSignal { si.data(), si.size() };
    auto p = DoubleSignal { pi.data(), pi.size() };
    auto x = OverlapSave { s, p };

    x.convolute();
    REQUIRE(approxEqual(x.extractResult(), std::array { 0.0, 1.0, 1.5, 2.5, 3.5, 3.5, 2.5, 2.0, 1.0 }));

    x.crossCorrelate();
    REQUIRE(approxEqual(x.extractResult(), std::array { 1.0, 2.0, 2.5, 3.5, 3.5, 2.5, 1.5, 1.0, 0.0 }));

    return EXIT_SUCCESS;
}