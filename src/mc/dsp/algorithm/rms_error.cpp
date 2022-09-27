#include "rms_error.hpp"

#include <mc/core/cmath.hpp>
#include <mc/core/numeric.hpp>
#include <mc/core/span.hpp>

namespace mc {

auto rmsError(float const* data, float const* rec, std::size_t n) -> float
{
    auto diffSquared = [](auto x, auto y) {
        auto const diff = x - y;
        return diff * diff;
    };

    auto sum = std::transform_reduce(data, data + n, rec, 0.0F, std::plus<>{}, diffSquared);
    return sqrt(sum / static_cast<float>(n - 1U));
}

}  // namespace mc
