#include "WindowFunction.hpp"

#include "lt/cassert.hpp"

namespace lt
{
namespace dsp
{
auto toString(WindowFunction wf) -> std::string
{
    switch (wf)
    {
        case WindowFunction::rectangular: return "Rectangular";
        case WindowFunction::triangular: return "Triangular";
        case WindowFunction::hann: return "Hann";
        case WindowFunction::hamming: return "Hamming";
        case WindowFunction::blackman: return "Blackman";
        case WindowFunction::blackmanHarris: return "Blackman-Harris";
    }

    LT_ASSERT(false);  // NOLINT
    return "";
}
}  // namespace dsp
}  // namespace lt
