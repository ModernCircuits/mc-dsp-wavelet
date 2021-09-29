#pragma once

#include "lt/preprocessor.hpp"
#include "lt/string.hpp"

enum struct SignalExtension {
    periodic,
    symmetric,
};

LT_NODISCARD inline auto toString(SignalExtension ext) -> std::string
{
    if (ext == SignalExtension::periodic) {
        return "periodic";
    }
    return "symmetric";
}
