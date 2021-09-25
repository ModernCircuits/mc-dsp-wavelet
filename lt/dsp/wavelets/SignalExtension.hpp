#pragma once

#include "lt/string.hpp"

enum struct SignalExtension {
    periodic,
    symmetric,
};

[[nodiscard]] inline auto toString(SignalExtension ext) -> std::string
{
    if (ext == SignalExtension::periodic) {
        return "periodic";
    }
    return "symmetric";
}
