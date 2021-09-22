#pragma once

#include <string>

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
