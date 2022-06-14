#pragma once

#include <memory>

template<typename T>
auto makeZeros(std::size_t length) -> std::unique_ptr<T[]>
{
    auto ptr = std::make_unique<T[]>(length);
    for (std::size_t i{0}; i < length; ++i) { ptr[i] = T{}; }
    return ptr;
}
