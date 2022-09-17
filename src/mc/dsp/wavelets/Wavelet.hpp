#pragma once

#include <mc/core/config.hpp>

#include <mc/core/memory.hpp>
#include <mc/core/span.hpp>
#include <mc/core/string.hpp>

namespace mc::dsp {
struct Wavelet
{
    explicit Wavelet(char const* wname);

    [[nodiscard]] auto size() const noexcept -> std::size_t;
    [[nodiscard]] auto name() const noexcept -> String const&;

    [[nodiscard]] auto lpd() const noexcept -> Span<float>;
    [[nodiscard]] auto hpd() const noexcept -> Span<float>;
    [[nodiscard]] auto lpr() const noexcept -> Span<float>;
    [[nodiscard]] auto hpr() const noexcept -> Span<float>;

private:
    String _name;

    // When all filters are of the same length.
    // [Matlab uses zero-padding to make all filters of the same length]
    std::size_t _size;
    UniquePtr<float[]> _params;
};

[[nodiscard]] auto summary(Wavelet const& obj) -> String;

}  // namespace mc::dsp

namespace mc {
template<typename T>
auto makeZeros(std::size_t length) -> UniquePtr<T[]>
{
    auto ptr = makeUnique<T[]>(length);
    for (std::size_t i{0}; i < length; ++i) { ptr[i] = T{}; }
    return ptr;
}
}  // namespace mc
