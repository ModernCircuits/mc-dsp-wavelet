#include "lt/bit.hpp"

#include "lt/testing/test.hpp"

template <typename T>
auto test() -> bool
{
    REQUIRE(lt::bit_ceil(T(0)) == T(1));
    REQUIRE(lt::bit_ceil(T(1)) == T(1));
    REQUIRE(lt::bit_ceil(T(2)) == T(2));
    REQUIRE(lt::bit_ceil(T(3)) == T(4));
    REQUIRE(lt::bit_ceil(T(4)) == T(4));
    REQUIRE(lt::bit_ceil(T(15)) == T(16));
    REQUIRE(lt::bit_ceil(T(21)) == T(32));
    REQUIRE(lt::bit_ceil(T(33)) == T(64));
    REQUIRE(lt::bit_ceil(T(55)) == T(64));
    return true;
}

auto main() -> int
{
    REQUIRE(test<unsigned char>());
    REQUIRE(test<unsigned short>());
    REQUIRE(test<unsigned int>());
    // REQUIRE(test<unsigned long>());
    // REQUIRE(test<unsigned long long>());
    return 0;
}