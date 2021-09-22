#pragma once

auto waveletFilterLength(char const* name) -> int;
auto waveletFilterCoefficients(char const* name, double* lp1, double* hp1, double* lp2, double* hp2) -> int;
