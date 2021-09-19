#pragma once

#include <sstream>
#include <string>

auto absmax(double* array, int N) -> double;
auto sum1(double const* array, int N) -> double;
auto sum2(double const* array, int N) -> double;
auto sum3(double const* array, int N) -> double;

// np.sum(w[2*m:(2*N+2*m)]*w[0:2*N])
auto sum4(double const* array, int N) -> double;
// np.sum(w[2 * m:(2 * N)] * w[0:2 * N - 2 * m])
auto sum5(double const* array, int N, int m) -> double;

auto RMS_Error(double const* data, double const* rec, int N) -> double;
auto REL_Error(double const* data, double const* rec, int N) -> double;

auto generate_rnd() -> double;