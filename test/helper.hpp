#pragma once

#include <sstream>
#include <string>

auto absmax(double* array, int N) -> double;
auto sum1(const double* array, int N) -> double;
auto sum2(const double* array, int N) -> double;
auto sum3(const double* array, int N) -> double;

// np.sum(w[2*m:(2*N+2*m)]*w[0:2*N])
auto sum4(const double* array, int N) -> double;
// np.sum(w[2 * m:(2 * N)] * w[0:2 * N - 2 * m])
auto sum5(const double* array, int N, int m) -> double;

auto RMS_Error(const double* data, const double* rec, int N) -> double;
auto REL_Error(const double* data, const double* rec, int N) -> double;

auto generate_rnd() -> double;