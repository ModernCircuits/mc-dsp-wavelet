#pragma once

#include <cmath>
#include <cstdlib>
#include <sstream>
#include <string>

namespace patch {
template <typename T>
auto to_string(const T& n) -> std::string
{
    std::ostringstream stm;
    stm << n;
    return stm.str();
}
} // namespace patch

inline auto absmax(double* array, int N) -> double
{
    double max;
    int i;

    max = 0.0;
    for (i = 0; i < N; ++i) {
        if (fabs(array[i]) >= max) {
            max = fabs(array[i]);
        }
    }

    return max;
}

inline auto sum1(const double* array, int N) -> double
{
    double sum;
    int i;

    sum = 0.0;
    for (i = 0; i < N; ++i) {
        sum += array[i];
    }
    return sum;
}
inline auto sum2(const double* array, int N) -> double
{
    double sum;
    int i;

    sum = 0.0;
    for (i = 0; i < N; i += 2) {
        sum += array[i];
    }
    return sum;
}
inline auto sum3(const double* array, int N) -> double
{
    double sum;
    int i;

    sum = 0.0;
    for (i = 1; i < N; i += 2) {
        sum += array[i];
    }
    return sum;
}
// np.sum(w[2*m:(2*N+2*m)]*w[0:2*N])
inline auto sum4(const double* array, int N) -> double
{
    double sum;
    int i;

    sum = 0.0;
    for (i = 0; i < N; i += 1) {
        sum += array[i] * array[i];
    }
    return sum;
}
// np.sum(w[2 * m:(2 * N)] * w[0:2 * N - 2 * m])
inline auto sum5(const double* array, int N, int m) -> double
{
    double sum;
    int i;

    sum = 0.0;
    for (i = 2 * m; i < N; i += 1) {
        sum += array[i] * array[i - 2 * m];
    }
    return sum;
}

inline auto RMS_Error(const double* data, const double* rec, int N) -> double
{
    int i;
    double sum = 0;
    for (i = 0; i < N; ++i) {
        sum += (data[i] - rec[i]) * (data[i] - rec[i]);
    }
    return sqrt(sum / ((double)N - 1));
}

inline auto REL_Error(const double* data, const double* rec, int N) -> double
{
    int i;
    double sum1 = 0;
    double sum2 = 0;
    for (i = 0; i < N; ++i) {
        sum1 += (data[i] - rec[i]) * (data[i] - rec[i]);
        sum2 += data[i] * data[i];
    }
    return sqrt(sum1) / sqrt(sum2);
}

inline auto generate_rnd() -> double
{
    return (double)(std::rand() % 100 + 1);
}