#include "helper.hpp"

#include <cmath>
#include <cstdlib>

auto absmax(double* array, int N) -> double
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

auto sum1(double const* array, int N) -> double
{
    double sum;
    int i;

    sum = 0.0;
    for (i = 0; i < N; ++i) {
        sum += array[i];
    }
    return sum;
}
auto sum2(double const* array, int N) -> double
{
    double sum;
    int i;

    sum = 0.0;
    for (i = 0; i < N; i += 2) {
        sum += array[i];
    }
    return sum;
}
auto sum3(double const* array, int N) -> double
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
auto sum4(double const* array, int N) -> double
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
auto sum5(double const* array, int N, int m) -> double
{
    double sum;
    int i;

    sum = 0.0;
    for (i = 2 * m; i < N; i += 1) {
        sum += array[i] * array[i - 2 * m];
    }
    return sum;
}

auto RMS_Error(double const* data, double const* rec, int N) -> double
{
    int i;
    double sum = 0;
    for (i = 0; i < N; ++i) {
        sum += (data[i] - rec[i]) * (data[i] - rec[i]);
    }
    return std::sqrt(sum / ((double)N - 1));
}

auto REL_Error(double const* data, double const* rec, int N) -> double
{
    int i;
    double sum1 = 0;
    double sum2 = 0;
    for (i = 0; i < N; ++i) {
        sum1 += (data[i] - rec[i]) * (data[i] - rec[i]);
        sum2 += data[i] * data[i];
    }
    return std::sqrt(sum1) / std::sqrt(sum2);
}

auto generate_rnd() -> double
{
    return (double)(std::rand() % 100 + 1);
}