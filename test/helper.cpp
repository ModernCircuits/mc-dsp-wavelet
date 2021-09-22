#include "helper.hpp"

#include <cmath>
#include <cstdlib>
#include <random>

auto absmax(double* array, int n) -> double
{
    double max;
    int i;

    max = 0.0;
    for (i = 0; i < n; ++i) {
        if (fabs(array[i]) >= max) {
            max = fabs(array[i]);
        }
    }

    return max;
}

auto sum1(double const* array, int n) -> double
{
    double sum;
    int i;

    sum = 0.0;
    for (i = 0; i < n; ++i) {
        sum += array[i];
    }
    return sum;
}
auto sum2(double const* array, int n) -> double
{
    double sum;
    int i;

    sum = 0.0;
    for (i = 0; i < n; i += 2) {
        sum += array[i];
    }
    return sum;
}
auto sum3(double const* array, int n) -> double
{
    double sum;
    int i;

    sum = 0.0;
    for (i = 1; i < n; i += 2) {
        sum += array[i];
    }
    return sum;
}
// np.sum(w[2*m:(2*N+2*m)]*w[0:2*N])
auto sum4(double const* array, int n) -> double
{
    double sum;
    int i;

    sum = 0.0;
    for (i = 0; i < n; i += 1) {
        sum += array[i] * array[i];
    }
    return sum;
}
// np.sum(w[2 * m:(2 * N)] * w[0:2 * N - 2 * m])
auto sum5(double const* array, int n, int m) -> double
{
    double sum;
    int i;

    sum = 0.0;
    for (i = 2 * m; i < n; i += 1) {
        sum += array[i] * array[i - 2 * m];
    }
    return sum;
}

auto rmsError(double const* data, double const* rec, int n) -> double
{
    int i;
    double sum = 0;
    for (i = 0; i < n; ++i) {
        sum += (data[i] - rec[i]) * (data[i] - rec[i]);
    }
    return std::sqrt(sum / ((double)n - 1));
}

auto relError(double const* data, double const* rec, int n) -> double
{
    int i;
    double sum1 = 0;
    double sum2 = 0;
    for (i = 0; i < n; ++i) {
        sum1 += (data[i] - rec[i]) * (data[i] - rec[i]);
        sum2 += data[i] * data[i];
    }
    return std::sqrt(sum1) / std::sqrt(sum2);
}

auto generateRnd() -> double
{
    auto rd = std::random_device {};
    auto gen = std::mt19937 { rd() };
    auto dis = std::uniform_real_distribution<double> { 1.0, 100.0 };
    return dis(gen);
}