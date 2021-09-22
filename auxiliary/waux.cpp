#include "waux.h"
#include "wauxlib.h"

#include <algorithm>
#include <memory>
#include <numeric>

auto mean(double const* vec, int n) -> double
{
    return std::accumulate(vec, vec + n, 0.0) / static_cast<double>(n);
}

auto var(double const* vec, int n) -> double
{
    double v;
    double temp;
    double m;
    v = 0.0;
    m = mean(vec, n);

    for (auto i = 0; i < n; ++i) {
        temp = vec[i] - m;
        v += temp * temp;
    }

    v = v / n;

    return v;
}

auto median(double* const x, int n) -> double
{
    std::sort(x, x + n, std::less<double> {});

    double sigma;
    if ((n % 2) == 0) {
        sigma = (x[n / 2 - 1] + x[n / 2]) / 2.0;
    } else {
        sigma = x[n / 2];
    }

    return sigma;
}

auto mad(double* x, int n) -> double
{
    double sigma;

    sigma = median(x, n);

    for (auto i = 0; i < n; ++i) {
        x[i] = (x[i] - sigma) > 0 ? (x[i] - sigma) : -(x[i] - sigma);
    }

    sigma = median(x, n);

    return sigma;
}

auto minindex(double const* arr, int n) -> int
{
    double min;
    int index;

    min = DBL_MAX;
    index = 0;
    for (auto i = 0; i < n; ++i) {
        if (arr[i] < min) {
            min = arr[i];
            index = i;
        }
    }

    return index;
}

void autocovar(double const* vec, int n, double* acov, int m)
{
    double temp1;
    double temp2;
    int t;

    if (m > n) {
        m = n - 1;
        printf("\n Lag is greater than the length N of the input vector. It is automatically set to length N - 1.\n");
        printf("\n The Output Vector only contains N calculated values.");
    } else if (m < 0) {
        m = 0;
    }

    auto const me = mean(vec, n);
    for (auto i = 0; i < m; i++) {
        acov[i] = 0.0;
        for (t = 0; t < n - i; t++) {
            temp1 = vec[t] - me;
            temp2 = vec[t + i] - me;
            acov[i] += temp1 * temp2;
        }
        acov[i] = acov[i] / n;
    }
}

void autocorr(double const* vec, int n, double* acorr, int m)
{
    double var;
    if (m > n) {
        m = n - 1;
        printf("\n Lag is greater than the length N of the input vector. It is automatically set to length N - 1.\n");
        printf("\n The Output Vector only contains N calculated values.");
    } else if (m < 0) {
        m = 0;
    }
    autocovar(vec, n, acorr, m);
    var = acorr[0];
    acorr[0] = 1.0;

    for (auto i = 1; i < m; i++) {
        acorr[i] = acorr[i] / var;
    }
}
