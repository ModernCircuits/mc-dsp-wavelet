#include "wavelib.h"

#include "readFileToVector.hpp"

#include <cmath>
#include <memory>

auto main() -> int
{
    int i;
    int N;
    int J;
    int subscale;
    int a0;
    int iter;
    int nd;
    int k;
    double dt;
    double dj;
    double s0;
    double param;
    double mn;
    double td;
    double tn;
    double den;
    double num;
    double recon_mean;
    double recon_var;
    cwt_object wt;

    // Set Morlet wavelet. Other options "paul" and "dog"
    auto const* wave = "morlet";
    auto const* type = "pow";

    N = 504;
    param = 6.0;
    subscale = 4;
    dt = 0.25;
    s0 = dt;
    dj = 1.0 / (double)subscale;
    J = 11 * subscale; // Total Number of scales
    a0 = 2; //power

    wt = cwt_init(wave, param, N, dt, J);

    auto const inp = readFileToVector("testData/sst_nino3.dat");
    auto oup = std::make_unique<double[]>(N);

    setCWTScales(wt, s0, dj, type, a0);

    cwt(wt, inp.data());

    std::printf("\n MEAN %g \n", wt->smean);

    mn = 0.0;

    for (i = 0; i < N; ++i) {
        mn += sqrt(wt->output[i].re * wt->output[i].re + wt->output[i].im * wt->output[i].im);
    }

    cwt_summary(wt);

    std::printf("\n abs mean %g \n", mn / N);

    std::printf("\n\n");
    std::printf("Let CWT w = w(j, n/2 - 1) where n = %d\n\n", N);
    nd = N / 2 - 1;

    std::printf("%-15s%-15s%-15s%-15s \n", "j", "Scale", "Period", "ABS(w)^2");
    for (k = 0; k < wt->J; ++k) {
        iter = nd + k * N;
        std::printf("%-15d%-15lf%-15lf%-15lf \n", k, wt->scale[k], wt->period[k],
            wt->output[iter].re * wt->output[iter].re + wt->output[iter].im * wt->output[iter].im);
    }

    icwt(wt, oup.get());

    num = den = recon_var = recon_mean = 0.0;
    std::printf("\n\n");
    std::printf("Signal Reconstruction\n");
    std::printf("%-15s%-15s%-15s \n", "i", "Input(i)", "Output(i)");

    for (i = N - 10; i < N; ++i) {
        std::printf("%-15d%-15lf%-15lf \n", i, inp[i], oup[i]);
    }

    for (i = 0; i < N; ++i) {
        //std::printf("%g %g \n", oup[i] ,inp[i] - wt->smean);
        td = inp[i];
        tn = oup[i] - td;
        num += (tn * tn);
        den += (td * td);
        recon_mean += oup[i];
    }

    recon_var = sqrt(num / N);
    recon_mean /= N;

    std::printf("\nRMS Error %g \n", sqrt(num) / sqrt(den));
    std::printf("\nVariance %g \n", recon_var);
    std::printf("\nMean %g \n", recon_mean);

    cwt_free(wt);
    return 0;
}
