#include "mc/dsp/wavelets/Denoise.hpp"

#include "mc/cmath.hpp"
#include "mc/cstdlib.hpp"
#include "mc/cstring.hpp"
#include "mc/format.hpp"
#include "mc/memory.hpp"
#include "mc/testing/test.hpp"

namespace dsp = mc::dsp;

auto main() -> int
{
    char const* wname  = "db5";
    char const* method = "dwt";
    char const* ext    = "sym";
    char const* thresh = "soft";
    char const* level  = "all";

    auto* ifp = std::fopen("testData/pieceregular1024.txt", "r");
    auto i    = 0;
    if (ifp == nullptr)
    {
        fmt::printf("Cannot Open File");
        std::exit(EXIT_FAILURE);
    }

    float temp[2400];
    while (std::feof(ifp) == 0)
    {
        std::fscanf(ifp, "%f \n", &temp[i]);
        i++;
    }
    std::fclose(ifp);

    auto n = i;
    auto j = 4;

    auto inp = std::make_unique<float[]>(n);
    auto oup = std::make_unique<float[]>(n);
    auto sig = std::make_unique<float[]>(n);

    for (i = 0; i < n; ++i) { sig[i] = temp[i]; }

    ifp = std::fopen("testData/PieceRegular10.txt", "r");
    i   = 0;
    if (ifp == nullptr)
    {
        fmt::printf("Cannot Open File");
        std::exit(EXIT_FAILURE);
    }

    while (std::feof(ifp) == 0)
    {
        std::fscanf(ifp, "%f \n", &temp[i]);
        i++;
    }

    std::fclose(ifp);

    for (i = 0; i < n; ++i) { inp[i] = temp[i]; }
    auto obj = dsp::DenoiseSet(n, j, wname);
    setDenoiseMethod(obj,
                     "visushrink");  // sureshrink is also the default. The other option with dwt and swt is visushrink.
    // modwt works only with modwtshrink method
    setDenoiseWTMethod(obj, method);           // Default is dwt. the other options are swt and modwt
    setDenoiseWTExtension(obj, ext);           // Default is sym. the other option is per
    setDenoiseParameters(obj, thresh, level);  // Default for thresh is soft. Other option is hard
    // Default for level is all. The other option is first

    denoise(obj, inp.get(), oup.get());

    // Alternative to denoise_set*
    // Just use visushrink, modwtshrink and sureshrink functions
    // visushrink(inp.get(),N,J,wname,method,ext,thresh,level,oup.get());
    // sureshrink(inp.get(),N,J,wname,method,ext,thresh,level,oup.get());
    // modwtshrink(sig.get(),N,J,wname,cmethod,ext,thresh,oup.get()); See modwtdenoisetest.c
    // ofp = std::fopen("testData/denoiseds.txt", "w");

    fmt::printf("Signal - Noisy Signal Stats \n");
    fmt::printf("RMSE %g\n", rmsError(sig.get(), inp.get(), n));
    fmt::printf("Corr Coeff %g\n", corrcoef(n, sig.get(), inp.get()));

    fmt::printf("Signal - DeNoised Signal Stats \n");
    fmt::printf("RMSE %g\n", rmsError(sig.get(), oup.get(), n));
    fmt::printf("Corr Coeff %g\n", corrcoef(n, sig.get(), oup.get()));

    return 0;
}
