#include <mc/dsp/wavelets/Denoise.hpp>

#include <mc/core/cmath.hpp>
#include <mc/core/cstdlib.hpp>
#include <mc/core/cstring.hpp>
#include <mc/core/format.hpp>
#include <mc/core/memory.hpp>
#include <mc/testing/test.hpp>

using namespace mc;

auto main() -> int
{
    char const* wname  = "db5";
    char const* method = "dwt";
    char const* ext    = "sym";
    char const* thresh = "soft";
    char const* level  = "all";

    auto* ifp = std::fopen("test_data/raw/pieceregular1024.txt", "r");
    auto i    = 0;
    if (ifp == nullptr) {
        print("Cannot Open File");
        std::exit(EXIT_FAILURE);
    }

    float temp[2400];
    while (std::feof(ifp) == 0) {
        std::fscanf(ifp, "%f \n", &temp[i]);
        i++;
    }
    std::fclose(ifp);

    auto n = i;
    auto j = 4;

    auto inp = makeUnique<float[]>(n);
    auto oup = makeUnique<float[]>(n);
    auto sig = makeUnique<float[]>(n);

    for (i = 0; i < n; ++i) { sig[i] = temp[i]; }

    ifp = std::fopen("test_data/raw/PieceRegular10.txt", "r");
    i   = 0;
    if (ifp == nullptr) {
        print("Cannot Open File");
        std::exit(EXIT_FAILURE);
    }

    while (std::feof(ifp) == 0) {
        std::fscanf(ifp, "%f \n", &temp[i]);
        i++;
    }

    std::fclose(ifp);

    for (i = 0; i < n; ++i) { inp[i] = temp[i]; }
    auto obj = dsp::DenoiseSet(n, j, wname);
    setDenoiseMethod(
        obj,
        "visushrink"
    );  // sureshrink is also the default. The other option with dwt and swt is visushrink.
    // modwt works only with modwtshrink method
    setDenoiseWTMethod(obj, method);  // Default is dwt. the other options are swt and modwt
    setDenoiseWTExtension(obj, ext);  // Default is sym. the other option is per
    setDenoiseParameters(obj, thresh, level);
    // Default for thresh is soft. Other option is hard
    // Default for level is all. The other option is first

    denoise(obj, inp.get(), oup.get());

    // Alternative to denoise_set*
    // Just use visushrink, modwtshrink and sureshrink functions
    // visushrink(inp.get(),N,J,wname,method,ext,thresh,level,oup.get());
    // sureshrink(inp.get(),N,J,wname,method,ext,thresh,level,oup.get());
    // modwtshrink(sig.get(),N,J,wname,cmethod,ext,thresh,oup.get()); See modwtdenoisetest.c
    // ofp = std::fopen("test_data/raw/denoiseds.txt", "w");

    print("Signal - Noisy Signal Stats \n");
    print("RMSE {}\n", rmsError(sig.get(), inp.get(), n));
    print("Corr Coeff {}\n", corrcoef(n, sig.get(), inp.get()));

    print("Signal - DeNoised Signal Stats \n");
    print("RMSE {}\n", rmsError(sig.get(), oup.get(), n));
    print("Corr Coeff {}\n", corrcoef(n, sig.get(), oup.get()));

    return 0;
}
