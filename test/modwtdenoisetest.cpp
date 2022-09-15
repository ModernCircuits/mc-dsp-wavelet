#include "mc/dsp/wavelets/Denoise.hpp"

#include "mc/cmath.hpp"
#include "mc/cstdlib.hpp"
#include "mc/cstring.hpp"
#include "mc/format.hpp"
#include "mc/memory.hpp"

#include "mc/testing/test.hpp"

namespace dsp = mc::dsp;

// modwtshrink can also be called from the denoise object.
// See denoisetest.cpp for more information.
auto main() -> int
{
    auto const* wname   = "db5";
    auto const* ext     = "per";
    auto const* thresh  = "soft";
    auto const* cmethod = "direct";

    auto const inp = readFileToVector<float>("test_data/raw/PieceRegular10.txt");
    auto sig       = readFileToVector<float>("test_data/raw/pieceregular1024.txt");

    auto const n = sig.size();
    auto const j = 4;
    auto out     = std::make_unique<float[]>(n);

    dsp::modwtshrink(sig.data(), n, j, wname, cmethod, ext, thresh, out.get());

    fmt::printf("Signal - Noisy Signal Stats \n");
    fmt::printf("RMSE %g\n", rmsError(sig.data(), inp.data(), n));
    fmt::printf("Corr Coeff %g\n", corrcoef(n, sig.data(), inp.data()));

    fmt::printf("Signal - DeNoised Signal Stats \n");
    fmt::printf("RMSE %g\n", rmsError(sig.data(), out.get(), n));
    fmt::printf("Corr Coeff %g\n", corrcoef(n, sig.data(), out.get()));

    return 0;
}
