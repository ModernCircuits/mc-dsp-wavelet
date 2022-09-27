// SPDX-License-Identifier: BSL-1.0

#include <mc/fft/algorithm.hpp>

#include <mc/wavelet/algorithm.hpp>

#include <mc/core/iterator.hpp>
#include <mc/core/memory.hpp>
#include <mc/testing/test.hpp>
#include <mc/wavelet.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace mc;

static constexpr auto epsilon = 1e-5F;

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define DWT_IDWT_ROUNDTRIP(waveletName)                                                    \
    TEST_CASE(                                                                             \
        "wavelet: WaveletPacketTransform(dwpt/idwpt) - " waveletName,                      \
        "[dsp][wavelet]"                                                                   \
    )                                                                                      \
    {                                                                                      \
        auto extension = GENERATE("sym", "per");                                           \
        auto entropy   = GENERATE("shannon", "logenergy");                                 \
        auto levels    = GENERATE(as<size_t>{}, 1, 2);                                     \
        auto n         = 8096;                                                             \
        auto inp       = generateRandomTestData(n);                                        \
        auto out       = makeUnique<float[]>(n);                                           \
        auto obj       = Wavelet{waveletName};                                             \
        auto wt        = WaveletPacketTransform(&obj, n, levels);                          \
        setDWPTExtension(wt, extension);                                                   \
        setDWPTEntropy(wt, entropy, 0);                                                    \
        dwpt(wt, data(inp));                                                               \
        idwpt(wt, out.get());                                                              \
        REQUIRE_THAT(                                                                      \
            rmsError(out.get(), data(inp), wt.signalLength()),                             \
            Catch::Matchers::WithinAbs(0.0F, epsilon)                                      \
        );                                                                                 \
    }

DWT_IDWT_ROUNDTRIP("db1")     // NOLINT
DWT_IDWT_ROUNDTRIP("db2")     // NOLINT
DWT_IDWT_ROUNDTRIP("db3")     // NOLINT
DWT_IDWT_ROUNDTRIP("db4")     // NOLINT
DWT_IDWT_ROUNDTRIP("db5")     // NOLINT
DWT_IDWT_ROUNDTRIP("db6")     // NOLINT
DWT_IDWT_ROUNDTRIP("db7")     // NOLINT
DWT_IDWT_ROUNDTRIP("db8")     // NOLINT
DWT_IDWT_ROUNDTRIP("db9")     // NOLINT
DWT_IDWT_ROUNDTRIP("db10")    // NOLINT
DWT_IDWT_ROUNDTRIP("db11")    // NOLINT
DWT_IDWT_ROUNDTRIP("db12")    // NOLINT
DWT_IDWT_ROUNDTRIP("db13")    // NOLINT
DWT_IDWT_ROUNDTRIP("db14")    // NOLINT
DWT_IDWT_ROUNDTRIP("db15")    // NOLINT
DWT_IDWT_ROUNDTRIP("db16")    // NOLINT
DWT_IDWT_ROUNDTRIP("db17")    // NOLINT
DWT_IDWT_ROUNDTRIP("db18")    // NOLINT
DWT_IDWT_ROUNDTRIP("db19")    // NOLINT
DWT_IDWT_ROUNDTRIP("db20")    // NOLINT
DWT_IDWT_ROUNDTRIP("db21")    // NOLINT
DWT_IDWT_ROUNDTRIP("db22")    // NOLINT
DWT_IDWT_ROUNDTRIP("db23")    // NOLINT
DWT_IDWT_ROUNDTRIP("db24")    // NOLINT
DWT_IDWT_ROUNDTRIP("db25")    // NOLINT
DWT_IDWT_ROUNDTRIP("db26")    // NOLINT
DWT_IDWT_ROUNDTRIP("db27")    // NOLINT
DWT_IDWT_ROUNDTRIP("db28")    // NOLINT
DWT_IDWT_ROUNDTRIP("db29")    // NOLINT
DWT_IDWT_ROUNDTRIP("db30")    // NOLINT
DWT_IDWT_ROUNDTRIP("db31")    // NOLINT
DWT_IDWT_ROUNDTRIP("db32")    // NOLINT
DWT_IDWT_ROUNDTRIP("db33")    // NOLINT
DWT_IDWT_ROUNDTRIP("db34")    // NOLINT
DWT_IDWT_ROUNDTRIP("db35")    // NOLINT
DWT_IDWT_ROUNDTRIP("coif1")   // NOLINT
DWT_IDWT_ROUNDTRIP("coif2")   // NOLINT
DWT_IDWT_ROUNDTRIP("coif3")   // NOLINT
DWT_IDWT_ROUNDTRIP("coif4")   // NOLINT
DWT_IDWT_ROUNDTRIP("coif5")   // NOLINT
DWT_IDWT_ROUNDTRIP("coif6")   // NOLINT
DWT_IDWT_ROUNDTRIP("coif7")   // NOLINT
DWT_IDWT_ROUNDTRIP("coif8")   // NOLINT
DWT_IDWT_ROUNDTRIP("coif9")   // NOLINT
DWT_IDWT_ROUNDTRIP("coif10")  // NOLINT
DWT_IDWT_ROUNDTRIP("coif11")  // NOLINT
DWT_IDWT_ROUNDTRIP("coif12")  // NOLINT
DWT_IDWT_ROUNDTRIP("coif13")  // NOLINT
DWT_IDWT_ROUNDTRIP("coif14")  // NOLINT
DWT_IDWT_ROUNDTRIP("coif15")  // NOLINT
DWT_IDWT_ROUNDTRIP("coif16")  // NOLINT
DWT_IDWT_ROUNDTRIP("coif17")  // NOLINT
DWT_IDWT_ROUNDTRIP("sym2")    // NOLINT
DWT_IDWT_ROUNDTRIP("sym3")    // NOLINT
DWT_IDWT_ROUNDTRIP("sym4")    // NOLINT
DWT_IDWT_ROUNDTRIP("sym5")    // NOLINT
DWT_IDWT_ROUNDTRIP("sym6")    // NOLINT
DWT_IDWT_ROUNDTRIP("sym7")    // NOLINT
DWT_IDWT_ROUNDTRIP("sym8")    // NOLINT
DWT_IDWT_ROUNDTRIP("sym9")    // NOLINT
DWT_IDWT_ROUNDTRIP("sym10")   // NOLINT
DWT_IDWT_ROUNDTRIP("sym11")   // NOLINT
DWT_IDWT_ROUNDTRIP("sym12")   // NOLINT
DWT_IDWT_ROUNDTRIP("sym13")   // NOLINT
DWT_IDWT_ROUNDTRIP("sym14")   // NOLINT
DWT_IDWT_ROUNDTRIP("sym15")   // NOLINT
DWT_IDWT_ROUNDTRIP("sym16")   // NOLINT
DWT_IDWT_ROUNDTRIP("sym17")   // NOLINT
DWT_IDWT_ROUNDTRIP("sym18")   // NOLINT
DWT_IDWT_ROUNDTRIP("sym19")   // NOLINT
DWT_IDWT_ROUNDTRIP("sym20")   // NOLINT
