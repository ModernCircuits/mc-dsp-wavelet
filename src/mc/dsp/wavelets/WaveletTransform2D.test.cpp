#include "mc/dsp/wavelets.hpp"

#include "mc/memory.hpp"
#include "mc/testing/test.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

namespace dsp = mc::dsp;

static constexpr auto const epsilon = 6e-7F;

#define DWT_IDWT_ROUNDTRIP(waveletName)                                                                                \
    TEST_CASE("dsp/wavelet: WaveletTransform2D(dwt/idwt) - " waveletName, "[dsp][wavelet]")                            \
    {                                                                                                                  \
        auto extension  = GENERATE(dsp::SignalExtension::periodic, dsp::SignalExtension::symmetric);                   \
        auto levels     = GENERATE(1);                                                                                 \
        auto const rows = 256;                                                                                         \
        auto const cols = 200;                                                                                         \
        auto const n    = rows * cols;                                                                                 \
        auto const inp  = generateRandomTestData(n);                                                                   \
        auto out        = makeZeros<float>(n);                                                                         \
        auto wavelet    = dsp::Wavelet{waveletName};                                                                   \
        auto wt         = dsp::WaveletTransform2D(wavelet, "dwt", rows, cols, levels);                                 \
        setDWT2Extension(wt, extension == dsp::SignalExtension::symmetric ? "sym" : "per");                            \
        auto wavecoeffs = dwt(wt, data(inp));                                                                          \
        idwt(wt, wavecoeffs.get(), out.get());                                                                         \
        REQUIRE(rmsError(out.get(), data(inp), n) <= epsilon);                                                         \
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
// DWT_IDWT_ROUNDTRIP("coif17") // NOLINT
DWT_IDWT_ROUNDTRIP("sym2")      // NOLINT
DWT_IDWT_ROUNDTRIP("sym3")      // NOLINT
DWT_IDWT_ROUNDTRIP("sym4")      // NOLINT
DWT_IDWT_ROUNDTRIP("sym5")      // NOLINT
DWT_IDWT_ROUNDTRIP("sym6")      // NOLINT
DWT_IDWT_ROUNDTRIP("sym7")      // NOLINT
DWT_IDWT_ROUNDTRIP("sym8")      // NOLINT
DWT_IDWT_ROUNDTRIP("sym9")      // NOLINT
DWT_IDWT_ROUNDTRIP("sym10")     // NOLINT
DWT_IDWT_ROUNDTRIP("sym11")     // NOLINT
DWT_IDWT_ROUNDTRIP("sym12")     // NOLINT
DWT_IDWT_ROUNDTRIP("sym13")     // NOLINT
DWT_IDWT_ROUNDTRIP("sym14")     // NOLINT
DWT_IDWT_ROUNDTRIP("sym15")     // NOLINT
DWT_IDWT_ROUNDTRIP("sym16")     // NOLINT
DWT_IDWT_ROUNDTRIP("sym17")     // NOLINT
DWT_IDWT_ROUNDTRIP("sym18")     // NOLINT
DWT_IDWT_ROUNDTRIP("sym19")     // NOLINT
DWT_IDWT_ROUNDTRIP("sym20")     // NOLINT
DWT_IDWT_ROUNDTRIP("bior1.1")   // NOLINT
DWT_IDWT_ROUNDTRIP("bior1.3")   // NOLINT
DWT_IDWT_ROUNDTRIP("bior1.5")   // NOLINT
DWT_IDWT_ROUNDTRIP("bior2.2")   // NOLINT
DWT_IDWT_ROUNDTRIP("bior2.4")   // NOLINT
DWT_IDWT_ROUNDTRIP("bior2.6")   // NOLINT
DWT_IDWT_ROUNDTRIP("bior2.8")   // NOLINT
DWT_IDWT_ROUNDTRIP("bior3.1")   // NOLINT
DWT_IDWT_ROUNDTRIP("bior3.3")   // NOLINT
DWT_IDWT_ROUNDTRIP("bior3.5")   // NOLINT
DWT_IDWT_ROUNDTRIP("bior3.7")   // NOLINT
DWT_IDWT_ROUNDTRIP("bior3.9")   // NOLINT
DWT_IDWT_ROUNDTRIP("bior4.4")   // NOLINT
DWT_IDWT_ROUNDTRIP("bior5.5")   // NOLINT
DWT_IDWT_ROUNDTRIP("bior6.8")   // NOLINT
DWT_IDWT_ROUNDTRIP("rbior1.1")  // NOLINT
DWT_IDWT_ROUNDTRIP("rbior1.3")  // NOLINT
DWT_IDWT_ROUNDTRIP("rbior1.5")  // NOLINT
DWT_IDWT_ROUNDTRIP("rbior2.2")  // NOLINT
DWT_IDWT_ROUNDTRIP("rbior2.4")  // NOLINT
DWT_IDWT_ROUNDTRIP("rbior2.6")  // NOLINT
DWT_IDWT_ROUNDTRIP("rbior2.8")  // NOLINT
DWT_IDWT_ROUNDTRIP("rbior3.1")  // NOLINT
DWT_IDWT_ROUNDTRIP("rbior3.3")  // NOLINT
DWT_IDWT_ROUNDTRIP("rbior3.5")  // NOLINT
DWT_IDWT_ROUNDTRIP("rbior3.7")  // NOLINT
DWT_IDWT_ROUNDTRIP("rbior3.9")  // NOLINT
DWT_IDWT_ROUNDTRIP("rbior4.4")  // NOLINT
DWT_IDWT_ROUNDTRIP("rbior5.5")  // NOLINT
DWT_IDWT_ROUNDTRIP("rbior6.8")  // NOLINT

#define MODWT_IMODWT_ROUNDTRIP(waveletName)                                                                            \
    TEST_CASE("dsp/wavelet: WaveletTransform2D(modwt/imodwt) - " waveletName, "[dsp][wavelet]")                        \
    {                                                                                                                  \
        auto levels     = GENERATE(1);                                                                                 \
        auto const rows = 256;                                                                                         \
        auto const cols = 256;                                                                                         \
        auto const n    = rows * cols;                                                                                 \
        auto inp        = generateRandomTestData(n);                                                                   \
        auto out        = makeZeros<float>(n);                                                                         \
        auto wavelet    = dsp::Wavelet{waveletName};                                                                   \
        auto wt         = dsp::WaveletTransform2D(wavelet, "modwt", rows, cols, levels);                               \
        setDWT2Extension(wt, "per");                                                                                   \
        auto wavecoeffs = modwt(wt, data(inp));                                                                        \
        imodwt(wt, wavecoeffs.get(), out.get());                                                                       \
        REQUIRE(rmsError(out.get(), data(inp), n) <= epsilon);                                                         \
    }

MODWT_IMODWT_ROUNDTRIP("db1")     // NOLINT
MODWT_IMODWT_ROUNDTRIP("db2")     // NOLINT
MODWT_IMODWT_ROUNDTRIP("db3")     // NOLINT
MODWT_IMODWT_ROUNDTRIP("db4")     // NOLINT
MODWT_IMODWT_ROUNDTRIP("db5")     // NOLINT
MODWT_IMODWT_ROUNDTRIP("db6")     // NOLINT
MODWT_IMODWT_ROUNDTRIP("db7")     // NOLINT
MODWT_IMODWT_ROUNDTRIP("db8")     // NOLINT
MODWT_IMODWT_ROUNDTRIP("db9")     // NOLINT
MODWT_IMODWT_ROUNDTRIP("db10")    // NOLINT
MODWT_IMODWT_ROUNDTRIP("db11")    // NOLINT
MODWT_IMODWT_ROUNDTRIP("db12")    // NOLINT
MODWT_IMODWT_ROUNDTRIP("db13")    // NOLINT
MODWT_IMODWT_ROUNDTRIP("db14")    // NOLINT
MODWT_IMODWT_ROUNDTRIP("db15")    // NOLINT
MODWT_IMODWT_ROUNDTRIP("db16")    // NOLINT
MODWT_IMODWT_ROUNDTRIP("db17")    // NOLINT
MODWT_IMODWT_ROUNDTRIP("db18")    // NOLINT
MODWT_IMODWT_ROUNDTRIP("db19")    // NOLINT
MODWT_IMODWT_ROUNDTRIP("db20")    // NOLINT
MODWT_IMODWT_ROUNDTRIP("db21")    // NOLINT
MODWT_IMODWT_ROUNDTRIP("db22")    // NOLINT
MODWT_IMODWT_ROUNDTRIP("db23")    // NOLINT
MODWT_IMODWT_ROUNDTRIP("db24")    // NOLINT
MODWT_IMODWT_ROUNDTRIP("db25")    // NOLINT
MODWT_IMODWT_ROUNDTRIP("db26")    // NOLINT
MODWT_IMODWT_ROUNDTRIP("db27")    // NOLINT
MODWT_IMODWT_ROUNDTRIP("db28")    // NOLINT
MODWT_IMODWT_ROUNDTRIP("db29")    // NOLINT
MODWT_IMODWT_ROUNDTRIP("db30")    // NOLINT
MODWT_IMODWT_ROUNDTRIP("db31")    // NOLINT
MODWT_IMODWT_ROUNDTRIP("db32")    // NOLINT
MODWT_IMODWT_ROUNDTRIP("db33")    // NOLINT
MODWT_IMODWT_ROUNDTRIP("db34")    // NOLINT
MODWT_IMODWT_ROUNDTRIP("db35")    // NOLINT
MODWT_IMODWT_ROUNDTRIP("coif1")   // NOLINT
MODWT_IMODWT_ROUNDTRIP("coif2")   // NOLINT
MODWT_IMODWT_ROUNDTRIP("coif3")   // NOLINT
MODWT_IMODWT_ROUNDTRIP("coif4")   // NOLINT
MODWT_IMODWT_ROUNDTRIP("coif5")   // NOLINT
MODWT_IMODWT_ROUNDTRIP("coif6")   // NOLINT
MODWT_IMODWT_ROUNDTRIP("coif7")   // NOLINT
MODWT_IMODWT_ROUNDTRIP("coif8")   // NOLINT
MODWT_IMODWT_ROUNDTRIP("coif9")   // NOLINT
MODWT_IMODWT_ROUNDTRIP("coif10")  // NOLINT
MODWT_IMODWT_ROUNDTRIP("coif11")  // NOLINT
MODWT_IMODWT_ROUNDTRIP("coif12")  // NOLINT
MODWT_IMODWT_ROUNDTRIP("coif13")  // NOLINT
MODWT_IMODWT_ROUNDTRIP("coif14")  // NOLINT
MODWT_IMODWT_ROUNDTRIP("coif15")  // NOLINT
MODWT_IMODWT_ROUNDTRIP("coif16")  // NOLINT
// MODWT_IMODWT_ROUNDTRIP("coif17") // NOLINT
MODWT_IMODWT_ROUNDTRIP("sym2")   // NOLINT
MODWT_IMODWT_ROUNDTRIP("sym3")   // NOLINT
MODWT_IMODWT_ROUNDTRIP("sym4")   // NOLINT
MODWT_IMODWT_ROUNDTRIP("sym5")   // NOLINT
MODWT_IMODWT_ROUNDTRIP("sym6")   // NOLINT
MODWT_IMODWT_ROUNDTRIP("sym7")   // NOLINT
MODWT_IMODWT_ROUNDTRIP("sym8")   // NOLINT
MODWT_IMODWT_ROUNDTRIP("sym9")   // NOLINT
MODWT_IMODWT_ROUNDTRIP("sym10")  // NOLINT
MODWT_IMODWT_ROUNDTRIP("sym11")  // NOLINT
MODWT_IMODWT_ROUNDTRIP("sym12")  // NOLINT
MODWT_IMODWT_ROUNDTRIP("sym13")  // NOLINT
MODWT_IMODWT_ROUNDTRIP("sym14")  // NOLINT
MODWT_IMODWT_ROUNDTRIP("sym15")  // NOLINT
MODWT_IMODWT_ROUNDTRIP("sym16")  // NOLINT
MODWT_IMODWT_ROUNDTRIP("sym17")  // NOLINT
MODWT_IMODWT_ROUNDTRIP("sym18")  // NOLINT
MODWT_IMODWT_ROUNDTRIP("sym19")  // NOLINT
MODWT_IMODWT_ROUNDTRIP("sym20")  // NOLINT

#define SWT2_ISWT2_ROUNDTRIP(waveletName)                                                                              \
    TEST_CASE("dsp/wavelet: WaveletTransform2D(swt2/iswt2) - " waveletName, "[dsp][wavelet]")                          \
    {                                                                                                                  \
        auto levels     = GENERATE(1, 2);                                                                              \
        auto const rows = 512;                                                                                         \
        auto const cols = 500;                                                                                         \
        auto const n    = rows * cols;                                                                                 \
        auto inp        = generateRandomTestData(n);                                                                   \
        auto out        = makeZeros<float>(n);                                                                         \
        auto wavelet    = dsp::Wavelet{waveletName};                                                                   \
        auto wt         = dsp::WaveletTransform2D(wavelet, "swt", rows, cols, levels);                                 \
        setDWT2Extension(wt, "per");                                                                                   \
        auto wavecoeffs = swt2(wt, data(inp));                                                                         \
        iswt2(wt, wavecoeffs.get(), out.get());                                                                        \
        REQUIRE(rmsError(out.get(), data(inp), n) <= epsilon);                                                         \
    }

SWT2_ISWT2_ROUNDTRIP("db1")     // NOLINT
SWT2_ISWT2_ROUNDTRIP("db2")     // NOLINT
SWT2_ISWT2_ROUNDTRIP("db3")     // NOLINT
SWT2_ISWT2_ROUNDTRIP("db4")     // NOLINT
SWT2_ISWT2_ROUNDTRIP("db5")     // NOLINT
SWT2_ISWT2_ROUNDTRIP("db6")     // NOLINT
SWT2_ISWT2_ROUNDTRIP("db7")     // NOLINT
SWT2_ISWT2_ROUNDTRIP("db8")     // NOLINT
SWT2_ISWT2_ROUNDTRIP("db9")     // NOLINT
SWT2_ISWT2_ROUNDTRIP("db10")    // NOLINT
SWT2_ISWT2_ROUNDTRIP("db11")    // NOLINT
SWT2_ISWT2_ROUNDTRIP("db12")    // NOLINT
SWT2_ISWT2_ROUNDTRIP("db13")    // NOLINT
SWT2_ISWT2_ROUNDTRIP("db14")    // NOLINT
SWT2_ISWT2_ROUNDTRIP("db15")    // NOLINT
SWT2_ISWT2_ROUNDTRIP("db16")    // NOLINT
SWT2_ISWT2_ROUNDTRIP("db17")    // NOLINT
SWT2_ISWT2_ROUNDTRIP("db18")    // NOLINT
SWT2_ISWT2_ROUNDTRIP("db19")    // NOLINT
SWT2_ISWT2_ROUNDTRIP("db20")    // NOLINT
SWT2_ISWT2_ROUNDTRIP("db21")    // NOLINT
SWT2_ISWT2_ROUNDTRIP("db22")    // NOLINT
SWT2_ISWT2_ROUNDTRIP("db23")    // NOLINT
SWT2_ISWT2_ROUNDTRIP("db24")    // NOLINT
SWT2_ISWT2_ROUNDTRIP("db25")    // NOLINT
SWT2_ISWT2_ROUNDTRIP("db26")    // NOLINT
SWT2_ISWT2_ROUNDTRIP("db27")    // NOLINT
SWT2_ISWT2_ROUNDTRIP("db28")    // NOLINT
SWT2_ISWT2_ROUNDTRIP("db29")    // NOLINT
SWT2_ISWT2_ROUNDTRIP("db30")    // NOLINT
SWT2_ISWT2_ROUNDTRIP("db31")    // NOLINT
SWT2_ISWT2_ROUNDTRIP("db32")    // NOLINT
SWT2_ISWT2_ROUNDTRIP("db33")    // NOLINT
SWT2_ISWT2_ROUNDTRIP("db34")    // NOLINT
SWT2_ISWT2_ROUNDTRIP("db35")    // NOLINT
SWT2_ISWT2_ROUNDTRIP("coif1")   // NOLINT
SWT2_ISWT2_ROUNDTRIP("coif2")   // NOLINT
SWT2_ISWT2_ROUNDTRIP("coif3")   // NOLINT
SWT2_ISWT2_ROUNDTRIP("coif4")   // NOLINT
SWT2_ISWT2_ROUNDTRIP("coif5")   // NOLINT
SWT2_ISWT2_ROUNDTRIP("coif6")   // NOLINT
SWT2_ISWT2_ROUNDTRIP("coif7")   // NOLINT
SWT2_ISWT2_ROUNDTRIP("coif8")   // NOLINT
SWT2_ISWT2_ROUNDTRIP("coif9")   // NOLINT
SWT2_ISWT2_ROUNDTRIP("coif10")  // NOLINT
SWT2_ISWT2_ROUNDTRIP("coif11")  // NOLINT
SWT2_ISWT2_ROUNDTRIP("coif12")  // NOLINT
SWT2_ISWT2_ROUNDTRIP("coif13")  // NOLINT
SWT2_ISWT2_ROUNDTRIP("coif14")  // NOLINT
SWT2_ISWT2_ROUNDTRIP("coif15")  // NOLINT
SWT2_ISWT2_ROUNDTRIP("coif16")  // NOLINT
// SWT2_ISWT2_ROUNDTRIP("coif17") // NOLINT
SWT2_ISWT2_ROUNDTRIP("sym2")      // NOLINT
SWT2_ISWT2_ROUNDTRIP("sym3")      // NOLINT
SWT2_ISWT2_ROUNDTRIP("sym4")      // NOLINT
SWT2_ISWT2_ROUNDTRIP("sym5")      // NOLINT
SWT2_ISWT2_ROUNDTRIP("sym6")      // NOLINT
SWT2_ISWT2_ROUNDTRIP("sym7")      // NOLINT
SWT2_ISWT2_ROUNDTRIP("sym8")      // NOLINT
SWT2_ISWT2_ROUNDTRIP("sym9")      // NOLINT
SWT2_ISWT2_ROUNDTRIP("sym10")     // NOLINT
SWT2_ISWT2_ROUNDTRIP("sym11")     // NOLINT
SWT2_ISWT2_ROUNDTRIP("sym12")     // NOLINT
SWT2_ISWT2_ROUNDTRIP("sym13")     // NOLINT
SWT2_ISWT2_ROUNDTRIP("sym14")     // NOLINT
SWT2_ISWT2_ROUNDTRIP("sym15")     // NOLINT
SWT2_ISWT2_ROUNDTRIP("sym16")     // NOLINT
SWT2_ISWT2_ROUNDTRIP("sym17")     // NOLINT
SWT2_ISWT2_ROUNDTRIP("sym18")     // NOLINT
SWT2_ISWT2_ROUNDTRIP("sym19")     // NOLINT
SWT2_ISWT2_ROUNDTRIP("sym20")     // NOLINT
SWT2_ISWT2_ROUNDTRIP("bior1.1")   // NOLINT
SWT2_ISWT2_ROUNDTRIP("bior1.3")   // NOLINT
SWT2_ISWT2_ROUNDTRIP("bior1.5")   // NOLINT
SWT2_ISWT2_ROUNDTRIP("bior2.2")   // NOLINT
SWT2_ISWT2_ROUNDTRIP("bior2.4")   // NOLINT
SWT2_ISWT2_ROUNDTRIP("bior2.6")   // NOLINT
SWT2_ISWT2_ROUNDTRIP("bior2.8")   // NOLINT
SWT2_ISWT2_ROUNDTRIP("bior3.1")   // NOLINT
SWT2_ISWT2_ROUNDTRIP("bior3.3")   // NOLINT
SWT2_ISWT2_ROUNDTRIP("bior3.5")   // NOLINT
SWT2_ISWT2_ROUNDTRIP("bior3.7")   // NOLINT
SWT2_ISWT2_ROUNDTRIP("bior3.9")   // NOLINT
SWT2_ISWT2_ROUNDTRIP("bior4.4")   // NOLINT
SWT2_ISWT2_ROUNDTRIP("bior5.5")   // NOLINT
SWT2_ISWT2_ROUNDTRIP("bior6.8")   // NOLINT
SWT2_ISWT2_ROUNDTRIP("rbior1.1")  // NOLINT
SWT2_ISWT2_ROUNDTRIP("rbior1.3")  // NOLINT
SWT2_ISWT2_ROUNDTRIP("rbior1.5")  // NOLINT
SWT2_ISWT2_ROUNDTRIP("rbior2.2")  // NOLINT
SWT2_ISWT2_ROUNDTRIP("rbior2.4")  // NOLINT
SWT2_ISWT2_ROUNDTRIP("rbior2.6")  // NOLINT
SWT2_ISWT2_ROUNDTRIP("rbior2.8")  // NOLINT
SWT2_ISWT2_ROUNDTRIP("rbior3.1")  // NOLINT
SWT2_ISWT2_ROUNDTRIP("rbior3.3")  // NOLINT
SWT2_ISWT2_ROUNDTRIP("rbior3.5")  // NOLINT
SWT2_ISWT2_ROUNDTRIP("rbior3.7")  // NOLINT
SWT2_ISWT2_ROUNDTRIP("rbior3.9")  // NOLINT
SWT2_ISWT2_ROUNDTRIP("rbior4.4")  // NOLINT
SWT2_ISWT2_ROUNDTRIP("rbior5.5")  // NOLINT
SWT2_ISWT2_ROUNDTRIP("rbior6.8")  // NOLINT