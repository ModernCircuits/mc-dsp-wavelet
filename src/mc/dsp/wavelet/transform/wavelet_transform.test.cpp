#include <mc/dsp/wavelet.hpp>

#include <mc/core/memory.hpp>
#include <mc/testing/test.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace mc;

static constexpr auto const epsilon = 6e-6F;

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define DWT_IDWT_ROUNDTRIP(waveletName)                                                    \
    TEST_CASE("dsp/wavelet: WaveletTransform(dwt/idwt) - " waveletName, "[dsp][wavelet]")  \
    {                                                                                      \
        using namespace mc::dsp;                                                           \
        auto method    = GENERATE(ConvolutionMethod::fft, ConvolutionMethod::direct);      \
        auto extension = GENERATE(SignalExtension::periodic, SignalExtension::symmetric);  \
        auto levels    = GENERATE(as<std::size_t>{}, 1, 2, 3);                             \
        auto const n   = 11'025;                                                           \
        auto const inp = generateRandomTestData(n);                                        \
        auto out       = makeZeros<float>(n);                                              \
        auto wavelet   = Wavelet{waveletName};                                             \
        auto wt        = WaveletTransform(wavelet, "dwt", n, levels);                      \
        wt.extension(extension);                                                           \
        wt.convMethod(method);                                                             \
        dwt(wt, data(inp));                                                                \
        idwt(wt, out.get());                                                               \
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

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define SWT_ISWT_ROUNDTRIP(waveletName)                                                    \
    TEST_CASE("dsp/wavelet: WaveletTransform(swt/iswt) - " waveletName, "[dsp][wavelet]")  \
    {                                                                                      \
        auto convolutionMethod                                                             \
            = GENERATE(dsp::ConvolutionMethod::fft, dsp::ConvolutionMethod::direct);       \
        auto extension = GENERATE(dsp::SignalExtension::periodic);                         \
        auto levels    = GENERATE(as<std::size_t>{}, 1, 2, 3);                             \
        auto const n   = 4096U;                                                            \
        auto const inp = generateRandomTestData(n);                                        \
        auto out       = makeZeros<float>(n);                                              \
        auto wavelet   = dsp::Wavelet{waveletName};                                        \
        auto wt        = dsp::WaveletTransform(wavelet, "swt", n, levels);                 \
        wt.extension(extension);                                                           \
        wt.convMethod(convolutionMethod);                                                  \
        swt(wt, data(inp));                                                                \
        iswt(wt, out.get());                                                               \
        REQUIRE_THAT(                                                                      \
            rmsError(out.get(), data(inp), wt.signalLength()),                             \
            Catch::Matchers::WithinAbs(0.0F, epsilon)                                      \
        );                                                                                 \
    }

SWT_ISWT_ROUNDTRIP("db1")     // NOLINT
SWT_ISWT_ROUNDTRIP("db2")     // NOLINT
SWT_ISWT_ROUNDTRIP("db3")     // NOLINT
SWT_ISWT_ROUNDTRIP("db4")     // NOLINT
SWT_ISWT_ROUNDTRIP("db5")     // NOLINT
SWT_ISWT_ROUNDTRIP("db6")     // NOLINT
SWT_ISWT_ROUNDTRIP("db7")     // NOLINT
SWT_ISWT_ROUNDTRIP("db8")     // NOLINT
SWT_ISWT_ROUNDTRIP("db9")     // NOLINT
SWT_ISWT_ROUNDTRIP("db10")    // NOLINT
SWT_ISWT_ROUNDTRIP("db11")    // NOLINT
SWT_ISWT_ROUNDTRIP("db12")    // NOLINT
SWT_ISWT_ROUNDTRIP("db13")    // NOLINT
SWT_ISWT_ROUNDTRIP("db14")    // NOLINT
SWT_ISWT_ROUNDTRIP("db15")    // NOLINT
SWT_ISWT_ROUNDTRIP("db16")    // NOLINT
SWT_ISWT_ROUNDTRIP("db17")    // NOLINT
SWT_ISWT_ROUNDTRIP("db18")    // NOLINT
SWT_ISWT_ROUNDTRIP("db19")    // NOLINT
SWT_ISWT_ROUNDTRIP("db20")    // NOLINT
SWT_ISWT_ROUNDTRIP("db21")    // NOLINT
SWT_ISWT_ROUNDTRIP("db22")    // NOLINT
SWT_ISWT_ROUNDTRIP("db23")    // NOLINT
SWT_ISWT_ROUNDTRIP("db24")    // NOLINT
SWT_ISWT_ROUNDTRIP("db25")    // NOLINT
SWT_ISWT_ROUNDTRIP("db26")    // NOLINT
SWT_ISWT_ROUNDTRIP("db27")    // NOLINT
SWT_ISWT_ROUNDTRIP("db28")    // NOLINT
SWT_ISWT_ROUNDTRIP("db29")    // NOLINT
SWT_ISWT_ROUNDTRIP("db30")    // NOLINT
SWT_ISWT_ROUNDTRIP("db31")    // NOLINT
SWT_ISWT_ROUNDTRIP("db32")    // NOLINT
SWT_ISWT_ROUNDTRIP("db33")    // NOLINT
SWT_ISWT_ROUNDTRIP("db34")    // NOLINT
SWT_ISWT_ROUNDTRIP("db35")    // NOLINT
SWT_ISWT_ROUNDTRIP("coif1")   // NOLINT
SWT_ISWT_ROUNDTRIP("coif2")   // NOLINT
SWT_ISWT_ROUNDTRIP("coif3")   // NOLINT
SWT_ISWT_ROUNDTRIP("coif4")   // NOLINT
SWT_ISWT_ROUNDTRIP("coif5")   // NOLINT
SWT_ISWT_ROUNDTRIP("coif6")   // NOLINT
SWT_ISWT_ROUNDTRIP("coif7")   // NOLINT
SWT_ISWT_ROUNDTRIP("coif8")   // NOLINT
SWT_ISWT_ROUNDTRIP("coif9")   // NOLINT
SWT_ISWT_ROUNDTRIP("coif10")  // NOLINT
SWT_ISWT_ROUNDTRIP("coif11")  // NOLINT
SWT_ISWT_ROUNDTRIP("coif12")  // NOLINT
SWT_ISWT_ROUNDTRIP("coif13")  // NOLINT
SWT_ISWT_ROUNDTRIP("coif14")  // NOLINT
SWT_ISWT_ROUNDTRIP("coif15")  // NOLINT
SWT_ISWT_ROUNDTRIP("coif16")  // NOLINT
SWT_ISWT_ROUNDTRIP("coif17")  // NOLINT
SWT_ISWT_ROUNDTRIP("sym2")    // NOLINT
SWT_ISWT_ROUNDTRIP("sym3")    // NOLINT
SWT_ISWT_ROUNDTRIP("sym4")    // NOLINT
SWT_ISWT_ROUNDTRIP("sym5")    // NOLINT
SWT_ISWT_ROUNDTRIP("sym6")    // NOLINT
SWT_ISWT_ROUNDTRIP("sym7")    // NOLINT
SWT_ISWT_ROUNDTRIP("sym8")    // NOLINT
SWT_ISWT_ROUNDTRIP("sym9")    // NOLINT
SWT_ISWT_ROUNDTRIP("sym10")   // NOLINT
SWT_ISWT_ROUNDTRIP("sym11")   // NOLINT
SWT_ISWT_ROUNDTRIP("sym12")   // NOLINT
SWT_ISWT_ROUNDTRIP("sym13")   // NOLINT
SWT_ISWT_ROUNDTRIP("sym14")   // NOLINT
SWT_ISWT_ROUNDTRIP("sym15")   // NOLINT
SWT_ISWT_ROUNDTRIP("sym16")   // NOLINT
SWT_ISWT_ROUNDTRIP("sym17")   // NOLINT
SWT_ISWT_ROUNDTRIP("sym18")   // NOLINT
SWT_ISWT_ROUNDTRIP("sym19")   // NOLINT
SWT_ISWT_ROUNDTRIP("sym20")   // NOLINT

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define MODWT_IMODWT_ROUNDTRIP(waveletName)                                                \
    TEST_CASE(                                                                             \
        "dsp/wavelet: WaveletTransform(modwt/imodwt) - " waveletName,                      \
        "[dsp][wavelet]"                                                                   \
    )                                                                                      \
    {                                                                                      \
        auto convolutionMethod                                                             \
            = GENERATE(dsp::ConvolutionMethod::fft, dsp::ConvolutionMethod::direct);       \
        auto extension = GENERATE(dsp::SignalExtension::periodic);                         \
        auto levels    = GENERATE(as<std::size_t>{}, 1, 2);                                \
        auto const n   = 4096U;                                                            \
        auto const inp = generateRandomTestData(n);                                        \
        auto out       = makeZeros<float>(n);                                              \
        auto wavelet   = dsp::Wavelet{waveletName};                                        \
        auto wt        = dsp::WaveletTransform(wavelet, "modwt", n, levels);               \
        wt.extension(extension);                                                           \
        wt.convMethod(convolutionMethod);                                                  \
        modwt(wt, data(inp));                                                              \
        imodwt(wt, out.get());                                                             \
        REQUIRE_THAT(                                                                      \
            rmsError(out.get(), data(inp), wt.signalLength()),                             \
            Catch::Matchers::WithinAbs(0.0F, epsilon)                                      \
        );                                                                                 \
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
MODWT_IMODWT_ROUNDTRIP("coif17")  // NOLINT
MODWT_IMODWT_ROUNDTRIP("sym2")    // NOLINT
MODWT_IMODWT_ROUNDTRIP("sym3")    // NOLINT
MODWT_IMODWT_ROUNDTRIP("sym4")    // NOLINT
MODWT_IMODWT_ROUNDTRIP("sym5")    // NOLINT
MODWT_IMODWT_ROUNDTRIP("sym6")    // NOLINT
MODWT_IMODWT_ROUNDTRIP("sym7")    // NOLINT
MODWT_IMODWT_ROUNDTRIP("sym8")    // NOLINT
MODWT_IMODWT_ROUNDTRIP("sym9")    // NOLINT
MODWT_IMODWT_ROUNDTRIP("sym10")   // NOLINT
MODWT_IMODWT_ROUNDTRIP("sym11")   // NOLINT
MODWT_IMODWT_ROUNDTRIP("sym12")   // NOLINT
MODWT_IMODWT_ROUNDTRIP("sym13")   // NOLINT
MODWT_IMODWT_ROUNDTRIP("sym14")   // NOLINT
MODWT_IMODWT_ROUNDTRIP("sym15")   // NOLINT
MODWT_IMODWT_ROUNDTRIP("sym16")   // NOLINT
MODWT_IMODWT_ROUNDTRIP("sym17")   // NOLINT
MODWT_IMODWT_ROUNDTRIP("sym18")   // NOLINT
MODWT_IMODWT_ROUNDTRIP("sym19")   // NOLINT
MODWT_IMODWT_ROUNDTRIP("sym20")   // NOLINT
