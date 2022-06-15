#include "mc/dsp/wavelets.hpp"

#include "mc/memory.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <random>

namespace dsp = mc::dsp;

static auto rmsError(float const* data, float const* rec, std::size_t n) -> float
{
    float sum = 0;
    for (std::size_t i = 0; i < n; ++i) { sum += (data[i] - rec[i]) * (data[i] - rec[i]); }
    return std::sqrt(sum / ((float)n - 1));
}

static auto generateRandomTestData(std::size_t n) -> std::vector<float>
{
    std::vector<float> data(n);
    auto rd  = std::random_device{};
    auto gen = std::mt19937{rd()};
    auto dis = std::uniform_real_distribution<float>{-1.0F, 1.0F};
    std::generate(data.begin(), data.end(), [&]() { return dis(gen); });
    return data;
}
static constexpr auto const epsilon = 6e-7f;

#define DWT_IDWT_ROUNDTRIP(waveletName)                                                                                \
    TEST_CASE("dsp/wavelet: WaveletTransform2D(dwt/idwt) - " waveletName, "[dsp][wavelet]")                            \
    {                                                                                                                  \
        auto extension  = GENERATE(dsp::SignalExtension::periodic, dsp::SignalExtension::symmetric);                   \
        auto levels     = GENERATE(1, 2, 3);                                                                           \
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

DWT_IDWT_ROUNDTRIP("db1")
DWT_IDWT_ROUNDTRIP("db2")
DWT_IDWT_ROUNDTRIP("db3")
DWT_IDWT_ROUNDTRIP("db4")
DWT_IDWT_ROUNDTRIP("db5")
DWT_IDWT_ROUNDTRIP("db6")
DWT_IDWT_ROUNDTRIP("db7")
DWT_IDWT_ROUNDTRIP("db8")
DWT_IDWT_ROUNDTRIP("db9")
DWT_IDWT_ROUNDTRIP("db10")
DWT_IDWT_ROUNDTRIP("db11")
DWT_IDWT_ROUNDTRIP("db12")
DWT_IDWT_ROUNDTRIP("db13")
DWT_IDWT_ROUNDTRIP("db14")
DWT_IDWT_ROUNDTRIP("db15")
DWT_IDWT_ROUNDTRIP("db16")
DWT_IDWT_ROUNDTRIP("db17")
DWT_IDWT_ROUNDTRIP("db18")
DWT_IDWT_ROUNDTRIP("db19")
DWT_IDWT_ROUNDTRIP("db20")
DWT_IDWT_ROUNDTRIP("db21")
DWT_IDWT_ROUNDTRIP("db22")
DWT_IDWT_ROUNDTRIP("db23")
DWT_IDWT_ROUNDTRIP("db24")
DWT_IDWT_ROUNDTRIP("db25")
DWT_IDWT_ROUNDTRIP("db26")
DWT_IDWT_ROUNDTRIP("db27")
DWT_IDWT_ROUNDTRIP("db28")
DWT_IDWT_ROUNDTRIP("db29")
DWT_IDWT_ROUNDTRIP("db30")
DWT_IDWT_ROUNDTRIP("db31")
DWT_IDWT_ROUNDTRIP("db32")
DWT_IDWT_ROUNDTRIP("db33")
DWT_IDWT_ROUNDTRIP("db34")
DWT_IDWT_ROUNDTRIP("db35")
DWT_IDWT_ROUNDTRIP("coif1")
DWT_IDWT_ROUNDTRIP("coif2")
DWT_IDWT_ROUNDTRIP("coif3")
DWT_IDWT_ROUNDTRIP("coif4")
DWT_IDWT_ROUNDTRIP("coif5")
DWT_IDWT_ROUNDTRIP("coif6")
DWT_IDWT_ROUNDTRIP("coif7")
DWT_IDWT_ROUNDTRIP("coif8")
DWT_IDWT_ROUNDTRIP("coif9")
DWT_IDWT_ROUNDTRIP("coif10")
DWT_IDWT_ROUNDTRIP("coif11")
DWT_IDWT_ROUNDTRIP("coif12")
DWT_IDWT_ROUNDTRIP("coif13")
DWT_IDWT_ROUNDTRIP("coif14")
DWT_IDWT_ROUNDTRIP("coif15")
DWT_IDWT_ROUNDTRIP("coif16")
// DWT_IDWT_ROUNDTRIP("coif17")
DWT_IDWT_ROUNDTRIP("sym2")
DWT_IDWT_ROUNDTRIP("sym3")
DWT_IDWT_ROUNDTRIP("sym4")
DWT_IDWT_ROUNDTRIP("sym5")
DWT_IDWT_ROUNDTRIP("sym6")
DWT_IDWT_ROUNDTRIP("sym7")
DWT_IDWT_ROUNDTRIP("sym8")
DWT_IDWT_ROUNDTRIP("sym9")
DWT_IDWT_ROUNDTRIP("sym10")
DWT_IDWT_ROUNDTRIP("sym11")
DWT_IDWT_ROUNDTRIP("sym12")
DWT_IDWT_ROUNDTRIP("sym13")
DWT_IDWT_ROUNDTRIP("sym14")
DWT_IDWT_ROUNDTRIP("sym15")
DWT_IDWT_ROUNDTRIP("sym16")
DWT_IDWT_ROUNDTRIP("sym17")
DWT_IDWT_ROUNDTRIP("sym18")
DWT_IDWT_ROUNDTRIP("sym19")
DWT_IDWT_ROUNDTRIP("sym20")
DWT_IDWT_ROUNDTRIP("bior1.1")
DWT_IDWT_ROUNDTRIP("bior1.3")
DWT_IDWT_ROUNDTRIP("bior1.5")
DWT_IDWT_ROUNDTRIP("bior2.2")
DWT_IDWT_ROUNDTRIP("bior2.4")
DWT_IDWT_ROUNDTRIP("bior2.6")
DWT_IDWT_ROUNDTRIP("bior2.8")
DWT_IDWT_ROUNDTRIP("bior3.1")
DWT_IDWT_ROUNDTRIP("bior3.3")
DWT_IDWT_ROUNDTRIP("bior3.5")
DWT_IDWT_ROUNDTRIP("bior3.7")
DWT_IDWT_ROUNDTRIP("bior3.9")
DWT_IDWT_ROUNDTRIP("bior4.4")
DWT_IDWT_ROUNDTRIP("bior5.5")
DWT_IDWT_ROUNDTRIP("bior6.8")
DWT_IDWT_ROUNDTRIP("rbior1.1")
DWT_IDWT_ROUNDTRIP("rbior1.3")
DWT_IDWT_ROUNDTRIP("rbior1.5")
DWT_IDWT_ROUNDTRIP("rbior2.2")
DWT_IDWT_ROUNDTRIP("rbior2.4")
DWT_IDWT_ROUNDTRIP("rbior2.6")
DWT_IDWT_ROUNDTRIP("rbior2.8")
DWT_IDWT_ROUNDTRIP("rbior3.1")
DWT_IDWT_ROUNDTRIP("rbior3.3")
DWT_IDWT_ROUNDTRIP("rbior3.5")
DWT_IDWT_ROUNDTRIP("rbior3.7")
DWT_IDWT_ROUNDTRIP("rbior3.9")
DWT_IDWT_ROUNDTRIP("rbior4.4")
DWT_IDWT_ROUNDTRIP("rbior5.5")
DWT_IDWT_ROUNDTRIP("rbior6.8")