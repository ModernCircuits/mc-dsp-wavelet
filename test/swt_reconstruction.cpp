#include "mc/dsp/wavelets.hpp"

#include "mc/cmath.hpp"
#include "mc/cstdlib.hpp"
#include "mc/cstring.hpp"
#include "mc/format.hpp"
#include "mc/random.hpp"
#include "mc/vector.hpp"

#include "mc/testing/test.hpp"

namespace dsp = mc::dsp;

auto swtReconstructionTest()
{

    int i         = 0;
    float epsilon = 1e-6;
    float err     = NAN;

    auto const n = 4000;

    // N = 256;

    auto inp = makeZeros<float>(n);
    auto out = makeZeros<float>(n);

    std::random_device rd{};
    auto gen = std::mt19937{rd()};
    auto dis = std::uniform_real_distribution<float>{0.0F, 1.0F};

    for (i = 0; i < n; ++i) { inp[i] = dis(gen); }
    std::vector<std::string> waveletNames;

    for (unsigned int j = 0; j < 15; j++) { waveletNames.push_back(std::string("db") + std::to_string(j + 1)); }
    for (unsigned int j = 0; j < 5; j++) { waveletNames.push_back(std::string("coif") + std::to_string(j + 1)); }
    for (unsigned int j = 1; j < 10; j++) { waveletNames.push_back(std::string("sym") + std::to_string(j + 1)); }

    waveletNames.emplace_back("bior1.1");
    waveletNames.emplace_back("bior1.3");
    waveletNames.emplace_back("bior1.5");
    waveletNames.emplace_back("bior2.2");
    waveletNames.emplace_back("bior2.4");
    waveletNames.emplace_back("bior2.6");
    waveletNames.emplace_back("bior2.8");
    waveletNames.emplace_back("bior3.1");
    waveletNames.emplace_back("bior3.3");
    waveletNames.emplace_back("bior3.5");
    waveletNames.emplace_back("bior3.7");
    waveletNames.emplace_back("bior3.9");
    waveletNames.emplace_back("bior4.4");
    waveletNames.emplace_back("bior5.5");
    waveletNames.emplace_back("bior6.8");

    waveletNames.emplace_back("rbior1.1");
    waveletNames.emplace_back("rbior1.3");
    waveletNames.emplace_back("rbior1.5");
    waveletNames.emplace_back("rbior2.2");
    waveletNames.emplace_back("rbior2.4");
    waveletNames.emplace_back("rbior2.6");
    waveletNames.emplace_back("rbior2.8");
    waveletNames.emplace_back("rbior3.1");
    waveletNames.emplace_back("rbior3.3");
    waveletNames.emplace_back("rbior3.5");
    waveletNames.emplace_back("rbior3.7");
    waveletNames.emplace_back("rbior3.9");
    waveletNames.emplace_back("rbior4.4");
    waveletNames.emplace_back("rbior5.5");
    waveletNames.emplace_back("rbior6.8");

    for (unsigned int directFft = 0; directFft < 2; directFft++)
    {
        for (unsigned int symPer = 0; symPer < 1; symPer++)
        {
            for (auto& name : waveletNames)
            {
                auto obj = dsp::Wavelet{name.c_str()};
                for (auto j = 1; j < 3; j++)
                {
                    auto wt = dsp::WaveletTransform(obj, "swt", n, j);

                    if (directFft == 0) { wt.convMethod(dsp::ConvolutionMethod::direct); }
                    else { wt.convMethod(dsp::ConvolutionMethod::fft); }

                    if (symPer == 0) { wt.extension(dsp::SignalExtension::periodic); }
                    else if (symPer == 1 && directFft == 1) { wt.extension(dsp::SignalExtension::symmetric); }
                    else { MC_REQUIRE(false); }

                    swt(wt, inp.get());
                    iswt(wt, out.get());

                    if (directFft == 0) { epsilon = 1e-5; }
                    else { epsilon = 1e-5; }

                    err = rmsError(out.get(), inp.get(), wt.signalLength());
                    MC_REQUIRE(err <= epsilon);
                }
            }
        }
    }
}

auto swT2ReconstructionTest()
{

    auto const rows = 512;
    auto const cols = 500;

    auto const n = rows * cols;

    auto inp = makeZeros<float>(n);
    auto out = makeZeros<float>(n);

    std::vector<std::string> waveletNames;

    for (unsigned int j = 0; j < 15; j++) { waveletNames.push_back(std::string("db") + std::to_string(j + 1)); }
    for (unsigned int j = 0; j < 5; j++) { waveletNames.push_back(std::string("coif") + std::to_string(j + 1)); }
    for (unsigned int j = 1; j < 10; j++) { waveletNames.push_back(std::string("sym") + std::to_string(j + 1)); }

    waveletNames.emplace_back("bior1.1");
    waveletNames.emplace_back("bior1.3");
    waveletNames.emplace_back("bior1.5");
    waveletNames.emplace_back("bior2.2");
    waveletNames.emplace_back("bior2.4");
    waveletNames.emplace_back("bior2.6");
    waveletNames.emplace_back("bior2.8");
    waveletNames.emplace_back("bior3.1");
    waveletNames.emplace_back("bior3.3");
    waveletNames.emplace_back("bior3.5");
    waveletNames.emplace_back("bior3.7");
    waveletNames.emplace_back("bior3.9");
    waveletNames.emplace_back("bior4.4");
    waveletNames.emplace_back("bior5.5");
    waveletNames.emplace_back("bior6.8");

    waveletNames.emplace_back("rbior1.1");
    waveletNames.emplace_back("rbior1.3");
    waveletNames.emplace_back("rbior1.5");
    waveletNames.emplace_back("rbior2.2");
    waveletNames.emplace_back("rbior2.4");
    waveletNames.emplace_back("rbior2.6");
    waveletNames.emplace_back("rbior2.8");
    waveletNames.emplace_back("rbior3.1");
    waveletNames.emplace_back("rbior3.3");
    waveletNames.emplace_back("rbior3.5");
    waveletNames.emplace_back("rbior3.7");
    waveletNames.emplace_back("rbior3.9");
    waveletNames.emplace_back("rbior4.4");
    waveletNames.emplace_back("rbior5.5");
    waveletNames.emplace_back("rbior6.8");

    for (auto i = 0; i < rows; ++i)
    {
        for (auto k = 0; k < cols; ++k)
        {
            // inp[i*cols + k] = i*cols + k;
            inp[i * cols + k] = generateRnd();
            out[i * cols + k] = 0.0F;
        }
    }

    for (unsigned int directFft = 0; directFft < 1; directFft++)
    {
        for (unsigned int symPer = 0; symPer < 1; symPer++)
        {
            for (auto& name : waveletNames)
            {
                auto obj = dsp::Wavelet{name.c_str()};
                for (auto j = 1; j < 3; j++)
                {
                    auto wt = dsp::WaveletTransform2D(obj, "swt", rows, cols, j);
                    if (symPer == 0)
                    {
                        setDWT2Extension(wt, "per");  // Options are "per"
                    }

                    auto wavecoeffs = swt2(wt, inp.get());
                    iswt2(wt, wavecoeffs.get(), out.get());

                    float epsilon{0.0F};
                    if (directFft == 0) { epsilon = 1e-4; }
                    else { epsilon = 1e-4; }

                    MC_REQUIRE(rmsError(out.get(), inp.get(), n) <= epsilon);
                }
            }
        }
    }
}

auto main() -> int
{
    fmt::print("Running SWT ReconstructionTests ... ");
    swtReconstructionTest();
    fmt::print("DONE \n");

    fmt::print("Running SWT2 ReconstructionTests ... ");
    swT2ReconstructionTest();
    fmt::print("DONE \n");

    return EXIT_SUCCESS;
}