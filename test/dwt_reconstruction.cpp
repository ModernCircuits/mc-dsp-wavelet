#include "lt/dsp/wavelets.hpp"

#include "lt/cmath.hpp"
#include "lt/cstdlib.hpp"
#include "lt/cstring.hpp"
#include "lt/format.hpp"
#include "lt/memory.hpp"
#include "lt/random.hpp"
#include "lt/vector.hpp"

#include "lt/testing/test.hpp"

#include "lt/sstream.hpp"

auto dwtReconstructionTest()
{
    auto const n = 22'050;

    auto inp = std::make_unique<float[]>(n);
    auto out = std::make_unique<float[]>(n);

    std::random_device rd {};
    auto gen = std::mt19937 { rd() };
    auto dis = std::uniform_real_distribution<float> { 0.0, 1.0 };

    for (auto i = 0; i < n; ++i) {
        inp[i] = dis(gen);
    }
    std::vector<std::string> waveletNames;

    for (unsigned int j = 0; j < 36; j++) {
        waveletNames.push_back(std::string("db") + std::to_string(j + 1));
    }
    for (unsigned int j = 0; j < 17; j++) {
        waveletNames.push_back(std::string("coif") + std::to_string(j + 1));
    }
    for (unsigned int j = 1; j < 20; j++) {
        waveletNames.push_back(std::string("sym") + std::to_string(j + 1));
    }

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

    for (unsigned int directFft = 0; directFft < 2; directFft++) {
        for (unsigned int symPer = 0; symPer < 2; symPer++) {
            for (auto& waveletName : waveletNames) {
                auto obj = Wavelet { waveletName.c_str() };
                for (auto j = 1; j < 3; j++) {
                    auto wt = WaveletTransform(obj, "dwt", n, j);
                    if (symPer == 0) {
                        wt.extension(SignalExtension::symmetric);
                    } else {
                        wt.extension(SignalExtension::periodic);
                    }
                    if (directFft == 0) {
                        wt.convMethod(ConvolutionMethod::direct);
                    } else {
                        wt.convMethod(ConvolutionMethod::fft);
                    }

                    dwt(wt, inp.get()); // Perform DWT

                    idwt(wt, out.get()); // Perform IDWT (if needed)
                    // Test Reconstruction

                    auto epsilon = 1e-5;
                    if (directFft == 0) {
                        epsilon = 1e-5;
                    } else {
                        epsilon = 1e-5;
                    }

                    REQUIRE(rmsError(out.get(), inp.get(), wt.signalLength()) <= epsilon);
                }
            }
        }
    }
}

auto dwT2ReconstructionTest()
{
    float epsilon = NAN;

    auto const rows = 256;
    auto const cols = 200;

    auto const n = rows * cols;

    auto inp = std::make_unique<float[]>(n);
    auto out = std::make_unique<float[]>(n);

    std::vector<std::string> waveletNames;

    for (unsigned int j = 0; j < 15; j++) {
        waveletNames.push_back(std::string("db") + std::to_string(j + 1));
    }
    for (unsigned int j = 0; j < 5; j++) {
        waveletNames.push_back(std::string("coif") + std::to_string(j + 1));
    }
    for (unsigned int j = 1; j < 10; j++) {
        waveletNames.push_back(std::string("sym") + std::to_string(j + 1));
    }

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

    for (auto i = 0; i < rows; ++i) {
        for (auto k = 0; k < cols; ++k) {
            // inp[i*cols + k] = i*cols + k;
            inp[i * cols + k] = generateRnd();
            out[i * cols + k] = 0.0;
        }
    }

    for (unsigned int directFft = 0; directFft < 1; directFft++) {
        for (unsigned int symPer = 0; symPer < 2; symPer++) {
            for (auto& waveletName : waveletNames) {
                auto obj = Wavelet { waveletName.c_str() };
                for (auto j = 1; j < 3; j++) {
                    // J = 3;

                    auto wt = WaveletTransform2D(obj, "dwt", rows, cols,
                        j); // Initialize the wavelet transform object
                    if (symPer == 0) {
                        setDWT2Extension(wt, "sym");
                    } else {
                        setDWT2Extension(wt, "per");
                    }

                    auto wavecoeffs = dwt(wt, inp.get()); // Perform DWT

                    idwt(wt, wavecoeffs.get(), out.get()); // Perform IDWT (if needed)
                    // Test Reconstruction

                    if (directFft == 0) {
                        epsilon = 1e-4;
                    } else {
                        epsilon = 1e-4;
                    }

                    REQUIRE(rmsError(out.get(), inp.get(), n) <= epsilon);
                }
            }
        }
    }
}

auto main() -> int
{
    fmt::printf("Running DWT ReconstructionTests ... ");
    dwtReconstructionTest();
    fmt::printf("DONE \n");
    fmt::printf("Running DWT2 ReconstructionTests ... ");
    dwT2ReconstructionTest();
    fmt::printf("DONE \n");
    return 0;
}