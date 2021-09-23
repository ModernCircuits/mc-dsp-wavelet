#include "lt/dsp/wavelets.hpp"

#include "testing.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <random>
#include <vector>

auto swtReconstructionTest()
{

    int i = 0;
    double epsilon = 1e-15;
    double err = NAN;

    auto const n = 4000;

    // N = 256;

    auto inp = makeZeros<double>(n);
    auto out = makeZeros<double>(n);

    auto rd = std::random_device {};
    auto gen = std::mt19937 { rd() };
    auto dis = std::uniform_real_distribution<double> { 0.0, 1.0 };

    for (i = 0; i < n; ++i) {
        inp[i] = dis(gen);
    }
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

    for (unsigned int directFft = 0; directFft < 2; directFft++) {
        for (unsigned int symPer = 0; symPer < 1; symPer++) {
            for (auto& name : waveletNames) {
                auto obj = Wavelet { name.c_str() };
                for (auto j = 1; j < 3; j++) {
                    auto wt = WaveletTransform(obj, "swt", n, j);

                    if (directFft == 0) {
                        wt.convMethod(ConvolutionMethod::direct);
                    } else {
                        wt.convMethod(ConvolutionMethod::fft);
                    }

                    if (symPer == 0) {
                        wt.extension(SignalExtension::periodic); // Options are "per" and "sym".
                        // Symmetric is the default option
                    } else if (symPer == 1 && directFft == 1) {
                        wt.extension(SignalExtension::symmetric);
                    } else {
                        break;
                    }

                    swt(wt, inp.get()); // Perform DWT

                    iswt(wt, out.get()); // Perform IDWT (if needed)
                    // Test Reconstruction

                    if (directFft == 0) {
                        epsilon = 1e-8;
                    } else {
                        epsilon = 1e-10;
                    }
                    // BOOST_CHECK_SMALL(RMS_Error(out, inp, wt.signalLength()), epsilon); //
                    // If Reconstruction succeeded then the output should be a small value.

                    // printf("%g ",RMS_Error(out, inp, wt.signalLength()));
                    err = rmsError(out.get(), inp.get(), wt.signalLength());
                    // printf("%d %d %g \n",direct_fft,sym_per,err);
                    if (err > epsilon) {
                        printf(
                            "\n ERROR : SWT Reconstruction Unit Test Failed. Exiting. \n");
                        exit(-1);
                    }
                }
            }
        }
    }
}

auto swT2ReconstructionTest()
{
    WaveletTransform2D* wt = nullptr;
    // int i;
    // int k;
    // int J;

    auto const rows = 512;
    auto const cols = 500;

    auto const n = rows * cols;

    auto inp = makeZeros<double>(n);
    auto out = makeZeros<double>(n);

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
        for (unsigned int symPer = 0; symPer < 1; symPer++) {
            for (auto& name : waveletNames) {
                auto obj = Wavelet { name.c_str() };
                for (auto j = 1; j < 3; j++) {
                    wt = wt2Init(obj, "swt", rows, cols,
                        j); // Initialize the wavelet transform object
                    if (symPer == 0) {
                        setDWT2Extension(wt, "per"); // Options are "per"
                    }

                    auto wavecoeffs = swt2(wt, inp.get()); // Perform DWT

                    iswt2(wt, wavecoeffs.get(), out.get()); // Perform IDWT (if needed)
                    // Test Reconstruction

                    double epsilon { 0.0 };
                    if (directFft == 0) {
                        epsilon = 1e-8;
                    } else {
                        epsilon = 1e-10;
                    }
                    // BOOST_CHECK_SMALL(RMS_Error(out, inp, wt.signalLength()), epsilon); //
                    // If Reconstruction succeeded then the output should be a small value.

                    // printf("%g ",RMS_Error(out, inp, wt.signalLength()));
                    if (rmsError(out.get(), inp.get(), n) > epsilon) {
                        printf(
                            "\n ERROR : SWT2 Reconstruction Unit Test Failed. Exiting. \n");
                        exit(-1);
                    }
                    wt2Free(wt);
                }
            }
        }
    }
}

auto main() -> int
{
    printf("Running SWT ReconstructionTests ... ");
    swtReconstructionTest();
    printf("DONE \n");
    printf("Running SWT2 ReconstructionTests ... ");
    swT2ReconstructionTest();
    printf("DONE \n");
    return 0;
}