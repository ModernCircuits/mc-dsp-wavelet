#include "lt/dsp/wavelets.hpp"

#include "testing.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include <vector>

auto dwtReconstructionTest()
{
    auto const n = 22'050;

    auto inp = std::make_unique<double[]>(n);
    auto out = std::make_unique<double[]>(n);

    auto rd = std::random_device {};
    auto gen = std::mt19937 { rd() };
    auto dis = std::uniform_real_distribution<double> { 0.0, 1.0 };

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

                    auto epsilon = 1e-15;
                    if (directFft == 0) {
                        epsilon = 1e-8;
                    } else {
                        epsilon = 1e-10;
                    }

                    // std::printf("%g ",RMS_Error(out.get(), inp.get(), wt.siglength));
                    if (rmsError(out.get(), inp.get(), wt.siglength) > epsilon) {
                        std::printf("\n ERROR : DWT Reconstruction Unit Test Failed. Exiting. \n");
                        std::exit(-1);
                    }
                }
            }
        }
    }
}

auto dwT2ReconstructionTest()
{
    WaveletTransform2D* wt;
    double epsilon;

    auto const rows = 256;
    auto const cols = 200;

    auto const n = rows * cols;

    auto inp = std::make_unique<double[]>(n);
    auto out = std::make_unique<double[]>(n);

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

                    wt = wt2Init(obj, "dwt", rows, cols,
                        j); // Initialize the wavelet transform object
                    if (symPer == 0) {
                        setDWT2Extension(wt,
                            "sym"); // Options are "per" and "sym".
                        // Symmetric is the default option
                    } else {
                        setDWT2Extension(wt, "per");
                    }

                    auto wavecoeffs = dwt(wt, inp.get()); // Perform DWT

                    idwt(wt, wavecoeffs.get(), out.get()); // Perform IDWT (if needed)
                    // Test Reconstruction

                    if (directFft == 0) {
                        epsilon = 1e-8;
                    } else {
                        epsilon = 1e-10;
                    }

                    if (rmsError(out.get(), inp.get(), n) > epsilon) {
                        std::printf("\n ERROR : DWT2 Reconstruction Unit Test Failed. Exiting. \n");
                        std::exit(-1);
                    }
                    wt2Free(wt);
                }
            }
        }
    }
}

auto main() -> int
{
    std::printf("Running DWT ReconstructionTests ... ");
    dwtReconstructionTest();
    std::printf("DONE \n");
    std::printf("Running DWT2 ReconstructionTests ... ");
    dwT2ReconstructionTest();
    std::printf("DONE \n");
    return 0;
}