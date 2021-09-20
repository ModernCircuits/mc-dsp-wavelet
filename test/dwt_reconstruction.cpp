#include "wavelib.h"

#include "helper.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <iostream>
#include <memory>
#include <sstream>
#include <vector>

void DWTReconstructionTest()
{
    wt_set* wt;
    int J;
    auto const N = 79926;

    auto inp = std::make_unique<double[]>(N);
    auto out = std::make_unique<double[]>(N);

    for (auto i = 0; i < N; ++i) {
        inp[i] = (rand() / (double)(RAND_MAX));
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

    for (unsigned int direct_fft = 0; direct_fft < 2; direct_fft++) {
        for (unsigned int sym_per = 0; sym_per < 2; sym_per++) {
            for (auto& waveletName : waveletNames) {
                char* name = new char[waveletName.size() + 1];
                memcpy(name, waveletName.c_str(), waveletName.size() + 1);
                auto* obj = wave_init(name); // Initialize the wavelet
                for (auto J = 1; J < 3; J++) {
                    auto* wt = wt_init(obj, "dwt", N, J);
                    if (sym_per == 0) {
                        setDWTExtension(wt, "sym");
                    } else {
                        setDWTExtension(wt, "per");
                    }
                    if (direct_fft == 0) {
                        setWTConv(wt, "direct");
                    } else {
                        setWTConv(wt, "fft");
                    }

                    dwt(wt, inp.get()); // Perform DWT

                    idwt(wt, out.get()); // Perform IDWT (if needed)
                    // Test Reconstruction

                    auto epsilon = 1e-15;
                    if (direct_fft == 0) {
                        epsilon = 1e-8;
                    } else {
                        epsilon = 1e-10;
                    }

                    // std::printf("%g ",RMS_Error(out.get(), inp.get(), wt->siglength));
                    if (RMS_Error(out.get(), inp.get(), wt->siglength) > epsilon) {
                        std::printf("\n ERROR : DWT Reconstruction Unit Test Failed. Exiting. \n");
                        std::exit(-1);
                    }
                    wt_free(wt);
                }
                wave_free(obj);
                delete[] name;
            }
        }
    }
}

void DWT2ReconstructionTest()
{
    wt2_set* wt;
    int i;
    int k;
    int J;
    int N;
    int rows;
    int cols;
    double* wavecoeffs;
    double epsilon;

    rows = 1024;
    cols = 1000;

    N = rows * cols;

    auto inp = std::make_unique<double[]>(N);
    auto out = std::make_unique<double[]>(N);

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
        for (k = 0; k < cols; ++k) {
            // inp[i*cols + k] = i*cols + k;
            inp[i * cols + k] = generate_rnd();
            out[i * cols + k] = 0.0;
        }
    }

    for (unsigned int direct_fft = 0; direct_fft < 1; direct_fft++) {
        for (unsigned int sym_per = 0; sym_per < 2; sym_per++) {
            for (auto& waveletName : waveletNames) {
                char* name = new char[waveletName.size() + 1];
                memcpy(name, waveletName.c_str(), waveletName.size() + 1);
                auto* obj = wave_init(name); // Initialize the wavelet
                for (J = 1; J < 3; J++) {
                    // J = 3;

                    wt = wt2_init(obj, "dwt", rows, cols,
                        J); // Initialize the wavelet transform object
                    if (sym_per == 0) {
                        setDWT2Extension(wt,
                            "sym"); // Options are "per" and "sym".
                        // Symmetric is the default option
                    } else {
                        setDWT2Extension(wt, "per");
                    }

                    wavecoeffs = dwt2(wt, inp.get()); // Perform DWT

                    idwt2(wt, wavecoeffs, out.get()); // Perform IDWT (if needed)
                    // Test Reconstruction

                    if (direct_fft == 0) {
                        epsilon = 1e-8;
                    } else {
                        epsilon = 1e-10;
                    }

                    if (RMS_Error(out.get(), inp.get(), N) > epsilon) {
                        std::printf("\n ERROR : DWT2 Reconstruction Unit Test Failed. Exiting. \n");
                        std::exit(-1);
                    }
                    wt2_free(wt);
                    free(wavecoeffs);
                }
                wave_free(obj);
                delete[] name;
            }
        }
    }
}

auto main() -> int
{
    std::printf("Running DWT ReconstructionTests ... ");
    DWTReconstructionTest();
    std::printf("DONE \n");
    std::printf("Running DWT2 ReconstructionTests ... ");
    DWT2ReconstructionTest();
    std::printf("DONE \n");
    return 0;
}