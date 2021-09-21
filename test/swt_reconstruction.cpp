#include "wavelib.h"

#include "helper.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <iostream>
#include <vector>

void SWTReconstructionTest()
{

    int i;
    double epsilon = 1e-15;
    double err;

    auto const N = 4000;

    // N = 256;

    auto inp = makeZeros<double>(N);
    auto out = makeZeros<double>(N);
    // wmean = mean(temp, N);

    for (i = 0; i < N; ++i) {
        inp[i] = (rand() / (double)(RAND_MAX));
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

    for (unsigned int direct_fft = 0; direct_fft < 2; direct_fft++) {
        for (unsigned int sym_per = 0; sym_per < 1; sym_per++) {
            for (auto& name : waveletNames) {
                auto obj = wavelet { name.c_str() };
                for (auto J = 1; J < 3; J++) {
                    auto wt = wavelet_transform(obj, "swt", N, J);

                    if (direct_fft == 0) {
                        wt.convolution_method("direct");
                    } else {
                        wt.convolution_method("fft");
                    }

                    if (sym_per == 0) {
                        wt.extension(signal_extension::periodic); // Options are "per" and "sym".
                        // Symmetric is the default option
                    } else if (sym_per == 1 && direct_fft == 1) {
                        wt.extension(signal_extension::symmetric);
                    } else {
                        break;
                    }

                    swt(&wt, inp.get()); // Perform DWT

                    iswt(&wt, out.get()); // Perform IDWT (if needed)
                    // Test Reconstruction

                    if (direct_fft == 0) {
                        epsilon = 1e-8;
                    } else {
                        epsilon = 1e-10;
                    }
                    // BOOST_CHECK_SMALL(RMS_Error(out, inp, wt.siglength), epsilon); //
                    // If Reconstruction succeeded then the output should be a small value.

                    // printf("%g ",RMS_Error(out, inp, wt.siglength));
                    err = RMS_Error(out.get(), inp.get(), wt.siglength);
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

void SWT2ReconstructionTest()
{
    wt2_set* wt;
    // int i;
    // int k;
    // int J;

    auto const rows = 512;
    auto const cols = 500;

    auto const N = rows * cols;

    auto inp = makeZeros<double>(N);
    auto out = makeZeros<double>(N);

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
            inp[i * cols + k] = generate_rnd();
            out[i * cols + k] = 0.0;
        }
    }

    for (unsigned int direct_fft = 0; direct_fft < 1; direct_fft++) {
        for (unsigned int sym_per = 0; sym_per < 1; sym_per++) {
            for (auto& name : waveletNames) {
                auto obj = wavelet { name.c_str() };
                for (auto J = 1; J < 3; J++) {
                    wt = wt2_init(obj, "swt", rows, cols,
                        J); // Initialize the wavelet transform object
                    if (sym_per == 0) {
                        setDWT2Extension(wt, "per"); // Options are "per"
                    }

                    auto wavecoeffs = swt2(wt, inp.get()); // Perform DWT

                    iswt2(wt, wavecoeffs.get(), out.get()); // Perform IDWT (if needed)
                    // Test Reconstruction

                    double epsilon { 0.0 };
                    if (direct_fft == 0) {
                        epsilon = 1e-8;
                    } else {
                        epsilon = 1e-10;
                    }
                    // BOOST_CHECK_SMALL(RMS_Error(out, inp, wt.siglength), epsilon); //
                    // If Reconstruction succeeded then the output should be a small value.

                    // printf("%g ",RMS_Error(out, inp, wt.siglength));
                    if (RMS_Error(out.get(), inp.get(), N) > epsilon) {
                        printf(
                            "\n ERROR : SWT2 Reconstruction Unit Test Failed. Exiting. \n");
                        exit(-1);
                    }
                    wt2_free(wt);
                }
            }
        }
    }
}

auto main() -> int
{
    printf("Running SWT ReconstructionTests ... ");
    SWTReconstructionTest();
    printf("DONE \n");
    printf("Running SWT2 ReconstructionTests ... ");
    SWT2ReconstructionTest();
    printf("DONE \n");
    return 0;
}