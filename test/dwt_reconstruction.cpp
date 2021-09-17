#include "wavelib.h"

#include "helper.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <iostream>
#include <sstream>
#include <vector>

void DWTReconstructionTest()
{

    wave_object obj;
    wt_object wt;
    double* inp;
    double* out;
    int N;
    int i;
    int J;
    double epsilon = 1e-15;

    N = 79926;

    // N = 256;

    inp = (double*)malloc(sizeof(double) * N);
    out = (double*)malloc(sizeof(double) * N);
    // wmean = mean(temp, N);

    for (i = 0; i < N; ++i) {
        inp[i] = (rand() / (double)(RAND_MAX));
    }
    std::vector<std::string> waveletNames;

    for (unsigned int j = 0; j < 36; j++) {
        waveletNames.push_back(std::string("db") + patch::to_string(j + 1));
    }
    for (unsigned int j = 0; j < 17; j++) {
        waveletNames.push_back(std::string("coif") + patch::to_string(j + 1));
    }
    for (unsigned int j = 1; j < 20; j++) {
        waveletNames.push_back(std::string("sym") + patch::to_string(j + 1));
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
                obj = wave_init(name); // Initialize the wavelet
                for (J = 1; J < 3; J++) {
                    // J = 3;

                    wt = wt_init(obj, (char*)"dwt", N,
                        J); // Initialize the wavelet transform object
                    if (sym_per == 0) {
                        setDWTExtension(wt,
                            (char*)"sym"); // Options are "per" and "sym".
                        // Symmetric is the default option
                    } else {
                        setDWTExtension(wt, (char*)"per");
                    }
                    if (direct_fft == 0) {
                        setWTConv(wt, (char*)"direct");
                    } else {
                        setWTConv(wt, (char*)"fft");
                    }

                    dwt(wt, inp); // Perform DWT

                    idwt(wt, out); // Perform IDWT (if needed)
                    // Test Reconstruction

                    if (direct_fft == 0) {
                        epsilon = 1e-8;
                    } else {
                        epsilon = 1e-10;
                    }
                    // BOOST_CHECK_SMALL(RMS_Error(out, inp, wt->siglength), epsilon); //
                    // If Reconstruction succeeded then the output should be a small value.

                    // printf("%g ",RMS_Error(out, inp, wt->siglength));
                    if (RMS_Error(out, inp, wt->siglength) > epsilon) {
                        printf(
                            "\n ERROR : DWT Reconstruction Unit Test Failed. Exiting. \n");
                        exit(-1);
                    }
                    wt_free(wt);
                }
                wave_free(obj);
                delete[] name;
            }
        }
    }

    free(out);
    free(inp);
}

void DWT2ReconstructionTest()
{
    wave_object obj;
    wt2_object wt;
    int i;
    int k;
    int J;
    int N;
    int rows;
    int cols;
    double* inp;
    double* wavecoeffs;
    double* out;
    double epsilon;

    rows = 1024;
    cols = 1000;

    N = rows * cols;

    inp = (double*)malloc(sizeof(double) * N);
    out = (double*)malloc(sizeof(double) * N);

    std::vector<std::string> waveletNames;

    for (unsigned int j = 0; j < 15; j++) {
        waveletNames.push_back(std::string("db") + patch::to_string(j + 1));
    }
    for (unsigned int j = 0; j < 5; j++) {
        waveletNames.push_back(std::string("coif") + patch::to_string(j + 1));
    }
    for (unsigned int j = 1; j < 10; j++) {
        waveletNames.push_back(std::string("sym") + patch::to_string(j + 1));
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

    for (i = 0; i < rows; ++i) {
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
                obj = wave_init(name); // Initialize the wavelet
                for (J = 1; J < 3; J++) {
                    // J = 3;

                    wt = wt2_init(obj, (char*)"dwt", rows, cols,
                        J); // Initialize the wavelet transform object
                    if (sym_per == 0) {
                        setDWT2Extension(wt,
                            (char*)"sym"); // Options are "per" and "sym".
                        // Symmetric is the default option
                    } else {
                        setDWT2Extension(wt, (char*)"per");
                    }

                    wavecoeffs = dwt2(wt, inp); // Perform DWT

                    idwt2(wt, wavecoeffs, out); // Perform IDWT (if needed)
                    // Test Reconstruction

                    if (direct_fft == 0) {
                        epsilon = 1e-8;
                    } else {
                        epsilon = 1e-10;
                    }
                    // BOOST_CHECK_SMALL(RMS_Error(out, inp, wt->siglength), epsilon); //
                    // If Reconstruction succeeded then the output should be a small value.

                    // printf("%g ",RMS_Error(out, inp, wt->siglength));
                    if (RMS_Error(out, inp, N) > epsilon) {
                        printf(
                            "\n ERROR : DWT2 Reconstruction Unit Test Failed. Exiting. \n");
                        exit(-1);
                    }
                    wt2_free(wt);
                    free(wavecoeffs);
                }
                wave_free(obj);
                delete[] name;
            }
        }
    }

    free(inp);
    free(out);
}

auto main() -> int
{
    printf("Running DWT ReconstructionTests ... ");
    DWTReconstructionTest();
    printf("DONE \n");
    printf("Running DWT2 ReconstructionTests ... ");
    DWT2ReconstructionTest();
    printf("DONE \n");
    return 0;
}