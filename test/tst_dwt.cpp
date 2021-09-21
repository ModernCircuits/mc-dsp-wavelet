/*
 * Copyright (c) 2016 Holger Nahrstaedt (TU Berlin)
 */
#include "wavelib.h"

#include "helper.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <iostream>
#include <memory>
#include <sstream>
#include <vector>

void MODWTReconstructionTest()
{
    auto const N = 4096;
    auto epsilon = 1e-15;

    auto out = std::make_unique<double[]>(N);
    auto inp = std::make_unique<double[]>(N);
    std::generate_n(inp.get(), N, [] { return (rand() / (double)(RAND_MAX)); });

    auto waveletNames = std::vector<std::string> {};

    for (std::size_t j = 0; j < 15; j++) {
        waveletNames.push_back(std::string("db") + std::to_string(j + 1));
    }
    for (std::size_t j = 0; j < 5; j++) {
        waveletNames.push_back(std::string("coif") + std::to_string(j + 1));
    }
    for (std::size_t j = 1; j < 10; j++) {
        waveletNames.push_back(std::string("sym") + std::to_string(j + 1));
    }

    for (std::size_t direct_fft = 0; direct_fft < 2; direct_fft++) {
        for (std::size_t sym_per = 0; sym_per < 1; sym_per++) {
            for (auto const& name : waveletNames) {
                auto obj = wavelet { name.c_str() };
                for (auto J = 1; J < 3; J++) {
                    auto wt = wavelet_transform(obj, "modwt", N, J);

                    if (direct_fft == 0) {
                        wt.convolution_method("direct");
                    } else {
                        wt.convolution_method("fft");
                    }

                    if (sym_per == 0) {
                        wt.dwt_extension("per");
                    } else if (sym_per == 1 && direct_fft == 1) {
                        wt.dwt_extension("sym");
                    } else {
                        break;
                    }

                    modwt(&wt, inp.get());
                    imodwt(&wt, out.get());

                    if (direct_fft == 0) {
                        epsilon = 1e-8;
                    } else {
                        epsilon = 1e-10;
                    }

                    auto const err = RMS_Error(out.get(), inp.get(), wt.siglength);
                    if (err > epsilon) {
                        printf("\n ERROR : DWT Reconstruction Unit Test Failed. Exiting. \n");
                        exit(-1);
                    }
                }
            }
        }
    }
}

void MODWT2ReconstructionTest()
{
    wt2_set* wt;
    int i;
    int k;
    int J;
    int N;
    int rows;
    int cols;

    double epsilon;

    rows = 512;
    cols = 500;

    N = rows * cols;

    auto inp = std::make_unique<double[]>(N);
    auto out = std::make_unique<double[]>(N);

    std::vector<std::string> waveletNames;

    for (std::size_t j = 0; j < 15; j++) {
        waveletNames.push_back(std::string("db") + std::to_string(j + 1));
    }
    for (std::size_t j = 0; j < 5; j++) {
        waveletNames.push_back(std::string("coif") + std::to_string(j + 1));
    }
    for (std::size_t j = 1; j < 10; j++) {
        waveletNames.push_back(std::string("sym") + std::to_string(j + 1));
    }

    for (i = 0; i < rows; ++i) {
        for (k = 0; k < cols; ++k) {
            // inp[i*cols + k] = i*cols + k;
            inp[i * cols + k] = generate_rnd();
            out[i * cols + k] = 0.0;
        }
    }

    for (std::size_t direct_fft = 0; direct_fft < 1; direct_fft++) {
        for (std::size_t sym_per = 0; sym_per < 1; sym_per++) {
            for (auto const& name : waveletNames) {
                auto obj = wavelet { name.c_str() };
                for (J = 1; J < 3; J++) {
                    wt = wt2_init(obj, "modwt", rows, cols, J);
                    if (sym_per == 0) {
                        setDWT2Extension(wt, "per");
                    }

                    auto wavecoeffs = modwt2(wt, inp.get());
                    imodwt2(wt, wavecoeffs.get(), out.get());

                    if (direct_fft == 0) {
                        epsilon = 1e-8;
                    } else {
                        epsilon = 1e-10;
                    }
                    if (RMS_Error(out.get(), inp.get(), N) > epsilon) {
                        printf("\n ERROR : MODWT2 Reconstruction Unit Test Failed. "
                               "Exiting. \n");
                        exit(-1);
                    }
                    wt2_free(wt);
                }
            }
        }
    }
}

void DWPTReconstructionTest()
{

    wpt_set* wt;

    int N;
    int i;
    int J;
    double epsilon = 1e-8;

    N = 8096;

    // N = 256;

    auto inp = std::make_unique<double[]>(N);
    auto out = std::make_unique<double[]>(N);
    // wmean = mean(temp, N);

    for (i = 0; i < N; ++i) {
        inp[i] = (rand() / (double)(RAND_MAX));
    }
    std::vector<std::string> waveletNames;

    for (std::size_t j = 0; j < 36; j++) {
        waveletNames.push_back(std::string("db") + std::to_string(j + 1));
    }
    for (std::size_t j = 0; j < 17; j++) {
        waveletNames.push_back(std::string("coif") + std::to_string(j + 1));
    }
    for (std::size_t j = 1; j < 20; j++) {
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

    for (std::size_t ent = 0; ent < 2; ent++) {
        for (std::size_t sym_per = 0; sym_per < 2; sym_per++) {
            for (auto const& name : waveletNames) {
                auto obj = wavelet { name.c_str() };
                for (J = 1; J < 3; J++) {
                    // J = 3;

                    wt = wpt_init(&obj, N, J); // Initialize the wavelet transform object
                    if (sym_per == 0) {
                        setDWPTExtension(wt,
                            "sym"); // Options are "per" and "sym".
                        // Symmetric is the default option
                    } else {
                        setDWPTExtension(wt, "per");
                    }

                    if (ent == 0) {
                        setDWPTEntropy(wt, "shannon", 0);
                    } else {
                        setDWPTEntropy(wt, "logenergy", 0);
                    }

                    dwpt(wt, inp.get()); // Perform DWT

                    idwpt(wt, out.get()); // Perform IDWT (if needed)
                        // Test Reconstruction

                    // BOOST_CHECK_SMALL(RMS_Error(out.get(), inp.get(), wt->siglength), epsilon); //
                    // If Reconstruction succeeded then the output should be a small value.

                    // printf("%s %g \n",name,RMS_Error(out.get(), inp.get(), wt->siglength));

                    if (RMS_Error(out.get(), inp.get(), wt->siglength) > epsilon) {
                        printf(
                            "\n ERROR : DWPT Reconstruction Unit Test Failed. Exiting. \n");
                        exit(-1);
                    }
                    wpt_free(wt);
                }
            }
        }
    }
}

void CWTReconstructionTest()
{
    int i;
    int N;
    int J;
    int subscale;
    int a0;
    double dt;
    double dj;
    double s0;
    double pi;
    double t;
    double epsilon;
    int it1;
    int it2;
    cwavelet_transform* wt;

    char const* wave[3] {
        "morl",
        "paul",
        "dog",
    };
    double param[30] = {
        4.5,
        5,
        5.5,
        6,
        6.5,
        8,
        10,
        13,
        17,
        20,
        4,
        5,
        7,
        8,
        10,
        12,
        13,
        14,
        17,
        20,
        2,
        4,
        6,
        8,
        10,
        12,
        14,
        16,
        18,
        20,
    };
    char const* type = "pow";

    epsilon = 0.01;
    N = 2048;
    dt = 0.000125;
    subscale = 20;
    dj = 1.0 / (double)subscale;
    s0 = dt / 32;
    J = 32 * subscale;
    a0 = 2; // power

    auto inp = std::make_unique<double[]>(N);
    auto oup = std::make_unique<double[]>(N);

    pi = 4.0 * atan(1.0);

    for (i = 0; i < N; ++i) {
        t = dt * i;
        inp[i] = sin(2 * pi * 500 * t) + sin(2 * pi * 1000 * t) + 0.1 * sin(2 * pi * 8 * t);
        if (i == 1200 || i == 1232) {
            inp[i] += 5.0;
        }
    }

    for (it1 = 0; it1 < 3; ++it1) {
        for (it2 = 0; it2 < 10; ++it2) {

            wt = cwt_init(wave[it1], param[it1 * 10 + it2], N, dt, J);

            setCWTScales(wt, s0, dj, type, a0);

            cwt(wt, inp.get());

            icwt(wt, oup.get());

            // printf("\nWavelet : %s Parameter %g Error %g \n",
            // wave[it1],param[it1*10+it2],REL_Error(inp.get(),oup.get(), wt->siglength));
            if (REL_Error(inp.get(), oup.get(), wt->siglength) > epsilon) {
                printf("\n ERROR : DWPT Reconstruction Unit Test Failed. Exiting. \n");
                exit(-1);
            }

            cwt_free(wt);
        }
    }
}

void DBCoefTests()
{
    constexpr auto epsilon = 1e-15;
    auto waveletNames = std::vector<std::string>(38);
    std::generate(begin(waveletNames), end(waveletNames), [i = 1]() mutable {
        return std::string("db") + std::to_string(i);
        ++i;
    });

    for (auto const& name : waveletNames) {
        auto obj = wavelet { name.c_str() };
        auto t1 = sum1(obj.lpr(), obj.lpr_len()) - std::sqrt(2.0);
        auto t2 = sum2(obj.lpr(), obj.lpr_len()) - 1.0 / std::sqrt(2.0);
        auto t3 = sum3(obj.lpr(), obj.lpr_len()) - 1.0 / std::sqrt(2.0);
        auto t4 = sum4(obj.lpr(), obj.lpr_len()) - 1.0;

        if (fabs(t1) > epsilon || fabs(t2) > epsilon || fabs(t3) > epsilon || fabs(t4) > epsilon) {
            printf("\n ERROR : DB Coefficients Unit Test Failed. Exiting. \n");
            exit(-1);
        }

        for (int m = 1; m < (obj.lpr_len() / 2) - 1; m++) {
            auto t5 = sum5(obj.lpr(), obj.lpr_len(), m);
            if (fabs(t5) > epsilon) {
                printf("\n ERROR : DB Coefficients Unit Test Failed. Exiting. \n");
                exit(-1);
            }
        }
    }
}

void CoifCoefTests()
{
    double epsilon = 1e-15;
    double t1;
    double t2;
    double t3;
    double t4;
    double t5;
    std::vector<std::string> waveletNames;
    waveletNames.resize(17);
    for (std::size_t i = 0; i < waveletNames.size(); i++) {
        waveletNames[i] = std::string("coif") + std::to_string(i + 1);
    }

    for (auto const& name : waveletNames) {
        auto obj = wavelet { name.c_str() };
        t1 = sum1(obj.lpr(), obj.lpr_len()) - std::sqrt(2.0);
        t2 = sum2(obj.lpr(), obj.lpr_len()) - 1.0 / std::sqrt(2.0);
        t3 = sum3(obj.lpr(), obj.lpr_len()) - 1.0 / std::sqrt(2.0);
        t4 = sum4(obj.lpr(), obj.lpr_len()) - 1.0;

        if (fabs(t1) > epsilon || fabs(t2) > epsilon || fabs(t3) > epsilon || fabs(t4) > epsilon) {
            printf("\n ERROR : Coif Coefficients Unit Test Failed. Exiting. \n");
            exit(-1);
        }

        for (int m = 1; m < (obj.lpr_len() / 2) - 1; m++) {
            t5 = sum5(obj.lpr(), obj.lpr_len(), m);
            if (fabs(t5) > epsilon) {
                printf("\n ERROR : Coif Coefficients Unit Test Failed. Exiting. \n");
                exit(-1);
            }
        }
    }
}

void SymCoefTests()
{
    double epsilon = 1e-10;
    double t1;
    double t2;
    double t3;
    double t4;
    double t5;
    std::vector<std::string> waveletNames;
    for (std::size_t i = 1; i < 20; i++) {
        waveletNames.push_back(std::string("sym") + std::to_string(i + 1));
    }

    for (auto const& name : waveletNames) {
        auto obj = wavelet { name.c_str() };
        t1 = sum1(obj.lpr(), obj.lpr_len()) - std::sqrt(2.0);
        t2 = sum2(obj.lpr(), obj.lpr_len()) - 1.0 / std::sqrt(2.0);
        t3 = sum3(obj.lpr(), obj.lpr_len()) - 1.0 / std::sqrt(2.0);
        t4 = sum4(obj.lpr(), obj.lpr_len()) - 1.0;

        if (fabs(t1) > epsilon || fabs(t2) > epsilon || fabs(t3) > epsilon || fabs(t4) > epsilon) {
            printf("\n ERROR : Sym Coefficients Unit Test Failed. Exiting. \n");
            exit(-1);
        }

        for (int m = 1; m < (obj.lpr_len() / 2) - 1; m++) {
            t5 = sum5(obj.lpr(), obj.lpr_len(), m);
            if (fabs(t5) > epsilon) {
                printf("\n ERROR : Sym Coefficients Unit Test Failed. Exiting. \n");
                exit(-1);
            }
        }
    }
}

void BiorCoefTests()
{
    double epsilon = 1e-10;
    double t1;
    double t2;
    double t3;
    double t4;
    double t5;
    double t6;
    std::vector<std::string> waveletNames;
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

    for (auto const& name : waveletNames) {
        auto obj = wavelet { name.c_str() };

        t1 = sum1(obj.lpr(), obj.lpr_len()) - std::sqrt(2.0);
        t2 = sum1(obj.lpd(), obj.lpd_len()) - std::sqrt(2.0);

        t3 = sum2(obj.lpr(), obj.lpr_len()) - 1.0 / std::sqrt(2.0);
        t4 = sum2(obj.lpd(), obj.lpd_len()) - 1.0 / std::sqrt(2.0);

        t5 = sum3(obj.lpr(), obj.lpr_len()) - 1.0 / std::sqrt(2.0);
        t6 = sum3(obj.lpd(), obj.lpd_len()) - 1.0 / std::sqrt(2.0);

        if (fabs(t1) > epsilon || fabs(t2) > epsilon || fabs(t3) > epsilon || fabs(t4) > epsilon || fabs(t5) > epsilon || fabs(t6) > epsilon) {
            printf("\n ERROR : Bior Coefficients Unit Test Failed. Exiting. \n");
            exit(-1);
        }
    }
}

void RBiorCoefTests()
{
    double epsilon = 1e-10;
    double t1;
    double t2;
    double t3;
    double t4;
    double t5;
    double t6;
    std::vector<std::string> waveletNames;
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

    for (auto const& name : waveletNames) {
        auto obj = wavelet { name.c_str() };

        t1 = sum1(obj.lpr(), obj.lpr_len()) - std::sqrt(2.0);
        t2 = sum1(obj.lpd(), obj.lpd_len()) - std::sqrt(2.0);

        t3 = sum2(obj.lpr(), obj.lpr_len()) - 1.0 / std::sqrt(2.0);
        t4 = sum2(obj.lpd(), obj.lpd_len()) - 1.0 / std::sqrt(2.0);

        t5 = sum3(obj.lpr(), obj.lpr_len()) - 1.0 / std::sqrt(2.0);
        t6 = sum3(obj.lpd(), obj.lpd_len()) - 1.0 / std::sqrt(2.0);

        if (fabs(t1) > epsilon || fabs(t2) > epsilon || fabs(t3) > epsilon || fabs(t4) > epsilon || fabs(t5) > epsilon || fabs(t6) > epsilon) {
            printf("\n ERROR : RBior Coefficients Unit Test Failed. Exiting. \n");
            exit(-1);
        }
    }
}

auto main() -> int
{
    printf("Running Unit Tests : \n \n");
    printf("Running DBCoefTests ... ");
    DBCoefTests();
    printf("DONE \n");
    printf("Running CoifCoefTests ... ");
    CoifCoefTests();
    printf("DONE \n");
    printf("Running SymCoefTests ... ");
    SymCoefTests();
    printf("DONE \n");
    printf("Running BiorCoefTests ... ");
    BiorCoefTests();
    printf("DONE \n");
    printf("Running RBiorCoefTests ... ");
    RBiorCoefTests();
    printf("DONE \n");
    printf("Running MODWT ReconstructionTests ... ");
    MODWTReconstructionTest();
    printf("DONE \n");
    printf("Running DWPT ReconstructionTests ... ");
    DWPTReconstructionTest();
    printf("DONE \n");
    printf("Running CWT ReconstructionTests ... ");
    CWTReconstructionTest();
    printf("DONE \n");
    printf("Running MODWT2 ReconstructionTests ... ");
    MODWT2ReconstructionTest();
    printf("DONE \n");
    printf("\n\nUnit Tests Successful\n\n");
    return 0;
}
