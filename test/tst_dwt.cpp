#include "lt/cmath.hpp"
#include "lt/dsp/wavelets.hpp"

#include "lt/testing/test.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include <vector>

auto modwtReconstructionTest()
{
    auto const n = 4096;
    auto epsilon = 1e-15;

    auto out = std::make_unique<double[]>(n);
    auto inp = std::make_unique<double[]>(n);

    auto rd = std::random_device {};
    auto gen = std::mt19937 { rd() };
    auto dis = std::uniform_real_distribution<double> { 0.0, 1.0 };
    std::generate_n(inp.get(), n, [&] { return dis(gen); });

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

    for (std::size_t directFft = 0; directFft < 2; directFft++) {
        for (std::size_t symPer = 0; symPer < 1; symPer++) {
            for (auto const& name : waveletNames) {
                auto obj = Wavelet { name.c_str() };
                for (auto j = 1; j < 3; j++) {
                    auto wt = WaveletTransform(obj, "modwt", n, j);

                    if (directFft == 0) {
                        wt.convMethod(ConvolutionMethod::direct);
                    } else {
                        wt.convMethod(ConvolutionMethod::fft);
                    }

                    if (symPer == 0) {
                        wt.extension(SignalExtension::periodic);
                    } else if (symPer == 1 && directFft == 1) {
                        wt.extension(SignalExtension::symmetric);
                    } else {
                        break;
                    }

                    modwt(wt, inp.get());
                    imodwt(wt, out.get());

                    if (directFft == 0) {
                        epsilon = 1e-8;
                    } else {
                        epsilon = 1e-10;
                    }

                    auto const err = rmsError(out.get(), inp.get(), wt.signalLength());
                    REQUIRE(err <= epsilon);
                }
            }
        }
    }
}

auto modwT2ReconstructionTest()
{
    int i = 0;
    int k = 0;
    int rows = 0;
    int cols = 0;

    double epsilon = NAN;

    rows = 512;
    cols = 500;

    auto n = rows * cols;

    auto inp = std::make_unique<double[]>(n);
    auto out = std::make_unique<double[]>(n);

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
            inp[i * cols + k] = generateRnd();
            out[i * cols + k] = 0.0;
        }
    }

    for (std::size_t directFft = 0; directFft < 1; directFft++) {
        for (std::size_t symPer = 0; symPer < 1; symPer++) {
            for (auto const& name : waveletNames) {
                auto obj = Wavelet { name.c_str() };
                for (auto j = 1; j < 3; j++) {
                    auto wt = WaveletTransform2D(obj, "modwt", rows, cols, j);
                    if (symPer == 0) {
                        setDWT2Extension(wt, "per");
                    }

                    auto wavecoeffs = modwt(wt, inp.get());
                    imodwt(wt, wavecoeffs.get(), out.get());

                    if (directFft == 0) {
                        epsilon = 1e-8;
                    } else {
                        epsilon = 1e-10;
                    }
                    REQUIRE(rmsError(out.get(), inp.get(), n) <= epsilon);
                }
            }
        }
    }
}

auto dwtReconstructionTest()
{

    int i = 0;
    double epsilon = 1e-8;

    auto n = 8096;

    auto inp = std::make_unique<double[]>(n);
    auto out = std::make_unique<double[]>(n);

    auto rd = std::random_device {};
    auto gen = std::mt19937 { rd() };
    auto dis = std::uniform_real_distribution<double> { 0.0, 1.0 };

    for (i = 0; i < n; ++i) {
        inp[i] = dis(gen);
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
        for (std::size_t symPer = 0; symPer < 2; symPer++) {
            for (auto const& name : waveletNames) {
                auto obj = Wavelet { name.c_str() };
                for (auto j = 1; j < 3; j++) {
                    // J = 3;

                    auto wt = WaveletPacketTransform(&obj, n, j); // Initialize the wavelet transform object
                    if (symPer == 0) {
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

                    dwt(wt, inp.get());
                    idwt(wt, out.get());

                    REQUIRE(rmsError(out.get(), inp.get(), wt.siglength) <= epsilon);
                }
            }
        }
    }
}

auto cwtReconstructionTest()
{
    int i = 0;
    int n = 0;
    int j = 0;
    int subscale = 0;
    int a0 = 0;
    double dt = NAN;
    double dj = NAN;
    double s0 = NAN;
    double pi = NAN;
    double t = NAN;
    double epsilon = NAN;
    int it1 = 0;
    int it2 = 0;

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
    n = 2048;
    dt = 0.000125;
    subscale = 20;
    dj = 1.0 / (double)subscale;
    s0 = dt / 32;
    j = 32 * subscale;
    a0 = 2; // power

    auto inp = std::make_unique<double[]>(n);
    auto oup = std::make_unique<double[]>(n);

    pi = 4.0 * atan(1.0);

    for (i = 0; i < n; ++i) {
        t = dt * i;
        inp[i] = sin(2 * pi * 500 * t) + sin(2 * pi * 1000 * t) + 0.1 * sin(2 * pi * 8 * t);
        if (i == 1200 || i == 1232) {
            inp[i] += 5.0;
        }
    }

    for (it1 = 0; it1 < 3; ++it1) {
        for (it2 = 0; it2 < 10; ++it2) {
            auto wt = ComplexWaveletTransform { wave[it1], param[it1 * 10 + it2], n, dt, j };
            wt.scales(s0, dj, type, a0);
            cwt(wt, inp.get());
            icwt(wt, oup.get());

            REQUIRE(relError(inp.get(), oup.get(), wt.signalLength) <= epsilon);
        }
    }
}

auto dbCoefTests()
{
    constexpr auto epsilon = 1e-15;
    auto waveletNames = std::vector<std::string>(38);
    std::generate(begin(waveletNames), end(waveletNames), [i = 1]() mutable {
        return std::string("db") + std::to_string(i);
        ++i;
    });

    for (auto const& name : waveletNames) {
        auto obj = Wavelet { name.c_str() };
        auto t1 = sum1(obj.lpr(), obj.lprLen()) - std::sqrt(2.0);
        auto t2 = sum2(obj.lpr(), obj.lprLen()) - 1.0 / std::sqrt(2.0);
        auto t3 = sum3(obj.lpr(), obj.lprLen()) - 1.0 / std::sqrt(2.0);
        auto t4 = sum4(obj.lpr(), obj.lprLen()) - 1.0;

        REQUIRE(fabs(t1) <= epsilon);
        REQUIRE(fabs(t2) <= epsilon);
        REQUIRE(fabs(t3) <= epsilon);
        REQUIRE(fabs(t4) <= epsilon);

        for (int m = 1; m < (obj.lprLen() / 2) - 1; m++) {
            auto t5 = sum5(obj.lpr(), obj.lprLen(), m);
            REQUIRE(fabs(t5) <= epsilon);
        }
    }
}

auto coifCoefTests()
{
    double epsilon = 1e-15;
    double t1 = NAN;
    double t2 = NAN;
    double t3 = NAN;
    double t4 = NAN;
    double t5 = NAN;
    std::vector<std::string> waveletNames;
    waveletNames.resize(17);
    for (std::size_t i = 0; i < waveletNames.size(); i++) {
        waveletNames[i] = std::string("coif") + std::to_string(i + 1);
    }

    for (auto const& name : waveletNames) {
        auto obj = Wavelet { name.c_str() };
        t1 = sum1(obj.lpr(), obj.lprLen()) - std::sqrt(2.0);
        t2 = sum2(obj.lpr(), obj.lprLen()) - 1.0 / std::sqrt(2.0);
        t3 = sum3(obj.lpr(), obj.lprLen()) - 1.0 / std::sqrt(2.0);
        t4 = sum4(obj.lpr(), obj.lprLen()) - 1.0;

        REQUIRE(fabs(t1) <= epsilon);
        REQUIRE(fabs(t2) <= epsilon);
        REQUIRE(fabs(t3) <= epsilon);
        REQUIRE(fabs(t4) <= epsilon);
        for (int m = 1; m < (obj.lprLen() / 2) - 1; m++) {
            t5 = sum5(obj.lpr(), obj.lprLen(), m);
            REQUIRE(fabs(t5) <= epsilon);
        }
    }
}

auto symCoefTests()
{
    double epsilon = 1e-10;
    double t1 = NAN;
    double t2 = NAN;
    double t3 = NAN;
    double t4 = NAN;
    double t5 = NAN;
    std::vector<std::string> waveletNames;
    for (std::size_t i = 1; i < 20; i++) {
        waveletNames.push_back(std::string("sym") + std::to_string(i + 1));
    }

    for (auto const& name : waveletNames) {
        auto obj = Wavelet { name.c_str() };
        t1 = sum1(obj.lpr(), obj.lprLen()) - std::sqrt(2.0);
        t2 = sum2(obj.lpr(), obj.lprLen()) - 1.0 / std::sqrt(2.0);
        t3 = sum3(obj.lpr(), obj.lprLen()) - 1.0 / std::sqrt(2.0);
        t4 = sum4(obj.lpr(), obj.lprLen()) - 1.0;

        REQUIRE(fabs(t1) <= epsilon);
        REQUIRE(fabs(t2) <= epsilon);
        REQUIRE(fabs(t3) <= epsilon);
        REQUIRE(fabs(t4) <= epsilon);

        for (int m = 1; m < (obj.lprLen() / 2) - 1; m++) {
            t5 = sum5(obj.lpr(), obj.lprLen(), m);
            REQUIRE(fabs(t5) <= epsilon);
        }
    }
}

auto biorCoefTests()
{
    constexpr double epsilon = 1e-10;
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
        auto obj = Wavelet { name.c_str() };

        auto const t1 = sum1(obj.lpr(), obj.lprLen()) - std::sqrt(2.0);
        auto const t2 = sum1(obj.lpd(), obj.lpdLen()) - std::sqrt(2.0);

        auto const t3 = sum2(obj.lpr(), obj.lprLen()) - 1.0 / std::sqrt(2.0);
        auto const t4 = sum2(obj.lpd(), obj.lpdLen()) - 1.0 / std::sqrt(2.0);

        auto const t5 = sum3(obj.lpr(), obj.lprLen()) - 1.0 / std::sqrt(2.0);
        auto const t6 = sum3(obj.lpd(), obj.lpdLen()) - 1.0 / std::sqrt(2.0);

        REQUIRE(fabs(t1) <= epsilon);
        REQUIRE(fabs(t2) <= epsilon);
        REQUIRE(fabs(t3) <= epsilon);
        REQUIRE(fabs(t4) <= epsilon);
        REQUIRE(fabs(t5) <= epsilon);
        REQUIRE(fabs(t6) <= epsilon);
    }
}

auto rBiorCoefTests()
{
    double epsilon = 1e-10;
    double t1 = NAN;
    double t2 = NAN;
    double t3 = NAN;
    double t4 = NAN;
    double t5 = NAN;
    double t6 = NAN;
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
        auto obj = Wavelet { name.c_str() };

        t1 = sum1(obj.lpr(), obj.lprLen()) - std::sqrt(2.0);
        t2 = sum1(obj.lpd(), obj.lpdLen()) - std::sqrt(2.0);

        t3 = sum2(obj.lpr(), obj.lprLen()) - 1.0 / std::sqrt(2.0);
        t4 = sum2(obj.lpd(), obj.lpdLen()) - 1.0 / std::sqrt(2.0);

        t5 = sum3(obj.lpr(), obj.lprLen()) - 1.0 / std::sqrt(2.0);
        t6 = sum3(obj.lpd(), obj.lpdLen()) - 1.0 / std::sqrt(2.0);

        REQUIRE(fabs(t1) <= epsilon);
        REQUIRE(fabs(t2) <= epsilon);
        REQUIRE(fabs(t3) <= epsilon);
        REQUIRE(fabs(t4) <= epsilon);
        REQUIRE(fabs(t5) <= epsilon);
        REQUIRE(fabs(t6) <= epsilon);
    }
}

auto main() -> int
{
    printf("Running Unit Tests : \n \n");
    printf("Running DBCoefTests ... ");
    dbCoefTests();
    printf("DONE \n");
    printf("Running CoifCoefTests ... ");
    coifCoefTests();
    printf("DONE \n");
    printf("Running SymCoefTests ... ");
    symCoefTests();
    printf("DONE \n");
    printf("Running BiorCoefTests ... ");
    biorCoefTests();
    printf("DONE \n");
    printf("Running RBiorCoefTests ... ");
    rBiorCoefTests();
    printf("DONE \n");
    printf("Running MODWT ReconstructionTests ... ");
    modwtReconstructionTest();
    printf("DONE \n");
    printf("Running DWPT ReconstructionTests ... ");
    dwtReconstructionTest();
    printf("DONE \n");
    printf("Running CWT ReconstructionTests ... ");
    cwtReconstructionTest();
    printf("DONE \n");
    printf("Running MODWT2 ReconstructionTests ... ");
    modwT2ReconstructionTest();
    printf("DONE \n");
    printf("\n\nUnit Tests Successful\n\n");
    return 0;
}
