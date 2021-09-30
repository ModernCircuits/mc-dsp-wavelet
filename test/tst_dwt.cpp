#include "lt/dsp/wavelets.hpp"

#include "lt/algorithm.hpp"
#include "lt/cmath.hpp"
#include "lt/cstdlib.hpp"
#include "lt/cstring.hpp"
#include "lt/format.hpp"
#include "lt/memory.hpp"
#include "lt/random.hpp"
#include "lt/vector.hpp"

#include "lt/testing/test.hpp"

#include "lt/sstream.hpp"

auto modwtReconstructionTest()
{
    auto const n = 4096;
    auto epsilon = 1e-6;

    auto out = std::make_unique<float[]>(n);
    auto inp = std::make_unique<float[]>(n);

    std::random_device rd {};
    auto gen = std::mt19937 { rd() };
    auto dis = std::uniform_real_distribution<float> { 0.0F, 1.0F };
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
                        epsilon = 1e-6;
                    } else {
                        epsilon = 1e-6;
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

    float epsilon = NAN;

    rows = 512;
    cols = 500;

    auto n = rows * cols;

    auto inp = std::make_unique<float[]>(n);
    auto out = std::make_unique<float[]>(n);

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
            out[i * cols + k] = 0.0F;
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

auto dwtReconstructionTest()
{

    int i = 0;
    float epsilon = 1e-5;

    auto n = 8096;

    auto inp = std::make_unique<float[]>(n);
    auto out = std::make_unique<float[]>(n);

    std::random_device rd {};
    auto gen = std::mt19937 { rd() };
    auto dis = std::uniform_real_distribution<float> { 0.0F, 1.0F };

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

                    REQUIRE(rmsError(out.get(), inp.get(), wt.signalLength()) <= epsilon);
                }
            }
        }
    }
}

auto dbCoefTests()
{
    constexpr auto epsilon = 1e-6;
    auto waveletNames = std::vector<std::string>(38);
    std::generate(begin(waveletNames), end(waveletNames), [i = 1]() mutable {
        return std::string("db") + std::to_string(i);
        ++i;
    });

    for (auto const& name : waveletNames) {
        auto obj = Wavelet { name.c_str() };
        auto t1 = sum1(obj.lpr().data(), obj.lpr().size()) - std::sqrt(2.0F);
        auto t2 = sum2(obj.lpr().data(), obj.lpr().size()) - 1.0F / std::sqrt(2.0F);
        auto t3 = sum3(obj.lpr().data(), obj.lpr().size()) - 1.0F / std::sqrt(2.0F);
        auto t4 = sum4(obj.lpr().data(), obj.lpr().size()) - 1.0F;

        REQUIRE(fabs(t1) <= epsilon);
        REQUIRE(fabs(t2) <= epsilon);
        REQUIRE(fabs(t3) <= epsilon);
        REQUIRE(fabs(t4) <= epsilon);

        for (std::size_t m = 1; m < (obj.lpr().size() / 2) - 1; m++) {
            auto t5 = sum5(obj.lpr().data(), obj.lpr().size(), m);
            REQUIRE(fabs(t5) <= epsilon);
        }
    }
}

auto coifCoefTests()
{
    float epsilon = 1e-6;
    float t1 = NAN;
    float t2 = NAN;
    float t3 = NAN;
    float t4 = NAN;
    float t5 = NAN;
    std::vector<std::string> waveletNames;
    waveletNames.resize(17);
    for (std::size_t i = 0; i < waveletNames.size(); i++) {
        waveletNames[i] = std::string("coif") + std::to_string(i + 1);
    }

    for (auto const& name : waveletNames) {
        auto obj = Wavelet { name.c_str() };
        t1 = sum1(obj.lpr().data(), obj.lpr().size()) - std::sqrt(2.0F);
        t2 = sum2(obj.lpr().data(), obj.lpr().size()) - 1.0F / std::sqrt(2.0F);
        t3 = sum3(obj.lpr().data(), obj.lpr().size()) - 1.0F / std::sqrt(2.0F);
        t4 = sum4(obj.lpr().data(), obj.lpr().size()) - 1.0F;

        REQUIRE(std::fabs(t1) <= epsilon);
        REQUIRE(std::fabs(t2) <= epsilon);
        REQUIRE(std::fabs(t3) <= epsilon);
        REQUIRE(std::fabs(t4) <= epsilon);
        for (std::size_t m = 1; m < (obj.lpr().size() / 2) - 1; m++) {
            t5 = sum5(obj.lpr().data(), obj.lpr().size(), m);
            REQUIRE(std::fabs(t5) <= epsilon);
        }
    }
}

auto symCoefTests()
{
    float epsilon = 1e-6;
    float t1 = NAN;
    float t2 = NAN;
    float t3 = NAN;
    float t4 = NAN;
    float t5 = NAN;
    std::vector<std::string> waveletNames;
    for (std::size_t i = 1; i < 20; i++) {
        waveletNames.push_back(std::string("sym") + std::to_string(i + 1));
    }

    for (auto const& name : waveletNames) {
        auto obj = Wavelet { name.c_str() };
        t1 = sum1(obj.lpr().data(), obj.lpr().size()) - std::sqrt(2.0F);
        t2 = sum2(obj.lpr().data(), obj.lpr().size()) - 1.0F / std::sqrt(2.0F);
        t3 = sum3(obj.lpr().data(), obj.lpr().size()) - 1.0F / std::sqrt(2.0F);
        t4 = sum4(obj.lpr().data(), obj.lpr().size()) - 1.0F;

        REQUIRE(std::fabs(t1) <= epsilon);
        REQUIRE(std::fabs(t2) <= epsilon);
        REQUIRE(std::fabs(t3) <= epsilon);
        REQUIRE(std::fabs(t4) <= epsilon);

        for (std::size_t m = 1; m < (obj.lpr().size() / 2) - 1; m++) {
            t5 = sum5(obj.lpr().data(), obj.lpr().size(), m);
            REQUIRE(std::fabs(t5) <= epsilon);
        }
    }
}

auto biorCoefTests()
{
    constexpr float epsilon = 1e-6;
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

        auto const t1 = sum1(obj.lpr().data(), obj.lpr().size()) - std::sqrt(2.0F);
        auto const t2 = sum1(obj.lpd().data(), obj.lpd().size()) - std::sqrt(2.0F);

        auto const t3 = sum2(obj.lpr().data(), obj.lpr().size()) - 1.0F / std::sqrt(2.0F);
        auto const t4 = sum2(obj.lpd().data(), obj.lpd().size()) - 1.0F / std::sqrt(2.0F);

        auto const t5 = sum3(obj.lpr().data(), obj.lpr().size()) - 1.0F / std::sqrt(2.0F);
        auto const t6 = sum3(obj.lpd().data(), obj.lpd().size()) - 1.0F / std::sqrt(2.0F);

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
    float epsilon = 1e-6;
    float t1 = NAN;
    float t2 = NAN;
    float t3 = NAN;
    float t4 = NAN;
    float t5 = NAN;
    float t6 = NAN;
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

        t1 = sum1(obj.lpr().data(), obj.lpr().size()) - std::sqrt(2.0F);
        t2 = sum1(obj.lpd().data(), obj.lpd().size()) - std::sqrt(2.0F);

        t3 = sum2(obj.lpr().data(), obj.lpr().size()) - 1.0F / std::sqrt(2.0F);
        t4 = sum2(obj.lpd().data(), obj.lpd().size()) - 1.0F / std::sqrt(2.0F);

        t5 = sum3(obj.lpr().data(), obj.lpr().size()) - 1.0F / std::sqrt(2.0F);
        t6 = sum3(obj.lpd().data(), obj.lpd().size()) - 1.0F / std::sqrt(2.0F);

        REQUIRE(std::fabs(t1) <= epsilon);
        REQUIRE(std::fabs(t2) <= epsilon);
        REQUIRE(std::fabs(t3) <= epsilon);
        REQUIRE(std::fabs(t4) <= epsilon);
        REQUIRE(std::fabs(t5) <= epsilon);
        REQUIRE(std::fabs(t6) <= epsilon);
    }
}

auto main() -> int
{
    fmt::printf("Running Unit Tests : \n \n");

    fmt::printf("Running DBCoefTests ... ");
    dbCoefTests();
    fmt::printf("DONE \n");

    fmt::printf("Running CoifCoefTests ... ");
    coifCoefTests();
    fmt::printf("DONE \n");

    fmt::printf("Running SymCoefTests ... ");
    symCoefTests();
    fmt::printf("DONE \n");

    fmt::printf("Running BiorCoefTests ... ");
    biorCoefTests();
    fmt::printf("DONE \n");

    fmt::printf("Running RBiorCoefTests ... ");
    rBiorCoefTests();
    fmt::printf("DONE \n");

    fmt::printf("Running MODWT ReconstructionTests ... ");
    modwtReconstructionTest();
    fmt::printf("DONE \n");

    fmt::printf("Running DWPT ReconstructionTests ... ");
    dwtReconstructionTest();
    fmt::printf("DONE \n");

    fmt::printf("Running MODWT2 ReconstructionTests ... ");
    modwT2ReconstructionTest();
    fmt::printf("DONE \n");

    fmt::printf("\n\nUnit Tests Successful\n\n");
    return 0;
}
