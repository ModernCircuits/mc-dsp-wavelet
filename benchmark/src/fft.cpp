// SPDX-License-Identifier: BSL-1.0

#include <mc/dsp/fft.hpp>

#include <mc/core/complex.hpp>
#include <mc/core/vector.hpp>
#include <mc/testing/test.hpp>

#include <fftw3.h>
#include <pffft.h>

#include <benchmark/benchmark.h>

using namespace mc;

static auto BM_RFFT(benchmark::State& state) -> void
{
    auto size        = static_cast<size_t>(state.range(0));
    auto out         = Vector<Complex<float>>(size);
    auto const input = generateRandomTestData(size);

    auto engine = mc::dsp::RFFT{(int)size};
    while (state.KeepRunning()) {
        rfft(engine, data(input), data(out));
        benchmark::DoNotOptimize(out.front());
        benchmark::DoNotOptimize(out.back());
    }
}

BENCHMARK(BM_RFFT)->Arg(128)->Arg(256)->Arg(512)->Arg(8192 * 4);

static auto BM_PFFFT(benchmark::State& state) -> void
{
    auto size        = static_cast<size_t>(state.range(0));
    auto out         = Vector<Complex<float>>(size);
    auto const input = generateRandomTestData(size);

    auto engine = pffft_new_setup(static_cast<int>(size), PFFFT_REAL);

    while (state.KeepRunning()) {
        pffft_transform_ordered(
            engine,
            data(input),
            (float*)data(out),
            nullptr,
            PFFFT_FORWARD
        );
        benchmark::DoNotOptimize(out.front());
        benchmark::DoNotOptimize(out.back());
    }
}

BENCHMARK(BM_PFFFT)->Arg(128)->Arg(256)->Arg(512)->Arg(8192 * 4);

static auto BM_FFTW(benchmark::State& state) -> void
{
    auto size    = static_cast<size_t>(state.range(0));
    auto in      = Vector<float>(size);
    auto out     = Vector<Complex<float>>(size);
    auto* output = (fftwf_complex*)data(out);
    auto flags   = FFTW_UNALIGNED | FFTW_ESTIMATE;

    auto engine = fftwf_plan_dft_r2c_1d((int)size, data(in), output, flags);

    in  = generateRandomTestData(size);
    out = Vector<Complex<float>>(size);

    while (state.KeepRunning()) {
        fftwf_execute_dft_r2c(engine, data(in), output);
        benchmark::DoNotOptimize(out.front());
        benchmark::DoNotOptimize(out.back());
    }

    fftwf_destroy_plan(engine);
}

BENCHMARK(BM_FFTW)->Arg(128)->Arg(256)->Arg(512)->Arg(8192 * 4);

BENCHMARK_MAIN();
