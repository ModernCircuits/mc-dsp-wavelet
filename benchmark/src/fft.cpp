// SPDX-License-Identifier: BSL-1.0

#include <mc/dsp/fft.hpp>

#include <mc/core/complex.hpp>
#include <mc/core/vector.hpp>
#include <mc/testing/test.hpp>

#include <benchmark/benchmark.h>

using namespace mc;

void BM_RFFT(benchmark::State& state)
{
    auto size        = static_cast<size_t>(state.range(0));
    auto fft         = mc::dsp::RFFT{(int)size, mc::dsp::FFTDirection::forward};
    auto const input = generateRandomTestData(size);

    auto output = Vector<Complex<float>>(size);
    while (state.KeepRunning()) {
        fft.performRealToComplex(data(input), data(output));
        benchmark::DoNotOptimize(output.front());
        benchmark::DoNotOptimize(output.back());
    }
}

BENCHMARK(BM_RFFT)->Arg(128)->Arg(256)->Arg(512)->Arg(1024)->Arg(2048)->Arg(4096);

BENCHMARK_MAIN();
