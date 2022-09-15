#include <mc/core/config.hpp>

#include <mc/dsp/widget/TempoDetect.hpp>

#include <mc/core/algorithm.hpp>
#include <mc/core/AudioFile.h>
#include <mc/core/cmath.hpp>
#include <mc/core/format.hpp>
#include <mc/core/iterator.hpp>
#include <mc/core/numeric.hpp>
#include <mc/core/span.hpp>
#include <mc/core/utility.hpp>

#include <catch2/catch_test_macros.hpp>

auto median(mc::Span<float> data) -> float
{
    if (mc::empty(data)) { return 0.0; }

    std::sort(std::begin(data), std::end(data));
    auto const size = mc::size(data);
    auto const mid  = size / 2;
    return size % 2 == 0 ? (data[mid] + data[mid - 1]) / 2 : data[mid];
}

auto mode(mc::Span<float> arr) -> float
{
    auto const n   = arr.size();
    float count    = 1;
    float countmax = 0;
    float current  = arr[0];
    float moda     = 0;
    for (std::size_t i = 1; i < n; i++) {
        if (arr[i] == current) {
            count++;
        } else if (count > countmax) {
            countmax = count;
            count    = 1;
            moda     = arr[i - 1];
        }
        current = arr[i];
    }
    return moda;
}

TEST_CASE("dsp/wavelet: TempoDetect", "[dsp][wavelet]")
{
    auto audioFile = AudioFile<float>{};
    audioFile.load("./test_data/wav/Drums.wav");

    auto const fs = static_cast<float>(audioFile.getSampleRate());
    auto channel
        = mc::Span<float>(audioFile.samples[0].data(), audioFile.samples[0].size());

    auto const windowSize     = static_cast<std::size_t>(std::floor(3.0 * fs));
    auto const maxWindowIndex = mc::size(channel) / windowSize;

    auto sampsNdx = 0U;
    auto bpms     = std::vector<float>{};
    for (auto windowNdx{0U}; windowNdx < maxWindowIndex; ++windowNdx) {
        if (sampsNdx + windowSize >= mc::size(channel)) { continue; }

        auto detector  = mc::dsp::TempoDetect{windowSize, 4};
        auto subBuffer = channel.subspan(sampsNdx, windowSize);
        auto bpm       = detector(subBuffer, fs);

        if (bpm == 0.0) { continue; }

        bpms.push_back(bpm);
        sampsNdx += windowSize;
    }

    auto outOfRange = [](auto b) { return b < 100.0 || b > 200.0; };
    bpms.erase(std::remove_if(begin(bpms), end(bpms), outOfRange), end(bpms));

    REQUIRE(std::round(median(bpms)) == 120.0);
    // REQUIRE(std::round(mode(bpms)) == 120.0);
}
