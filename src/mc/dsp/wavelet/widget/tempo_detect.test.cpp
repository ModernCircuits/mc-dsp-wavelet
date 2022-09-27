// SPDX-License-Identifier: BSL-1.0

#include <mc/core/config.hpp>

#include <mc/dsp/algorithm/median.hpp>
#include <mc/dsp/algorithm/mode.hpp>
#include <mc/dsp/wavelet/widget/tempo_detect.hpp>

#include <mc/audio/audio_file.hpp>
#include <mc/core/algorithm.hpp>
#include <mc/core/cmath.hpp>
#include <mc/core/iterator.hpp>
#include <mc/core/numeric.hpp>
#include <mc/core/span.hpp>
#include <mc/core/utility.hpp>

#include <catch2/catch_test_macros.hpp>

TEST_CASE("dsp/wavelet: TempoDetect", "[dsp][wavelet]")
{
    auto audioFile = mc::AudioFile<float>{};
    audioFile.load("./test_data/wav/Drums.wav");

    auto const fs = static_cast<float>(audioFile.getSampleRate());
    auto channel  = mc::Span<float>(audioFile.samples[0]);

    auto const windowSize     = static_cast<size_t>(std::floor(3.0 * fs));
    auto const maxWindowIndex = mc::size(channel) / windowSize;

    auto sampsNdx = size_t{0};
    auto bpms     = mc::Vector<float>{};
    for (auto windowNdx{0U}; windowNdx < maxWindowIndex; ++windowNdx) {
        if (sampsNdx + windowSize >= mc::size(channel)) { continue; }

        auto detector  = mc::TempoDetect{windowSize, 4};
        auto subBuffer = channel.subspan(sampsNdx, windowSize);
        auto bpm       = detector(subBuffer, fs);

        if (bpm == 0.0) { continue; }

        bpms.push_back(bpm);
        sampsNdx += windowSize;
    }

    auto outOfRange = [](auto b) { return b < 100.0 || b > 200.0; };
    bpms.erase(ranges::remove_if(bpms, outOfRange), end(bpms));

    ranges::sort(bpms);
    REQUIRE(std::round(mc::median<float>(bpms)) == 120.0F);
    // REQUIRE(std::round(mode(bpms)) == 120.0);
}
