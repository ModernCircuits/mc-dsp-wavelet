// SPDX-License-Identifier: BSL-1.0

#include <mc/core/config.hpp>

#include <mc/wavelet/algorithm/median.hpp>
#include <mc/wavelet/algorithm/mode.hpp>
#include <mc/wavelet/widget/tempo_detect.hpp>

#include <mc/core/_utility/pair.hpp>
#include <mc/core/algorithm.hpp>
#include <mc/core/cmath.hpp>
#include <mc/core/iterator.hpp>
#include <mc/core/numeric.hpp>
#include <mc/core/span.hpp>
#include <mc/core/utility.hpp>

#include <catch2/catch_test_macros.hpp>

#include "dr_wav.h"

using namespace mc;

auto readWavFile(char const* path) -> Pair<Vector<Vector<float>>, unsigned>
{
    auto ch     = unsigned{0};
    auto sr     = unsigned{0};
    auto frames = drwav_uint64{};
    auto* s     = drwav_open_file_and_read_pcm_frames_f32(path, &ch, &sr, &frames, nullptr);
    if (s == nullptr) { return {}; }

    auto buffer = Vector<Vector<float>>(ch);
    for (auto i{0U}; i < ch; ++i) { buffer[i].resize(frames); }

    for (auto f{0U}; f < frames; ++f) {
        for (auto i{0U}; i < ch; ++i) { buffer[i][f] = s[(f * ch) + i]; }
    }

    return {buffer, sr};
}

TEST_CASE("wavelet: TempoDetect", "[dsp][wavelet]")
{
    auto const [audioFile, sr] = readWavFile("./test_data/wav/Drums.wav");

    auto const fs = static_cast<float>(sr);
    auto channel  = Span<float const>(audioFile[0]);

    auto const windowSize     = static_cast<size_t>(std::floor(3.0 * fs));
    auto const maxWindowIndex = size(channel) / windowSize;

    auto sampsNdx = size_t{0};
    auto bpms     = Vector<float>{};
    for (auto windowNdx{0U}; windowNdx < maxWindowIndex; ++windowNdx) {
        if (sampsNdx + windowSize >= size(channel)) { continue; }

        auto detector  = TempoDetect{windowSize, 4};
        auto subBuffer = channel.subspan(sampsNdx, windowSize);
        auto bpm       = detector(subBuffer, fs);

        if (bpm == 0.0) { continue; }

        bpms.push_back(bpm);
        sampsNdx += windowSize;
    }

    auto outOfRange = [](auto b) { return b < 100.0 || b > 200.0; };
    bpms.erase(ranges::remove_if(bpms, outOfRange), end(bpms));

    ranges::sort(bpms);
    REQUIRE(std::round(median<float>(bpms)) == 120.0F);
    // REQUIRE(std::round(mode(bpms)) == 120.0);
}
