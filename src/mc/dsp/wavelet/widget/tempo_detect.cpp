// SPDX-License-Identifier: BSL-1.0

#include "tempo_detect.hpp"

#include <mc/dsp/algorithm.hpp>
#include <mc/dsp/convolution.hpp>

#include <mc/core/algorithm.hpp>
#include <mc/core/cmath.hpp>
#include <mc/core/numeric.hpp>
#include <mc/core/utility.hpp>

namespace mc {

TempoDetect::TempoDetect(std::size_t n, std::size_t levels)
    : _wave{"db4"}
    , _wt{_wave, "dwt", n, levels}
{
    _wt.extension(SignalExtension::symmetric);
    _wt.convMethod(ConvolutionMethod::direct);
}

auto TempoDetect::operator()(Span<float> input, float sampleRate) -> float
{

    auto const levels        = _wt.levels();
    auto const maxDecimation = std::pow(2.0F, static_cast<float>(levels - 1));
    auto const minNdx        = std::floor(60.0F / 220.0F * (sampleRate / maxDecimation));
    auto const maxNdx        = std::floor(60.0F / 40.0F * (sampleRate / maxDecimation));

    dwt(_wt, input.data());

    auto cD             = _wt.detail(1);
    auto const cA       = _wt.approx();
    auto const cDMinlen = static_cast<size_t>((float)(mc::size(cD)) / maxDecimation + 1.0F);

    auto cDSum = Vector<float>(cDMinlen);

    for (auto loop{0}; loop < levels; ++loop) {
        cD = _wt.detail(loop + 1);

        // 2) Filter
        // cD = signal.lfilter([0.01], [1 - 0.99], cD)

        // 4) Subtract out the mean.
        // cD = cD - np.mean(cD)
        auto const m = mean(cD);
        ranges::transform(cD, begin(cD), [m](auto v) { return v - m; });

        // 5) Decimate for reconstruction later.
        // cD = abs(cD[:: (2 ** (levels - loop - 1))])
        cD = cD.subspan(0, static_cast<std::size_t>(std::pow(2, levels - loop - 1)));
        ranges::transform(cD, begin(cD), [](auto v) { return std::abs(v); });

        // 6) Recombine the signal before ACF
        //    Essentially, each level the detail coefs (i.e. the HPF values)
        //    are concatenated to the beginning of the array
        // cD_sum = cD[0: math.floor(cD_minlen)] + cD_sum
        cDSum.insert(begin(cDSum), std::begin(cD), std::next(std::begin(cD), cDMinlen));
    }

    // if [b for b in cA if b != 0.0] == []:
    //     return no_audio_data()
    if (ranges::none_of(cA, [](auto s) { return s != 0.0F; })) { return 0.0F; }

    // # Adding in the approximate data as well...
    // cA = signal.lfilter([0.01], [1 - 0.99], cA)
    // cA = abs(cA)
    // cA = cA - np.mean(cA)
    // cD_sum = cA[0: math.floor(cD_minlen)] + cD_sum
    ranges::transform(cA, std::begin(cA), [](auto v) { return std::abs(v); });
    auto const m = mean(cA);
    ranges::transform(cA, std::begin(cA), [m](auto v) { return v - m; });
    cDSum.insert(begin(cDSum), std::begin(cA), std::next(std::begin(cA), cDMinlen));

    // # ACF
    // correl = np.correlate(cD_sum, cD_sum, "full")
    // cD_sum.clear();
    std::copy(
        std::begin(_wt.output()),
        std::begin(_wt.output()) + _wt.outlength,
        std::back_inserter(cDSum)
    );
    // std::transform(begin(cD_sumf_), end(cD_sumf_), begin(cD_sumf_), [](auto v) { return
    // std::abs(v); }); auto const m = mean(begin(cD_sumf_), end(cD_sumf_));
    // std::transform(begin(cD_sumf_), end(cD_sumf_), begin(cD_sumf_), [m](auto v) { return
    // v - m; });

    auto s = FloatSignal(cDSum.data(), cDSum.size());
    auto x = OverlapSaveConvolver(s, s);
    x.crossCorrelate();
    auto correl = x.extractResult();

    auto midpoint       = static_cast<std::size_t>(std::floor(correl.size() / 2.0F));
    auto correlMidpoint = Span<float>{correl.data(), correl.size()}.subspan(midpoint);
    auto searchRange    = correlMidpoint.subspan((size_t)minNdx, (size_t)(maxNdx - minNdx));
    auto const peakNdx  = indexOfPeak(searchRange);

    auto const peakNdxAdjusted = peakNdx + minNdx;
    auto const bpm             = 60.0F / peakNdxAdjusted * (sampleRate / maxDecimation);
    return bpm;
}

}  // namespace mc
