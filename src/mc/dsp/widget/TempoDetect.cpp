#include "TempoDetect.hpp"

#include <mc/core/config.hpp>

#include <mc/dsp/convolution.hpp>

#include <mc/core/algorithm.hpp>
#include <mc/core/cmath.hpp>
#include <mc/core/format.hpp>
#include <mc/core/numeric.hpp>
#include <mc/core/utility.hpp>

namespace mc::dsp {

template<typename It>
auto mean(It f, It l) -> float
{
    auto const sum = std::accumulate(f, l, 0.0F);
    return sum / static_cast<float>(std::distance(f, l));
}

auto peakDetect(Span<float> data) -> std::size_t
{
    auto peaks = std::minmax_element(data.begin(), data.end());
    if (std::fabs(*peaks.first) >= std::fabs(*peaks.second)) {
        return std::distance(data.begin(), peaks.first);
    }
    return std::distance(data.begin(), peaks.second);
}

TempoDetect::TempoDetect(std::size_t n, std::size_t levels)
    : wave_{"db4"}
    , wt_{wave_, "dwt", n, levels}
{
    wt_.extension(dsp::SignalExtension::symmetric);
    wt_.convMethod(dsp::ConvolutionMethod::direct);
}

auto TempoDetect::operator()(Span<float> input, float sampleRate) -> float
{

    auto const levels        = wt_.levels();
    auto const maxDecimation = std::pow(2.0F, static_cast<float>(levels - 1));
    auto const minNdx        = std::floor(60.0F / 220.0F * (sampleRate / maxDecimation));
    auto const maxNdx        = std::floor(60.0F / 40.0F * (sampleRate / maxDecimation));

    auto cA = Span<float>{};
    auto cD = Span<float>{};

    auto cDMinlen = 0.0F;
    auto cDSum    = Vector<float>{};
    dwt(wt_, input.data());
    for (auto loop{0}; loop < levels; ++loop) {
        if (loop == 0) {
            // dwt(wt_, input.data());
            cA       = wt_.approx();
            cD       = wt_.detail(loop + 1);
            cDMinlen = static_cast<float>(mc::size(cD)) / maxDecimation + 1.0F;
            cDSum.resize(static_cast<std::size_t>(std::floor(cDMinlen)));
            std::fill(begin(cDSum), end(cDSum), 0.0F);
        } else {
            // dwt(wt_, cA.data());
            cA = wt_.approx();
            cD = wt_.detail(loop + 1);
        }

        // 2) Filter
        // cD = signal.lfilter([0.01], [1 - 0.99], cD)

        // 4) Subtract out the mean.
        // cD = cD - np.mean(cD)
        auto const m = mean(std::begin(cD), std::end(cD));
        std::transform(std::begin(cD), std::end(cD), std::begin(cD), [m](auto v) {
            return v - m;
        });

        // 5) Decimate for reconstruction later.
        // cD = abs(cD[:: (2 ** (levels - loop - 1))])
        cD = cD.subspan(0, static_cast<std::size_t>(std::pow(2, levels - loop - 1)));
        std::transform(std::begin(cD), std::end(cD), std::begin(cD), [](auto v) {
            return std::fabs(v);
        });

        // 6) Recombine the signal before ACF
        //    Essentially, each level the detail coefs (i.e. the HPF values)
        //    are concatenated to the beginning of the array
        // cD_sum = cD[0: math.floor(cD_minlen)] + cD_sum
        cDSum.insert(
            begin(cDSum),
            std::begin(cD),
            std::next(std::begin(cD), static_cast<size_t>(std::floor(cDMinlen)))
        );
    }

    // if [b for b in cA if b != 0.0] == []:
    //     return no_audio_data()
    if (std::none_of(std::begin(cA), std::end(cA), [](auto s) { return s != 0.0F; })) {
        return 0.0F;
    }

    // # Adding in the approximate data as well...
    // cA = signal.lfilter([0.01], [1 - 0.99], cA)
    // cA = abs(cA)
    // cA = cA - np.mean(cA)
    // cD_sum = cA[0: math.floor(cD_minlen)] + cD_sum
    std::transform(std::begin(cA), std::end(cA), std::begin(cA), [](auto v) {
        return std::fabs(v);
    });
    auto const m = mean(std::begin(cA), std::end(cA));
    std::transform(std::begin(cA), std::end(cA), std::begin(cA), [m](auto v) {
        return v - m;
    });
    cDSum.insert(
        begin(cDSum),
        std::begin(cA),
        std::next(std::begin(cA), static_cast<size_t>(std::floor(cDMinlen)))
    );

    // # ACF
    // correl = np.correlate(cD_sum, cD_sum, "full")
    // cD_sum.clear();
    std::copy(
        std::begin(wt_.output()),
        std::begin(wt_.output()) + wt_.outlength,
        std::back_inserter(cDSum)
    );
    // std::transform(begin(cD_sumf_), end(cD_sumf_), begin(cD_sumf_), [](auto v) { return
    // std::fabs(v); }); auto const m = mean(begin(cD_sumf_), end(cD_sumf_));
    // std::transform(begin(cD_sumf_), end(cD_sumf_), begin(cD_sumf_), [m](auto v) { return
    // v - m; });

    auto s = dsp::FloatSignal(cDSum.data(), cDSum.size());
    auto x = dsp::OverlapSaveConvolver(s, s);
    x.crossCorrelate();
    auto correl = x.extractResult();

    auto midpoint          = static_cast<std::size_t>(std::floor(correl.size() / 2.0F));
    auto correlMidpointTmp = Span<float>{correl.data(), correl.size()}.subspan(midpoint);
    auto const peakNdx     = peakDetect(correlMidpointTmp.subspan(minNdx, maxNdx - minNdx));

    auto const peakNdxAdjusted = peakNdx + minNdx;
    auto const bpm             = 60.0F / peakNdxAdjusted * (sampleRate / maxDecimation);
    return bpm;
}

}  // namespace mc::dsp
