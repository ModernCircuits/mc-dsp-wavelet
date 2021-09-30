#include "lt/core/AudioFile.h"
#include "lt/dsp/convolution.hpp"
#include "lt/dsp/wavelets.hpp"

#include "lt/algorithm.hpp"
#include "lt/cmath.hpp"
#include "lt/format.hpp"
#include "lt/numeric.hpp"
#include "lt/preprocessor.hpp"
#include "lt/utility.hpp"

#include "lt/testing/test.hpp"
#include "readFileToVector.hpp"

template <typename It>
auto mean(It f, It l) -> float
{
    auto const sum = std::accumulate(f, l, 0.0F);
    return sum / static_cast<float>(std::distance(f, l));
}

auto peakDetect(lt::span<float> data) -> std::size_t
{
    auto peaks = std::minmax_element(data.begin(), data.end());
    if (std::fabs(*peaks.first) >= std::fabs(*peaks.second)) {
        return std::distance(data.begin(), peaks.first);
    }
    return std::distance(data.begin(), peaks.second);
}

struct BpmDetect {
    BpmDetect(std::size_t n, std::size_t levels)
        : wave_ { "db4" }
        , wt_ { wave_, "dwt", n, levels }
    {
        wt_.extension(SignalExtension::symmetric);
        wt_.convMethod(ConvolutionMethod::direct);
    }

    LT_NODISCARD auto perform(lt::span<float> input, float sampleRate) -> float
    {

        auto const levels = 4;
        auto const maxDecimation = std::pow(2.0F, levels - 1);
        auto const minNdx = std::floor(60.0F / 220.0F * (sampleRate / maxDecimation));
        auto const maxNdx = std::floor(60.0F / 40.0F * (sampleRate / maxDecimation));

        auto cA = lt::span<float> {};
        auto cD = lt::span<float> {};

        auto cDMinlen = 0.0F;
        auto cDSum = std::vector<float> {};
        dwt(wt_, input.data());
        for (auto loop { 0 }; loop < levels; ++loop) {
            if (loop == 0) {
                // dwt(wt_, input.data());
                cA = wt_.approx();
                cD = wt_.detail(loop + 1);
                cDMinlen = static_cast<float>(lt::size(cD)) / maxDecimation + 1.0F;
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
            std::transform(std::begin(cD), std::end(cD), std::begin(cD), [m](auto v) { return v - m; });

            // 5) Decimate for reconstruction later.
            // cD = abs(cD[:: (2 ** (levels - loop - 1))])
            cD = cD.subspan(0, static_cast<std::size_t>(std::pow(2, levels - loop - 1)));
            std::transform(std::begin(cD), std::end(cD), std::begin(cD), [](auto v) { return std::fabs(v); });

            // 6) Recombine the signal before ACF
            //    Essentially, each level the detail coefs (i.e. the HPF values)
            //    are concatenated to the beginning of the array
            // cD_sum = cD[0: math.floor(cD_minlen)] + cD_sum
            cDSum.insert(begin(cDSum), std::begin(cD), std::next(std::begin(cD), static_cast<size_t>(std::floor(cDMinlen))));
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
        std::transform(std::begin(cA), std::end(cA), std::begin(cA), [](auto v) { return std::fabs(v); });
        auto const m = mean(std::begin(cA), std::end(cA));
        std::transform(std::begin(cA), std::end(cA), std::begin(cA), [m](auto v) { return v - m; });
        cDSum.insert(begin(cDSum), std::begin(cA), std::next(std::begin(cA), static_cast<size_t>(std::floor(cDMinlen))));

        // # ACF
        // correl = np.correlate(cD_sum, cD_sum, "full")
        // cD_sum.clear();
        std::copy(std::begin(wt_.output()), std::begin(wt_.output()) + wt_.outlength, std::back_inserter(cDSum));
        // std::transform(begin(cD_sumf_), end(cD_sumf_), begin(cD_sumf_), [](auto v) { return std::fabs(v); });
        // auto const m = mean(begin(cD_sumf_), end(cD_sumf_));
        // std::transform(begin(cD_sumf_), end(cD_sumf_), begin(cD_sumf_), [m](auto v) { return v - m; });

        auto s = FloatSignal(cDSum.data(), cDSum.size());
        auto x = OverlapSaveConvolver(s, s);
        x.crossCorrelate();
        auto correl = x.extractResult();

        auto midpoint = static_cast<std::size_t>(std::floor(correl.size() / 2.0F));
        auto correlMidpointTmp = lt::span<float> { correl.data(), correl.size() }.subspan(midpoint);
        auto const peakNdx = peakDetect(correlMidpointTmp.subspan(minNdx, maxNdx - minNdx));

        auto const peakNdxAdjusted = peakNdx + minNdx;
        auto const bpm = 60.0F / peakNdxAdjusted * (sampleRate / maxDecimation);
        return bpm;
    }

private:
    Wavelet wave_;
    WaveletTransform wt_;
};

auto median(lt::span<float> data) -> float
{
    if (lt::empty(data)) {
        return 0.0;
    }

    std::sort(std::begin(data), std::end(data));
    auto const size = lt::size(data);
    auto const mid = size / 2;
    return size % 2 == 0 ? (data[mid] + data[mid - 1]) / 2 : data[mid];
}

auto mode(lt::span<float> arr) -> float
{
    auto const n = arr.size();
    float count = 1;
    float countmax = 0;
    float current = arr[0];
    float moda = 0;
    for (std::size_t i = 1; i < n; i++) {
        if (arr[i] == current) {
            count++;
        } else if (count > countmax) {
            countmax = count;
            count = 1;
            moda = arr[i - 1];
        }
        current = arr[i];
    }
    return moda;
}

auto main(int argc, char** argv) -> int
{
    if (argc != 2) {
        fmt::print("no file path provided");
        return EXIT_FAILURE;
    }

    AudioFile<float> audioFile;
    audioFile.load(argv[1]);
    audioFile.printSummary();

    auto const fs = static_cast<float>(audioFile.getSampleRate());
    // auto const numSamples = static_cast<size_t>(audioFile.getNumSamplesPerChannel());

    auto channel = lt::span<float>(audioFile.samples[0].data(), audioFile.samples[0].size());
    // channel = channel.subspan(static_cast<std::size_t>(lt::size(audioFile.samples[0]) * (0.05 / 100.0)));
    // channel = channel.last(static_cast<std::size_t>(lt::size(audioFile.samples[0]) * (0.05 / 100.0)));

    auto const windowSize = static_cast<std::size_t>(std::floor(3.0 * fs));
    auto const maxWindowIndex = lt::size(channel) / windowSize;

    auto sampsNdx = 0U;
    auto bpms = std::vector<float> {};
    for (auto windowNdx { 0U }; windowNdx < maxWindowIndex; ++windowNdx) {
        if (sampsNdx + windowSize >= lt::size(channel)) {
            fmt::print("ERROR");
            continue;
        }

        auto detector = BpmDetect { windowSize, 4 };
        auto subBuffer = channel.subspan(sampsNdx, windowSize);
        auto bpm = detector.perform(subBuffer, fs);

        fmt::printf("BPM:%.1f\n", std::round(bpm));
        if (bpm == 0.0) {
            continue;
        }

        bpms.push_back(bpm);
        sampsNdx += windowSize;
    }

    auto outOfRange = [](auto b) { return b < 100.0 || b > 200.0; };
    bpms.erase(std::remove_if(begin(bpms), end(bpms), outOfRange), end(bpms));

    fmt::printf("Detected BPM (median): %.1f\n", std::round(median(bpms)));
    fmt::printf("Detected BPM (mode): %.1f\n", std::round(mode(bpms)));
    return 0;
}