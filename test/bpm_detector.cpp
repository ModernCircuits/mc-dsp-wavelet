#include "wavelib.h"

#include "AudioFile.h"
#include "helper.hpp"
#include "overlap_save_convolver.hpp"
#include "readFileToVector.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <utility>

auto approx_coeffs(wavelet_transform const& wt) -> std::vector<double>
{
    auto size = static_cast<std::size_t>(wt.length[0]);
    auto result = std::vector<double>(size);
    std::copy(wt.output, wt.output + size, begin(result));
    return result;
}

auto detail_coeffs(wavelet_transform const& wt) -> std::vector<double>
{
    auto const* first = wt.output + wt.length[0];
    auto size = static_cast<std::size_t>(wt.length[1]);
    auto result = std::vector<double>(size);
    std::copy(first, first + size, begin(result));
    return result;
}

// t = wt->length[0];
// for (auto i = 0; i < J; ++i) {
//     printf("Level %d Access : output[%d] Length : %d \n", J - i, t, wt->length[i + 1]);
//     t += wt->length[i + 1];
// }

template <typename It>
auto mean(It f, It l) -> double
{
    auto const sum = std::accumulate(f, l, 0.0);
    return sum / static_cast<double>(std::distance(f, l));
}

auto peak_detect(lt::span<float> data) -> std::size_t
{
    auto peaks = std::minmax_element(data.begin(), data.end());
    if (fabs(*peaks.first) >= fabs(*peaks.second)) {
        return std::distance(data.begin(), peaks.first);
    }
    return std::distance(data.begin(), peaks.second);
}

struct bpm_detect {
    bpm_detect(std::size_t N, std::size_t levels)
        : wave_ { "db4" }
        , wt_ { wave_, "dwt", static_cast<int>(N), static_cast<int>(levels) }
    {
        wt_.extension(signal_extension::symmetric);
        wt_.convolution_method("fft");
    }

    [[nodiscard]] auto perform(lt::span<double> input, double sampleRate) -> double
    {

        auto const levels = 4;
        auto const maxDecimation = std::pow(2.0, levels - 1);
        auto const minNdx = std::floor(60.0 / 220.0 * (sampleRate / maxDecimation));
        auto const maxNdx = std::floor(60.0 / 40.0 * (sampleRate / maxDecimation));

        auto cA = lt::span<double> {};
        auto cD = lt::span<double> {};

        auto cD_minlen = 0.0;
        cD_sum_.clear();
        dwt(&wt_, input.data());
        for (auto loop { 0 }; loop < levels; ++loop) {
            cD_.clear();
            cA_.clear();
            if (loop == 0) {
                // dwt(&wt_, input.data());
                cA = wt_.approx();
                cD = wt_.detail(loop + 1);
                cD_minlen = static_cast<double>(std::size(cD)) / maxDecimation + 1.0;
                cD_sum_.resize(static_cast<std::size_t>(std::floor(cD_minlen)));
                std::fill(begin(cD_sum_), end(cD_sum_), 0.0);
            } else {
                // dwt(&wt_, cA.data());
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
            cD_sum_.insert(begin(cD_sum_), std::begin(cD), std::next(std::begin(cD), static_cast<size_t>(std::floor(cD_minlen))));
        }

        // if [b for b in cA if b != 0.0] == []:
        //     return no_audio_data()
        if (std::none_of(std::begin(cA), std::end(cA), [](auto s) { return s != 0.0; })) {
            return 0.0;
        }

        // # Adding in the approximate data as well...
        // cA = signal.lfilter([0.01], [1 - 0.99], cA)
        // cA = abs(cA)
        // cA = cA - np.mean(cA)
        // cD_sum = cA[0: math.floor(cD_minlen)] + cD_sum
        std::transform(std::begin(cA), std::end(cA), std::begin(cA), [](auto v) { return std::fabs(v); });
        auto const m = mean(std::begin(cA), std::end(cA));
        std::transform(std::begin(cA), std::end(cA), std::begin(cA), [m](auto v) { return v - m; });
        cD_sum_.insert(begin(cD_sum_), std::begin(cA), std::next(std::begin(cA), static_cast<size_t>(std::floor(cD_minlen))));

        // # ACF
        // correl = np.correlate(cD_sum, cD_sum, "full")
        cD_sumf_.clear();
        std::copy(wt_.output, wt_.output + wt_.outlength, std::back_inserter(cD_sumf_));
        // std::transform(begin(cD_sumf_), end(cD_sumf_), begin(cD_sumf_), [](auto v) { return std::fabs(v); });
        // auto const m = mean(begin(cD_sumf_), end(cD_sumf_));
        // std::transform(begin(cD_sumf_), end(cD_sumf_), begin(cD_sumf_), [m](auto v) { return v - m; });

        auto s = FloatSignal(cD_sumf_.data(), cD_sumf_.size());
        auto x = OverlapSaveConvolver(s, s);
        x.executeXcorr();
        auto correl = x.extractResult();

        auto midpoint = static_cast<std::size_t>(std::floor(correl.size() / 2.0));
        auto correl_midpoint_tmp = lt::span<float> { correl.data(), correl.size() }.subspan(midpoint);
        auto const peak_ndx = peak_detect(correl_midpoint_tmp.subspan(minNdx, maxNdx - minNdx));

        auto const peak_ndx_adjusted = peak_ndx + minNdx;
        auto const bpm = 60.0 / peak_ndx_adjusted * (sampleRate / maxDecimation);
        return bpm;
    }

private:
    wavelet wave_;
    wavelet_transform wt_;

    std::vector<double> cA_ {};
    std::vector<double> cD_ {};
    std::vector<double> correl_ {};
    std::vector<double> cD_sum_ {};
    std::vector<float> cD_sumf_ {};
};

auto main(int argc, char** argv) -> int
{
    if (argc != 2) {
        std::puts("no file path provided");
        return EXIT_FAILURE;
    }

    AudioFile<double> audioFile;
    audioFile.load(argv[1]);
    audioFile.printSummary();

    auto const fs = static_cast<double>(audioFile.getSampleRate());
    // auto const numSamples = static_cast<size_t>(audioFile.getNumSamplesPerChannel());

    auto channel = lt::span<double>(audioFile.samples[0].data(), audioFile.samples[0].size());
    // channel = channel.subspan(static_cast<std::size_t>(std::size(audioFile.samples[0]) * (0.05 / 100.0)));
    // channel = channel.last(static_cast<std::size_t>(std::size(audioFile.samples[0]) * (0.05 / 100.0)));

    auto const windowSize = static_cast<std::size_t>(std::floor(3.0 * fs));
    auto const maxWindowIndex = std::size(channel) / windowSize;

    auto samps_ndx = 0U;
    auto bpms = makeZeros<double>(maxWindowIndex);
    for (auto window_ndx { 0U }; window_ndx < maxWindowIndex; ++window_ndx) {
        if (samps_ndx + windowSize >= std::size(channel)) {
            std::puts("ERROR");
            continue;
        }

        auto detector = bpm_detect { windowSize, 4 };
        auto subBuffer = channel.subspan(samps_ndx, windowSize);
        auto bpm = detector.perform(subBuffer, fs);

        printf("BPM:%f\n", bpm);
        if (bpm == 0.0) {
            continue;
        }

        bpms[window_ndx] = bpm;
        samps_ndx += windowSize;
    }

    return 0;
}