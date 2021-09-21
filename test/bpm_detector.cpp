#include "wavelib.h"

#include "AudioFile.h"
#include "helper.hpp"
#include "overlap_save_convolver.hpp"
#include "readFileToVector.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <utility>

auto approx_coeffs(wt_set const& wt) -> std::vector<double>
{
    auto size = static_cast<std::size_t>(wt.length[0]);
    auto result = std::vector<double>(size);
    std::copy(wt.output, wt.output + size, begin(result));
    return result;
}

auto detail_coeffs(wt_set const& wt) -> std::vector<double>
{
    auto const* first = wt.output + wt.length[0];
    auto size = static_cast<std::size_t>(wt.length[1]);
    auto result = std::vector<double>(size);
    std::copy(first, first + size, begin(result));
    return result;
}

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
        , wt_ { wt_init(wave_, "dwt", N, levels) }
    {
        setDWTExtension(wt_.get(), "sym");
        setWTConv(wt_.get(), "fft");
    }

    [[nodiscard]] auto perform(lt::span<double> input, double sampleRate) -> double
    {

        auto const levels = 4;
        auto const maxDecimation = std::pow(2.0, levels - 1);
        auto const minNdx = std::floor(60.0 / 220.0 * (sampleRate / maxDecimation));
        auto const maxNdx = std::floor(60.0 / 40.0 * (sampleRate / maxDecimation));

        auto cD_minlen = 0.0;
        for (auto loop { 0 }; loop < levels; ++loop) {
            cD_.clear();
            if (loop == 0) {
                dwt(wt_.get(), input.data());
                cA_ = approx_coeffs(*wt_);
                cD_ = detail_coeffs(*wt_);
                cD_minlen = static_cast<double>(size(cD_)) / maxDecimation + 1.0;
                cD_sum_.resize(static_cast<std::size_t>(std::floor(cD_minlen)));
                std::fill(begin(cD_sum_), end(cD_sum_), 0.0);
            } else {
                dwt(wt_.get(), cA_.data());
            }

            // 2) Filter
            // cD = signal.lfilter([0.01], [1 - 0.99], cD)

            // 4) Subtract out the mean.
            // cD = cD - np.mean(cD)
            auto const m = mean(begin(cD_), end(cD_));
            std::transform(begin(cD_), end(cD_), begin(cD_), [m](auto v) { return v - m; });

            // 5) Decimate for reconstruction later.
            // cD = abs(cD[:: (2 ** (levels - loop - 1))])
            cD_.resize(static_cast<std::size_t>(std::pow(2, levels - loop - 1)));
            std::transform(begin(cD_), end(cD_), begin(cD_), [](auto v) { return std::fabs(v); });

            // 6) Recombine the signal before ACF
            //    Essentially, each level the detail coefs (i.e. the HPF values)
            //    are concatenated to the beginning of the array
            // cD_sum = cD[0: math.floor(cD_minlen)] + cD_sum
            cD_sum_.insert(begin(cD_sum_), begin(cD_), std::next(begin(cD_), static_cast<size_t>(std::floor(cD_minlen))));
        }

        // if [b for b in cA if b != 0.0] == []:
        //     return no_audio_data()
        if (std::none_of(begin(cA_), end(cA_), [](auto s) { return s != 0.0; })) {
            return 0.0;
        }

        // # Adding in the approximate data as well...
        // cA = signal.lfilter([0.01], [1 - 0.99], cA)
        // cA = abs(cA)
        // cA = cA - np.mean(cA)
        // cD_sum = cA[0: math.floor(cD_minlen)] + cD_sum
        std::transform(begin(cA_), end(cA_), begin(cA_), [](auto v) { return std::fabs(v); });
        auto const m = mean(begin(cA_), end(cA_));
        std::transform(begin(cA_), end(cA_), begin(cA_), [m](auto v) { return v - m; });
        cD_sum_.insert(begin(cD_sum_), begin(cA_), std::next(begin(cA_), static_cast<size_t>(std::floor(cD_minlen))));

        // # ACF
        // correl = np.correlate(cD_sum, cD_sum, "full")
        auto cD_sumf = std::vector<float> {};
        std::copy(begin(cD_sum_), end(cD_sum_), std::back_inserter(cD_sumf));
        FloatSignal s(cD_sumf.data(), cD_sumf.size());
        OverlapSaveConvolver x(s, s);
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
    std::unique_ptr<wt_set> wt_;

    std::vector<double> cA_ {};
    std::vector<double> cD_ {};
    std::vector<double> correl_ {};
    std::vector<double> cD_sum_ {};
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
    auto detector = bpm_detect { std::size(channel), 4 };

    printf("windowSize: %zu\n", windowSize);
    printf("maxWindowIndex: %zu\n", maxWindowIndex);

    auto samps_ndx = 0U;
    auto bpms = makeZeros<double>(maxWindowIndex);
    for (auto window_ndx { 0U }; window_ndx < maxWindowIndex; ++window_ndx) {
        if (samps_ndx + windowSize >= std::size(channel)) {
            std::puts("ERROR");
            continue;
        }

        auto subBuffer = channel.subspan(samps_ndx, windowSize);
        printf("size(subBuffer): %zu\n", std::size(subBuffer));

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