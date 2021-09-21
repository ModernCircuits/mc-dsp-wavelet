import argparse
import array
import math
from pathlib import Path
import wave

import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import pywt
from scipy import signal, stats


def read_wav(filename):
    # open file, get metadata for audio
    try:
        wf = wave.open(filename, "rb")
    except IOError as e:
        print(e)
        return

    # typ = choose_type( wf.getsampwidth() ) # TODO: implement choose_type
    nsamps = wf.getnframes()
    assert nsamps > 0

    fs = wf.getframerate()
    assert fs > 0

    # Read entire file and make into an array
    samps = list(array.array("i", wf.readframes(nsamps)))

    try:
        assert nsamps == len(samps)
    except AssertionError:
        print(nsamps, "not equal to", len(samps))

    return samps, fs


# print an error when no data can be found
def no_audio_data():
    print("No audio data for sample, skipping...")
    return None, None


# simple peak detection
def peak_detect(data):
    max_val = np.amax(abs(data))
    peak_ndx = np.where(data == max_val)
    if len(peak_ndx[0]) == 0:  # if nothing found then the max must be negative
        peak_ndx = np.where(data == -max_val)
    return peak_ndx


def bpm_detector(data, fs):
    cA = []
    cD = []
    correl = []
    cD_sum = []
    levels = 4
    max_decimation = 2 ** (levels - 1)
    min_ndx = math.floor(60.0 / 220 * (fs / max_decimation))
    max_ndx = math.floor(60.0 / 40 * (fs / max_decimation))

    for loop in range(0, levels):
        cD = []
        # 1) DWT
        if loop == 0:
            [cA, cD] = pywt.dwt(data, "db4")
            cD_minlen = len(cD) / max_decimation + 1
            cD_sum = np.zeros(math.floor(cD_minlen))
        else:
            [cA, cD] = pywt.dwt(cA, "db4")

        # 2) Filter
        cD = signal.lfilter([0.01], [1 - 0.99], cD)

        # 4) Subtract out the mean.

        # 5) Decimate for reconstruction later.
        cD = abs(cD[:: (2 ** (levels - loop - 1))])
        cD = cD - np.mean(cD)

        # 6) Recombine the signal before ACF
        #    Essentially, each level the detail coefs (i.e. the HPF values) are concatenated to the beginning of the array
        cD_sum = cD[0: math.floor(cD_minlen)] + cD_sum

    if [b for b in cA if b != 0.0] == []:
        return no_audio_data()

    # Adding in the approximate data as well...
    cA = signal.lfilter([0.01], [1 - 0.99], cA)
    cA = abs(cA)
    cA = cA - np.mean(cA)
    cD_sum = cA[0: math.floor(cD_minlen)] + cD_sum

    # ACF
    correl = np.correlate(cD_sum, cD_sum, "full")

    midpoint = math.floor(len(correl) / 2)
    correl_midpoint_tmp = correl[midpoint:]
    peak_ndx = peak_detect(correl_midpoint_tmp[min_ndx:max_ndx])
    if len(peak_ndx) > 1:
        return no_audio_data()

    peak_ndx_adjusted = peak_ndx[0] + min_ndx
    bpm = 60.0 / peak_ndx_adjusted * (fs / max_decimation)
    print(bpm)
    return bpm, correl


def running_median(data):
    medians = [0]
    for i in range(1, len(data)):
        medians.append(np.median(data[:i]))
    return medians


def running_mode(data):
    modes = [0]
    for i in range(1, len(data)):
        modes.append(stats.mode(data[:i], axis=None)[0])
    return modes


def parse_beat_fraction(beat):
    parts = beat.split("/")
    beat = [int(part) for part in parts]
    return beat[0], beat[1]


def main():
    parser = argparse.ArgumentParser(
        description="Process .wav file to determine the Beats Per Minute.")
    parser.add_argument("--output", required=True,
                        help="Output directory for plots")
    parser.add_argument("--filename", required=True,
                        help=".wav file for processing")
    parser.add_argument(
        "--window",
        type=float,
        default=3,
        help="Size of the the window (seconds) that will be scanned to determine the bpm. Typically less than 10 seconds. [3]",
    )
    parser.add_argument(
        "--min",
        type=float,
        default=60,
        help="Minimum bpm count to consider in detection",
    )
    parser.add_argument(
        "--max",
        type=float,
        default=200,
        help="Maximum bpm count to consider in detection",
    )
    parser.add_argument(
        "--skip-start",
        type=float,
        default=0,
        help="Skip x percent of the samples. From the start",
    )
    parser.add_argument(
        "--skip-end",
        type=float,
        default=0,
        help="Skip x percent of the samples. From the end",
    )
    parser.add_argument(
        "--beat",
        "-B",
        default='4/4',
        help="Beats (4/4, 5/4, etc.).",
    )

    args = parser.parse_args()
    samps, fs = read_wav(args.filename)
    data = []
    # correl = []
    bpm = 0
    samps = samps[int(len(samps) * (args.skip_start/100.0)):]
    samps = samps[:len(samps)-int(len(samps) * (args.skip_end/100.0))]
    nsamps = len(samps)
    window_samps = int(args.window * fs)
    samps_ndx = 0  # First sample in window_ndx
    max_window_ndx = math.floor(nsamps / window_samps)
    bpms = np.zeros(max_window_ndx)

    # Iterate through all windows
    for window_ndx in range(0, max_window_ndx):
        data = samps[samps_ndx: samps_ndx + window_samps]
        if not ((len(data) % window_samps) == 0):
            raise AssertionError(str(len(data)))

        bpm, correl_temp = bpm_detector(data, fs)
        if bpm is None:
            continue
        bpms[window_ndx] = bpm
        # correl = correl_temp

        # Iterate at the end of the loop
        samps_ndx = samps_ndx + window_samps

    beat_den, beat_num = parse_beat_fraction(args.beat)
    bpms = bpms / beat_num * beat_den

    filtered_bpms = bpms[(bpms >= args.min) & (bpms <= args.max)]

    bpm_median = round(np.median(filtered_bpms), 1)
    bpm_mode = round(stats.mode(filtered_bpms, axis=None)[0][0], 1)
    print(
        f"Completed!  Estimated Beats Per Minute: {bpm_median} / {bpm_mode}")

    medians = running_median(filtered_bpms)
    n = range(0, len(medians))
    modes = running_mode(filtered_bpms)
    # plt.bar(n, medians)

    path = Path(args.filename)
    fig, axs = plt.subplots(4)
    scale = 6
    fig.set_size_inches(3*scale, 2*scale)
    fig.suptitle(
        f'File: {path.name} -- Range: {args.min}-{args.max} -- Time: {args.skip_start}% - {100-args.skip_end}% -- Signature: {args.beat} -- ME: {bpm_median} bpm -- MO: {bpm_mode} bpm')

    signal = np.array(samps)
    downsampled_signal = np.delete(signal, np.arange(0, signal.size, 128))
    axs[0].set_title('signal')
    axs[0].plot(downsampled_signal)
    axs[1].set_title('bpms')
    axs[1].bar(n, filtered_bpms)
    axs[2].set_title('prediction (median)')
    axs[2].plot(n, medians)
    axs[3].set_title('prediction (mode)')
    axs[3].plot(n, modes)

    plt.savefig(f'{args.output}/{path.name}.pdf', dpi=200)


if __name__ == "__main__":
    main()
