import array
import sys
import wave

import matplotlib.pyplot as plt
import numpy as np
from scipy.fft import fft, ifft
from scipy import signal


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


data, fs = read_wav(sys.argv[1])
for s in data:
    print(f'{s}.0')

# cA = signal.lfilter([0.01], [1 - 0.99], data)

# dataY = np.fft.fft(data)
# cAY = np.fft.fft(cA)
# freq = np.fft.fftfreq(len(data), 1000.0/fs)

# fig, axs = plt.subplots(4)
# axs[0].plot(data, color='g')
# axs[1].plot(cA, color='r')
# axs[2].plot(freq, np.abs(dataY), color='g')
# axs[3].plot(freq, np.abs(cAY), color='r')

# plt.show(block=True)
