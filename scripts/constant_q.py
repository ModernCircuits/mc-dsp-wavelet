import sys

import librosa
import librosa.display

import matplotlib.pyplot as plt
import numpy as np

path = librosa.ex('trumpet')
if len(sys.argv) == 2:
    path = sys.argv[1]

y, sr = librosa.load(path, sr=44100/2)

C = np.abs(librosa.cqt(y, sr=sr))

fig, ax = plt.subplots()

img = librosa.display.specshow(librosa.amplitude_to_db(
    C, ref=np.max),         sr=sr, x_axis='time', y_axis='cqt_hz', ax=ax)

ax.set_title('Constant-Q power spectrum')

fig.colorbar(img, ax=ax, format="%+2.0f dB")
plt.show(block=True)
