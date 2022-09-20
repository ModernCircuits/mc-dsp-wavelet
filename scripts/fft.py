import numpy as np


data = [0.0, 1.0, 0.3, 0.0]
print(data)
y = np.fft.fft(data)
print(y)
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
