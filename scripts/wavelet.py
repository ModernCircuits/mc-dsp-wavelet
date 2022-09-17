import pywt


def format_array(arr):
    return "[%s]" % ", ".join(["%.14f" % x for x in arr])


wavelet = pywt.Wavelet('db2')
print(wavelet)
print(format_array(wavelet.dec_lo), format_array(wavelet.dec_hi))
print(format_array(wavelet.rec_lo), format_array(wavelet.rec_hi))
