import allantools.noise as noise
import numpy as np
import matplotlib.pyplot as plt


def many_psds(k=2, fs=1.0, b0=1.0, N=1024):
    """ compute average of many PSDs """
    psd = []
    for j in range(k):
        print(j)
        x = noise.white(num_points=2*4096, b0=b0, fs=fs)
        f, tmp = noise.numpy_psd(x, fs)
        if j == 0:
            psd = tmp
        else:
            psd = psd + tmp
    return f, psd/k


fs = 256  # sample rate
b0 = 0.0002
N = 2*4096
x = noise.white(num_points=2*4096, b0=b0, fs=fs)
t = np.linspace(0, (1.0/fs)*N, len(x))

plt.figure()
plt.plot(t, x)
plt.xlabel('Time / s')
plt.ylabel('Amplitude / V')
print(x)
f, psd = many_psds(k=50, fs=fs, b0=b0, N=N)

fxx, Pxx_den = noise.scipy_psd(x, fs)

plt.figure()
plt.semilogy(f, psd, label='numpy.fft()')
plt.semilogy(fxx, Pxx_den, label='scipy.signal.welch()')
plt.semilogy(f, [b0]*len(f), label='b_0 = %.3g' % b0)
plt.legend(framealpha=0.5)
plt.xlabel('Frequeny / Hz')
plt.ylabel('one-sided PSD / V^2/Hz')
plt.show()
