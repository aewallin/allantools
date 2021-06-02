import allantools
import allantools.noise as noise
import numpy as np
import matplotlib.pyplot as plt

import math

# ieee1139-table for white FM
# produces synthetic dataset with given PSD and compares against
# the predicted S_y, S_fi, S_x, and ADEV given in the table
# AW 2015-07-29

# from the ieee1139 table
# PSD_y(f)    = h0 * f^0                   fractional frequency PSD
# PSD_fi(f)   = h0 * vo^2 * f^-2           phase (radians) PSD
# PSD_x(f)    = h0 * (2 pi f)^-2           phase (time) PSD
# ADEV_y(tau) = sqrt{ 1/2 * h0 * tau^-1 }  Allan deviation


fs = 100
h0 = 2e-16
N = 16*4096
v0 = 1e6  # nominal oscillator frequency

y = noise.white(num_points=N, b0=h0, fs=fs)  # fractional frequency
x = allantools.frequency2phase(y, fs)  # phase in seconds
fi = [2*math.pi*v0*xx for xx in x]  # phase in radians
t = np.linspace(0, (1.0/fs)*N, len(y))

plt.figure()
plt.plot(t, y)
plt.xlabel('Time / s')
plt.ylabel('Fractional frequency')
f_y, psd_y = noise.numpy_psd(y, fs)
f_fi, psd_fi = noise.numpy_psd(fi, fs)
f_x, psd_x = noise.numpy_psd(x, fs)

fxx, Pxx_den = noise.scipy_psd(y, fs)
f_fi2, psd_fi2 = noise.scipy_psd(fi, fs)
f_x2, psd_x2 = noise.scipy_psd(x, fs)

plt.figure()
plt.semilogy(f_y, psd_y, label='numpy.fft()')
plt.semilogy(fxx, Pxx_den, label='scipy.signal.welch()')
plt.semilogy(f_y, [h0]*len(f_y), label='h_0 = %.3g' % h0)
plt.legend(framealpha=0.5)
plt.title('PSD of fractional frequency')
plt.xlabel('Frequeny / Hz')
plt.ylabel('one-sided PSD / S_y(f)')

plt.figure()
plt.loglog(f_fi, psd_fi, label='numpy.fft()')
plt.loglog(f_fi2, psd_fi2, label='scipy.signal.welch()')
plt.loglog(f_fi[1:], [h0*v0*v0/(ff*ff)
                      for ff in f_fi[1:]], label='h_0 * v0^2 * f^-2')
plt.legend(framealpha=0.5)
plt.title('PSD of phase (radians)')
plt.xlabel('Frequeny / Hz')
plt.ylabel('one-sided PSD / S_fi(f)')

plt.figure()
plt.loglog(f_x, psd_x, label='numpy.fft()')
plt.loglog(f_x2, psd_x2, label='scipy.signal.welch()')
plt.loglog(f_x[1:], [h0/(2*math.pi*ff)**2 for ff in f_x[1:]],
           label='h0/(2*pi*f)**2')
plt.legend(framealpha=0.5)
plt.title('PSD of phase (time)')
plt.xlabel('Frequeny / Hz')
plt.ylabel('one-sided PSD / S_x(f)')

plt.figure()
taus = [tt for tt in np.logspace(-2.2, 4, 100)]
(taus_y, devs_y, errs_y, ns_y) = allantools.oadev(
    y, data_type='freq', rate=fs, taus=taus)
(taus_x, devs_x, errs_x, ns_x) = allantools.oadev(
    x, data_type='phase', rate=fs, taus=taus)
plt.loglog(taus_y, devs_y, 'o', label='ADEV from y')
plt.loglog(taus_x, devs_x, '*', label='ADEV from x')

adev_y = [math.sqrt(0.5*h0*(1.0/tt)) for tt in taus]
plt.loglog(taus, adev_y, label='sqrt( 0.5*h0*tau^-1 )')
plt.xlim((8e-3, 1e3))
plt.legend(framealpha=0.6)
plt.title('Allan deviation')
plt.xlabel('Tau / s')
plt.ylabel('Allan deviation')
plt.show()
