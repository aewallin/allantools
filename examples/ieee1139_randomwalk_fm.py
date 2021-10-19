import allantools
import allantools.noise as noise
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

import math

# ieee1139-table for random walk (Brownian) FM
# produces synthetic dataset with given PSD and compares against
# the predicted S_y, S_fi, S_x, and ADEV given in the table
# AW 2015-08-06

# from the ieee1139 table
# PSD_y(f)    = h2 * f^-2                     fractional frequency PSD
# PSD_fi(f)   = h2 * vo^2 * f^-4              phase (radians) PSD
# PSD_x(f)    = h2 * (2 pi)^-2 * f^-4         phase (time) PSD
# ADEV_y(tau) = sqrt{ 2*pi^2/3 * h2 * tau }   Allan deviation


# sampling frequency in Hz (code should work for any non-zero value here)
fs = 12.8
h2 = 2e-20  # PSD f^-2 coefficient
N = 10*4096  # number of samples
v0 = 1.2345e6  # nominal oscillator frequency

y = noise.brown(num_points=N, b_minus2=h2, fs=fs)  # fractional frequency
x = allantools.frequency2phase(y, fs)  # phase in seconds
fi = [2*math.pi*v0*xx for xx in x]  # phase in radians
t = np.linspace(0, (1.0/fs)*N, len(y))  # time-series time axis

# time-series figure
#plt.figure()
fig, ax1 = plt.subplots()
ax1.plot(t, y, 'b', label='y')
ax1.set_ylabel('Fractional frequency y', color='b')
ax2 = ax1.twinx()
ax2.plot(t, x[1:], 'r', label='x')
ax2.set_ylabel('Phase x / s', color='r')
plt.legend()
plt.title('Time-series of y and x')
plt.xlabel('Time / s')
#plt.ylabel('Fractional frequency')

# note: calculating the PSD of an 1/f^4 signal using fft seems to be difficult
# the welch method returns a correct 1/f^4 shaped PSD, but fft often does not
# things that may help
# - using a longer time-series (increase N above)
# - detrend using signal.detrend()
# - calculate PSD for a short window of data, and average over windows
#   (this is done in the Welch method)
# - read the Welch code to see what is going on
f_y,  psd_y = noise.numpy_psd(y, fs)
f_fi, psd_fi = noise.numpy_psd(signal.detrend(fi[:int(len(fi)/20)]), fs)
f_x,  psd_x = noise.numpy_psd(x[:int(len(x)/20)], fs)

fxx, Pxx_den = noise.scipy_psd(y,  fs)
f_fi2, psd_fi2 = noise.scipy_psd(fi, fs)
f_x2, psd_x2 = noise.scipy_psd(x,  fs)


# Fractional frequency PSD
plt.figure()
plt.loglog(f_y, psd_y, label='numpy.fft()')
plt.loglog(fxx, Pxx_den, label='scipy.signal.welch()')
plt.loglog(f_y[1:], [h2/(ff*ff) for ff in f_y[1:]],
           label='h_2/f^2 with h_2 = %.3g' % h2)
plt.legend(framealpha=0.5)
plt.title('PSD of fractional frequency')
plt.grid()
plt.xlabel('Frequeny / Hz')
plt.ylabel('one-sided PSD / S_y(f)')

# Phase (radians) PSD
plt.figure()
plt.loglog(f_fi, psd_fi, label='numpy.fft()')
plt.loglog(f_fi2, psd_fi2, label='scipy.signal.welch()')
plt.loglog(f_fi[1:], [h2*v0*v0/(ff**4.0)
                      for ff in f_fi[1:]], label='h_2 * v0^2 * f^-4')
plt.legend(framealpha=0.5)
plt.title('PSD of phase (radians)')
plt.xlabel('Frequeny / Hz')
plt.ylabel('one-sided PSD / S_fi(f)')
plt.grid()


# Phase (time) PSD
plt.figure()
plt.loglog(f_x, psd_x, label='numpy.fft()')
plt.loglog(f_x2, psd_x2, label='scipy.signal.welch()')
plt.loglog(f_x[1:], [h2/((2*math.pi)**2 * ff**4.0)
                     for ff in f_x[1:]], label='h2 * (2 pi)^-2 * f^-4')
plt.legend(framealpha=0.5)
plt.title('PSD of phase (time)')
plt.xlabel('Frequeny / Hz')
plt.ylabel('one-sided PSD / S_x(f)')
plt.grid()

plt.figure()
taus = [tt for tt in np.logspace(-2.2, 4, 100)]
(taus_y, devs_y, errs_y, ns_y) = allantools.oadev(
    y, rate=fs, data_type='freq',  taus=taus)
(taus_x, devs_x, errs_x, ns_x) = allantools.oadev(x, rate=fs, taus=taus)
plt.loglog(taus_y, devs_y, 'o', label='ADEV from y')
plt.loglog(taus_x, devs_x, '*', label='ADEV from x')

# sqrt{ 2*pi^2/3 * h2 * tau }
adev_y = [math.sqrt(((2*math.pi**2)/3)*h2*tt) for tt in taus]
plt.loglog(taus, adev_y, label='sqrt{ 2*pi^2/3 * h2 * tau }')
# plt.xlim((8e-3,1e3))
plt.legend(framealpha=0.6)
plt.title('Allan deviation')
plt.xlabel('Tau / s')
plt.ylabel('Allan deviation')
plt.grid()

plt.show()
