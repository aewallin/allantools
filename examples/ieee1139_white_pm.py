import allantools
import allantools.noise as noise
import numpy as np
import matplotlib.pyplot as plt

import math

# ieee1139-table for white PM
# produces synthetic dataset with given PSD and compares against
# the predicted S_y, S_fi, S_x, and ADEV given in the table
# AW 2015-07-29

# from the ieee1139 table
# PSD_y(f)    = h2 * f^2                   fractional frequency PSD
# PSD_fi(f)   = h2 * v0^2 * f^0            phase (radians) PSD
# PSD_x(f)    = h2 * (2 pi )^-2            phase (time) PSD
# ADEV_y(tau) = sqrt{ 3*fh / (4*pi^2) * h2 * tau^-2 }  Allan deviation
# fh is the upper limit for the noise process, otherwise we would have infinite power...

fs=100
h2=2e-16
N=16*4096
v0 = 1e6 # nominal oscillator frequency

fi = noise.white(N=N,b0=h2*v0*v0,fs=fs)  # phase in radians
x = [fifi/(2*math.pi*v0) for fifi in fi] # phase in seconds
y = (fs)*np.diff(x)                      # fractional frequency
t = np.linspace(0, (1.0/fs)*N, len(y))

plt.figure()
plt.plot(t,y)
plt.xlabel('Time / s')
plt.ylabel('Fractional frequency')
plt.title('Fractional frequency time-series')

f_y, psd_y = noise.numpy_psd(y,fs)
f_fi, psd_fi = noise.numpy_psd(fi,fs)
f_x, psd_x = noise.numpy_psd(x,fs)

fxx, Pxx_den = noise.scipy_psd(y, fs)
f_fi2, psd_fi2 = noise.scipy_psd(fi, fs)
f_x2, psd_x2 = noise.scipy_psd(x, fs)

plt.figure()
plt.loglog(f_y,psd_y,label='numpy.fft()')
plt.loglog(fxx,Pxx_den,label='scipy.signal.welch()')
plt.loglog(f_y,[h2*ff*ff for ff in f_y],label='h2 * f^2 ' )
plt.legend(framealpha=0.5)
plt.title('PSD of fractional frequency')
plt.xlabel('Frequeny / Hz')
plt.ylabel('one-sided PSD / S_y(f)')

plt.figure()
plt.loglog(f_fi,psd_fi,label='numpy.fft()')
plt.loglog(f_fi2,psd_fi2,label='scipy.signal.welch()')
plt.loglog(f_fi,[h2*v0*v0]*len(f_fi),label='h2 * v0^2'  )
plt.legend(framealpha=0.5)
plt.title('PSD of phase (radians)')
plt.xlabel('Frequeny / Hz')
plt.ylabel('one-sided PSD / S_fi(f)')

plt.figure()
plt.loglog(f_x,psd_x,label='numpy.fft()')
plt.loglog(f_x2,psd_x2,label='scipy.signal.welch()')
plt.loglog(f_x,[h2/(2*math.pi)**2 ]*len(f_x),label='h2/(4*pi^2)'  )
plt.legend(framealpha=0.5)
plt.title('PSD of phase (time)')
plt.xlabel('Frequeny / Hz')
plt.ylabel('one-sided PSD / S_x(f)')

plt.figure()
taus=np.logspace(-2.2,4,100)
(taus_y, devs_y, errs_y, ns_y) = allantools.oadev(y, fs, taus)
(taus_x, devs_x, errs_x, ns_x) = allantools.oadev_phase(x, fs, taus)
plt.loglog(taus_y,devs_y,'o',label='ADEV from y')
plt.loglog(taus_x,devs_x,'*',label='ADEV from x')


fh = 0.5*fs # this restricts the noise power to below fh
adev_y = [math.sqrt( 3*fh/(4*math.pi**2) * h2*(1.0/tt**2) ) for tt in taus]
plt.loglog(taus,adev_y,label='sqrt( 3*fh/(4*pi^2) * h2*tau^-2 )')
plt.xlim((8e-3,1e3))
plt.legend(framealpha=0.6)
plt.title('Allan deviation')
plt.xlabel('Tau / s')
plt.ylabel('Allan deviation')
plt.show()
