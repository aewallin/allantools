#!/usr/bin/env python

# import numpy as np
import matplotlib.pyplot as plt
import allantools as at
import math
# note that latex requires e.g. a texlive installation and dvipng
plt.rc('text', usetex=True)  # for latex
plt.rc('font', family='serif')

"""

colorednoise example with phase-PSD, frequency-PSD, and ADEV/MDEV
for different noise types.

Anders Wallin 2016-04-09

This example uses the Noise() class and Kasdin&Walter
updated 2018-01
"""


def main():

    nr = 2**14  # number of datapoints in time-series
    adev0 = 1.0e-11

    #
    # for each noise type generated Noise class
    # compute
    #   - phase PSD and coefficient g_b
    #   - frequency PSD and coefficient h_a
    #   - ADEV
    #   - MDEV
    #
    # discrete variance for noiseGen()
    qd0 = qd1 = qd2 = qd3 = qd4 = pow(adev0, 2)
    tau0 = 1.0  # sample interval
    sample_rate = 1.0/tau0

    x0 = at.Noise(nr, qd0, 0)   # white phase noise (WPM)
    x0.generateNoise()
    x1 = at.Noise(nr, qd1, -1)  # flicker phase noise (FPM)
    x1.generateNoise()
    x2 = at.Noise(nr, qd2, -2)  # white frequency noise (WFM)
    x2.generateNoise()
    x3 = at.Noise(nr, qd3, -3)  # flicker frequency noise (FFM)
    x3.generateNoise()
    x4 = at.Noise(nr, qd4, -4)  # random walk frequency noise (RWFM)
    x4.generateNoise()

    # compute frequency time-series
    y0 = at.phase2frequency(x0.time_series, sample_rate)
    y1 = at.phase2frequency(x1.time_series, sample_rate)
    y2 = at.phase2frequency(x2.time_series, sample_rate)
    y3 = at.phase2frequency(x3.time_series, sample_rate)
    y4 = at.phase2frequency(x4.time_series, sample_rate)

    # compute phase PSD
    (f0, psd0) = at.noise.scipy_psd(
        x0.time_series, f_sample=sample_rate, nr_segments=4)
    (f1, psd1) = at.noise.scipy_psd(
        x1.time_series, f_sample=sample_rate, nr_segments=4)
    (f2, psd2) = at.noise.scipy_psd(
        x2.time_series, f_sample=sample_rate, nr_segments=4)
    (f3, psd3) = at.noise.scipy_psd(
        x3.time_series, f_sample=sample_rate, nr_segments=4)
    (f4, psd4) = at.noise.scipy_psd(
        x4.time_series, f_sample=sample_rate, nr_segments=4)

    # compute phase PSD prefactor g_b
    g0 = x0.phase_psd_from_qd(tau0)
    g1 = x1.phase_psd_from_qd(tau0)
    g2 = x2.phase_psd_from_qd(tau0)
    g3 = x3.phase_psd_from_qd(tau0)
    g4 = x4.phase_psd_from_qd(tau0)

    # compute frequency PSD
    (ff0, fpsd0) = at.noise.scipy_psd(y0, f_sample=sample_rate, nr_segments=4)
    (ff1, fpsd1) = at.noise.scipy_psd(y1, f_sample=sample_rate, nr_segments=4)
    (ff2, fpsd2) = at.noise.scipy_psd(y2, f_sample=sample_rate, nr_segments=4)
    (ff3, fpsd3) = at.noise.scipy_psd(y3, f_sample=sample_rate, nr_segments=4)
    (ff4, fpsd4) = at.noise.scipy_psd(y4, f_sample=sample_rate, nr_segments=4)

    # compute frequency PSD prefactor h_a
    a0 = x0.frequency_psd_from_qd(tau0)
    a1 = x1.frequency_psd_from_qd(tau0)
    a2 = x2.frequency_psd_from_qd(tau0)
    a3 = x3.frequency_psd_from_qd(tau0)
    a4 = x4.frequency_psd_from_qd(tau0)

    # compute ADEV
    (t0, d0, e, n) = at.oadev(x0.time_series, rate=sample_rate)
    (t1, d1, e, n) = at.oadev(x1.time_series, rate=sample_rate)
    (t2, d2, e, n) = at.oadev(x2.time_series, rate=sample_rate)
    (t3, d3, e, n) = at.oadev(x3.time_series, rate=sample_rate)
    (t4, d4, e, n) = at.oadev(x4.time_series, rate=sample_rate)

    # compute MDEV
    (mt0, md0, e, n) = at.mdev(x0.time_series, rate=sample_rate)
    (mt1, md1, e, n) = at.mdev(x1.time_series, rate=sample_rate)
    (mt2, md2, e, n) = at.mdev(x2.time_series, rate=sample_rate)
    (mt3, md3, e, n) = at.mdev(x3.time_series, rate=sample_rate)
    (mt4, md4, e, n) = at.mdev(x4.time_series, rate=sample_rate)

    plt.figure()
    # Phase PSD figure
    plt.subplot(2, 2, 1)
    plt.loglog(f0, [g0*pow(xx, 0.0)
                    for xx in f0], '--', label=r'$g_0f^0$', color='black')
    plt.loglog(f0, psd0, '.', color='black')

    plt.loglog(f1[1:], [g1*pow(xx, -1.0) for xx in f1[1:]],
               '--', label=r'$g_{-1}f^{-1}$', color='red')
    plt.loglog(f1, psd1, '.', color='red')

    plt.loglog(f2[1:], [g2*pow(xx, -2.0) for xx in f2[1:]],
               '--', label=r'$g_{-2}f^{-2}$', color='green')
    plt.loglog(f2, psd2, '.',  color='green')

    plt.loglog(f3[1:], [g3*pow(xx, -3.0) for xx in f3[1:]],
               '--', label=r'$g_{-3}f^{-3}$', color='pink')
    plt.loglog(f3, psd3, '.', color='pink')

    plt.loglog(f4[1:], [g4*pow(xx, -4.0) for xx in f4[1:]],
               '--', label=r'$g_{-4}f^{-4}$', color='blue')
    plt.loglog(f4, psd4, '.', color='blue')
    plt.grid()
    plt.legend(framealpha=0.9)
    plt.title(r'Phase Power Spectral Density')
    plt.xlabel(r'Frequency (Hz)')
    plt.ylabel(r' $S_x(f)$ $(s^2/Hz)$')

    # frequency PSD figure
    plt.subplot(2, 2, 2)
    plt.loglog(ff0[1:], [a0*pow(xx, 2) for xx in ff0[1:]],
               '--', label=r'$h_{2}f^{2}$', color='black')
    plt.loglog(ff0, fpsd0, '.',  color='black')

    plt.loglog(ff1[1:], [a1*pow(xx, 1) for xx in ff1[1:]],
               '--', label=r'$h_{1}f^{1}$', color='red')
    plt.loglog(ff1, fpsd1, '.',  color='red')

    plt.loglog(ff2[1:], [a2*pow(xx, 0) for xx in ff2[1:]],
               '--', label=r'$h_{0}f^{0}$', color='green')
    plt.loglog(ff2, fpsd2, '.', color='green')

    plt.loglog(ff3[1:], [a3*pow(xx, -1) for xx in ff3[1:]],
               '--', label=r'$h_{-1}f^{-1}$', color='pink')
    plt.loglog(ff3, fpsd3, '.', color='pink')

    plt.loglog(ff4[1:], [a4*pow(xx, -2) for xx in ff4[1:]],
               '--', label=r'$h_{-2}f^{-2}$', color='blue')
    plt.loglog(ff4, fpsd4, '.', color='blue')
    plt.grid()
    plt.legend(framealpha=0.9)
    plt.title(r'Frequency Power Spectral Density')
    plt.xlabel(r'Frequency (Hz)')
    plt.ylabel(r' $S_y(f)$ $(1/Hz)$')

    # ADEV figure
    plt.subplot(2, 2, 3)

    plt.loglog(t0, [x0.adev_from_qd(tau0, xx)/xx for xx in t0],
               '--', label=r'$\propto\sqrt{h_{2}}\tau^{-1}$', color='black')
    plt.loglog(t0, d0, 'o', color='black')

    plt.loglog(t1, [x1.adev_from_qd(tau0, xx)/xx for xx in t1],
               '--', label=r'$\propto\sqrt{h_{1}}\tau^{-1}$', color='red')
    plt.loglog(t1, d1, 'o', color='red')

    plt.loglog(t2, [x2.adev_from_qd(tau0, xx)/math.sqrt(xx) for xx in t2],
               '--', label=r'$\propto\sqrt{h_{0}}\tau^{-1/2}$', color='green')
    plt.loglog(t2, d2, 'o', color='green')

    plt.loglog(t3, [x3.adev_from_qd(tau0, xx)*1 for xx in t3],
               '--', label=r'$\propto\sqrt{h_{-1}}\tau^0$', color='pink')
    plt.loglog(t3, d3, 'o', color='pink')

    plt.loglog(t4, [x4.adev_from_qd(tau0, xx)*math.sqrt(xx) for xx in t4],
               '--', label=r'$\propto\sqrt{h_{-2}}\tau^{+1/2}$', color='blue')
    plt.loglog(t4, d4, 'o', color='blue')

    plt.legend(framealpha=0.9, loc='lower left')
    plt.grid()
    plt.title(r'Allan Deviation')
    plt.xlabel(r'$\tau$ (s)')
    plt.ylabel(r'ADEV')

    # MDEV
    plt.subplot(2, 2, 4)

    plt.loglog(t0, [x0.mdev_from_qd(tau0, xx)/pow(xx, 3.0/2.0) for xx in t0],
               '--', label=r'$\propto\sqrt{h_{2}}\tau^{-3/2}$', color='black')
    plt.loglog(mt0, md0, 'o', color='black')

    plt.loglog(t1, [x1.mdev_from_qd(tau0, xx)/xx for xx in t1],
               '--', label=r'$\propto\sqrt{h_{1}}\tau^{-1}$', color='red')
    plt.loglog(mt1, md1, 'o', color='red')

    plt.loglog(t2, [x2.mdev_from_qd(tau0, xx)/math.sqrt(xx) for xx in t2],
               '--', label=r'$\propto\sqrt{h_{0}}\tau^{-1/2}$', color='green')
    plt.loglog(mt2, md2, 'o', color='green')

    plt.loglog(t3, [x3.mdev_from_qd(tau0, xx)**1 for xx in t3],
               '--', label=r'$\propto\sqrt{h_{-1}}\tau^0$', color='pink')
    plt.loglog(mt3, md3, 'o', color='pink')

    plt.loglog(t4, [x4.mdev_from_qd(tau0, xx)*math.sqrt(xx) for xx in t4],
               '--', label=r'$\propto\sqrt{h_{-2}}\tau^{+1/2}$', color='blue')
    plt.loglog(mt4, md4, 'o', color='blue')

    plt.legend(framealpha=0.9, loc='lower left')
    plt.grid()
    plt.title(r'Modified Allan Deviation')
    plt.xlabel(r'$\tau$ (s)')
    plt.ylabel(r'MDEV')
    plt.show()


if __name__ == "__main__":
    main()
