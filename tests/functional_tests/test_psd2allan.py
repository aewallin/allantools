#!/usr/bin/python

import sys
sys.path.append("..")

from allantools import noise
import numpy as np
import matplotlib
matplotlib.use('Agg') # automated tests run without a display
import matplotlib.pyplot as plt

import allantools as at
from  scipy.signal import welch

def test_psd2allan_figure():
    f= np.arange(1e4+1) # generate f-vector 0...10^4 Hz in 1 Hz steps -> Nyquist freq is 5 kHz
    S_y0=1e-24 # some arbitrarily chosen noise level
    S_y_WFM= S_y0*np.ones(np.size(f)) # generate white frequency noise S_y(f)
    S_y_WPM= S_y_WFM*f**2 # white phase noise S_y(f)

    plt.rc('text', usetex=True)

    plt.close('all')
    plt.figure(1)
    plt.loglog(f[f>0],S_y_WFM[f>0])
    plt.loglog(f[f>0],S_y_WPM[f>0])

    y_WFM_ind= noise.white(num_points= int(2e4), b0=S_y0, fs=2e4)
    y_WPM_ind= noise.violet(num_points= int(2e4), b2=S_y0, fs=2e4)
    f, S_y_WFM_ind= welch(y_WFM_ind,fs=2e4, nperseg=y_WFM_ind.size, window='hamming')
    f, S_y_WPM_ind= welch(y_WPM_ind,fs=2e4, nperseg=y_WPM_ind.size, window='hamming')

    plt.loglog(f,S_y_WFM_ind)
    plt.loglog(f,S_y_WPM_ind)
    plt.xlabel('$f$ [Hz]')
    plt.ylabel('$S_y$')
    plt.legend(('WFM direct', 'WPM direct', 'WFM indirect', 'WPM indirect'))

    #(tau, sigma_WFM)= at.psd2allan(S_y_WFM, f, kind= 'a')
    #(tau, sigma_WPM)= at.psd2allan(S_y_WPM, f, kind= 'a')
    #(tau, modsigma_WFM)= at.psd2allan(S_y_WFM, f, kind= 'm')
    #(tau, modsigma_WPM)= at.psd2allan(S_y_WPM, f, kind= 'm')
    (tau, sigma_WFM)= at.psd2allan(S_y_WFM, kind= 'a')
    (tau, sigma_WPM)= at.psd2allan(S_y_WPM, kind= 'a')
    (tau, modsigma_WFM)= at.psd2allan(S_y_WFM, kind= 'm')
    (tau, modsigma_WPM)= at.psd2allan(S_y_WPM, kind= 'm')

    plt.figure(2)
    plt.loglog(tau, sigma_WFM)
    plt.loglog(tau, sigma_WPM)
    plt.loglog(tau, modsigma_WFM)
    plt.loglog(tau, modsigma_WPM)

    (tau, sigma_WFM_ind)= at.psd2allan(S_y_WFM_ind, kind= 'a')
    (tau, sigma_WPM_ind)= at.psd2allan(S_y_WPM_ind, kind= 'a')
    (tau, modsigma_WFM_ind)= at.psd2allan(S_y_WFM_ind, kind= 'm')
    (tau, modsigma_WPM_ind)= at.psd2allan(S_y_WPM_ind, kind= 'm')

    plt.loglog(tau, sigma_WFM_ind, ':C0')
    plt.loglog(tau, sigma_WPM_ind, ':C1')
    plt.loglog(tau, modsigma_WFM_ind, ':C2')
    plt.loglog(tau, modsigma_WPM_ind, ':C3')

    f_h= 1e4
    sigma_WFM_theo= np.sqrt(S_y0 / tau / 2.0)
    sigma_WPM_theo= np.sqrt(3.0 * f_h * S_y0 / tau**2 / (2.0 * np.pi)**2)
    modsigma_WFM_theo= np.sqrt(S_y0 / tau / 4.0)
    modsigma_WPM_theo= np.sqrt(3.0 * S_y0 / tau**3 / (2.0 * np.pi)**2 / 2)

    plt.loglog(tau, sigma_WFM_theo, '.C0')
    plt.loglog(tau, sigma_WPM_theo, '.C1')
    plt.loglog(tau, modsigma_WFM_theo, '.C2')
    plt.loglog(tau, modsigma_WPM_theo, '.C3')
    plt.xlabel(r'$\tau$ [s]')
    plt.ylabel(r'$\sigma_y$')
    plt.legend(('psd2allan direct, ADEV, WFM', 'psd2allan direct, ADEV, WPM',
                'psd2allan direct, modADEV, WFM', 'psd2allan direct, modADEV, WPM',
                'psd2allan indirect, ADEV, WFM', 'psd2allan  indirect, ADEV, WPM',
                'psd2allan  indirect, modADEV, WFM', 'psd2allan  indirect, modADEV, WPM',
                'theoretical, ADEV, WFM', 'theoretical, ADEV, WPM',
                'theoretical, modADEV, WFM', 'theoretical, modADEV, WPM'))

    # Resume:
    # psd2allan works as expected, but artefacts due to the limited resolution of
    # FFT can easily lead to artefacts at large tau if the input S_y has been
    # computed using FFT (e.g. welch, periodogram, ...)

if __name__=="__main__":
    test_psd2allan_figure()
