"""
Allan deviation tools, Functions for generating noise
=====================================================

See: http://en.wikipedia.org/wiki/Colors_of_noise

**Author:** Anders Wallin (anders.e.e.wallin "at" gmail.com)

Version history
---------------

**2019.07**
- Homogenized output types and parameter names
- PEP8 + docstrings update

**v1.00**
- initial version Anders Wallin, 2014 January


This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GGNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import math
import numpy as np
import scipy.signal  # for welch PSD


def numpy_psd(x, f_sample=1.0):
    """ calculate power spectral density of input signal x
        x = signal
        f_sample = sampling frequency in Hz. i.e. 1/fs is the time-interval
             in seconds between datapoints
        scale fft so that output corresponds to 1-sided PSD
        output has units of [X^2/Hz] where X is the unit of x
    """
    psd_of_x = ((2.0 / (float(len(x)) * f_sample))
                * np.abs(np.fft.rfft(x))**2)
    f_axis = np.linspace(0, f_sample/2.0, len(psd_of_x))  # frequency axis
    return f_axis, psd_of_x


def scipy_psd(x, f_sample=1.0, nr_segments=4):
    """ PSD routine from scipy
        we can compare our own numpy result against this one
    """
    f_axis, psd_of_x = scipy.signal.welch(x,
                                          f_sample,
                                          nperseg=len(x)/nr_segments)
    return f_axis, psd_of_x


def white(num_points=1024, b0=1.0, fs=1.0):
    """ White noise generator

        Generate time series with white noise that has constant PSD = b0,
        up to the nyquist frequency fs/2.

        The PSD is at 'height' b0 and extends from 0 Hz up to the nyquist
        frequency fs/2 (prefactor math.sqrt(b0*fs/2.0))

        Parameters
        ----------
        num_points: int, optional
            number of samples
        b0: float, optional
            desired power-spectral density in [X^2/Hz] where X is the unit of x
        fs: float, optional
            sampling frequency, i.e. 1/fs is the time-interval between
            datapoints

        Returns
        -------
        White noise sample: numpy.array
    """
    return math.sqrt(b0*fs/2.0)*np.random.randn(num_points)


def brown(num_points=1024, b_minus2=1.0, fs=1.0):
    """ Brownian or random walk (diffusion) noise with 1/f^2 PSD

        Not really a color... rather Brownian or random-walk.
        Obtained by integrating white-noise.

        Parameters
        ----------
        num_points: int, optional
            number of samples
        b_minus2: float, optional
            desired power-spectral density is b2*f^-2
        fs: float, optional
            sampling frequency, i.e. 1/fs is the time-interval between
            datapoints

        Returns
        -------
        Random walk sample: numpy.array
    """
    return (1.0/float(fs))*np.cumsum(
        white(num_points,
              b0=b_minus2*(4.0*math.pi*math.pi),
              fs=fs))


def violet(num_points=1024, b2=1, fs=1):
    """ Violet noise with f^2 PSD

        Obtained by differentiating white noise

        Parameters
        ----------
        num_points: int, optional
            number of samples
        b2: float, optional
            desired power-spectral density is b2*f^2
        fs: float, optional
            sampling frequency, i.e. 1/fs is the time-interval between
            datapoints

        Returns
        -------
        Violet noise sample: numpy.array
    """
    # diff() reduces number of points by one.
    return (float(fs))*np.diff(
        white(num_points+1, b0=b2/(2.0*math.pi)**2, fs=fs))


def pink(num_points=1024, depth=80):
    """ Pink noise (approximation) with 1/f PSD

        Fills a sample with results from a pink noise generator
        from http://pydoc.net/Python/lmj.sound/0.1.1/lmj.sound.noise/,
        based on the Voss-McCartney algorithm, discussion and code examples at
        http://www.firstpr.com.au/dsp/pink-noise/

        Parameters
        ----------
        num_points: int, optional
            number of samples
        depth: int, optional
            number of iteration for each point. High numbers are slower but
            generates a more correct spectrum on low-frequencies end.

        Returns
        -------
        Pink noise sample: numpy.array
    """
    # FIXME: couldn't we implement here the normalization as for the other
    # noise types using a noise power law coefficient b_minus1?
    a = []
    s = iterpink(depth)
    for n in range(num_points):  # FIXME: num_points is unused here.
        a.append(next(s))
    return np.array(a)


def iterpink(depth=20):
    """Generate a sequence of samples of pink noise.

    pink noise generator
    from http://pydoc.net/Python/lmj.sound/0.1.1/lmj.sound.noise/

    Based on the Voss-McCartney algorithm, discussion and code examples at
    http://www.firstpr.com.au/dsp/pink-noise/

    depth: Use this many samples of white noise to calculate the output. A
      higher  number is slower to run, but renders low frequencies with more
      correct power spectra.

    Generates a never-ending sequence of floating-point values. Any continuous
    set of these samples will tend to have a 1/f power spectrum.
    """
    values = np.random.randn(depth)
    smooth = np.random.randn(depth)
    source = np.random.randn(depth)
    sumvals = values.sum()
    i = 0
    while True:
        yield sumvals + smooth[i]

        # advance the index by 1. if the index wraps, generate noise to use in
        # the calculations, but do not update any of the pink noise values.
        i += 1
        if i == depth:
            i = 0
            smooth = np.random.randn(depth)
            source = np.random.randn(depth)
            continue

        # count trailing zeros in i
        c = 0
        while not (i >> c) & 1:
            c += 1

        # replace value c with a new source element
        sumvals += source[i] - values[c]
        values[c] = source[i]


def timmer_koenig_from_psd(f_nodes, h, alpha, duration, timestep, output='phase', seed=None):
    """
    Generate time-domain clock noise from a piecewise power-law one-sided PSD :math:`S_y(f)`
    using a Timmer–Koenig-style Fourier synthesis.

     Frequency grid (old convention used in your validated script):
        n  = int(duration/dt)
        f1 = 1 / ((n - 1) * dt)
        fn = 1 / (2 * dt)
        f_k = linspace(f1, fn, n/2 + 1)

    Sample one-sided fractional-frequency PSD:
        :math:`S_y(f) = h_i f^{\\alpha_i}`

    Convert to phase/time PSD:
        :math:`S_x(f) = S_y(f) / (2 \\pi f)^2`

    Generate complex spectrum coefficients (Timmer–Koenig style):
        :math:`X_k = \\sqrt{S_x(f_k)}/2 * (N(0,1) + i N(0,1))`

    Impose Hermitian symmetry and inverse FFT:
        x = ifft(X) * sqrt((n-1)/dt)

    Parameters
    ----------
    f_nodes : array_like
        Break frequencies (Hz) defining PSD intervals, ascending.
        Length is typically nseg-1.
    h : array_like
        PSD coefficients for each interval, length nseg.
        Defines Sy(f) = h_i * f**alpha_i in each interval.
    alpha : array_like
        PSD slopes for each interval, length nseg.
    duration : float
        Total duration of generated time series (seconds).
    timestep : float
        Sampling interval dt (seconds). Sampling rate is fs = 1/dt.
    output : {'phase', 'freq'}, optional
        'phase' returns phase/time fluctuation x(t) [s].
        'freq' returns fractional frequency y(t) = dx/dt [-].
    seed : int or None, optional
        Random seed for reproducible synthesis.

    Returns
    -------
    x_or_y : ndarray
        Time series of length n. If output='phase', returns x(t) [s].
        If output='freq', returns y(t) [-].

    References
    ----------
    J. Timmer and M. Koenig, "On generating power law noise",
    Astronomy and Astrophysics, vol. 300, pp. 707–710, 1995.
    """
    f_nodes = np.asarray(f_nodes, dtype=float)
    h = np.asarray(h, dtype=float)
    alpha = np.asarray(alpha, dtype=float)

    if duration <= 0 or timestep <= 0:
        raise ValueError("duration and timestep must be positive.")
    if h.size != alpha.size:
        raise ValueError("h and alpha must have the same length.")
    if f_nodes.size not in (0, h.size - 1):
        # Typical output of your reconstruction: len(f_nodes)=len(h)-1
        raise ValueError("Expected len(f_nodes) == len(h)-1 (or 0 for single segment).")

    # Seed behavior: match numpy global RNG behavior used in your original code
    if seed is not None:
        np.random.seed(int(seed))

    # Extrapolate behavior beyond the last node (matches your original)
    h_ext = np.concatenate([h, [h[-1]]])
    alpha_ext = np.concatenate([alpha, [alpha[-1]]])

    n = int(duration / timestep)
    if n < 4:
        raise ValueError("duration/timestep too small; need at least 4 samples.")

    f1 = 1.0 / ((n - 1) * timestep)
    fn = 1.0 / (2.0 * timestep)
    frequencies = np.linspace(f1, fn, n // 2 + 1)

    # Sample Sy on this grid using the same interval logic as your original
    Sy = np.zeros_like(frequencies)
    for i, f in enumerate(frequencies):
        if f_nodes.size > 0 and f < f_nodes[0]:
            Sy[i] = h_ext[0] * f ** alpha_ext[0]
        elif f_nodes.size > 0 and f > f_nodes[-1]:
            Sy[i] = h_ext[-1] * f ** alpha_ext[-1]
        else:
            if f_nodes.size == 0:
                Sy[i] = h_ext[0] * f ** alpha_ext[0]
            else:
                for j in range(1, f_nodes.size):
                    if f_nodes[j - 1] <= f < f_nodes[j]:
                        Sy[i] = h_ext[j] * f ** alpha_ext[j]
                        break

    # Convert Sy -> Sx
    Sx = Sy / (2.0 * np.pi * frequencies) ** 2
    Sx[0] = 0.0

    # Build full complex spectrum and enforce Hermitian symmetry
    X = np.zeros(n, dtype=complex)
    for j in range(1, n // 2 + 1):
        X[j] = (np.sqrt(Sx[j]) / 2.0) * (np.random.normal(0, 1) + 1j * np.random.normal(0, 1))

    X[n // 2 + 1:] = np.conj(X[1:n // 2][::-1])

    x = np.fft.ifft(X) * np.sqrt((n - 1) / timestep)
    x = x.real

    if output == 'phase':
        return x
    elif output == 'freq':
        return np.concatenate(([0.0], np.diff(x) / timestep))
    else:
        raise ValueError("output must be 'phase' or 'freq'.")
        
