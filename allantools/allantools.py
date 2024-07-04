"""
Allan deviation tools
=====================

**Author:** Anders Wallin (anders.e.e.wallin "at" gmail.com)

Version history
---------------

**2024.06**
- fix import scipy.integrate.simpson

**2024.04**
- ITU PRC, PRTC, ePRTC masks for TDEV and MTIE in new file mask.py
- psd2allan() - convert PSD to ADEV/MDEV
- GCODEV, Grosslambert deviation, an improved three-cornered-hat analysis
- PDEV, Parabolic deviation


**2019.09** 2019 September
- packaging changes, for conda package
  (see https://anaconda.org/conda-forge/allantools)

**2019.07** 2019 August 3
- move edf-functions and noise-ID functions to ci.py
- mtotdev() htotdev() speed improvements
- save dataset results to text file
- real-time adev/mdev/hdev, in new file realtime.py
- travis testing on Linux, OSX, and Windows

**2018.03** 2018 March 27
- change license to LGPL v3 or later
- lag-1 autocorrelation noise identification function
- B1 noise identification
- R(n) noise identification
- Noise() class using Kasdin & Walter algorithm
- work on Greenhall's EDF and confidence intervals
- tests for confidence intervals

**2016.11** 2016 November 18
- Dataset class
- plotting with a Plot class
- confidence intervals based on Greenhall's EDF algorithm
- testing on multiple python versions with tox
- continuous integration with https://travis-ci.org/aewallin/allantools
- test coverage report on
  https://coveralls.io/github/aewallin/allantools?branch=master

**2016.4** 2016 April 8
- convert tests to use pytest
- split tests into individual pytests, make them all pass
- accept a numpy.array as taus parameter.
- Switch to new signature (https://github.com/aewallin/allantools/issues/29)
- remove old and slow pure-python implementation

**2016.3** 2016 March
- improve documentation and add __version__
- added Theo1 deviation theo1()
- added Hadamard Total Deviatio htotdev()
- added Modified Total Deviation mtotdev(), and Time Total Deviation ttotdev()
  http://www.anderswallin.net/2016/03/modified-total-deviation-in-allantools/
- automatic tau-lists:  taus=[ "all" | "octave" | "decade" ]
- merge adev() and adev_phase() into one, requiring phase= or
  frequency= argument
- add GPS dataset as example and test

**2016.2** 2016 February
- update release on PyPi https://pypi.python.org/pypi/AllanTools
- pytest and coverage
- setuptools
- change version number to year.month

**v1.2.1** 2015 July
- Python 3 compatibility using 2to3 tool, by kuzavas
- IPython notebook examples
- sphinx documentation, auto-built on readthedocs

**v1.2** 2014 November, Cantwell G. Carson conrtibuted:
- A gap-robust version of ADEV based on the paper by Sesia et al.
   gradev_phase() and gradev()
- Improved uncertainty estimates: uncertainty_estimate()
  This introduces a new dependency: scipy.stats.chi2()

**v1.1** 2014 August
- Danny Price converted the library to use numpy.
- many functions in allantools are now 100x faster than before.
- see http://www.anderswallin.net/2014/08/faster-allantools-with-numpy/

**v1.01** 2014 August
- PEP8 compliance improvements by Danny Price.

**v1.00** 2014 January, first version of allantools.
- see http://www.anderswallin.net/2014/01/allantools/

License
-------

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""

import os
import json
import numpy as np
from scipy import interpolate      # used in psd2allan()
from scipy.integrate import simpson  # used in psd2allan()

from . import ci  # edf, confidence intervals

# Get version number from json metadata
pkginfo_path = os.path.join(os.path.dirname(__file__),
                            'allantools_info.json')
with open(pkginfo_path) as fp:
    pkginfo = json.load(fp)
__version__ = pkginfo["version"]


def tdev(data, rate=1.0, data_type="phase", taus=None):
    """ Time deviation.

    Based on modified Allan variance.

    .. math::

        \\sigma^2_{TDEV}( \\tau ) = { \\tau^2 \\over 3 }
        \\sigma^2_{MDEV}( \\tau )

    Note that TDEV has a unit of seconds.

    Parameters
    ----------
    data: np.array
        Input data. Provide either phase or frequency (fractional,
        adimensional).
    rate: float
        The sampling rate for data, in Hz. Defaults to 1.0
    data_type: {'phase', 'freq'}
        Data type, i.e. phase or frequency. Defaults to "phase".
    taus: np.array
        Array of tau values, in seconds, for which to compute statistic.
        Optionally set taus=["all"|"octave"|"decade"] for automatic
        tau-list generation.

    Returns
    -------
    (taus, tdev, tdev_error, ns): tuple
          Tuple of values
    taus: np.array
        Tau values for which td computed
    tdev: np.array
        Computed time deviations (in seconds) for each tau value
    tdev_errors: np.array
        Time deviation errors
    ns: np.array
        Values of N used in mdev_phase()

    References
    ----------
    * http://en.wikipedia.org/wiki/Time_deviation
    * NIST [SP1065]_ eqn (15), page 18.
    """
    phase = input_to_phase(data, rate, data_type)
    (taus, md, mde, ns) = mdev(phase, rate=rate, taus=taus)
    td = taus * md / np.sqrt(3.0)
    tde = td / np.sqrt(ns)
    return taus, td, tde, ns


def mdev(data, rate=1.0, data_type="phase", taus=None):
    """  Modified Allan deviation.
    Used to distinguish between White and Flicker Phase Modulation.

    .. math::

        \\sigma^2_{MDEV}(m\\tau_0) = { 1 \\over 2 (m \\tau_0 )^2 (N-3m+1) }
        \\sum_{j=1}^{N-3m+1} \\left[
        \\sum_{i=j}^{j+m-1} {x}_{i+2m} - 2x_{i+m} + x_{i} \\right]^2

    Parameters
    ----------
    data: np.array
        Input data. Provide either phase or frequency (fractional,
        adimensional).
    rate: float
        The sampling rate for data, in Hz. Defaults to 1.0
    data_type: {'phase', 'freq'}
        Data type, i.e. phase or frequency. Defaults to "phase".
    taus: np.array
        Array of tau values, in seconds, for which to compute statistic.
        Optionally set taus=["all"|"octave"|"decade"] for automatic
        tau-list generation.

    Returns
    -------
    (taus2, md, mde, ns): tuple
          Tuple of values
    taus2: np.array
        Tau values for which td computed
    md: np.array
        Computed mdev for each tau value
    mde: np.array
        mdev errors
    ns: np.array
        Values of N used in each mdev calculation

    References
    ----------
    * NIST [SP1065]_ eqn (14), page 17.
    * http://www.leapsecond.com/tools/adev_lib.c
    """
    phase = input_to_phase(data, rate, data_type)
    (phase, ms, taus_used) = tau_generator(phase, rate, taus=taus)
    data, taus = np.array(phase), np.array(taus)

    md = np.zeros_like(ms)
    mderr = np.zeros_like(ms)
    ns = np.zeros_like(ms)

    # this is a 'loop-unrolled' algorithm following
    # http://www.leapsecond.com/tools/adev_lib.c
    for idx, m in enumerate(ms):
        m = int(m)  # without this we get: VisibleDeprecationWarning:
        # using a non-integer number instead of an integer
        # will result in an error in the future
        tau = taus_used[idx]

        # First loop sum
        d0 = phase[0:m]
        d1 = phase[m:2*m]
        d2 = phase[2*m:3*m]
        e = min(len(d0), len(d1), len(d2))
        v = np.sum(d2[:e] - 2*d1[:e] + d0[:e])
        s = v * v

        # Second part of sum
        d3 = phase[3*m:]
        d2 = phase[2*m:]
        d1 = phase[1*m:]
        d0 = phase[0:]

        e = min(len(d0), len(d1), len(d2), len(d3))
        n = e + 1

        v_arr = v + np.cumsum(d3[:e] - 3 * d2[:e] + 3 * d1[:e] - d0[:e])

        s = s + np.sum(v_arr * v_arr)
        s /= 2.0 * m * m * tau * tau * n
        s = np.sqrt(s)

        md[idx] = s
        mderr[idx] = (s / np.sqrt(n))
        ns[idx] = n

    return remove_small_ns(taus_used, md, mderr, ns)


def adev(data, rate=1.0, data_type="phase", taus=None):
    """ Allan deviation.
        Classic - use only if required - relatively poor confidence.

    .. math::

        \\sigma^2_{ADEV}(\\tau) = { 1 \\over 2 \\tau^2 }
        \\langle ( {x}_{n+2} - 2x_{n+1} + x_{n} )^2 \\rangle
        = { 1 \\over 2 (N-2) \\tau^2 }
        \\sum_{n=1}^{N-2} ( {x}_{n+2} - 2x_{n+1} + x_{n} )^2

    where :math:`x_n` is the time-series of phase observations, spaced
    by the measurement interval :math:`\\tau`, and with length :math:`N`.

    Or alternatively calculated from a time-series of fractional frequency:

    .. math::

        \\sigma^{2}_{ADEV}(\\tau) =  { 1 \\over 2 }
        \\langle ( \\bar{y}_{n+1} - \\bar{y}_n )^2 \\rangle

    where :math:`\\bar{y}_n` is the time-series of fractional frequency
    at averaging time :math:`\\tau`


    Parameters
    ----------
    data: np.array
        Input data. Provide either phase or frequency (fractional,
        adimensional).
    rate: float
        The sampling rate for data, in Hz. Defaults to 1.0
    data_type: {'phase', 'freq'}
        Data type, i.e. phase or frequency. Defaults to "phase".
    taus: np.array
        Array of tau values, in seconds, for which to compute statistic.
        Optionally set taus=["all"|"octave"|"decade"] for automatic
        tau-list generation.

    Returns
    -------
    (taus2, ad, ade, ns): tuple
          Tuple of values
    taus2: np.array
        Tau values for which td computed
    ad: np.array
        Computed adev for each tau value
    ade: np.array
        adev errors
    ns: np.array
        Values of N used in each adev calculation

    References
    ----------
    * NIST [SP1065]_ eqn (6) and (7), pages 14 and 15.
    * [wikipedia_adev]_
    """
    phase = input_to_phase(data, rate, data_type)
    (phase, m, taus_used) = tau_generator(phase, rate, taus)

    ad = np.zeros_like(taus_used)
    ade = np.zeros_like(taus_used)
    adn = np.zeros_like(taus_used)

    for idx, mj in enumerate(m):  # loop through each tau value m(j)
        (ad[idx], ade[idx], adn[idx]) = calc_adev_phase(phase, rate, mj, mj)

    return remove_small_ns(taus_used, ad, ade, adn)


def calc_adev_phase(phase, rate, mj, stride):
    """  Main algorithm for adev() (stride=mj) and oadev() (stride=1)

    Parameters
    ----------
    phase: np.array
        Phase data in seconds.
    rate: float
        The sampling rate for phase or frequency, in Hz
    mj: int
        averaging factor, we evaluate at tau = m*tau0
    stride: int
        Size of stride

    Returns
    -------
    (dev, deverr, n): tuple
        Array of computed values.

    Notes
    -----
    stride = mj for nonoverlapping Allan deviation
    stride = 1 for overlapping Allan deviation

    References
    ----------
    * http://en.wikipedia.org/wiki/Allan_variance
    * http://www.leapsecond.com/tools/adev_lib.c
    * NIST [SP1065]_ eqn (7) and (11) page 16
    """
    mj = int(mj)
    stride = int(stride)
    d2 = phase[2 * mj::stride]
    d1 = phase[1 * mj::stride]
    d0 = phase[::stride]

    n = min(len(d0), len(d1), len(d2))

    if n == 0:
        RuntimeWarning("Data array length is too small: %i" % len(phase))
        n = 1

    v_arr = d2[:n] - 2 * d1[:n] + d0[:n]
    s = np.sum(v_arr * v_arr)

    dev = np.sqrt(s / (2.0*n)) / mj*rate
    deverr = dev / np.sqrt(n)

    return dev, deverr, n


def pdev(data, rate=1.0, data_type="phase", taus=None):
    """ Parabolic deviation.

    Use for evaluating uncertainty of omega-average of frequency.

    .. math::

        \\sigma^2_{PDEV}(m\\tau_0) = { 72 \\over (N-2m) (m \\tau_0 )^2  }
        \\sum_{i=0}^{N-2m-1} \\left[
        \\sum_{k=0}^{m-1} \\left( { m-1 \\over 2} - k \\right) {x}_{i+k} - x_{i+k+m}  \\right]^2

    for :math:`m>1` and for an averaging-factor of :math:`m=1` PDEV equals ADEV/MDEV: :math:`\\sigma_{PDEV}(\\tau_0)=\\sigma_{ADEV}(\\tau_0)`.

    Parameters
    ----------
    data: np.array
        Input data. Provide either phase or frequency (fractional,
        adimensional).
    rate: float
        The sampling rate for data, in Hz. Defaults to 1.0
    data_type: {'phase', 'freq'}
        Data type, i.e. phase or frequency. Defaults to "phase".
    taus: np.array
        Array of tau values, in seconds, for which to compute statistic.
        Optionally set taus=["all"|"octave"|"decade"] for automatic
        tau-list generation.

    Returns
    -------
    (taus2, ad, ade, ns): tuple
          Tuple of values
    taus2: np.array
        Tau values for which td computed
    ad: np.array
        Computed adev for each tau value
    ade: np.array
        adev errors
    ns: np.array
        Values of N used in each adev calculation

    References
    ----------
    * [Vernotte2020]_
    * [Vernotte2015]_
    """
    phase = input_to_phase(data, rate, data_type)
    (phase, m, taus_used) = tau_generator(phase, rate, taus)

    ad = np.zeros_like(taus_used)
    ade = np.zeros_like(taus_used)
    adn = np.zeros_like(taus_used)

    for idx, mj in enumerate(m):  # loop through each tau value m(j)
        (ad[idx], ade[idx], adn[idx]) = calc_pdev_phase(phase, rate, mj)

    return remove_small_ns(taus_used, ad, ade, adn)


def calc_pdev_phase(phase, rate, mj):
    """  Parabolic deviation

    Parameters
    ----------
    phase: np.array
        Phase data in seconds.
    rate: float
        The sampling rate for phase or frequency, in Hz
    mj: int
        M index value for stride

    Returns
    -------
    (dev, deverr, n): tuple
        Array of computed values.


    References
    ----------
    * [Vernotte2020]_
    * [Vernotte2015]_
    """
    mj = int(mj)
    stride = int(1)
    tau0 = 1.0/rate
    if mj == 1:  # same as OADEV
        d2 = phase[2 * mj::stride]
        d1 = phase[1 * mj::stride]
        d0 = phase[::stride]
        n = min(len(d0), len(d1), len(d2))

        if n == 0:
            RuntimeWarning("Data array length is too small: %i" % len(phase))
            n = 1

        v_arr = d2[:n] - 2 * d1[:n] + d0[:n]
        s = np.sum(v_arr * v_arr)

        dev = np.sqrt(s / (2.0*n)) / mj*rate
        deverr = dev / np.sqrt(n)
    else:
        N = len(phase)  # number of frequency samples
        M = N-2*mj  # Vernotte2020 has the correct(?) M = N - 2m
        # Vernotte2015 has M = N-2m+2 which seems wrong, we get index out-of-bounds in the sum

        if M < 1:
            return 0, 0, 0
        Msum = 0
        Mi = 0
        for i in range(0, M):  # 0..M-1
            asum = 0
            krange = np.linspace(0, mj-1, mj)
            asum = sum(((mj-1.0)/2.0 - krange) * (phase[i:i+mj] - phase[i+mj:i+2*mj]))
            Msum = Msum + pow(asum, 2)
            Mi = Mi + 1
        dev = np.sqrt(72*Msum / ((M)*pow(mj, 4)*pow(mj*tau0, 2)))
        deverr = dev / np.sqrt(M)
        n = M

    return dev, deverr, n


def calc_gcodev_phase(phase_1, phase_2, rate, mj, stride):
    """
    Main algorithm for the Groslambert codeviation (see arXiv:1904.05849)

    Parameters
    ----------
    phase_1 : np.array
        Phase data of oscillator 1.
    phase_2 : np.array
        Phase data of oscillator 2.
    mj: int
        M index value for stride.
    stride: int
        Size of stride.

    Returns
    -------
    (dev, deverr, n): tuple
        Array of computed values.

    stride = mj for nonoverlapping Allan deviation
    stride = 1 for overlapping Allan deviation (used for GCODEV by default)

    """

    mj = int(mj)
    stride = int(stride)
    d2_1 = phase_1[2 * mj::stride]
    d1_1 = phase_1[1 * mj::stride]
    d0_1 = phase_1[::stride]

    d2_2 = phase_2[2 * mj::stride]
    d1_2 = phase_2[1 * mj::stride]
    d0_2 = phase_2[::stride]

    n_1 = min(len(d0_1), len(d1_1), len(d2_1))
    n_2 = min(len(d0_2), len(d1_2), len(d2_2))

    if n_1 == 0:
        RuntimeWarning("Data array length is too small: %i" % len(phase_1))
        n_1 = 1
        n_2 = 1

    v_arr_1 = d2_1[:n_1] - 2 * d1_1[:n_1] + d0_1[:n_1]
    v_arr_2 = d2_2[:n_2] - 2 * d1_2[:n_2] + d0_2[:n_2]
    s = np.sum(v_arr_1 * v_arr_2)

    # Result can be negative
    if(s >= 0):
        dev = np.sqrt(s / (2.0*n_1)) / mj*rate
    else:
        dev = np.sqrt(np.abs(s) / (2.0*n_1)) / mj*rate

    deverr = dev / np.sqrt(n_1)

    return dev, deverr, n_1


def gcodev(data_1, data_2, rate=1.0, data_type="phase", taus=None):
    """ Groslambert codeviation  a.k.a. Allan Covariance

    Similarly to the three-cornered hat method, we consider three uncorrelated
    oscillators A, B, C. The Groslambert codeviation estimates the noise of
    one oscillator (e.g. B), given two synchronous measurements AB and BC.
    Unlike three-cordenred hat, Gcodev is not affected by the (uncorrelated) noise of
    the measurement devices (time-interval or frequency counter) used for
    the measurements AB and BC.

    Parameters
    ----------
    data_1: np.array
        Oscillator 1 input data. Provide either phase or frequency
    data_2: np.array
        Oscillator 2 input data. Provide either phase or frequency
    rate: float
        The sampling rate for data, in Hz. Defaults to 1.0
    data_type: {'phase', 'freq'}
        Data type, i.e. phase or frequency. Defaults to "phase".
    taus: np.array
        Array of tau values, in seconds, for which to compute statistic.
        Optionally set taus=["all"|"octave"|"decade"] for automatic
        tau-list generation.

    Returns
    -------
    (taus, gd): tuple
        Tuple of values
    taus: np.array
        Tau values for which gcodev computed
    gd: np.array
        Computed gcodev for each tau value

    References
    ----------
    * [Vernotte2016]_
    * [Lantz2019]_
    """
    phase_1 = input_to_phase(data_1, rate, data_type)
    phase_2 = input_to_phase(data_2, rate, data_type)
    (phase_1, m, taus_used) = tau_generator(phase_1, rate, taus)
    (phase_2, m, taus_used) = tau_generator(phase_2, rate, taus)

    gd = np.zeros_like(taus_used)
    gde = np.zeros_like(taus_used)
    gdn = np.zeros_like(taus_used)

    for idx, mj in enumerate(m):  # stride=1 for overlapping ADEV
        (gd[idx], gde[idx], gdn[idx]) = calc_gcodev_phase(phase_1,
                                                          phase_2,
                                                          rate,
                                                          mj,
                                                          stride=1)

    return remove_small_ns(taus_used, gd, gde, gdn)


def oadev(data, rate=1.0, data_type="phase", taus=None):
    """ Overlapping Allan deviation.
    General purpose - most widely used - first choice.

    .. math::

        \\sigma^2_{OADEV}(m\\tau_0) = { 1 \\over 2 (m \\tau_0 )^2 (N-2m) }
        \\sum_{n=1}^{N-2m} ( {x}_{n+2m} - 2x_{n+1m} + x_{n} )^2

    where :math:`\\sigma_{OADEV}(m\\tau_0)` is the overlapping Allan
    deviation at an averaging time of :math:`\\tau=m\\tau_0`, and
    :math:`x_n` is the time-series of phase observations, spaced by the
    measurement interval :math:`\\tau_0`, with length :math:`N`.

    Parameters
    ----------
    data: np.array
        Input data. Provide either phase or frequency (fractional,
        adimensional).
    rate: float
        The sampling rate for data, in Hz. Defaults to 1.0
    data_type: {'phase', 'freq'}
        Data type, i.e. phase or frequency. Defaults to "phase".
    taus: np.array
        Array of tau values, in seconds, for which to compute statistic.
        Optionally set taus=["all"|"octave"|"decade"] for automatic
        tau-list generation.

    Returns
    -------
    (taus2, ad, ade, ns): tuple
          Tuple of values
    taus2: np.array
        Tau values for which td computed
    ad: np.array
        Computed oadev for each tau value
    ade: np.array
        oadev errors
    ns: np.array
        Values of N used in each oadev calculation

    References
    ----------
    * NIST [SP1065]_ eqn (11), page 16.
    """
    phase = input_to_phase(data, rate, data_type)
    (phase, m, taus_used) = tau_generator(phase, rate, taus)
    ad = np.zeros_like(taus_used)
    ade = np.zeros_like(taus_used)
    adn = np.zeros_like(taus_used)

    for idx, mj in enumerate(m):  # stride=1 for overlapping ADEV
        (ad[idx], ade[idx], adn[idx]) = calc_adev_phase(phase, rate, mj, 1)

    return remove_small_ns(taus_used, ad, ade, adn)


def ohdev(data, rate=1.0, data_type="phase", taus=None):
    """ Overlapping Hadamard deviation.
        Better confidence than normal Hadamard.

    .. math::

        \\sigma^2_{OHDEV}(m\\tau_0) = { 1 \\over 6 (m \\tau_0 )^2 (N-3m) }
        \\sum_{i=1}^{N-3m} ( {x}_{i+3m} - 3x_{i+2m} + 3x_{i+m} - x_{i} )^2

    where :math:`x_i` is the time-series of phase observations, spaced
    by the measurement interval :math:`\\tau_0`, and with length :math:`N`.

    Parameters
    ----------
    data: np.array
        Input data. Provide either phase or frequency (fractional,
        adimensional).
    rate: float
        The sampling rate for data, in Hz. Defaults to 1.0
    data_type: {'phase', 'freq'}
        Data type, i.e. phase or frequency. Defaults to "phase".
    taus: np.array
        Array of tau values, in seconds, for which to compute statistic.
        Optionally set taus=["all"|"octave"|"decade"] for automatic
        tau-list generation.

    Returns
    -------
    (taus2, hd, hde, ns): tuple
          Tuple of values
    taus2: np.array
        Tau values for which td computed
    hd: np.array
        Computed hdev for each tau value
    hde: np.array
        hdev errors
    ns: np.array
        Values of N used in each hdev calculation

    References
    ----------
    * NIST [SP1065]_ eqn (20), page 21
    """
    phase = input_to_phase(data, rate, data_type)
    (phase, m, taus_used) = tau_generator(phase, rate, taus)
    hdevs = np.zeros_like(taus_used)
    hdeverrs = np.zeros_like(taus_used)
    ns = np.zeros_like(taus_used)

    for idx, mj in enumerate(m):
        (hdevs[idx],
         hdeverrs[idx],
         ns[idx]) = calc_hdev_phase(phase, rate, mj, 1)

    return remove_small_ns(taus_used, hdevs, hdeverrs, ns)


def hdev(data, rate=1.0, data_type="phase", taus=None):
    """ Hadamard deviation.
        Rejects frequency drift, and handles divergent noise.

    .. math::

        \\sigma^2_{HDEV}( \\tau ) = { 1 \\over 6 \\tau^2 (N-3) }
        \\sum_{i=1}^{N-3} ( {x}_{i+3} - 3x_{i+2} + 3x_{i+1} - x_{i} )^2

    where :math:`x_i` is the time-series of phase observations, spaced
    by the measurement interval :math:`\\tau`, and with length :math:`N`.

    Parameters
    ----------
    data : np.array
        Input data. Provide either phase or frequency (fractional,
        adimensional).
    rate : float
        The sampling rate for data, in Hz. Defaults to 1.0
    data_type : string, {'phase', 'freq'}
        Data type, i.e. phase or frequency. Defaults to "phase".
    taus : np.array
        Array of tau values, in seconds, for which to compute statistic.
        Optionally set taus=["all"|"octave"|"decade"] for automatic
        tau-list generation.

    References
    ----------
    * NIST [SP1065]_ eqn (17) and (18), page 20
    """
    phase = input_to_phase(data, rate, data_type)
    (phase, m, taus_used) = tau_generator(phase, rate, taus)
    hdevs = np.zeros_like(taus_used)
    hdeverrs = np.zeros_like(taus_used)
    ns = np.zeros_like(taus_used)

    for idx, mj in enumerate(m):
        (hdevs[idx],
         hdeverrs[idx],
         ns[idx]) = calc_hdev_phase(phase, rate, mj, mj)  # stride = mj

    return remove_small_ns(taus_used, hdevs, hdeverrs, ns)


def calc_hdev_phase(phase, rate, mj, stride):
    """ main calculation function for HDEV and OHDEV

                         1        N-3
         s2y(t) = --------------- sum [x(i+3) - 3x(i+2) + 3x(i+1) - x(i) ]^2
                  6*tau^2 (N-3m)  i=1

        N=M+1 phase measurements
        m is averaging factor


    Parameters
    ----------
    phase: np.array
        Phase data in seconds.
    rate: float
        The sampling rate for phase or frequency, in Hz
    mj: int
        M index value for stride
    stride: int
        Size of stride

    Returns
    -------
    (dev, deverr, n): tuple
        Array of computed values.

    References
    ----------
    * http://www.leapsecond.com/tools/adev_lib.c
    * NIST [SP1065]_ eqn (18) and (20) pages 20 and 21
    """
    tau0 = 1.0 / float(rate)
    mj = int(mj)
    stride = int(stride)
    d3 = phase[3 * mj::stride]
    d2 = phase[2 * mj::stride]
    d1 = phase[1 * mj::stride]
    d0 = phase[::stride]

    n = min(len(d0), len(d1), len(d2), len(d3))

    v_arr = d3[:n] - 3 * d2[:n] + 3 * d1[:n] - d0[:n]

    s = np.sum(v_arr * v_arr)

    if n == 0:
        n = 1

    h = np.sqrt(s / 6.0 / float(n)) / float(tau0 * mj)
    e = h / np.sqrt(n)
    return h, e, n


def totdev(data, rate=1.0, data_type="phase", taus=None):
    """ Total deviation.
        Better confidence at long averages for Allan deviation.

    .. math::

        \\sigma^2_{TOTDEV}( m\\tau_0 ) = { 1 \\over 2 (m\\tau_0)^2 (N-2) }
            \\sum_{i=2}^{N-1} ( {x}^*_{i-m} - 2x^*_{i} + x^*_{i+m} )^2


    Where :math:`x^*_i` is a new time-series of length :math:`3N-4`
    derived from the original phase time-series :math:`x_n` of
    length :math:`N` by reflection at both ends.
    The original data :math:`x_n` is in the center of :math:`x^*`:

    .. math::

        x^*_{1-j} = 2x_1 - x_{1+j}  \\quad \\text{for} \\quad  j=1..N-2

        x^*_i   = x_i        \\quad    \\text{for} \\quad  i=1..N

        x^*_{N+j} = 2x_N - x_{N-j}  \\quad \\text{for} \\quad  j=1..N-2

    FIXME: bias correction http://www.wriley.com/CI2.pdf page 5

    Parameters
    ----------
    phase: np.array
        Phase data in seconds. Provide either phase or frequency.
    frequency: np.array
        Fractional frequency data (nondimensional). Provide either
        frequency or phase.
    rate: float
        The sampling rate for phase or frequency, in Hz
    taus: np.array
        Array of tau values for which to compute measurement

    References
    ----------
    * NIST [SP1065]_ eqn (25) page 23
    """
    phase = input_to_phase(data, rate, data_type)
    (phase, m, taus_used) = tau_generator(phase, rate, taus)
    N = len(phase)

    # totdev requires a new dataset
    # Begin by adding reflected data before dataset
    x1 = 2.0 * phase[0] * np.ones((N - 2,))
    x1 = x1 - phase[1:-1]
    x1 = x1[::-1]

    # Reflected data at end of dataset
    x2 = 2.0 * phase[-1] * np.ones((N - 2,))
    x2 = x2 - phase[1:-1][::-1]

    # check length of new dataset
    assert len(x1)+len(phase)+len(x2) == 3*N - 4
    # Combine into a single array
    x = np.zeros((3*N - 4))
    x[0:N-2] = x1
    x[N-2:2*(N-2)+2] = phase  # original data in the middle
    x[2*(N-2)+2:] = x2

    devs = np.zeros_like(taus_used)
    deverrs = np.zeros_like(taus_used)
    ns = np.zeros_like(taus_used)

    mid = len(x1)

    for idx, mj in enumerate(m):
        mj = int(mj)
        d0 = x[mid + 1:]
        d1 = x[mid + mj + 1:]
        d1n = x[mid - mj + 1:]
        e = min(len(d0), len(d1), len(d1n))

        v_arr = d1n[:e] - 2.0 * d0[:e] + d1[:e]
        dev = np.sum(v_arr[:mid] * v_arr[:mid])

        dev /= float(2 * pow(mj / rate, 2) * (N - 2))
        dev = np.sqrt(dev)
        devs[idx] = dev
        deverrs[idx] = dev / np.sqrt(mid)
        ns[idx] = mid

    return remove_small_ns(taus_used, devs, deverrs, ns)


def ttotdev(data, rate=1.0, data_type="phase", taus=None):
    """ Time Total Deviation

    Modified total variance scaled by :math:`\\tau^2 / 3`

    .. math::

        \\sigma^2_{TTOTDEV}( \\tau ) = { \\tau^2 \\over 3 }
        \\sigma^2_{MTOTDEV}( \\tau )

    Note that [SP1065]_ erroneously has tau-cubed here (!).

    References
    ----------
    * NIST [SP1065]_ eqn (28) page 26.
    """

    (taus, mtotdevs, mde, ns) = mtotdev(data, data_type=data_type,
                                        rate=rate, taus=taus)
    td = taus*mtotdevs / np.sqrt(3.0)
    tde = td / np.sqrt(ns)
    return taus, td, tde, ns


def mtotdev(data, rate=1.0, data_type="phase", taus=None):
    """ Modified Total deviation.

    Better confidence at long averages for modified Allan

    FIXME: bias-correction http://www.wriley.com/CI2.pdf page 6

    The variance is scaled up (divided by this number) based on the
    noise-type identified.

    +------------+------------------+
    | noise type | bias correction  |
    +============+==================+
    | WPM        |  0.94            |
    +------------+------------------+
    | FPM        | 0.83             |
    +------------+------------------+
    | WFM        |   0.73           |
    +------------+------------------+
    | FFM        |  0.70            |
    +------------+------------------+
    | RWFM       |    0.69          |
    +------------+------------------+

    Parameters
    ----------
    data: np.array
        Input data. Provide either phase or frequency (fractional,
        adimensional).
    rate: float
        The sampling rate for data, in Hz. Defaults to 1.0
    data_type: {'phase', 'freq'}
        Data type, i.e. phase or frequency. Defaults to "phase".
    taus: np.array
        Array of tau values, in seconds, for which to compute statistic.
        Optionally set taus=["all"|"octave"|"decade"] for automatic
        tau-list generation.

    References
    ----------
    * NIST [SP1065]_ eqn (27) page 25
    """
    phase = input_to_phase(data, rate, data_type)
    (phase, ms, taus_used) = tau_generator(phase, rate, taus,
                                           maximum_m=float(len(phase))/3.0)
    devs = np.zeros_like(taus_used)
    deverrs = np.zeros_like(taus_used)
    ns = np.zeros_like(taus_used)

    for idx, mj in enumerate(ms):
        devs[idx], deverrs[idx], ns[idx] = calc_mtotdev_phase(phase, rate, mj)

    return remove_small_ns(taus_used, devs, deverrs, ns)


def calc_mtotdev_phase(phase, rate, m):
    """ PRELIMINARY - REQUIRES FURTHER TESTING.
        calculation of mtotdev for one averaging factor m
        tau = m*tau0

    Computed from a set of N - 3m + 1 subsequences of 3m points.
    1. A linear trend (frequency offset) is removed from the subsequence
       by averaging the first and last halves of the subsequence and
       dividing by half the interval.
    2. The offset-removed subsequence is extended at both ends
       by uninverted, even reflection.

    References
    ----------
    * [Howe1999]_
    * NIST [SP1065]_ Eqn (27), page 25.
    """
    tau0 = 1.0/rate
    N = len(phase)  # phase data, N points
    m = int(m)

    n = 0      # number of terms in the sum, for error estimation
    dev = 0.0  # the deviation we are computing
    # print('calc_mtotdev N=%d m=%d' % (N,m) )
    for i in range(0, N-3*m+1):
        # subsequence of length 3m, from the original phase data
        xs = phase[i:i+3*m]
        assert len(xs) == 3*m
        # Step 1.
        # remove linear trend. by averaging first/last half,
        # computing slope, and subtracting
        half1_idx = int(np.floor(3*m/2.0))
        half2_idx = int(np.ceil(3*m/2.0))
        # m
        # 1    0:1   2:2
        mean1 = np.mean(xs[:half1_idx])
        mean2 = np.mean(xs[half2_idx:])

        if int(3*m) % 2 == 1:  # m is odd
            # 3m = 2k+1 is odd, with the averages at both ends over k points
            # the distance between the averages is then k+1 = (3m-1)/2 +1
            slope = (mean2-mean1) / ((0.5*(3*m-1)+1)*tau0)
        else:  # m is even
            # 3m = 2k is even, so distance between averages is k=m/2
            slope = (mean2-mean1) / (0.5*3*m*tau0)

        # remove the linear trend
        x0 = [x - slope*idx*tau0 for (idx, x) in enumerate(xs)]
        x0_flip = x0[::-1]  # left-right flipped version of array

        # Step 2.
        # extend sequence, by uninverted even reflection
        # extended sequence xstar, of length 9m,
        xstar = np.concatenate((x0_flip, x0, x0_flip))
        assert len(xstar) == 9*m

        # now compute mdev on these 9m points
        # 6m unique groups of m-point averages,
        # use all possible overlapping second differences
        # one term in the 6m sum:  [ x_i - 2 x_i+m + x_i+2m ]^2
        squaresum = 0.0

        # below we want the following sums
        # (averages, see squaresum below where we divide by m)
        # xmean1=np.sum(xstar[j     :   j+m])
        # xmean2=np.sum(xstar[j+m   : j+2*m])
        # xmean3=np.sum(xstar[j+2*m : j+3*m])
        # for speed these are not computed with np.sum or np.mean in each loop
        # instead they are initialized at m=0, and then just updated
        for j in range(0, 6*m):  # summation of the 6m terms.
            # faster inner sum, based on Stable32 MTC.c code
            if j == 0:
                # intialize the sum
                xmean1 = np.sum(xstar[0:m])
                xmean2 = np.sum(xstar[m:2*m])
                xmean3 = np.sum(xstar[2*m:3*m])
            else:
                # j>=1, subtract old point, add new point
                xmean1 = xmean1 - xstar[j-1] + xstar[j+m-1]
                xmean2 = xmean2 - xstar[m+j-1] + xstar[j+2*m-1]
                xmean3 = xmean3 - xstar[2*m+j-1] + xstar[j+3*m-1]

            squaresum += pow((xmean1 - 2.0*xmean2 + xmean3)/float(m), 2)

        squaresum = (1.0/(6.0*m)) * squaresum
        dev += squaresum
        n = n+1

    # scaling in front of double-sum
    assert n == N-3*m+1  # sanity check on the number of terms n
    dev = dev * 1.0 / (2.0*pow(m*tau0, 2)*(N-3*m+1))
    dev = np.sqrt(dev)
    error = dev / np.sqrt(n)
    return (dev, error, n)


def htotdev(data, rate=1.0, data_type="phase", taus=None):
    """ Hadamard Total deviation.

    Better confidence at long averages for Hadamard deviation

    PRELIMINARY - REQUIRES FURTHER TESTING.

    Computed for N fractional frequency points y_i with sampling
    period tau0, analyzed at tau = m*tau0
    1. remove linear trend by averaging first and last half,
    and dividing by interval
    2. extend sequence by uninverted even reflection
    3. compute Hadamard for extended, length 9m, sequence.

    FIXME: bias corrections from http://www.wriley.com/CI2.pdf
    W FM    0.995      alpha= 0
    F FM    0.851      alpha=-1
    RW FM   0.771      alpha=-2
    FW FM   0.717      alpha=-3
    RR FM   0.679      alpha=-4

    Parameters
    ----------
    data: np.array
        Input data. Provide either phase or frequency (fractional,
        adimensional).
    rate: float
        The sampling rate for data, in Hz. Defaults to 1.0
    data_type: {'phase', 'freq'}
        Data type, i.e. phase or frequency. Defaults to "phase".
    taus: np.array
        Array of tau values, in seconds, for which to compute statistic.
        Optionally set taus=["all"|"octave"|"decade"] for automatic
        tau-list generation.

    """
    if data_type == "phase":
        phase = data
        freq = phase2frequency(phase, rate)
    elif data_type == "freq":
        phase = frequency2phase(data, rate)
        freq = data
    else:
        raise Exception("unknown data_type: " + data_type)

    rate = float(rate)
    (freq, ms, taus_used) = tau_generator(freq, rate, taus,
                                          maximum_m=float(len(freq))/3.0)
    phase = np.array(phase)
    freq = np.array(freq)
    devs = np.zeros_like(taus_used)
    deverrs = np.zeros_like(taus_used)
    ns = np.zeros_like(taus_used)

    # NOTE at mj==1 we use ohdev(), based on comment from here:
    # http://www.wriley.com/paper4ht.htm
    # "For best consistency, the overlapping Hadamard variance is used
    # instead of the Hadamard total variance at m=1"
    # FIXME: this uses both freq and phase datasets,
    # which uses double the memory really needed...
    for idx, mj in enumerate(ms):
        if int(mj) == 1:
            (devs[idx],
             deverrs[idx],
             ns[idx]) = calc_hdev_phase(phase, rate, mj, 1)
        else:
            (devs[idx],
             deverrs[idx],
             ns[idx]) = calc_htotdev_freq(freq, mj)

    return remove_small_ns(taus_used, devs, deverrs, ns)


def calc_htotdev_freq(freq, m):
    """ calculation of htotdev for one averaging factor m
        tau = m*tau0

    PRELIMINARY - REQUIRES FURTHER TESTING.

    Parameters
    ----------
    frequency: np.array
        Fractional frequency data (nondimensional).
    m: int
        Averaging factor. tau = m*tau0, where tau0=1/rate.
    """

    N = int(len(freq))  # frequency data, N points
    m = int(m)
    n = 0      # number of terms in the sum, for error estimation
    dev = 0.0  # the deviation we are computing
    for i in range(0, N-3*m+1):
        # subsequence of length 3m, from the original phase data
        xs = freq[i:i+3*m]
        assert len(xs) == 3*m
        # remove linear trend. by averaging first/last half,
        # computing slope, and subtracting
        half1_idx = int(np.floor(3*m/2.0))
        half2_idx = int(np.ceil(3*m/2.0))
        # m
        # 1    0:1   2:2
        mean1 = np.mean(xs[:half1_idx])
        mean2 = np.mean(xs[half2_idx:])

        if int(3*m) % 2 == 1:  # m is odd
            # 3m = 2k+1 is odd, with the averages at both ends over k points
            # the distance between the averages is then k+1 = (3m-1)/2 +1
            slope = (mean2-mean1) / ((0.5*(3*m-1)+1))
        else:  # m is even
            # 3m = 2k is even, so distance between averages is k=3m/2
            slope = (mean2-mean1) / (0.5*3*m)

        # remove the linear trend
        x0 = [x - slope*(idx-np.floor(3*m/2)) for (idx, x) in enumerate(xs)]
        x0_flip = x0[::-1]  # left-right flipped version of array
        # extended sequence, to length 9m, by uninverted even reflection
        xstar = np.concatenate((x0_flip, x0, x0_flip))
        assert len(xstar) == 9*m

        # now compute totdev on these 9m points
        # 6m unique groups of m-point averages,
        # all possible overlapping second differences
        # one term in the 6m sum:  [ x_i - 2 x_i+m + x_i+2m ]^2
        squaresum = 0.0
        k = 0
        for j in range(0, 6*int(m)):  # summation of the 6m terms.
            # old naive code
            # xmean1 = np.mean(xstar[j+0*m : j+1*m])
            # xmean2 = np.mean(xstar[j+1*m : j+2*m])
            # xmean3 = np.mean(xstar[j+2*m : j+3*m])
            # squaresum += pow(xmean1 - 2.0*xmean2 + xmean3, 2)
            # new faster way of doing the sums
            if j == 0:
                # intialize the sum
                xmean1 = np.sum(xstar[0:m])
                xmean2 = np.sum(xstar[m:2*m])
                xmean3 = np.sum(xstar[2*m:3*m])
            else:
                # j>=1, subtract old point, add new point
                xmean1 = xmean1 - xstar[j-1] + xstar[j+m-1]
                xmean2 = xmean2 - xstar[m+j-1] + xstar[j+2*m-1]
                xmean3 = xmean3 - xstar[2*m+j-1] + xstar[j+3*m-1]

            squaresum += pow((xmean1 - 2.0*xmean2 + xmean3)/float(m), 2)

            k = k+1
        assert k == 6*m  # check number of terms in the sum
        squaresum = (1.0/(6.0*k)) * squaresum
        dev += squaresum
        n = n+1

    # scaling in front of double-sum
    assert n == N-3*m+1  # sanity check on the number of terms n
    dev = dev * 1.0/(N-3*m+1)
    dev = np.sqrt(dev)
    error = dev / np.sqrt(n)
    return (dev, error, n)


def theo1(data, rate=1.0, data_type="phase", taus=None):
    """ Theo1 is a two-sample variance with improved confidence and extended averaging factor range.

    PRELIMINARY - REQUIRES FURTHER TESTING.

    .. math::

        \\sigma^2_{THEO1}(m\\tau_0) = { 1 \\over  (m \\tau_0 )^2 (N-m) }
            \\sum_{i=1}^{N-m}   \\sum_{\\delta=0}^{m/2-1}
            {1\\over m/2-\\delta}\\lbrace
                ({x}_{i} - x_{i-\\delta +m/2}) +
                (x_{i+m}- x_{i+\\delta +m/2}) \\rbrace^2


    Where :math:`10<=m<=N-1` is even.

    FIXME: bias correction

    Parameters
    ----------
    data: np.array
        Input data. Provide either phase or frequency (fractional,
        adimensional).
    rate: float
        The sampling rate for data, in Hz. Defaults to 1.0
    data_type: {'phase', 'freq'}
        Data type, i.e. phase or frequency. Defaults to "phase".
    taus: np.array
        Array of tau values, in seconds, for which to compute statistic.
        Optionally set taus=["all"|"octave"|"decade"] for automatic
        tau-list generation.

    References
    ----------
    * NIST [SP1065]_ eq (30) page 29
    """
    phase = input_to_phase(data, rate, data_type)

    tau0 = 1.0/rate
    (phase, ms, taus_used) = tau_generator(phase, rate, taus, even=True)

    devs = np.zeros_like(taus_used)
    deverrs = np.zeros_like(taus_used)
    ns = np.zeros_like(taus_used)

    N = len(phase)
    for idx, m in enumerate(ms):
        m = int(m)  # to avoid: VisibleDeprecationWarning: using a
        # non-integer number instead of an integer will
        # result in an error in the future
        assert m % 2 == 0  # m must be even
        dev = 0
        n = 0
        for i in range(int(N-m)):
            s = 0
            for d in range(int(m/2)):  # inner sum
                pre = 1.0 / (float(m)/2 - float(d))
                s += pre*pow(phase[i]-phase[i-d+int(m/2)] +
                             phase[i+m]-phase[i+d+int(m/2)], 2)
                n = n+1
            dev += s
        assert n == (N-m)*m/2  # N-m outer sums, m/2 inner sums
        dev = dev/(0.75*(N-m)*pow(m*tau0, 2))
        # factor 0.75 used here? http://tf.nist.gov/general/pdf/1990.pdf
        # but not here? http://tf.nist.gov/timefreq/general/pdf/2220.pdf
        # (page 29)
        devs[idx] = np.sqrt(dev)
        deverrs[idx] = devs[idx] / np.sqrt(N-m)
        ns[idx] = n

    return remove_small_ns(taus_used, devs, deverrs, ns)


def tierms(data, rate=1.0, data_type="phase", taus=None):
    """ Time Interval Error RMS.

    Parameters
    ----------
    data: np.array
        Input data. Provide either phase or frequency (fractional,
        adimensional).
    rate: float
        The sampling rate for data, in Hz. Defaults to 1.0
    data_type: {'phase', 'freq'}
        Data type, i.e. phase or frequency. Defaults to "phase".
    taus: np.array
        Array of tau values, in seconds, for which to compute statistic.
        Optionally set taus=["all"|"octave"|"decade"] for automatic
        tau-list generation.

    """
    phase = input_to_phase(data, rate, data_type)
    (data, m, taus_used) = tau_generator(phase, rate, taus)

    count = len(phase)

    devs = np.zeros_like(taus_used)
    deverrs = np.zeros_like(taus_used)
    ns = np.zeros_like(taus_used)

    for idx, mj in enumerate(m):
        mj = int(mj)

        # This seems like an unusual way to
        phases = np.column_stack((phase[:-mj], phase[mj:]))
        p_max = np.max(phases, axis=1)
        p_min = np.min(phases, axis=1)
        phases = p_max - p_min
        tie = np.sqrt(np.mean(phases * phases))

        ncount = count - mj

        devs[idx] = tie
        deverrs[idx] = 0 / np.sqrt(ncount)  # TODO! I THINK THIS IS WRONG!
        ns[idx] = ncount

    return remove_small_ns(taus_used, devs, deverrs, ns)


def mtie_rolling_window(a, window):
    """
    Make an ndarray with a rolling window of the last dimension, from
    http://mail.scipy.org/pipermail/numpy-discussion/2011-January/054401.html

    Parameters
    ----------
    a : array_like
        Array to add rolling window to
    window : int
        Size of rolling window

    Returns
    -------
    Array that is a view of the original array with a added dimension
    of size window.

    Note
    ----
    This may consume large amounts of memory. See discussion:
    https://mail.python.org/pipermail/numpy-discussion/2011-January/054364.html
    https://mail.python.org/pipermail/numpy-discussion/2011-January/054370.html

    """
    if window < 1:
        raise ValueError("`window` must be at least 1.")
    if window > a.shape[-1]:
        raise ValueError("`window` is too long.")
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)


def mtie(data, rate=1.0, data_type="phase", taus=None):
    """ Maximum Time Interval Error.

    Parameters
    ----------
    data: np.array
        Input data. Provide either phase or frequency (fractional,
        adimensional).
    rate: float
        The sampling rate for data, in Hz. Defaults to 1.0
    data_type: {'phase', 'freq'}
        Data type, i.e. phase or frequency. Defaults to "phase".
    taus: np.array
        Array of tau values, in seconds, for which to compute statistic.
        Optionally set taus=["all"|"octave"|"decade"] for automatic
        tau-list generation.

    Notes
    -----
    this seems to correspond to Stable32 setting "Fast(u)"
    Stable32 also has "Decade" and "Octave" modes where the
    dataset is extended somehow?
    """
    phase = input_to_phase(data, rate, data_type)
    (phase, m, taus_used) = tau_generator(phase, rate, taus)
    devs = np.zeros_like(taus_used)
    deverrs = np.zeros_like(taus_used)
    ns = np.zeros_like(taus_used)

    for idx, mj in enumerate(m):
        try:
            # the older algorithm uses a lot of memory
            # but can be used for short datasets.
            rw = mtie_rolling_window(phase, int(mj + 1))
            win_max = np.max(rw, axis=1)
            win_min = np.min(rw, axis=1)
            tie = win_max - win_min
            dev = np.max(tie)
        except ValueError:
            if int(mj + 1) < 1:
                raise ValueError("`window` must be at least 1.")
            if int(mj + 1) > phase.shape[-1]:
                raise ValueError("`window` is too long.")

            mj = int(mj)
            currMax = np.max(phase[0:mj])
            currMin = np.min(phase[0:mj])
            dev = currMax - currMin
            for winStartIdx in range(1, int(phase.shape[0] - mj)):
                winEndIdx = mj + winStartIdx
                if currMax == phase[winStartIdx - 1]:
                    currMax = np.max(phase[winStartIdx:winEndIdx])
                elif currMax < phase[winEndIdx]:
                    currMax = phase[winEndIdx]

                if currMin == phase[winStartIdx - 1]:
                    currMin = np.min(phase[winStartIdx:winEndIdx])
                elif currMin > phase[winEndIdx]:
                    currMin = phase[winEndIdx]

                if dev < currMax - currMin:
                    dev = currMax - currMin

        ncount = phase.shape[0] - mj
        devs[idx] = dev
        deverrs[idx] = dev / np.sqrt(ncount)
        ns[idx] = ncount

    return remove_small_ns(taus_used, devs, deverrs, ns)


def mtie_phase_fast(phase, rate=1.0, data_type="phase", taus=None):
    """ fast binary decomposition algorithm for MTIE

    References
    ----------
    * [Bregni2001]
    """
    #
    # !!!!!!!
    # FIXME: mtie_phase_fast() is incomplete.
    # !!!!!!!
    #
    rate = float(rate)
    phase = np.asarray(phase)
    k_max = int(np.floor(np.log2(len(phase))))
    phase = phase[0:pow(2, k_max)]  # truncate data to 2**k_max datapoints
    assert len(phase) == pow(2, k_max)
    taus = [pow(2, k) for k in range(k_max)]

    print("taus N=", len(taus), " ", taus)
    devs = np.zeros(len(taus))
    deverrs = np.zeros(len(taus))
    ns = np.zeros(len(taus))
    taus_used = np.array(taus)  # [(1.0/rate)*t for t in taus]
    # matrices to store results
    mtie_max = np.zeros((len(phase)-1, k_max))
    mtie_min = np.zeros((len(phase)-1, k_max))
    for kidx in range(k_max):
        k = kidx+1
        imax = len(phase)-pow(2, k)+1
        # print k, imax
        tie = np.zeros(imax)
        ns[kidx] = imax
        # print np.max( tie )
        for i in range(imax):
            if k == 1:
                mtie_max[i, kidx] = max(phase[i], phase[i+1])
                mtie_min[i, kidx] = min(phase[i], phase[i+1])
            else:
                p = int(pow(2, k-1))
                mtie_max[i, kidx] = max(mtie_max[i, kidx-1],
                                        mtie_max[i+p, kidx-1])
                mtie_min[i, kidx] = min(mtie_min[i, kidx-1],
                                        mtie_min[i+p, kidx-1])

        # for i in range(imax):
            tie[i] = mtie_max[i, kidx] - mtie_min[i, kidx]
            # print tie[i]
        devs[kidx] = np.amax(tie)  # maximum along axis

    devs = np.array(devs)
    print("devs N=", len(devs), " ", devs)
    print("taus N=", len(taus_used), " ", taus_used)
    return remove_small_ns(taus_used, devs, deverrs, ns)


########################################################################
#
#  gap resistant Allan deviation
#


def gradev(data, rate=1.0, data_type="phase", taus=None, ci=0.9, noisetype='wp'):
    """ Gap resistant overlapping Allan deviation

    Parameters
    ----------
    data: np.array
        Input data. Provide either phase or frequency (fractional,
        adimensional). Warning : phase data works better (frequency data is
        first trantformed into phase using numpy.cumsum() function, which can
        lead to poor results).
    rate: float
        The sampling rate for data, in Hz. Defaults to 1.0
    data_type: {'phase', 'freq'}
        Data type, i.e. phase or frequency. Defaults to "phase".
    taus: np.array
        Array of tau values, in seconds, for which to compute statistic.
        Optionally set taus=["all"|"octave"|"decade"] for automatic
        tau-list generation.
    ci: float
        the total confidence interval desired, i.e. if ci = 0.9, the bounds
        will be at 0.05 and 0.95.
    noisetype: string
        the type of noise desired:
        'wp' returns white phase noise.
        'wf' returns white frequency noise.
        'fp' returns flicker phase noise.
        'ff' returns flicker frequency noise.
        'rf' returns random walk frequency noise.
        If the input is not recognized, it defaults to idealized, uncorrelated
        noise with (N-1) degrees of freedom.

    Returns
    -------
    taus: np.array
        list of tau vales in seconds
    adev: np.array
        deviations
    [err_l, err_h] : list of len()==2, np.array
        the upper and lower bounds of the confidence interval taken as
        distances from the the estimated two sample variance.
    ns: np.array
        numper of terms n in the adev estimate.

    """
    if data_type == "freq":
        print("Warning : phase data is preferred as input to gradev()")
    phase = input_to_phase(data, rate, data_type)
    (data, m, taus_used) = tau_generator(phase, rate, taus)

    ad = np.zeros_like(taus_used)
    ade_l = np.zeros_like(taus_used)
    ade_h = np.zeros_like(taus_used)
    adn = np.zeros_like(taus_used)

    for idx, mj in enumerate(m):
        (dev, deverr, n) = calc_gradev_phase(data,
                                             rate,
                                             mj,
                                             1,
                                             ci,
                                             noisetype)
        # stride=1 for overlapping ADEV
        ad[idx] = dev
        ade_l[idx] = deverr[0]
        ade_h[idx] = deverr[1]
        adn[idx] = n

    # Note that errors are split in 2 arrays
    return remove_small_ns(taus_used, ad, [ade_l, ade_h], adn)


def calc_gradev_phase(data, rate, mj, stride, confidence, noisetype):
    """ see http://www.leapsecond.com/tools/adev_lib.c
        stride = mj for nonoverlapping allan deviation
        stride = 1 for overlapping allan deviation

        see http://en.wikipedia.org/wiki/Allan_variance
             1       1
         s2y(t) = --------- sum [x(i+2) - 2x(i+1) + x(i) ]^2
                  2*tau^2
    """

    d2 = data[2 * int(mj)::int(stride)]
    d1 = data[1 * int(mj)::int(stride)]
    d0 = data[::int(stride)]

    n = min(len(d0), len(d1), len(d2))

    v_arr = d2[:n] - 2 * d1[:n] + d0[:n]

    # only average for non-nans
    n = len(np.where(np.isnan(v_arr) == False)[0])  # noqa

    if n == 0:
        RuntimeWarning("Data array length is too small: %i" % len(data))
        n = 1

    N = len(np.where(np.isnan(data) == False)[0])  # noqa

    # a summation robust to nans
    s = np.nansum(v_arr * v_arr)

    dev = np.sqrt(s / (2.0 * n)) / mj * rate
    # deverr = dev / np.sqrt(n) # old simple errorbars
    if noisetype == 'wp':
        alpha = 2
    elif noisetype == 'wf':
        alpha = 0
    elif noisetype == 'fp':
        alpha = -2
    else:
        alpha = None

    if n > 1:
        edf = ci.edf_simple(N, mj, alpha)
        deverr = ci.confidence_interval(dev, confidence, edf)
    else:
        deverr = [0, 0]

    return dev, deverr, n


def psd2allan(S_y, f=1.0, kind='adev', base=2):
    """ Convert a given (one-sided) power spectral density :math:`S_y(f)` to Allan
        deviation or modified Allan deviation

    For ergodic noise, the Allan variance or modified Allan variance
    is related to the power spectral density :math:`S_y` of the fractional
    frequency deviation:

    .. math::

        \\sigma^2_y(\\tau) = 2 \\int_0^\\infty S_y(f)
        \\left| \\sin(\\pi f \\tau)^{(k+1)} \\over (\\pi f \\tau)^k \\right|^2 df,


    where :math:`f` is the Fourier frequency and :math:`\\tau` the averaging
    time. The exponent :math:`k` is 1 for the Allan variance and 2 for the
    modified Allan variance.

    psd2allan() implements the integral by discrete numerical integration via
    a sum.

    Parameters
    ----------
    S_y: np.array
        Single-sided power spectral density (PSD) of fractional frequency
        deviation S_y in 1/Hz^2. First element is S_y(f=0).
    f: np.array or scalar numeric (float or int)
        if np.array: Spectral frequency vector in Hz
        if numeric scalar: Spectral frequency step in Hz
        default: Spectral frequency step 1 Hz
    kind: {'adev', 'mdev'}
        Which kind of Allan deviation to compute. Defaults to 'adev'
    base: float
        Base for logarithmic spacing of tau values. E.g. base= 10: decade,
        base= 2: octave, base <= 1: all

    Returns
    -------
    (taus_used, ad): tuple
          Tuple of 2 values
    taus_used: np.array
        tau values for which ad computed
    ad: np.array
        Computed Allan deviation of requested kind for each tau value

    References
    ----------
    * NIST [SP1065]_ eqs (65-66), page 73.
    * [Benkler2015]_ eqn (23).
    """
    # determine taus from df
    # first oversample S_y by a factor of 10 in order to avoid numerical
    # problem at tau > 1/2df
    if isinstance(S_y, np.ndarray):
        if isinstance(f, np.ndarray):  # f is frequency vector
            df = f[1]-f[0]
        elif np.isscalar(f):
            # assume that f is the frequency step, not frequency vector
            df = f
        else:
            raise ValueError(np.ndarray, float, int)
    else:
        raise ValueError(np.ndarray)  # raise error
    oversamplingfactor = 4
    df0 = oversamplingfactor * df
    f0 = np.arange(S_y.size * df0, step=df0)
    f = np.arange(df, (S_y.size - 1) * df0 + df, df)
    S_y = interpolate.interp1d(f0, S_y, kind='cubic')(f)
    f = f / oversamplingfactor

    tau0 = 1/np.max(f)  # minimum tau derived from the given frequency vector
    n = 1/df/tau0/2
    if base > 1:
        m = np.unique(np.round(np.append(base**np.arange(
            np.floor(np.log(n)/np.log(base))), n)))
    else:
        m = np.arange(1, n)
    taus_used = m*tau0

    # TODO: In principle, this approach can be extended to the other kinds of
    # Allan deviations, we just need to determine the respective transfer
    # function in the frequency domain.

    if kind[0].lower() == 'a':   # for ADEV
        exponent = 1.0
    elif kind[0].lower() == 'm':  # for modADEV
        exponent = 2.0

    integrand = np.array([
        S_y *
        np.abs(np.sin(np.pi * f * taus_used[idx])**(exponent + 1.0)
               / (np.pi * f * taus_used[idx])**exponent)**2.0
        for idx, mj in enumerate(m)])
    integrand = np.insert(integrand, 0, 0.0, axis=1)
    f = np.insert(f, 0, 0.0)
    ad = np.sqrt(2.0 * simpson(integrand, f))
    return taus_used, ad


########################################################################
#
#  Various helper functions and utilities
#


def input_to_phase(data, rate, data_type):
    """ Take either phase or frequency as input and return phase
    """
    if data_type == "phase":
        return data
    elif data_type == "freq":
        return frequency2phase(data, rate)
    else:
        raise Exception("unknown data_type: " + data_type)


def tau_generator(data, rate, taus=None, v=False, even=False, maximum_m=-1):
    """ pre-processing of the tau-list given by the user (Helper function)

    Does sanity checks, sorts data, removes duplicates and invalid values.
    Generates a tau-list based on keywords 'all', 'decade', 'octave'.
    Uses 'octave' by default if no taus= argument is given.

    Parameters
    ----------
    data: np.array
        data array
    rate: float
        Sample rate of data in Hz. Time interval between measurements
        is 1/rate seconds.
    taus: np.array
        Array of tau values for which to compute measurement.
        Alternatively one of the keywords: "all", "octave", "decade".
        Defaults to "octave" if omitted.

        +----------+--------------------------------+
        | keyword  |   averaging-factors            |
        +==========+================================+
        | "all"    |  1, 2, 3, 4, ..., len(data)    |
        +----------+--------------------------------+
        | "octave" |  1, 2, 4, 8, 16, 32, ...       |
        +----------+--------------------------------+
        | "decade" |  1, 2, 4, 10, 20, 40, 100, ... |
        +----------+--------------------------------+
        | "log10"  |  approx. 10 points per decade  |
        +----------+--------------------------------+
    v: bool
        verbose output if True
    even: bool
        require even m, where tau=m*tau0, for Theo1 statistic
    maximum_m: int
        limit m, where tau=m*tau0, to this value.
        used by mtotdev() and htotdev() to limit maximum tau.

    Returns
    -------
    (data, m, taus): tuple
        List of computed values
    data: np.array
        Data
    m: np.array
        Tau in units of data points
    taus: np.array
        Cleaned up list of tau values
    """

    if rate == 0:
        raise RuntimeError("Warning! rate==0")

    if taus is None:  # empty or no tau-list supplied
        taus = "octave"  # default to octave
    elif isinstance(taus, list) and taus == []:  # empty list
        taus = "octave"

    # numpy array or non-empty list detected first
    if isinstance(taus, np.ndarray) or isinstance(taus, list) and len(taus):
        pass
    elif taus == "all":  # was 'is'
        taus = (1.0/rate)*np.linspace(1.0, len(data), len(data))
    elif taus == "octave":
        maxn = np.floor(np.log2(len(data)))
        taus = (1.0/rate)*np.logspace(0, int(maxn), int(maxn+1), base=2.0)
    elif taus == "log10":
        maxn = np.log10(len(data))
        taus = (1.0/rate)*np.logspace(0, maxn, int(10*maxn), base=10.0)
        if v:
            print("tau_generator: maxn %.1f" % maxn)
            print("tau_generator: taus=" % taus)
    elif taus == "decade":  # 1, 2, 4, 10, 20, 40, spacing similar to Stable32
        maxn = np.floor(np.log10(len(data)))
        taus = []
        for k in range(int(maxn+1)):
            taus.append(1.0*(1.0/rate)*pow(10.0, k))
            taus.append(2.0*(1.0/rate)*pow(10.0, k))
            taus.append(4.0*(1.0/rate)*pow(10.0, k))

    data, taus = np.array(data), np.array(taus)
    rate = float(rate)
    m = []  # integer averaging factor. tau = m*tau0

    if maximum_m == -1:  # if no limit given
        maximum_m = len(data)
    # FIXME: should we use a "stop-ratio" like Stable32
    # found in Table III, page 9 of
    # "Evolution of frequency stability analysis software"
    # max(AF) = len(phase)/stop_ratio, where
    # function  stop_ratio
    # adev      5
    # oadev     4
    # mdev      4
    # tdev      4
    # hdev      5
    # ohdev     4
    # totdev    2
    # tierms    4
    # htotdev   3
    # mtie      2
    # theo1     1
    # theoH     1
    # mtotdev   2
    # ttotdev   2

    m = np.round(taus * rate)
    taus_valid1 = m < len(data)
    taus_valid2 = m > 0
    taus_valid3 = m <= maximum_m
    taus_valid = taus_valid1 & taus_valid2 & taus_valid3
    m = m[taus_valid]
    m = m[m != 0]       # m is tau in units of datapoints
    m = np.unique(m)    # remove duplicates and sort

    if v:
        print("tau_generator: ", m)

    if len(m) == 0:
        print("Warning: sanity-check on tau failed!")
        print("   len(data)=", len(data), " rate=", rate, "taus= ", taus)

    taus2 = m / float(rate)

    if even:  # used by Theo1
        m_even_mask = ((m % 2) == 0)
        m = m[m_even_mask]
        taus2 = taus2[m_even_mask]

    return data, m, taus2


def tau_reduction(ms, rate, n_per_decade):
    """Reduce the number of taus to maximum of n per decade (Helper function)

    takes in a tau list and reduces the number of taus to a maximum amount per
    decade. This is only useful if more than the "decade" and "octave" but
    less than the "all" taus are wanted. E.g. to show certain features of
    the data one might want 100 points per decade.

    NOTE: The algorithm is slightly inaccurate for ms under n_per_decade, and
    will also remove some points in this range, which is usually fine.

    Typical use would be something like:

    .. code-block:: python

        (data, m, taus) = tau_generator(data, rate, taus="all")
        (m, taus) = tau_reduction(m, rate, n_per_decade)

    Parameters
    ----------
    ms: array of integers
        List of m values (assumed to be an "all" list) to remove points from.
    rate: float
        Sample rate of data in Hz. Time interval between measurements
        is 1/rate seconds. Used to convert to taus.
    n_per_decade: int
        Number of ms/taus to keep per decade.

    Returns
    -------
    m: np.array
        Reduced list of m values
    taus: np.array
        Reduced list of tau values

    """
    ms = np.int64(ms)
    keep = np.bool_(np.rint(n_per_decade*np.log10(ms[1:])) -
                    np.rint(n_per_decade*np.log10(ms[:-1])))
    # Adjust ms size to fit above-defined mask
    ms = ms[:-1]
    assert len(ms) == len(keep)
    ms = ms[keep]
    taus = ms/float(rate)

    return ms, taus


def remove_small_ns(taus, devs, deverrs, ns):
    """ Remove results with small number of samples.

    If n is small (==1), reject the result

    Parameters
    ----------
    taus: array
        List of tau values for which deviation were computed
    devs: array
        List of deviations
    deverrs: array or list of arrays
        List of estimated errors (possibly a list containing two arrays :
        upper and lower values)
    ns: array
        Number of samples for each point

    Returns
    -------
    (taus, devs, deverrs, ns): tuple
        Identical to input, except that values with low ns have been removed.

    """
    ns_big_enough = ns > 1

    o_taus = taus[ns_big_enough]
    o_devs = devs[ns_big_enough]
    o_ns = ns[ns_big_enough]
    if isinstance(deverrs, list):
        assert len(deverrs) < 3
        o_deverrs = [deverrs[0][ns_big_enough], deverrs[1][ns_big_enough]]
    else:
        o_deverrs = deverrs[ns_big_enough]
    if len(o_devs) == 0:
        print("remove_small_ns() nothing remains!?")
        raise UserWarning

    return o_taus, o_devs, o_deverrs, o_ns


def trim_data(x):
    """Trim leading and trailing NaNs from dataset

    This is done by browsing the array from each end and store the index of the
    first non-NaN in each case, the return the appropriate slice of the array
    """
    # Find indices for first and last valid data
    first = 0
    while np.isnan(x[first]):
        first += 1
    last = len(x)
    while np.isnan(x[last - 1]):
        last -= 1
    return x[first:last]


def three_cornered_hat_phase(phasedata_ab, phasedata_bc, phasedata_ca, rate, taus, function):
    """Three Cornered Hat Method

    Given three clocks A, B, C, we seek to find their variances
    :math:`\\sigma^2_A`, :math:`\\sigma^2_B`, :math:`\\sigma^2_C`.
    We measure three phase differences, assuming no correlation between
    the clocks, the measurements have variances:

    .. math::

        \\sigma^2_{AB} = \\sigma^2_{A} + \\sigma^2_{B}

        \\sigma^2_{BC} = \\sigma^2_{B} + \\sigma^2_{C}

        \\sigma^2_{CA} = \\sigma^2_{C} + \\sigma^2_{A}

    Which allows solving for the variance of one clock as:

    .. math::

        \\sigma^2_{A}  = {1 \\over 2} ( \\sigma^2_{AB} +
        \\sigma^2_{CA} - \\sigma^2_{BC} )

    and similarly cyclic permutations for :math:`\\sigma^2_B` and
    :math:`\\sigma^2_C`

    Parameters
    ----------
    phasedata_ab: np.array
        phase measurements between clock A and B, in seconds
    phasedata_bc: np.array
        phase measurements between clock B and C, in seconds
    phasedata_ca: np.array
        phase measurements between clock C and A, in seconds
    rate: float
        The sampling rate for phase, in Hz
    taus: np.array
        The tau values for deviations, in seconds
    function: allantools deviation function
        The type of statistic to compute, e.g. allantools.oadev

    Returns
    -------
    tau_ab: np.array
        Tau values corresponding to output deviations
    dev_a: np.array
        List of computed values for clock A

    References
    ----------
    * http://www.wriley.com/3-CornHat.htm

    """
    (tau_ab, dev_ab, err_ab, ns_ab) = function(phasedata_ab,
                                               data_type='phase',
                                               rate=rate, taus=taus)
    (tau_bc, dev_bc, err_bc, ns_bc) = function(phasedata_bc,
                                               data_type='phase',
                                               rate=rate, taus=taus)
    (tau_ca, dev_ca, err_ca, ns_ca) = function(phasedata_ca,
                                               data_type='phase',
                                               rate=rate, taus=taus)

    var_ab = dev_ab * dev_ab
    var_bc = dev_bc * dev_bc
    var_ca = dev_ca * dev_ca
    assert len(var_ab) == len(var_bc) == len(var_ca)
    var_a = 0.5 * (var_ab + var_ca - var_bc)

    var_a[var_a < 0] = 0  # don't return imaginary deviations (?)
    dev_a = np.sqrt(var_a)
    err_a = [d/np.sqrt(nn) for (d, nn) in zip(dev_a, ns_ab)]

    return tau_ab, dev_a, err_a, ns_ab

########################################################################
#
# simple conversions between frequency, phase(seconds), phase(radians)
#


def frequency2phase(freqdata, rate):
    """ integrate fractional frequency data and output phase data

    Parameters
    ----------
    freqdata: np.array
        Data array of fractional frequency measurements (nondimensional)
    rate: float
        The sampling rate for phase or frequency, in Hz

    Returns
    -------
    phasedata: np.array
        Time integral of fractional frequency data, i.e. phase (time) data
        in units of seconds.
        For phase in units of radians, see phase2radians()
    """
    dt = 1.0 / float(rate)
    # Protect against NaN values in input array (issue #60)
    # Reintroduces data trimming as in commit 503cb82
    freqdata = trim_data(freqdata)
    # Erik Benkler (PTB): Subtract mean value before cumsum in order to
    # avoid precision issues when we have small frequency fluctuations on
    # a large average frequency
    freqdata = freqdata - np.nanmean(freqdata)
    phasedata = np.cumsum(freqdata) * dt
    phasedata = np.insert(phasedata, 0, 0)  # FIXME: why do we do this?
    # so that phase starts at zero and len(phase)=len(freq)+1 ??
    return phasedata


def phase2radians(phasedata, v0):
    """ Convert phase in seconds to phase in radians

    Parameters
    ----------
    phasedata: np.array
        Data array of phase in seconds
    v0: float
        Nominal oscillator frequency in Hz

    Returns
    -------
    fi: np.array
        phase data in radians
    """
    fi = [2*np.pi*v0*xx for xx in phasedata]
    return np.array(fi)


def phase2frequency(phase, rate):
    """ Convert phase in seconds to fractional frequency

    Parameters
    ----------
    phase: np.array
        Data array of phase in seconds, length N
    rate: float
        The sampling rate for phase, in Hz

    Returns
    -------
    y: np.array
        Data array of fractional frequency, length N-1
    """
    y = rate*np.diff(phase)
    return y


def frequency2fractional(frequency, mean_frequency=-1):
    """ Convert frequency in Hz to fractional frequency

    Parameters
    ----------
    frequency: np.array
        Data array of frequency in Hz
    mean_frequency: float
        (optional) The nominal mean frequency, in Hz
        if omitted, defaults to mean frequency=np.mean(frequency)

    Returns
    -------
    y: np.array
        Data array of fractional frequency
    """
    if mean_frequency == -1:
        mu = np.mean(frequency)
    else:
        mu = mean_frequency
    y = [(x-mu)/mu for x in frequency]
    return y

# end of file allantools.py
