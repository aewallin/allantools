"""
Allan deviation tools
=====================

**Author:** Anders Wallin (anders.e.e.wallin "at" gmail.com)

Version history
---------------

**2016.3** 2016 March
- improve documentation and add __version__
- added Theo1 deviation theo1()
- added Hadamard Total Deviatio htotdev()
- added Modified Total Deviation mtotdev(), and Time Total Deviation ttotdev()
  see http://www.anderswallin.net/2016/03/modified-total-deviation-in-allantools/
- automatic tau-lists:  taus=[ "all" | "octave" | "decade" ]
- merge adev() and adev_phase() into one, requiring phase= or frequency= argument
- add GPS dataset as example and test

**2016.2** 2016 February
- update release on PyPi
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
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""

import os
import json
import numpy as np
import scipy.stats # used in uncertainty_estimate()

# Get version number from json metadata
pkginfo_path = os.path.join(os.path.dirname(__file__),
                            'allantools_info.json')
pkginfo = json.load(open(pkginfo_path))
__version__ = pkginfo["version"]


def tdev(phase=None, frequency=None, rate=1.0, taus=[]):
    """ Time deviation.
        Based on modified Allan variance.

    Parameters
    ----------
    phase: np.array
        Phase data in seconds. Provide either phase or frequency.
    frequency: np.array
        Fractional frequency data (nondimensional). Provide either frequency or phase.
    rate: float
        The sampling rate for phase or frequency, in Hz
    taus: np.array
        Array of tau values, in seconds, for which to compute statistic.
        Optionally set taus=[autotau.alltau|autotau.octave|autotau.decade] for automatic
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

    Notes
    -----
    http://en.wikipedia.org/wiki/Time_deviation
    """

    if phase == None:
        phase = frequency2phase(frequency, rate)

    (taus, md, mde, ns) = mdev(phase=phase, rate=rate, taus=taus)
    td = taus * md / np.sqrt(3.0)
    tde = td / np.sqrt(ns)
    return taus, td, tde, ns

def mdev(phase=None, frequency=None, rate=1.0, taus=[]):
    """  Modified Allan deviation.
         Used to distinguish W and F PM.

    NIST SP 1065 eqn (14) and (15), page 17

    Parameters
    ----------
    phase: np.array
        Phase data in seconds. Provide either phase or frequency.
    frequency: np.array
        Fractional frequency data (nondimensional). Provide either frequency or phase.
    rate: float
        The sampling rate for phase or frequency, in Hz
    taus: np.array
        Array of tau values for which to compute measurement

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

    Notes
    -----
    see http://www.leapsecond.com/tools/adev_lib.c

    NIST SP 1065 eqn (14)
    """

    if phase == None:
        phase = frequency2phase(frequency, rate)

    (phase, ms, taus_used) = tau_generator(phase, rate, taus=taus)    
    data, taus = np.array(phase), np.array(taus)
    
    md    = np.zeros_like(ms)
    mderr = np.zeros_like(ms)
    ns    = np.zeros_like(ms)

    # this is a 'loop-unrolled' algorithm following
    # http://www.leapsecond.com/tools/adev_lib.c
    for idx, m in enumerate(ms):
        m = int(m) # without this we get: VisibleDeprecationWarning: using a non-integer number instead of an integer will result in an error in the future
        tau = taus_used[idx]

        # First loop sum
        d0 = phase[0:m]
        d1 = phase[m:2*m]
        d2 = phase[2*m:3*m]
        e = min(len(d0), len(d1), len(d2))
        v = np.sum(d2[:e] - 2* d1[:e] + d0[:e])
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

def adev(phase=None, frequency=None, rate=1.0, taus=[]):
    """ Allan deviation.
        Classic - use only if required - relatively poor confidence.

    .. math::

        \\sigma^2_x(\\tau) =

    .. math::

        \\sigma^{2}_y(\\tau) =  { 1 \\over 2 } \\langle ( \\bar{y}_{n+1} - \\bar{y}_n )^2 \\rangle

    NIST SP 1065 eqn (6) and (7), pages 14 and 15

    Parameters
    ----------
    phase: np.array
        Phase data in seconds. Provide either phase or frequency.
    frequency: np.array
        Fractional frequency data (nondimensional). Provide either frequency or phase.
    rate: float
        The sampling rate for phase or frequency, in Hz
    taus: np.array
        Array of tau values for which to compute measurement


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

    """
    if phase == None:
        phase = frequency2phase(frequency, rate)
        
    (phase, m, taus_used) = tau_generator(phase, rate, taus)

    ad  = np.zeros_like(taus_used)
    ade = np.zeros_like(taus_used)
    adn = np.zeros_like(taus_used)

    for idx, mj in enumerate(m):  # loop through each tau value m(j)
        (ad[idx], ade[idx], adn[idx]) = calc_adev_phase(phase, rate, mj, mj)

    return remove_small_ns(taus_used, ad, ade, adn)  # tau, adev, adeverror, naverages


def calc_adev_phase(phase, rate, mj, stride):
    """  Main algorithm for adev() (stride=mj) and oadev() (stride=1)

        see http://www.leapsecond.com/tools/adev_lib.c
        stride = mj for nonoverlapping allan deviation

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

    Notes
    -----
    stride = mj for nonoverlapping Allan deviation
    stride = 1 for overlapping Allan deviation

    References
    ----------
    * http://en.wikipedia.org/wiki/Allan_variance
    * http://www.leapsecond.com/tools/adev_lib.c
    NIST SP 1065, eqn (7) and (11) page 16
    """
    mj = int(mj)
    d2 = phase[2 * mj::stride]
    d1 = phase[1 * mj::stride]
    d0 = phase[::stride]

    n = min(len(d0), len(d1), len(d2))

    if n == 0:
        RuntimeWarning("Data array length is too small: %i" % len(phase))
        n = 1

    v_arr = d2[:n] - 2 * d1[:n] + d0[:n]
    s = np.sum(v_arr * v_arr)

    dev = np.sqrt(s / (2.0 * n)) / mj  * rate
    deverr = dev / np.sqrt(n)

    return dev, deverr, n

def oadev(phase=None, frequency=None, rate=1.0, taus=[]):
    """ overlapping Allan deviation.
        General purpose - most widely used - first choice



    Parameters
    ----------
    phase: np.array
        Phase data in seconds. Provide either phase or frequency.
    frequency: np.array
        Fractional frequency data (nondimensional). Provide either frequency or phase.
    rate: float
        The sampling rate for phase or frequency, in Hz
    taus: np.array
        Array of tau values for which to compute measurement


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

    """
    if phase == None:
        phase = frequency2phase(frequency, rate)

    (phase, m, taus_used) = tau_generator(phase, rate, taus)
    ad  = np.zeros_like(taus_used)
    ade = np.zeros_like(taus_used)
    adn = np.zeros_like(taus_used)

    for idx, mj in enumerate(m):
        (ad[idx], ade[idx], adn[idx]) = calc_adev_phase(phase, rate, mj, 1)  # stride=1 for overlapping ADEV

    return remove_small_ns(taus_used, ad, ade, adn)  # tau, adev, adeverror, naverages

def ohdev(phase=None, frequency=None, rate=1.0, taus=[]):
    """ Overlapping Hadamard deviation.
        Better confidence than normal Hadamard.

    Parameters
    ----------
    phase: np.array
        Phase data in seconds. Provide either phase or frequency.
    frequency: np.array
        Fractional frequency data (nondimensional). Provide either frequency or phase.
    rate: float
        The sampling rate for phase or frequency, in Hz
    taus: np.array
        Array of tau values for which to compute measurement

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

    """
    if phase == None:
        phase = frequency2phase(frequency, rate)

    rate = float(rate)
    (phase, m, taus_used) = tau_generator(phase, rate, taus)
    hdevs = np.zeros_like(taus_used)
    hdeverrs = np.zeros_like(taus_used)
    ns = np.zeros_like(taus_used)

    for idx, mj in enumerate(m):
        (hdevs[idx], hdeverrs[idx], ns[idx]) = calc_hdev_phase(phase, rate, mj, 1)  # stride = 1

    return remove_small_ns(taus_used, hdevs, hdeverrs, ns)

def hdev(phase=None, frequency=None, rate=1.0, taus=[]):
    """ Hadamard deviation.
        Rejects frequency drift, and handles divergent noise.

    NIST SP 1065 eqn (17) and (18), page 20

    Parameters
    ----------
    phase: np.array
        Phase data in seconds. Provide either phase or frequency.
    frequency: np.array
        Fractional frequency data (nondimensional). Provide either frequency or phase.
    rate: float
        The sampling rate for phase or frequency, in Hz
    taus: np.array
        Array of tau values for which to compute measurement
    """
    if phase == None:
        phase = frequency2phase(frequency, rate)

    rate = float(rate)
    (phase, m, taus_used) = tau_generator(phase, rate, taus)
    hdevs = np.zeros_like(taus_used)
    hdeverrs = np.zeros_like(taus_used)
    ns = np.zeros_like(taus_used)

    for idx, mj in enumerate(m):
        hdevs[idx], hdeverrs[idx], ns[idx] = calc_hdev_phase(phase, rate, mj, mj)  # stride = mj

    return remove_small_ns(taus_used, hdevs, hdeverrs, ns)


def calc_hdev_phase(phase, rate, mj, stride):
    """ main calculation fungtion for HDEV and OHDEV

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

    Notes
    -----
    http://www.leapsecond.com/tools/adev_lib.c
                         1        N-3
         s2y(t) = --------------- sum [x(i+3) - 3x(i+2) + 3x(i+1) - x(i) ]^2
                  6*tau^2 (N-3m)  i=1

        N=M+1 phase measurements
        m is averaging factor

    NIST SP 1065 eqn (18) and (20) pages 20 and 21
    """

    tau0 = 1.0 / float(rate)
    mj = int(mj)

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

def totdev(phase=None, frequency=None, rate=1.0, taus=[]):
    """ Total deviation.
        Better confidence at long averages for Allan.

        FIXME: bias correction http://www.wriley.com/CI2.pdf page 5

    Parameters
    ----------
    phase: np.array
        Phase data in seconds. Provide either phase or frequency.
    frequency: np.array
        Fractional frequency data (nondimensional). Provide either frequency or phase.
    rate: float
        The sampling rate for phase or frequency, in Hz
    taus: np.array
        Array of tau values for which to compute measurement


    References
    ----------
    David A. Howe,
    *The total deviation approach to long-term characterization
    of frequency stability*,
    IEEE tr. UFFC vol 47 no 5 (2000)

    Notes
    -----
    ..
                     1        N-1
    totvar(t) = ------------  sum   [ x*(i-m) - 2x*(i) + x*(i+m) ]**2
                2 t**2 (N-2)  i=2

    where x* is a new dataset with 'reflected' data at start/end.

    the original data x is in the center of x*:
    x*(1-j) = 2x(1) - x(1+j)  for j=1..N-2
    x*(i)   = x(i)            for i=1..N
    x*(N+j) = 2x(N) - x(N-j)  for j=1..N-2
    x* has length 3N-4
    tau = m*tau0

    NIST SP 1065 eqn (25) page 23

    """
    if phase == None:
        phase = frequency2phase(frequency, rate)

    rate = float(rate)
    (phase, m, taus_used) = tau_generator(phase, rate, taus)
    n = len(phase)

    # totdev requires a new dataset
    # Begin by adding reflected data before dataset
    x1 = 2.0 * phase[0] * np.ones((n - 2,))
    x1 = x1 - phase[1:-1]
    x1 = x1[::-1]

    # Reflected data at end of dataset
    x2 = 2.0 * phase[-1] * np.ones((n - 2,))
    x2 = x2 - phase[1:-1][::-1]

    # Combine into a single array
    x = np.zeros((3*n - 4))
    x[0:n-2] = x1
    x[n-2:2*(n-2)+2] = phase # original data in the middle
    x[2*(n-2)+2:] = x2

    devs = np.zeros_like(taus_used)
    deverrs = np.zeros_like(taus_used)
    ns = np.zeros_like(taus_used)

    mid = len(x1)

    for idx, mj in enumerate(m):

        d0 = x[mid + 1:]
        d1 = x[mid  + mj + 1:]
        d1n = x[mid - mj + 1:]
        e = min(len(d0), len(d1), len(d1n))

        v_arr = d1n[:e] - 2.0 * d0[:e] + d1[:e]
        dev = np.sum(v_arr[:mid] * v_arr[:mid])

        dev /= float(2 * pow(mj / rate, 2) * (n - 2))
        dev = np.sqrt(dev)
        devs[idx] = dev
        deverrs[idx] = dev / np.sqrt(mid)
        ns[idx] = mid

    return remove_small_ns(taus_used, devs, deverrs, ns)

def ttotdev(phase=None, frequency=None, rate=1.0, taus=[]):
    """ Time Total Deviation
        modified total variance scaled by tau^2 / 3
        NIST SP 1065 eqn (28) page 26  <--- formula should have tau squared !?!
    """

    (taus, mtotdevs, mde, ns) = mtotdev(phase=phase, frequency=frequency, rate=rate, taus=taus)
    td = taus*mtotdevs / np.sqrt(3.0)
    tde = td / np.sqrt(ns)
    return taus, td, tde, ns

def mtotdev(phase=None, frequency=None, rate=1.0, taus=[]):
    """ PRELIMINARY - REQUIRES FURTHER TESTING.
        Modified Total deviation.
        Better confidence at long averages for modified Allan

        FIXME: bias-correction http://www.wriley.com/CI2.pdf page 6
        The variance is scaled up (divided by this number) based on the noise-type identified.
        WPM 0.94
        FPM 0.83
        WFM 0.73
        FFM 0.70
        RWFM 0.69

    Parameters
    ----------
    phase: np.array
        Phase data in seconds. Provide either phase or frequency.
    frequency: np.array
        Fractional frequency data (nondimensional). Provide either frequency or phase.
    rate: float
        The sampling rate for phase or frequency, in Hz
    taus: np.array
        Array of tau values for which to compute measurement

    NIST SP 1065 eqn (27) page 25

    """
    if phase == None:
        phase = frequency2phase(frequency, rate)

    rate = float(rate)
    (phase, ms, taus_used) = tau_generator(phase, rate, taus)

    devs    = np.zeros_like(taus_used)
    deverrs = np.zeros_like(taus_used)
    ns      = np.zeros_like(taus_used)

    for idx, mj in enumerate(ms):
        devs[idx], deverrs[idx], ns[idx] = calc_mtotdev_phase(phase, rate, mj)

    return remove_small_ns(taus_used, devs, deverrs, ns)

def calc_mtotdev_phase(phase, rate, m):
    """ PRELIMINARY - REQUIRES FURTHER TESTING.
        calculation of mtotdev for one averaging factor m
        tau = m*tau0

    """
    tau0 = 1.0/rate
    N = int(len(phase)) # phase data, N points

    n=0    # number of terms in the sum, for error estimation
    dev=0.0 # the deviation we are computing
    err=0.0 # the error in the deviation
    #print('calc_mtotdev N=%d m=%d' % (N,m) )
    for i in range(0,N-3*int(m)+1):
        xs = phase[i:i+3*m] # subsequence of length 3m, from the original phase data
        assert( len(xs) == 3*m )
        # remove linear trend. by averaging first/last half, computing slope, and subtracting
        half1_idx =  np.floor(3*m/2.0)
        half2_idx =  np.ceil(3*m/2.0)
        # m
        # 1    0:1   2:2
        mean1 = np.mean( xs[:half1_idx] )
        mean2 = np.mean( xs[half2_idx:] )

        if int(3*m)%2==1: # m is odd
            # 3m = 2k+1 is odd, with the averages at both ends over k points
            # the distance between the averages is then k+1 = (3m-1)/2 +1
            slope = (mean2-mean1) / ( (0.5*(3*m-1)+1)*tau0)
        else: # m is even
            # 3m = 2k is even, so distance between averages is k=m/2
            slope = (mean2-mean1) / (0.5*3*m*tau0)

        x0 = [x - slope*idx*tau0 for (idx,x) in enumerate(xs)]  # remove the linear trend
        x0_flip = x0[::-1] # left-right flipped version of array
        # extended sequence of length 9m, by uninverted even reflection
        xstar = np.concatenate( (x0_flip,x0,x0_flip ))
        assert( len(xstar)==9*m )

        # now compute totdev on these 9m points
        # 6m unique groups of m-point averages, all possible overlapping second differences
        # one term in the 6m sum:  [ x_i - 2 x_i+m + x_i+2m ]^2
        squaresum=0.0
        #print('m=%d 9m=%d maxj+3*m=%d' %( m, len(xstar), 6*int(m)+3*int(m)) )
        for j in range(0,6*int(m)): # summation of the 6m terms.
            xmean1 = np.mean( xstar[j     : j+m  ] )
            xmean2 = np.mean( xstar[j+m   : j+2*m] )
            xmean3 = np.mean( xstar[j+2*m : j+3*m] )
            assert( len(xstar[j : j+m]) == m )
            assert( len(xstar[j+m : j+2*m]) == m )
            assert( len(xstar[j+2*m : j+3*m]) == m )
            squaresum += pow(xmean1-2.0*xmean2+xmean3, 2)

        squaresum = (1.0/(6.0*m)) * squaresum
        dev += squaresum
        n=n+1

    # scaling in front of double-sum
    if not ( n == N-3*int(m)+1 ):
        print('Error: not n == N-3*int(m)+1 !!')
        print(' n= %d  N= %d  m= %d '%(n,N,m))
        print(' N-3m+1 = %d'%(N-3*m+1))
    assert( n == N-3*int(m)+1 ) # sanity check on the number of terms n
    dev = dev* 1.0/ ( 2.0*pow(m*tau0,2)*(N-3*m+1) )
    dev = np.sqrt(dev)
    error = dev / np.sqrt(n)
    return (dev,error,n)


def htotdev(phase=None, frequency=None, rate=1.0, taus=[]):
    """ PRELIMINARY - REQUIRES FURTHER TESTING.
        Hadamard Total deviation.
        Better confidence at long averages for Hadamard deviation

        FIXME: bias corrections from http://www.wriley.com/CI2.pdf
        W FM    0.995      alpha= 0
        F FM    0.851      alpha=-1
        RW FM   0.771      alpha=-2
        FW FM   0.717      alpha=-3
        RR FM   0.679      alpha=-4
        
    Parameters
    ----------
    phase: np.array
        Phase data in seconds. Provide either phase or frequency.
    frequency: np.array
        Fractional frequency data (nondimensional). Provide either frequency or phase.
    rate: float
        The sampling rate for phase or frequency, in Hz
    taus: np.array
        Array of tau values for which to compute measurement
    
    """
    if phase == None:
        phase = frequency2phase(frequency, rate)
        
    rate = float(rate)
    (phase, ms, taus_used) = tau_generator(phase, rate, taus)

    devs    = np.zeros_like(taus_used)
    deverrs = np.zeros_like(taus_used)
    ns      = np.zeros_like(taus_used)
    
    freq = phase2frequency(phase, rate)
    
    # NOTE at mj==1 we use ohdev(), based on comment from here:
    # http://www.wriley.com/paper4ht.htm
    # "For best consistency, the overlapping Hadamard variance is used instead of the Hadamard total variance at m=1"
    for idx, mj in enumerate(ms):
        if mj==1:
            (devs[idx], deverrs[idx], ns[idx]) = calc_hdev_phase(phase, rate, mj, 1)
        else:
            devs[idx], deverrs[idx], ns[idx] = calc_htotdev_freq(freq, rate, mj)

    return remove_small_ns(taus_used, devs, deverrs, ns)

def calc_htotdev_freq(freq, rate, m):
    """ PRELIMINARY - REQUIRES FURTHER TESTING.
        calculation of htotdev for one averaging factor m
        tau = m*tau0
        
        Parameters
        ----------
        frequency: np.array
            Fractional frequency data (nondimensional).
        rate: float
            The sampling rate for frequency, in Hz
        m: int
            Averaging factor. tau = m*tau0, where tau0=1/rate.
    """

    N = int(len(freq)) # frequency data, N points
    m = int(m)
    n=0    # number of terms in the sum, for error estimation
    dev=0.0 # the deviation we are computing
    err=0.0 # the error in the deviation    
    for i in range(0,N-3*int(m)+1):
        xs = freq[i:i+3*m] # subsequence of length 3m, from the original phase data
        assert( len(xs) == 3*m )
        # remove linear trend. by averaging first/last half, computing slope, and subtracting
        half1_idx = int( np.floor(3*m/2.0) )
        half2_idx = int( np.ceil(3*m/2.0)  )
        # m
        # 1    0:1   2:2
        mean1 = np.mean( xs[:half1_idx] ) 
        mean2 = np.mean( xs[half2_idx:] )
        
        if int(3*m)%2==1: # m is odd
            # 3m = 2k+1 is odd, with the averages at both ends over k points
            # the distance between the averages is then k+1 = (3m-1)/2 +1
            slope = (mean2-mean1) / ( (0.5*(3*m-1)+1))
        else: # m is even
            # 3m = 2k is even, so distance between averages is k=3m/2
            slope = (mean2-mean1) / (0.5*3*m)
              
        x0 = [x - slope*(idx-np.floor(3*m/2)) for (idx,x) in enumerate(xs)]  # remove the linear trend
        x0_flip = x0[::-1] # left-right flipped version of array
        # extended sequence of length 9m, by uninverted even reflection
        xstar = np.concatenate( (x0_flip,x0,x0_flip ))
        assert( len(xstar)==9*m )
        
        # now compute totdev on these 9m points
        # 6m unique groups of m-point averages, all possible overlapping second differences
        # one term in the 6m sum:  [ x_i - 2 x_i+m + x_i+2m ]^2
        squaresum=0.0
        k=0
        for j in range(0,6*int(m)): # summation of the 6m terms.
            xmean1 = np.mean( xstar[j+0*m : j+1*m] )
            xmean2 = np.mean( xstar[j+1*m : j+2*m] )
            xmean3 = np.mean( xstar[j+2*m : j+3*m] )
            squaresum += pow(xmean1-2.0*xmean2+xmean3, 2)
            k=k+1
        assert( k == 6*int(m) )
        
        squaresum = (1.0/(6.0*k)) * squaresum
        dev += squaresum
        n=n+1
    
    # scaling in front of double-sum
    assert( n == N-3*int(m)+1 ) # sanity check on the number of terms n
    dev = dev* 1.0/ ( N-3*m+1 ) 
    dev = np.sqrt(dev)
    error = dev / np.sqrt(n)
    return (dev,error,n)

def theo1(phase=None, frequency=None, rate=1.0, taus=[]):
    """ PRELIMINARY - REQUIRES FURTHER TESTING.
        Theo1
        Theo1 is a two-sample variance with improved confidence and 
        extended averaging factor range. 

        NIST SP 1065 eq (30) page 29
        
        see tests/phasedat/phase_dat_test.py for initial test
        
        TODO: bias correction
        
    Parameters
    ----------
    phase: np.array
        Phase data in seconds. Provide either phase or frequency.
    frequency: np.array
        Fractional frequency data (nondimensional). Provide either frequency or phase.
    rate: float
        The sampling rate for phase or frequency, in Hz
    taus: np.array
        Array of tau values for which to compute measurement
    
    """
    if phase == None:
        phase = frequency2phase(frequency, rate)
        
    rate = float(rate)
    tau0 = 1.0/rate
    (phase, ms, taus_used) = tau_generator(phase, rate, taus, even=True)

    devs    = np.zeros_like(taus_used)
    deverrs = np.zeros_like(taus_used)
    ns      = np.zeros_like(taus_used)
    
    N=len(phase)
    for idx, m in enumerate(ms):
        m = int(m) # to avoid: VisibleDeprecationWarning: using a non-integer number instead of an integer will result in an error in the future
        assert( m % 2 == 0 ) # m must be even
        dev=0
        n=0
        for i in range( int(N-m) ):
            s=0
            for d in range( int(m)/2 ): # inner sum
                pre = 1.0 / (float(m)/2 - float(d))
                s += pre*pow( phase[i]-phase[i-d+m/2] + phase[i+m]-phase[i+d+m/2] , 2)
                n=n+1
            dev += s
        assert( n == (N-m)*m/2 ) # N-m outer sums, m/2 inner sums
        dev = dev/( 0.75*(N-m)*pow(m*tau0,2) )
        # factor 0.75 used here? http://tf.nist.gov/general/pdf/1990.pdf
        # but not here? http://tf.nist.gov/timefreq/general/pdf/2220.pdf page 29
        devs[idx] = np.sqrt( dev )
        deverrs[idx] = devs[idx] / np.sqrt(N-m)
        ns[idx] = n
        
    return remove_small_ns(taus_used, devs, deverrs, ns)
    

def tierms(phase=None, frequency=None, rate=1.0, taus=[]):
    """ Time Interval Error RMS.

    Parameters
    ----------
    phase: np.array
        Phase data in seconds. Provide either phase or frequency.
    frequency: np.array
        Fractional frequency data (nondimensional). Provide either frequency or phase.
    rate: float
        The sampling rate for phase or frequency, in Hz
    taus: np.array
        Array of tau values for which to compute measurement

    """
    if phase == None:
        phase = frequency2phase(frequency, rate)

    rate = float(rate)
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
        deverrs[idx] = 0 / np.sqrt(ncount) # TODO! I THINK THIS IS WRONG!
        ns[idx] = ncount

    return remove_small_ns(taus_used, devs, deverrs, ns)

def mtie_rolling_window(a, window):
    """
    Make an ndarray with a rolling window of the last dimension
    from http://mail.scipy.org/pipermail/numpy-discussion/2011-January/054401.html

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

    """
    if window < 1:
        raise ValueError("`window` must be at least 1.")
    if window > a.shape[-1]:
        raise ValueError("`window` is too long.")
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)

def mtie(phase=None, frequency=None, rate=1.0, taus=[]):
    """ Maximum Time Interval Error.

    Parameters
    ----------
    phase: np.array
        Phase data in seconds. Provide either phase or frequency.
    frequency: np.array
        Fractional frequency data (nondimensional). Provide either frequency or phase.
    rate: float
        The sampling rate for phase or frequency, in Hz

    Notes
    -----
    this seems to correspond to Stable32 setting "Fast(u)"
    Stable32 also has "Decade" and "Octave" modes where the
    dataset is extended somehow?
    """

    if phase == None:
        phase = frequency2phase(frequency, rate)

    rate = float(rate)
    (phase, m, taus_used) = tau_generator(phase, rate, taus)
    devs = np.zeros_like(taus_used)
    deverrs = np.zeros_like(taus_used)
    ns = np.zeros_like(taus_used)

    for idx, mj in enumerate(m):
        rw = mtie_rolling_window(phase, mj + 1)
        win_max = np.max(rw, axis=1)
        win_min = np.min(rw, axis=1)
        tie = win_max - win_min
        dev = np.max(tie)
        ncount = phase.shape[0] - mj
        devs[idx] = dev
        deverrs[idx] = dev / np.sqrt(ncount)
        ns[idx] = ncount

    return remove_small_ns(taus_used, devs, deverrs, ns)

#
# !!!!!!!
# FIXME: mtie_phase_fast() is incomplete.
# !!!!!!!
#
def mtie_phase_fast(phase, rate, taus):
    """ fast binary decomposition algorithm for MTIE

        See: STEFANO BREGNI "Fast Algorithms for TVAR and MTIE Computation in
        Characterization of Network Synchronization Performance"
    """
    rate = float(rate)
    phase = np.asarray(phase)
    k_max = int( np.floor( np.log2( len(phase) ) ) )
    phase = phase[0:pow(2,k_max)] # truncate data to 2**k_max datapoints
    assert( len(phase) == pow(2,k_max) )
    k = 1
    taus = []
    while k<=k_max:
        tau = pow(2,k)
        taus.append(tau)
        #print tau
        k += 1
    print("taus ", taus)
    devs = np.zeros(len(taus))
    deverrs = np.zeros(len(taus))
    ns = np.zeros(len(taus))
    # matrices to store results
    mtie_max = np.zeros( (len(phase)-1,k_max) )
    mtie_min = np.zeros( (len(phase)-1,k_max) )
    for kidx in range(k_max):
        k=kidx+1
        imax = len(phase)-pow(2,k)+1
        #print k, imax
        tie = np.zeros(imax)
        #print np.max( tie )
        for i in range(imax):
            if k==1:
                mtie_max[i,kidx] = max( phase[i], phase[i+1] )
                mtie_min[i,kidx] = min( phase[i], phase[i+1] )
            else:
                p = int( pow(2,k-1) )
                mtie_max[i,kidx] = max( mtie_max[i,kidx-1], mtie_max[i+p,kidx-1] )
                mtie_min[i,kidx] = min( mtie_min[i,kidx-1], mtie_min[i+p,kidx-1] )

        #for i in range(imax):
            tie[i] = mtie_max[i,kidx] - mtie_min[i,kidx]
            #print tie[i]
        devs[kidx] = np.amax( tie )
        #print "maximum %2.4f" % devs[kidx]
        #print np.amax( tie )
    #for tau in taus:
    #for

    print(devs)
    #print k_max
    #devs =

########################################################################
#
#  gap resistant Allan deviation
#

def gradev(phase=None, frequency=None, rate=1.0, taus=[], ci=0.9, noisetype='wp'):
    """ gap resistant overlapping Allan deviation

    """
    if phase == None:
        frequency= trim_data(frequency)
        phase = frequency2phase(frequency, rate)
        
    (data, m, taus_used) = tau_generator(phase, rate, taus)

    ad  = np.zeros_like(taus_used)
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
        ad[idx]  = dev
        ade_l[idx] = deverr[0]
        ade_h[idx] = deverr[1]
        adn[idx] = n

    return remove_small_ns(taus_used, ad, ade_l, ade_h, adn)

def calc_gradev_phase(data, rate, mj, stride, ci, noisetype):
    """ see http://www.leapsecond.com/tools/adev_lib.c
        stride = mj for nonoverlapping allan deviation
        stride = 1 for overlapping allan deviation

        see http://en.wikipedia.org/wiki/Allan_variance
             1       1
         s2y(t) = --------- sum [x(i+2) - 2x(i+1) + x(i) ]^2
                  2*tau^2
    """

    d2 = data[2 * mj::stride]
    d1 = data[1 * mj::stride]
    d0 = data[::stride]

    n = min(len(d0), len(d1), len(d2))

    if n == 0:
        RuntimeWarning("Data array length is too small: %i" % len(data))
        n = 1

    v_arr = d2[:n] - 2 * d1[:n] + d0[:n]

    n = len(np.where(np.isnan(v_arr)==False)[0]) # only average for non-nans
    N = len(np.where(np.isnan(data)==False)[0])

    s = np.nansum(v_arr * v_arr)   #  a summation robust to nans

    dev = np.sqrt(s / (2.0 * n)) / mj  * rate
    #deverr = dev / np.sqrt(n) # old simple errorbars
    if n > 1:
        deverr = uncertainty_estimate(N,
                                      mj,
                                      dev,
                                      ci,
                                      noisetype)
    else:
        deverr = [0,0]

    return dev, deverr, n


########################################################################
#
#  Various helper functions and utilities
#


def tau_generator(data, rate, taus=[], v=False, even=False):
    """ pre-processing of the tau-list given by the user (Helper function)

    Does sanity checks, sorts data, removes duplicates and invalid values.
    Generates a tau-list based on keywords 'all', 'decade', 'octave'.
    Uses 'octave' by default if no taus= argument is given.

    Parameters
    ----------
    data: np.array
        data array
    rate: float
        Sample rate of data in Hz. Time interval between measurements is 1/rate seconds.
    taus: np.array
        Array of tau values for which to compute measurement.
        Alternatively one of the keywords: "all", "octave", "decade".
        Defaults to 'octave' if omitted.
    v:
        verbose output if True
    even:
        require even m, where tau=m*tau0, for Theo1 statistic

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

    if not taus: # empty or no tau-list supplied
        taus='octave' # default to octave
    
    if taus == "all":
        taus = (1.0/rate)*np.linspace(1.0,len(data),len(data))
    elif taus == "octave":
        maxn = np.floor( np.log2( len(data) ) )
        taus = (1.0/rate)*np.logspace(0,maxn,maxn+1,base=2.0)
    elif taus == "decade":
        maxn = np.floor( np.log10( len(data) ) )
        taus = (1.0/rate)*np.logspace(0,maxn,maxn+1,base=10.0)
    
    data, taus = np.array(data), np.array(taus)
    rate = float(rate)
    m = [] # integer averaging factor. tau = m*tau0
    
    taus_valid1 = taus < (1 / float(rate)) * float(len(data))
    taus_valid2 = taus > 0
    taus_valid  = taus_valid1 & taus_valid2
    m = np.floor(taus[taus_valid] * rate)
    m = m[m != 0]       # m is tau in units of datapoints
    m = np.unique(m)    # remove duplicates and sort

    if v:
        print("tau_generator: ", m)

    if len(m) == 0:
        print("Warning: sanity-check on tau failed!")
        print("   len(data)=", len(data), " rate=", rate, "taus= ", taus)

    taus2 = m / float(rate)
    
    if even:
        m_even = ((m % 2) == 0)
        m = m[m_even]
        taus2 = taus2[m_even]
        
    return data, m, taus2

def remove_small_ns(*args):
    if len(args) == 4:
        return remove_small_ns_4(*args)
    elif len(args) == 5:
        return remove_small_ns_5(*args)
    else:
        assert(0)

# 4-parameter version of this function
def remove_small_ns_4(taus, devs, deverrs, ns):
    """ if n is small (==1), reject the result

        n is the number of averages in the deviation estimate at a
        certain tau value.
    """

    ns_big_enough = ns > 1

    o_taus = taus[ns_big_enough]
    o_dev  = devs[ns_big_enough]
    o_err  = deverrs[ns_big_enough]
    o_n    = ns[ns_big_enough]

    return o_taus, o_dev, o_err, o_n

# FIXME: 5-parameter version of the exact same function as above.
# try to merge both into one function.
# This function handles low-side and high-side errors separately
def remove_small_ns_5(taus, devs, deverrs_l, deverrs_h, ns):
    """ if n is small (==1), reject the result """
    ns_big_enough = ns > 1
    o_taus = taus[ns_big_enough]
    o_dev  = devs[ns_big_enough]
    o_err_l  = deverrs_l[ns_big_enough]
    o_err_h  = deverrs_h[ns_big_enough]
    o_n    = ns[ns_big_enough]
    return o_taus, o_dev, o_err_l, o_err_h, o_n

def trim_data(x):
    """
    Trim leading and trailing NaNs from dataset:
    This is done by creating a boolean mask that is False for all leading and
    trailing NaN values. The mask is True for all values from the first non-NaN
    to the last non-NaN. The mask is then applied to the data set and returned.
    """
    # create boolean mask default True:
    mask = np.ones(len(x), dtype=bool)
    # Set indices of mask to False for leading NaNs
    for i in range(0,len(x)):
        if np.isnan(x[i]) == True:
            mask[i] = False
        else:
            break

    # Set indices of mask to False for trailing NaNs
    for j in range(len(x)-1,0,-1):
        if np.isnan(x[j]) == True:
            mask[j] = False
        else:
            break

    return x[mask]

def uncertainty_estimate(N, m, s, ci=0.9, noisetype='wp'):
    """Determine the uncertainty of a given two-sample variance estimate for
    a given number of samples (N), sampling distance (m), variance estimate
    (s), confindence interval (ci), and type of noise (noisetype).

    Parameters
    ----------
    N : int
        the number of samples used in a two sample variance estimate
    m : int
        the length of the window, tau = m * tau0
    s : float
        the estimated two-sample variance
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

    Notes
    -----
       S. Stein, Frequency and Time - Their Measurement and
       Characterization. Precision Frequency Control Vol 2, 1985, pp 191-416.
       http://tf.boulder.nist.gov/general/pdf/666.pdf

    Returns
    -------
    [err_l, err_h] : list, float
        the upper and lower bounds of the confidence interval taken as
        distances from the the estimated two sample variance.

    """

    ci_l = min(np.abs(ci),np.abs((ci-1))) / 2
    ci_h = 1 - ci_l

    if noisetype in set(['wp', 'fp', 'wf', 'ff', 'rf']):

        if noisetype =='wp':
            df = (N + 1) * (N - 2*m) / (2 * (N - m))

        if noisetype =='wf':
            df = (((3 * (N - 1) / (2 * m)) - (2 * (N - 2) / N)) *
                  ((4*m**2) / ((4*m**2) + 5)))

        if noisetype =='fp':
            a = (N - 1)/(2 * m)
            b = (2 * m + 1) * (N - 1) / 4
            df = np.exp(np.sqrt(np.log(a) * np.log(b)))

        if noisetype =='ff':
            if m == 1:
                df = 2 * (N - 2) /(2.3 * N - 4.9)
            if m >= 2:
                df = 5 * N**2 / (4 * m * (N + (3 * m)))

        if noisetype == 'rf':
            a = (N - 2) / (m * (N - 3)**2)
            b = (N - 1)**2
            c = 3 * m * (N - 1)
            d = 4 * m **2
            df = a * (b - c + d)

    else:
        df = (N - 1)
        print("Noise type not recognized. Defaulting to N - 1 degrees of freedom.")

    chi2_l = scipy.stats.chi2.ppf(ci_l,df)
    chi2_h = scipy.stats.chi2.ppf(ci_h,df)

    err_h = np.abs(df * s / chi2_l - s)
    err_l = np.abs(df * s / chi2_h - s)

#    print N, m, s, df, chi2_l, err_h

    return [err_l, err_h]

def three_cornered_hat_phase(phasedata_ab, phasedata_bc,
                             phasedata_ca, rate, taus, function):
    """
    Three Cornered Hat Method

    We have three clocks with unknown variances sa^2, sb^2, sc^3
    Three pairwise measurements give variances:
    sab^2, sbc^2, sca^2
    Assuming covariances are zero (clocks are independent), we get:
    sa^2 = 0.5*( sab^2 + sca^2 - sbc^2 )
    (and cyclic permutations for sb and sc)

    References
    ----------
    http://www.wriley.com/3-CornHat.htm
    """
    # Until MTIE stuff is ported, need this fix:
    npa = np.array
    phasedata_ab, phasedata_bc, phasedata_ca = npa(phasedata_ab), npa(phasedata_bc), npa(phasedata_ca)
    taus = npa(taus)

    (tau_ab, dev_ab, err_ab, ns_ab) = function(phase=phasedata_ab, rate=rate, taus=taus)
    (tau_bc, dev_bc, err_bc, ns_bc) = function(phase=phasedata_bc, rate=rate, taus=taus)
    (tau_ca, dev_ca, err_ca, ns_ca) = function(phase=phasedata_ca, rate=rate, taus=taus)

    (tau_ab, dev_ab, err_ab, ns_ab) = npa(tau_ab), npa(dev_ab), npa(err_ab), npa(ns_ab)
    (tau_bc, dev_bc, err_bc, ns_bc) = npa(tau_bc), npa(dev_bc), npa(err_bc), npa(ns_bc)
    (tau_ca, dev_ca, err_ca, ns_ca) = npa(tau_ca), npa(dev_ca), npa(err_ca), npa(ns_ca)

    var_ab = dev_ab * dev_ab
    var_bc = dev_bc * dev_bc
    var_ca = dev_ca * dev_ca
    assert len(var_ab) == len(var_bc) == len(var_ca)
    var_a = 0.5 * (var_ab + var_ca - var_bc)

    dev_a = np.sqrt(var_a)
    dev_a[var_a < 0] = 0

    return tau_ab, dev_a

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
    phasedata = np.cumsum(freqdata) * dt
    phasedata = np.insert(phasedata, 0, 0) # FIXME: why do we do this? so that phase starts at zero and len(phase)=len(freq)+1 ??
    return phasedata

def phase2radians(phasedata,v0):
    """ Convert phase in seconds to phase in radians

    Parameters
    ----------
    phasedata: np.array
        Data array of phase in seconds
    v0: float
        Nominal oscillator frequency in Hz

    Returns
    -------
    fi:
        phase data in radians
    """
    fi = [2*math.pi*v0*xx for xx in phasedata]
    return fi

def phase2frequency(phase, rate):
    """ Convert phase in seconds to fractional frequency

    Parameters
    ----------
    phase: np.array
        Data array of phase in seconds
    rate: float
        The sampling rate for phase, in Hz

    Returns
    -------
    y:
        Data array of fractional frequency
    """
    y = rate*np.diff(phase)
    return y

if __name__ == "__main__":
    print("Nothing to see here.")
