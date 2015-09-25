"""
Allan deviation tools
=====================

**Author:** Anders Wallin (anders.e.e.wallin "at" gmail.com)

Version history
---------------

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
- http://www.anderswallin.net/2014/08/faster-allantools-with-numpy/

**v1.01** 2014 August
- PEP8 compliance improvements by Danny Price.

**v1.00** 2014 January, first version of allantools.
- http://www.anderswallin.net/2014/01/allantools/

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

import numpy as np
import scipy.stats # used in uncertainty_estimate()

def tdev_phase(phase, rate, taus):
    """ Time deviation of phase data
    
    TDEV = (tau / sqrt(3) ) * MDEV

    NIST SP 1065 eqn (15)
    
    Parameters
    ----------
    phase: np.array
        list of phase measurements in seconds
    rate: float
        Sample rate of data, i.e. interval between measurements is 1/rate (Hz)
    taus: np.array
        Array of tau values for which to compute measurement

    Returns
    -------
    (taus2, td, tde, ns): tuple  
          Tuple of values
    taus2: np.array
        Tau values for which td computed
    td: np.array
        Computed time deviations for each tau value
    tde: np.array
        Time deviation errors
    ns: np.array
        Values of N used in mdev_phase()
    """
    (taus2, md, mde, ns) = mdev_phase(phase, rate, taus)

    td = taus2 * md / np.sqrt(3.0)
    tde = td / np.sqrt(ns)
    return taus2, td, tde, ns


def tdev(data, rate, taus):
    """ Time deviation of fractional frequency data
    
    Parameters
    ----------
    data: np.array
        Data arra
    rate: float
        The sampling rate, in Hz
    taus: np.array
        Array of tau values for which to compute allan variance
    
    Returns
    -------
    (taus2, td, tde, ns): tuple  
          Tuple of values
    taus2: np.array
        Tau values for which td computed
    td: np.array
        Computed time deviations for each tau value
    tde: np.array
        Time deviation errors
    ns: np.array
        Values of N used in mdev_phase()
    
    Notes
    -----
    http://en.wikipedia.org/wiki/Time_deviation 
    """
    phase = frequency2phase(data, rate)
    return tdev_phase(phase, rate, taus)


def mdev_phase(data, rate, taus):
    """  Modified Allan deviation of phase data

    Parameters
    ----------
    phase: np.array
        list of phase measurements in seconds
    rate: float
        Sample rate of data, i.e. interval between measurements is 1/rate (Hz)
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
    data, taus = np.array(data), np.array(taus)
    (data, ms, taus_used) = tau_m(data, rate, taus)
    # taus = []
    md    = np.zeros_like(ms)
    mderr = np.zeros_like(ms)
    ns    = np.zeros_like(ms)

    # this is a 'loop-unrolled' algorithm following
    # http://www.leapsecond.com/tools/adev_lib.c 
    for idx, m in enumerate(ms):
        tau = taus_used[idx]

        # First loop sum
        d0 = data[0:m]
        d1 = data[m:2*m]
        d2 = data[2*m:3*m]
        e = min(len(d0), len(d1), len(d2))
        v = np.sum(d2[:e] - 2* d1[:e] + d0[:e])
        s = v * v

        # Second part of sum
        d3 = data[3*m:]
        d2 = data[2*m:]
        d1 = data[1*m:]
        d0 = data[0:]

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


def mdev(freqdata, rate, taus):
    """ modified Allan deviation, fractional frequency data 
    
    Parameters
    ----------
    freqdata: np.array
        Data array of fractional frequency measurements (nondimensional)
    rate: float
        Sample rate of data, i.e. interval between measurements is 1/rate (Hz)
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
        
    """
    phase = frequency2phase(freqdata, rate)
    return mdev_phase(phase, rate, taus)



def adev(data, rate, taus):
    """ Allan deviation for fractional frequency data

    Parameters
    ----------
    data: np.array
        Data time-series of evenly spaced fractional frequency measurements
    rate: float
        The samples/s in the time-series
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
    phase = frequency2phase(data, rate)
    return adev_phase(phase, rate, taus)


def adev_phase(data, rate, taus):
    """ Allan deviation for phase data 
    
    Parameters
    ----------
    data: np.array
        list of phase measurements in seconds
    rate: float
        Sample rate of data, i.e. interval between measurements is 1/rate (Hz)
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
    (data, m, taus_used) = tau_m(data, rate, taus)

    ad  = np.zeros_like(taus_used)
    ade = np.zeros_like(taus_used)
    adn = np.zeros_like(taus_used)

    for idx, mj in enumerate(m):  # loop through each tau value m(j)
        (ad[idx], ade[idx], adn[idx]) = adev_phase_calc(data, rate, mj, mj)

    return remove_small_ns(taus_used, ad, ade, adn)  # tau, adev, adeverror, naverages


def adev_phase_calc(data, rate, mj, stride):
    """ see http://www.leapsecond.com/tools/adev_lib.c
        stride = mj for nonoverlapping allan deviation
    Parameters
    ----------
    data: np.array
        list of phase measurements in seconds
    rate: float
        Sample rate of data, i.e. interval between measurements is 1/rate (Hz)
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
    stride = mj for nonoverlapping allan deviation
    stride = 1 for overlapping allan deviation

    References
    ----------        
    * http://en.wikipedia.org/wiki/Allan_variance
    * http://www.leapsecond.com/tools/adev_lib.c
    NIST SP 1065, eqn (7) and (11)
    """

    d2 = data[2 * mj::stride]
    d1 = data[1 * mj::stride]
    d0 = data[::stride]

    n = min(len(d0), len(d1), len(d2))

    if n == 0:
        RuntimeWarning("Data array length is too small: %i" % len(data))
        n = 1

    v_arr = d2[:n] - 2 * d1[:n] + d0[:n]
    s = np.sum(v_arr * v_arr)

    dev = np.sqrt(s / (2.0 * n)) / mj  * rate
    deverr = dev / np.sqrt(n)

    return dev, deverr, n

def oadev_phase(data, rate, taus):
    """ overlapping Allan deviation of phase data 
    
    Parameters
    ----------
    data: np.array
        list of phase measurements in seconds
    rate: float
        Sample rate of data, i.e. interval between measurements is 1/rate (Hz)
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
    (data, m, taus_used) = tau_m(data, rate, taus)
    ad  = np.zeros_like(taus_used)
    ade = np.zeros_like(taus_used)
    adn = np.zeros_like(taus_used)

    for idx, mj in enumerate(m):
        (ad[idx], ade[idx], adn[idx]) = adev_phase_calc(data, rate, mj, 1)  # stride=1 for overlapping ADEV

    return remove_small_ns(taus_used, ad, ade, adn)  # tau, adev, adeverror, naverages


def oadev(freqdata, rate, taus):
    """ overlapping Allan deviation for fractional frequency data 
    
    Parameters
    ----------
    freqdata: np.array
        Data array of fractional frequency measurements (nondimensional)
    rate: float
        Sample rate of data, i.e. interval between measurements is 1/rate (Hz)
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
    phase = frequency2phase(freqdata, rate)
    return oadev_phase(phase, rate, taus)

def ohdev(freqdata, rate, taus):
    """ Overlapping Hadamard deviation, fractional frequency data 
    
    Parameters
    ----------
    freqdata: np.array
        Data array of fractional frequency measurements (nondimensional)
    rate: float
        Sample rate of data, i.e. interval between measurements is 1/rate (Hz)
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
    phase = frequency2phase(freqdata, rate)
    return ohdev_phase(phase, rate, taus)

def ohdev_phase(data, rate, taus):
    """ Overlapping Hadamard deviation of phase data 
    
    Parameters
    ----------
    data: np.array
        list of phase measurements in seconds
    rate: float
        Sample rate of data, i.e. interval between measurements is 1/rate (Hz)
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
    rate = float(rate)
    (data, m, taus_used) = tau_m(data, rate, taus)
    hdevs = np.zeros_like(taus_used)
    hdeverrs = np.zeros_like(taus_used)
    ns = np.zeros_like(taus_used)

    for idx, mj in enumerate(m):
        (hdevs[idx], hdeverrs[idx], ns[idx]) = calc_hdev_phase(data, rate, mj, 1)  # stride = 1

    return remove_small_ns(taus_used, hdevs, hdeverrs, ns)

def hdev(freqdata, rate, taus):
    """ Hadamard deviation, fractional frequency data 
    
    Parameters
    ----------
    freqdata: np.array
        Data array of fractional frequency measurements (nondimensional)
    rate: float
        Sample rate of data, i.e. interval between measurements is 1/rate (Hz)
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
    phase = frequency2phase(freqdata, rate)
    return hdev_phase(phase, rate, taus)


def hdev_phase(data, rate, taus):
    """ Hadamard deviation, phase data 
    
    Parameters
    ----------
    data: np.array
        list of phase measurements in seconds
    rate: float
        Sample rate of data, i.e. interval between measurements is 1/rate (Hz)
    taus: np.array
        Array of tau values for which to compute measurement
    """
    rate = float(rate)
    (data, m, taus_used) = tau_m(data, rate, taus)
    hdevs = np.zeros_like(taus_used)
    hdeverrs = np.zeros_like(taus_used)
    ns = np.zeros_like(taus_used)

    for idx, mj in enumerate(m):
        hdevs[idx], hdeverrs[idx], ns[idx] = calc_hdev_phase(data, rate, mj, mj)  # stride = mj

    return remove_small_ns(taus_used, hdevs, hdeverrs, ns)


def calc_hdev_phase(data, rate, mj, stride):
    """ calc_hdev_phase helper function for hdev() ad hdev_phase()

    Parameters
    ----------
    data: np.array
        list of phase measurements in seconds
    rate: float
        Sample rate of data, i.e. interval between measurements is 1/rate (Hz)
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

    NIST SP 1065 eqn (18) and (20)
    """

    tau0 = 1.0 / float(rate)

    d3 = data[3 * mj::stride]
    d2 = data[2 * mj::stride]
    d1 = data[1 * mj::stride]
    d0 = data[::stride]

    n = min(len(d0), len(d1), len(d2), len(d3))

    v_arr = d3[:n] - 3 * d2[:n] + 3 * d1[:n] - d0[:n]

    s = np.sum(v_arr * v_arr)

    if n == 0:
        n = 1

    h = np.sqrt(s / 6.0 / float(n)) / float(tau0 * mj)
    e = h / np.sqrt(n)
    return h, e, n


def totdev(freqdata, rate, taus):
    """ Total deviation, fractional frequency data """
    phasedata = frequency2phase(freqdata, rate)
    return totdev_phase(phasedata, rate, taus)


def totdev_phase(data, rate, taus):
    """ Total deviation, phase data.
    
    Parameters
    ----------
    data: np.array
        Data array of fractional frequency measurements (nondimensional)
    rate: float
        Sample rate of data, i.e. interval between measurements is 1/rate (Hz)
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
    totvar(t) = ------------  sum   [ x*(i-m) - 2x*(i)+x*(i+m) ]**2
                2 t**2 (N-2)  i=2
                
    where x* is a new dataset with 'reflected' data at start/end.
    
    the original data x is in the center of x*:
    x*(1-j) = 2x(1) - x(1+j)  for j=1..N-2
    x*(i)   = x(i)            for i=1..N
    x*(N+j) = 2x(N) - x(N-j)  for j=1..N-2
    x* has length 3N-4
    tau = m*tau0
    
    NIST SP 1065 eqn (25)
    
    """

    rate = float(rate)
    (data, m, taus_used) = tau_m(data, rate, taus)
    n = len(data)

    # totdev requires a new dataset
    # Begin by adding reflected data before dataset
    x1 = 2.0 * data[0] * np.ones((n - 2,))
    x1 = x1 - data[1:-1]
    x1 = x1[::-1]

    # Reflected data at end of dataset
    x2 = 2.0 * data[-1] * np.ones((n - 2,))
    x2 = x2 - data[1:-1][::-1]

    # Combine into a single array
    x = np.zeros((3*n - 4))
    x[0:n-2] = x1
    x[n-2:2*(n-2)+2] = data
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


def tierms(freqdata, rate, taus):
    """ Time Interval Error RMS, for fractional frequency data """
    phasedata = frequency2phase(freqdata, rate)
    return tierms_phase(phasedata, rate, taus)


def tierms_phase(phase, rate, taus):
    """ Time Interval Error RMS, for phase data 
    
    Parameters
    ----------
    data: np.array
        Data array of fractional frequency measurements (nondimensional)
    rate: float
        Sample rate of data, i.e. interval between measurements is 1/rate (Hz)
    taus: np.array
        Array of tau values for which to compute measurement
    
    """
    rate = float(rate)
    (data, m, taus_used) = tau_m(phase, rate, taus)

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


def mtie(freqdata, rate, taus):
    """ Maximum Time Interval Error, fractional frequency data """
    phasedata = frequency2phase(freqdata, rate)
    return mtie_phase(phasedata, rate, taus)

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
    
def mtie_phase(phase, rate, taus):
    """ Maximum Time Interval Error, phase data
    
    this seems to correspond to Stable32 setting "Fast(u)"
    Stable32 also has "Decade" and "Octave" modes where the dataset is extended somehow?
    """
    rate = float(rate)
    (phase, m, taus_used) = tau_m(phase, rate, taus)
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

def gradev_phase(data, rate, taus, ci=0.9, noisetype='wp'):
    """ gap resistant overlapping Allan deviation of phase data 
    """
    (data, m, taus_used) = tau_m(data, rate, taus)
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
                                             noisetype)  # stride=1 for overlapping ADEV
        ad[idx]  = dev
        ade_l[idx] = deverr[0]
        ade_h[idx] = deverr[1]    
        adn[idx] = n

    return remove_small_ns(taus_used, ad, ade_l, ade_h, adn)

def gradev(freqdata, rate, taus, ci=0.9, noisetype='wp'):
    """ gap-resistant overlapping Allan deviation for fractional frequency data 
    """
    freqdata = trim_data(freqdata)     
    phase = frequency2phase(freqdata, rate)
    return gradev_phase(phase, rate, taus, ci=ci, noisetype=noisetype)

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


def tau_m(data, rate, taus, v=False):
    """ pre-processing of the tau-list given by the user (Helper function)

    Does sanity checks, sorts data, removes duplicates and invalid values.

    Parameters
    ----------
    data: np.array
        data array
    rate: float
        Sample rate of data, i.e. interval between measurements is 1/rate (Hz)
    taus: np.array
        Array of tau values for which to compute measurement

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
    data, taus = np.array(data), np.array(taus)

    if rate == 0:
        raise RuntimeError("Warning! rate==0")
    rate = float(rate)
    # n = len(data) # not used
    m = []

    taus_valid1 = taus < (1 / float(rate)) * float(len(data))
    taus_valid2 = taus > 0
    taus_valid  = taus_valid1 & taus_valid2
    m = np.floor(taus[taus_valid] * rate)
    m = m[m != 0]       # m is tau in units of datapoints
    m = np.unique(m)    # remove duplicates and sort

    if v:
        print("tau_m: ", m)

    if len(m) == 0:
        print("Warning: sanity-check on tau failed!")
        print("   len(data)=", len(data), " rate=", rate, "taus= ", taus)

    taus2 = m / float(rate)
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
    
    Returns
    -------
    [err_l, err_h] : list, float
        the upper and lower bounds of the confidence interval taken as
        distances from the the estimated two sample variance. 

    References
    ----------
       S. Stein, Frequency and Time - Their Measurement and 
       Characterization. Precision Frequency Control Vol 2, 1985, pp 191-416.
       http://tf.boulder.nist.gov/general/pdf/666.pdf
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

def three_cornered_hat_phase(phasedata_ab, phasedata_bc, phasedata_ca, rate, taus, function):
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

    (tau_ab, dev_ab, err_ab, ns_ab) = function(phasedata_ab, rate, taus)
    (tau_bc, dev_bc, err_bc, ns_bc) = function(phasedata_bc, rate, taus)
    (tau_ca, dev_ca, err_ca, ns_ca) = function(phasedata_ca, rate, taus)

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
        Sample rate of data, i.e. interval between measurements is 1/rate (Hz)
            
    Returns
    -------
    phasedata: np.array
        Time integral of fractional frequency data, i.e. phase (time) data in units of seconds.
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

def phase2frequency(phasedata, rate):
    """ Convert phase in seconds fractional frequency
    
    Parameters
    ----------
    phasedata: np.array
        Data array of phase in seconds
    rate: float
        Sample rate in Hz
            
    Returns
    -------
    y: 
        Data array of fractional frequency
    """
    y = rate*np.diff(phasedata)
    return y

if __name__ == "__main__":
    print("Nothing to see here.")
    
    # code to test mtie_phase_fast, incomplete!
    Nmax = pow(2,8)
    phase = 0.2345* np.random.randn(Nmax)
    taus = []
    rate = 1.0
    mtie_phase_fast(phase, rate, taus)
    # then try using old function
    (o_taus, o_dev, o_err, o_n)=mtie_phase(phase, rate, [1,3,7,16,32,64,128,255])
    print(o_taus)
    print(o_dev)
    
