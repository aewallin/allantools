"""
Allan deviation tools
=====================

**Author:** Anders Wallin (anders.e.e.wallin "at" gmail.com)

Version history
---------------

**v1.10** 2014 August, using numpy which is 100x faster than pure python

**v1.01** 2014 August, PEP8 compliance improvements by telegraphic

**v1.00** 2014 January

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
import allantools_pure_python as alpp
try:
    import bottleneck as bn
    USE_BOTTLENECK = True
except ImportError:
    USE_BOTTLENECK = False

def tdev_phase(phase, rate, taus):
    """ Time deviation of phase data

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
    """
    data, taus = np.array(data), np.array(taus)
    (data, ms, taus_used) = tau_m(data, rate, taus)
    # taus = []
    md    = np.zeros_like(ms)
    mderr = np.zeros_like(ms)
    ns    = np.zeros_like(ms)

    idx = 0
    for m in ms:
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

        idx += 1

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
        print "tau_m: ", m

    if len(m) == 0:
        print "Warning: sanity-check on tau failed!"
        print "   len(data)=", len(data), " rate=", rate, "taus= ", taus

    taus2 = m / float(rate)
    return data, m, taus2


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
    idx = 0

    for mj in m:  # loop through each tau value m(j)
        (dev, deverr, n) = calc_adev_phase(data, rate, mj, mj)
        ad[idx] = dev
        ade[idx] = deverr
        adn[idx] = n
        idx += 1

    return remove_small_ns(taus_used, ad, ade, adn)  # tau, adev, adeverror, naverages


def calc_adev_phase(data, rate, mj, stride):
    """" Calculate allan deviation (helper function, use adev_phase)
    
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


def remove_small_ns(taus, devs, deverrs, ns):
    """ if n is small (==1), reject the result """

    ns_big_enough = ns > 1

    o_taus = taus[ns_big_enough]
    o_dev  = devs[ns_big_enough]
    o_err  = deverrs[ns_big_enough]
    o_n    = ns[ns_big_enough]

    return o_taus, o_dev, o_err, o_n


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
    idx = 0

    for mj in m:
        (dev, deverr, n) = calc_adev_phase(data, rate, mj, 1)  # stride=1 for overlapping ADEV
        ad[idx]  = dev
        ade[idx] = deverr
        adn[idx] = n
        idx += 1

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
        Integrated fractional frequency data, i.e. phase data.
    """
    dt = 1.0 / float(rate)
    phasedata = np.cumsum(freqdata) * dt
    phasedata = np.insert(phasedata, 0, 0)
    return phasedata


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
    idx = 0

    for mj in m:
        (h, e, n) = calc_hdev_phase(data, rate, mj, 1)  # stride = 1
        hdevs[idx] = h
        hdeverrs[idx] = e
        ns[idx] = n
        idx += 1
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

    idx = 0
    for mj in m:
        h, e, n = calc_hdev_phase(data, rate, mj, mj)  # stride = mj
        hdevs[idx] = h
        hdeverrs[idx] = e
        ns[idx] = n
        idx += 1

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
    
    """

    rate = float(rate)
    (data, m, taus_used) = tau_m(data, rate, taus)
    n = len(data)

    # totdev requires a new dataset
    # Begin bu adding reflected data before dataset
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

    idx = 0
    for mj in m:

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

        idx += 1

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

    idx = 0
    for mj in m:
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

        idx += 1

    return remove_small_ns(taus_used, devs, deverrs, ns)


def mtie(freqdata, rate, taus):
    """ Maximum Time Interval Error, fractional frequency data """
    phasedata = frequency2phase(freqdata, rate)
    return mtie_phase(phasedata, rate, taus)

def mtie_phase(phase, rate, taus):
    """ maximum time interval error
    
    Parameters
    ----------
    phase: np.array
        Data array of fractional frequency measurements (nondimensional)
    rate: float
        Sample rate of data, i.e. interval between measurements is 1/rate (Hz)
    taus: np.array
        Array of tau values for which to compute measurement
    
    Notes
    -----
    * The fast version of MTIE_PHASE requires bottleneck to be installed. As this isn't
    the most widespread code, fall back to the pure python version if not installed.
    
    * this seems to correspond to Stable32 setting "Fast(u)"
    Stable32 also has "Decade" and "Octave" modes where the dataset is extended somehow?
    rate = float(rate)
    """
    if USE_BOTTLENECK:
        return mtie_phase_fast(phase, rate, taus)
    else:
        return alpp.mtie_phase(phase, rate, taus)


def mtie_phase_fast(phase, rate, taus):
    """ Bottleneck accelerated version of mtie_phase """

    (phase, m, taus_used) = tau_m(phase, rate, taus)

    devs = np.zeros_like(taus_used)
    deverrs = np.zeros_like(taus_used)
    ns = np.zeros_like(taus_used)

    idx = 0
    for mj in m:
        win_size = mj + 1
        win_max = bn.move_nanmax(phase, window=win_size)
        win_min = bn.move_nanmin(phase, window=win_size)
        tie = win_max - win_min
        dev = np.nanmax(tie)

        ncount = phase.shape[0] - mj

        devs[idx] = dev
        deverrs[idx] = dev / np.sqrt(ncount)
        ns[idx] = ncount
        idx += 1

    return remove_small_ns(taus_used, devs, deverrs, ns)


def three_cornered_hat_phase(phasedata_ab, phasedata_bc, phasedata_ca, rate, taus, function):
    """ Three Cornered Hat Method

    Three clocks with unknown variances sa^2, sb^2, sc^3
    Three pairwise measurements give variances:
    sab^2, sbc^2, sca^2
    Assuming covariances are zero, we get:
    sa^2 = 0.5*( sab^2 + sca^2 - sbc^2 )
    (and cyclic permutations for sb and sc) """

    (tau_ab, dev_ab, err_ab, ns_ab) = function(phasedata_ab, rate, taus)
    (tau_bc, dev_bc, err_bc, ns_bc) = function(phasedata_bc, rate, taus)
    (tau_ca, dev_ca, err_ca, ns_ca) = function(phasedata_ca, rate, taus)

    var_ab = dev_ab * dev_ab
    var_bc = dev_bc * dev_bc
    var_ca = dev_ca * dev_ca
    assert len(var_ab) == len(var_bc) == len(var_ca)
    var_a = 0.5 * (var_ab + var_ca - var_bc)

    dev_a = np.sqrt(var_a)
    dev_a[var_a < 0] = 0

    return tau_ab, dev_a


if __name__ == "__main__":
    print "Nothing to see here."

