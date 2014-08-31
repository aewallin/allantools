"""
Allan deviation tools
Anders Wallin (anders.e.e.wallin "at" gmail.com)
v1.0 2014 January
v1.01 2014 August (PEP8 compliance improvements by telegraphic)

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


Implemented functions:
adev() and adev_phase()     Allan deviation
oadev() and oadev_phase()   Overlapping Allan deviation
mdev() and mdev_phase()     Modified Allan deviation
tdev() and tdev_phase()     Time deviation ( modified variance scaled by (t^2/3) )
hdev() and hdev_phase()     Hadamard deviation
ohdev() and ohdev_phase()   Overlapping Hadamard deviation
totdev() and totdev_phase() Total deviation
mtie() and mtie_phase()     Maximum Time Interval Error
tierms() and tierms_phase() Time Interval Error RMS

to do:
Modified Total
The modified total variance, MTOT, is total version of the modified Allan variance.
It is defined for phase data as:
                         1           N-3m+1  1  N+3m-1
Mod s^2 total(t) = ----------------- sum     -- sum      [0zi*(m)]^2
                    2m^2t0^2(N-3m+1) n=1     6m i=n-3m
where the 0zi*(m) terms are the phase averages from a triply-extended
sequence created by uninverted even reflection at each end,
and the prefix 0 denotes that the linear trend has been removed.

Time Total (modified total variance scaled by (t^2/3) )
Hadamard Total


References
http://www.wriley.com/paper4ht.htm
http://en.wikipedia.org/wiki/Allan_variance
for code see e.g.:
http://www.mathworks.com/matlabcentral/fileexchange/26659-allan-v3-0
http://www.mathworks.com/matlabcentral/fileexchange/26637-allanmodified
http://www.leapsecond.com/tools/adev_lib.c

Inputs (phase data):
    phase = list of phase measurements in seconds, e.g. from a time-interval-counter
    rate  = sample rate of data, i.e. interval between phase measurements is 1/rate
    taus  = list of tau-values for ADEV computation
Inputs (frequency data):
    data = list of fractional frequency measurements (nondimensional), e.g. from a frequency-counter
    rate  = sample rate of data, i.e. gate time of a zero-dead-time counter is 1/rate
    taus  = list of tau-values for ADEV computation
Output (tau_out, adev, adeverr, n)
    tau_out = list of tau-values for which deviations were computed
    adev    = list of ADEV (or another statistic) deviations
    adeverr = list of estimated errors of allan deviations
    n       = list of number of pairs in allan computation. standard error is adeverr = adev/sqrt(n)
"""

import numpy as np

def rolling_window(a, window):
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
    of size w.

    """
    if window < 1:
        raise ValueError, "`window` must be at least 1."
    if window > a.shape[-1]:
        raise ValueError, "`window` is too long."
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)


def tdev_phase(phase, rate, taus):
    """ Time deviation of phase data

    Parameters
    ----------
    phase: np.array
        Phase data
    rate: float
        The sampling rate, in Hz
    taus: np.array
        Array of tau values for which to compute allan variance

    Returns
    -------
    (taus2, td, tde, ns): tuple
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
    http://en.wikipedia.org/wiki/Time_deviation """
    phase = frequency2phase(data, rate)
    return tdev_phase(phase, rate, taus)


def mdev_phase(data, rate, taus):
    """  Modified Allan deviation of phase data

                   N-3m+1     j+m-1
     Mod s2y(t) = ------------------  sum      { sum   [x(i+2m) - 2x(i+m) + x(i) ]  }**2
                  2m**2 t**2 (N-3m+1) j=1        i=j

     see http://www.leapsecond.com/tools/adev_lib.c """
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
    """ modified Allan deviation, fractional frequency data """
    phase = frequency2phase(freqdata, rate)
    return mdev_phase(phase, rate, taus)


def tau_m(data, rate, taus, v=False):
    """ pre-processing of the tau-list given by the user """
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
    data is a time-series of evenly spaced fractional frequency measurements
    rate is the samples/s in the time-series
    taus is a list of tau-values for which we compute ADEV """
    phase = frequency2phase(data, rate)
    return adev_phase(phase, rate, taus)


def adev_phase(data, rate, taus):
    """ Allan deviation for phase data """
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
    """ overlapping Allan deviation of phase data """
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
    """ overlapping Allan deviation for fractional frequency data """
    phase = frequency2phase(freqdata, rate)
    return oadev_phase(phase, rate, taus)


def frequency2phase(freqdata, rate):
    """ integrate fractional frequency data and output phase data """
    dt = 1.0 / float(rate)
    phasedata = np.cumsum(freqdata) * dt
    phasedata = np.insert(phasedata, 0, 0)
    return phasedata


def ohdev(freqdata, rate, taus):
    """ Overlapping Hadamard deviation, fractional frequency data """
    phase = frequency2phase(freqdata, rate)
    return ohdev_phase(phase, rate, taus)


def ohdev_phase(data, rate, taus):
    """ Overlapping Hadamard deviation of phase data """
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
    """ Hadamard deviation, fractional frequency data """
    phase = frequency2phase(freqdata, rate)
    return hdev_phase(phase, rate, taus)


def hdev_phase(data, rate, taus):
    """ Hadamard deviation, phase data """
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
    """ http://www.leapsecond.com/tools/adev_lib.c """

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
    """
    Total deviation, phase data.
    
    See:
    David A. Howe,
    The total deviation approach to long-term characterization
    of frequency stability,
    IEEE tr. UFFC vol 47 no 5 (2000)

                     1         N-1
    totvar(t) = ------------  sum   [ x*(i-m) - 2x*(i)+x*(i+m) ]**2
                2 t**2 (N-2)  i=2
    where x* is a new dataset with 'reflected' data at start/end """

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
    """ Time Interval Error RMS, for phase data """
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
    """ Maximum Time Interval Error, phase data
    
    this seems to correspond to Stable32 setting "Fast(u)"
    Stable32 also has "Decade" and "Octave" modes where the dataset is extended somehow?

    rate = float(rate) """
    (phase, m, taus_used) = tau_m(phase, rate, taus)
    devs = np.zeros_like(taus_used)
    deverrs = np.zeros_like(taus_used)
    ns = np.zeros_like(taus_used)

    idx = 0
    for mj in m:
        rw = rolling_window(phase, mj + 1)
        win_max = np.max(rw, axis=1)
        win_min = np.min(rw, axis=1)
        tie = win_max - win_min
        dev = np.max(tie)
        ncount = phase.shape[0] - mj
        devs[idx] = dev
        deverrs[idx] = dev / np.sqrt(ncount)
        ns[idx] = ncount
        idx += 1

    return remove_small_ns(taus_used, devs, deverrs, ns)


def three_cornered_hat_phase(phasedata_ab, phasedata_bc, phasedata_ca, rate, taus, function):
    """ Three Cornered Hat Method Three clocks with unknown variances sa^2, sb^2, sc^3
    Three pairwise measurements give variances:
    sab^2, sbc^2, sca^2
    Assuming covariances are zero (clocks are independent), we get:
    sa^2 = 0.5*( sab^2 + sca^2 - sbc^2 )
    (and cyclic permutations for sb and sc) """
    
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


if __name__ == "__main__":
    print "Nothing to see here."

