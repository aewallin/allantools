"""
This is the old pure python versio of allantools.
Please see https://github.com/aewallin/allantools for the latest version.


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

def tdev_phase(phase, rate, taus):
    """ Time deviation of phase data

    Parameters
    ----------
    phase: list
        Phase data
    rate: float
        The sampling rate, in Hz
    taus: list
        Array of tau values for which to compute allan variance

    Returns
    -------
    (taus2, td, tde, ns): tuple
        taus2: list
            Tau values for which td computed
        td: list
            Computed time deviations for each tau value
        tde: list
            Time deviation errors
        ns: list
            Values of N used in mdev_phase()
    """
    (taus2, md, mde, ns) = mdev_phase(phase, rate, taus)
    td = [t * m / np.sqrt(3.0) for (t, m) in zip(taus2, md)]
    tde = [x / np.sqrt(n) for (x, n) in zip(td, ns)]
    return taus2, td, tde, ns

def tdev(data, rate, taus):
    """ Time deviation of fractional frequency data
    http://en.wikipedia.org/wiki/Time_deviation """
    phase = frequency2phase(data, rate)
    return tdev_phase(phase, rate, taus)

def mdev_phase(data, rate, taus):
    """# Modified Allan deviation of phase data

     1             N-3m+1     j+m-1
     Mod s2y(t) = ------------------  sum      { sum   [x(i+2m) - 2x(i+m) + x(i) ]  }**2
                  2m**2 t**2 (N-3m+1) j=1        i=j

     see http://www.leapsecond.com/tools/adev_lib.c """

    (ms, taus_used) = tau_m(data, rate, taus)
    # taus = []
    md = []
    mderr = []
    ns = []
    for m in ms:
        s = 0
        v = 0
        n = 0
        tau = m / float(rate)
        i = 0
        while (i + 2 * m) < len(data) and i < m:
            v = v + data[i + 2 * m] - 2 * data[i + m] + data[i]
            i += 1
        s += v * v
        n += 1

        i = 0
        while (i + 3 * m) < len(data):
            v = v + data[i + 3 * m] - 3 * data[i + 2 * m] + 3 * data[i + m] - data[i]
            s += v * v
            n += 1
            i += 1
        s /= float(2.0 * m * m * tau * tau * n)
        # assert( n == len(data)-3*m+1 ) # n is the normalization (N-3m+1) before the sums
        s = np.sqrt(s)
        md.append(s)
        mderr.append(s / np.sqrt(n))
        ns.append(n)
    return remove_small_ns(taus_used, md, mderr, ns)

def mdev(freqdata, rate, taus):
    """ # modified Allan deviation, fractional frequency data """
    phase = frequency2phase(freqdata, rate)
    return mdev_phase(phase, rate, taus)

def tau_m(data, rate, taus, v=False):
    """ pre-processing of the tau-list given by the user """
    if rate == 0:
        print("Warning! rate==0")
    rate = float(rate)
    # n = len(data) # not used
    m = []
    for tau in taus:
        if 0 < tau < (1 / float(rate)) * float(len(data)):  # tau should be in [0, len(data)/rate]
            mvalue = int(np.floor(float(tau * rate)))
            if mvalue != 0:
                m.append(mvalue)  # m is tau in units of datapoints
    m = list(set(m))  # this removes duplicates
    m.sort()  # sort from small tau to large tau
    if v:
        print("tau_m: ", m)
    if len(m) == 0:
        print("Warning: sanity-check on tau failed!")
        print("   len(data)=", len(data), " rate=", rate, "taus= ", taus)
    taus2 = [x / float(rate) for x in m]
    return m, taus2

def adev(data, rate, taus):
    """Allan deviation
    data is a time-series of fractional frequency
    rate is the samples/s in the time-series
    taus is a list of tau-values for which we compute ADEV """
    phase = frequency2phase(data, rate)
    return adev_phase(phase, rate, taus)

def adev_phase(data, rate, taus):
    (m, taus_used) = tau_m(data, rate, taus)
    ad = []
    ade = []
    adn = []
    for mj in m:  # loop through each tau value m(j)
        (dev, deverr, n) = calc_adev_phase(data, rate, mj, mj)
        ad.append(dev)
        ade.append(deverr)
        adn.append(n)
    return remove_small_ns(taus_used, ad, ade, adn)  # tau, adev, adeverror, naverages

def calc_adev_phase(data, rate, mj, stride):
    s = 0
    n = 0
    count = len(data)
    i = 0
    while (i + 2 * mj) < count:
        v = data[i + 2 * mj] - 2 * data[i + mj] + data[i]
        s += v * v
        i = i + stride
        n += 1
    s /= float(2.0)
    dev = 0
    deverr = 0
    if not n == 0:
        dev = np.sqrt(s / float(n)) / float(mj * (1 / float(rate)))
        deverr = dev / np.sqrt(n)
    return dev, deverr, n

def remove_small_ns(taus, devs, deverrs, ns):
    """ if n is small (==1), reject the result """
    o_taus = []
    o_dev = []
    o_err = []
    o_n = []
    for (t, d, e, n) in zip(taus, devs, deverrs, ns):
        if n > 1:
            o_taus.append(t)
            o_dev.append(d)
            o_err.append(e)
            o_n.append(n)
    return o_taus, o_dev, o_err, o_n

def oadev_phase(data, rate, taus):
    """ overlapping Allan deviation of phase data """
    (m, taus_used) = tau_m(data, rate, taus)
    ad = []
    ade = []
    adn = []
    for mj in m:
        (dev, deverr, n) = calc_adev_phase(data, rate, mj, 1)  # stride=1 for overlapping ADEV
        ad.append(dev)
        ade.append(deverr)
        adn.append(n)
    return remove_small_ns(taus_used, ad, ade, adn)  # tau, adev, adeverror, naverages

def oadev(freqdata, rate, taus):
    """ overlapping Allan deviation """
    phase = frequency2phase(freqdata, rate)
    return oadev_phase(phase, rate, taus)

def frequency2phase(freqdata, rate):
    """ integrate fractional frequency data and output phase data """
    phase = [0] * (len(freqdata) + 1)  # initialize to zero
    dt = 1 / float(rate)
    for i in range(len(freqdata) + 1):
        if i == 0:
            phase[i] = 0
        else:
            phase[i] = phase[i - 1] + freqdata[i - 1] * dt
    return phase

def ohdev(freqdata, rate, taus):
    """ Overlapping Hadamard deviation """
    phase = frequency2phase(freqdata, rate)
    return ohdev_phase(phase, rate, taus)

def ohdev_phase(data, rate, taus):
    """ Overlapping Hadamard deviation of phase data """
    rate = float(rate)
    (m, taus_used) = tau_m(data, rate, taus)
    hdevs = []
    hdeverrs = []
    ns = []
    for mj in m:
        (h, e, n) = calc_hdev_phase(data, rate, mj, 1)  # stride = 1
        hdevs.append(h)
        hdeverrs.append(e)
        ns.append(n)
    return remove_small_ns(taus_used, hdevs, hdeverrs, ns)

def hdev(freqdata, rate, taus):
    """ Hadamard deviation """
    phase = frequency2phase(freqdata, rate)
    return hdev_phase(phase, rate, taus)

def hdev_phase(data, rate, taus):
    """ Hadamard deviation of phase data """
    rate = float(rate)
    (m, taus_used) = tau_m(data, rate, taus)
    hdevs = []
    hdeverrs = []
    ns = []
    for mj in m:
        h, e, n = calc_hdev_phase(data, rate, mj, mj)  # stride = mj
        hdevs.append(h)
        hdeverrs.append(e)
        ns.append(n)

    return remove_small_ns(taus_used, hdevs, hdeverrs, ns)

def calc_hdev_phase(data, rate, mj, stride):
    """ http://www.leapsecond.com/tools/adev_lib.c """
    s = 0
    n = 0
    i = 0
    tau0 = 1 / float(rate)

    while (i + 3 * mj) < len(data):
        v = data[i + 3 * mj] - 3 * data[i + 2 * mj] + 3 * data[i + mj] - data[i]
        s += v * v
        n += 1
        i = i + stride

    s /= 6.0

    if n == 0:
        n = 1
    h = np.sqrt(s / float(n)) / float(tau0 * mj)
    e = h / np.sqrt(n)
    return h, e, n

def totdev(freqdata, rate, taus):
    phasedata = frequency2phase(freqdata, rate)
    return totdev_phase(phasedata, rate, taus)

def totdev_phase(data, rate, taus):
    """
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
    (m, taus_used) = tau_m(data, rate, taus)
    n = len(data)

    # totdev requires a new dataset
    x1 = []
    for j in range(n - 2):
        x1.append(float(2.0) * data[0] - data[j + 1])
    x1.reverse()

    x2 = []
    for j in range(n - 2):
        x2.append(float(2.0) * data[len(data) - 1] - data[len(data) - 1 - (j + 1)])

    x = []
    for d in x1:  # reflected data at start
        x.append(d)
    for d in data:  # original data
        x.append(d)
    for d in x2:  # reflected data at end
        x.append(d)

    # original dataset is now in the middle of the new dataset
    try:
        assert data[0] == x[len(x1)]
    except AssertionError:
        raise RuntimeError("dataset reflection error: Data[0] != x[len(x1)]")
    devs = []
    deverrs = []
    ns = []
    for mj in m:
        dev = 0
        ncount = 0
        for i in range(1, len(data) - 1):
            i += len(x1)  # i=0 corresponds to start of original data
            dev += pow(x[i - mj] - 2 * x[i] + x[i + mj], 2)
            ncount += 1
        dev /= float(2 * pow(mj / rate, 2) * (n - 2))
        dev = np.sqrt(dev)
        devs.append(dev)
        deverrs.append(dev / np.sqrt(ncount))
        ns.append(ncount)

    return remove_small_ns(taus_used, devs, deverrs, ns)

def tierms(freqdata, rate, taus):
    phasedata = frequency2phase(freqdata, rate)
    return tierms_phase(phasedata, rate, taus)

def tierms_phase(phase, rate, taus):
    """ TIE rms """
    rate = float(rate)
    (m, taus_used) = tau_m(phase, rate, taus)
    count = len(phase)
    devs = []
    deverrs = []
    ns = []
    for mj in m:
        dev = 0
        ncount = 0
        tie = []
        for i in range(count - mj):
            phases = [phase[i], phase[i + mj]]  # pair of phases at distance mj from eachother
            tie.append(max(phases) - min(phases))  # phase error
            ncount += 1
        # RMS of tie vector
        tie = [pow(x, 2) for x in tie]  # square
        tie = np.mean(tie)  # mean
        tie = np.sqrt(tie)  # root
        devs.append(tie)
        deverrs.append(dev / np.sqrt(ncount))
        ns.append(ncount)

    return remove_small_ns(taus_used, devs, deverrs, ns)

def mtie(freqdata, rate, taus):
    """ maximum time interval error, for fractional frequency data """
    phasedata = frequency2phase(freqdata, rate)
    return mtie_phase_purepy(phasedata, rate, taus)

def mtie_phase_numpy(phase, rate, taus):
    """ maximum time interval error, for phase data
    this seems to correspond to Stable32 setting "Fast(u)"
    Stable32 also has "Decade" and "Octave" modes where the dataset is extended somehow?
    """
    rate = float(rate)
    (m, taus_used) = tau_m(phase, rate, taus)
    devs = np.zeros_like(taus_used)
    deverrs = np.zeros_like(taus_used)
    ns = np.zeros_like(taus_used)
    data = np.array(phase)
    idx = 0
    for mj in m:
        dev = 0
        ncount = 0
        win_max = 0
        win_min = 0
        # we start by finding the max/min in the first window
        # when we slide the window forward by one step, one of these things happen:
        # - neither the new sample that enters the window, nor the old sample that leaves the window
        #    changes the max/min of the window
        # - the new sample that enters the window is the new max/min of the window
        # - the old sample that leaves the window was the max/min of the window
        #    this is the most computationally expensive case, as we now have to search
        #    the entire window for a new max/min
        for i in range(0, len(data) - mj):  # slide the start of the window over the dataset
            if i == 0:  # initialize window max/min
                win_max = np.max(data[0:mj + 1])  # max and min in the first window
                win_min = np.min(data[0:mj + 1])
            else:
                newsample = data[i + mj]  # the new sample that enters the window
                if newsample > win_max:
                    win_max = newsample
                elif newsample < win_min:
                    win_min = newsample
                oldsample = data[i - 1]  # the old sample we throw away
                if win_max == oldsample:
                    win_max = np.max(data[i:i + mj + 1])  # must search for a new maximum
                if win_min == oldsample:
                    win_min = np.min(data[i:i + mj + 1])  # must search for new minimum sample
            tie = win_max - win_min  # largest error in this window
            if tie > dev:
                dev = tie # store max of all windows as the MTIE result
            ncount += 1

        devs[idx] = dev
        deverrs[idx] = dev / np.sqrt(ncount)
        ns[idx] = ncount
        idx += 1

    return remove_small_ns(taus_used, devs, deverrs, ns)

def mtie_phase_purepy(phase, rate, taus):
    """ maximum time interval error, for phase data
    this seems to correspond to Stable32 setting "Fast(u)"
    Stable32 also has "Decade" and "Octave" modes where the dataset is extended somehow?
    """
    rate = float(rate)
    (m, taus_used) = tau_m(phase, rate, taus)
    devs = []
    deverrs = []
    ns = []

    for mj in m:
        dev = 0
        ncount = 0
        win_max = 0
        win_min = 0
        # we start by finding the max/min in the first window
        # when we slide the window forward by one step, one of these things happen:
        # - neither the new sample that enters the window, nor the old sample that leaves the window
        #    changes the max/min of the window
        # - the new sample that enters the window is the new max/min of the window
        # - the old sample that leaves the window was the max/min of the window
        #    this is the most computationally expensive case, as we now have to search
        #    the entire window for a new max/min
        for i in range(0, len(phase) - mj):  # slide the start of the window over the dataset
            if i == 0:  # initialize window max/min
                win_max = max(phase[0:mj + 1])  # max and min in the first window
                win_min = min(phase[0:mj + 1])
            else:
                newsample = phase[i + mj]  # the new sample that enters the window
                if newsample > win_max:
                    win_max = newsample
                elif newsample < win_min:
                    win_min = newsample
                oldsample = phase[i - 1]  # the old sample we throw away
                if win_max == oldsample:
                    win_max = max(phase[i:i + mj + 1])  # must search for a new maximum
                if win_min == oldsample:
                    win_min = min(phase[i:i + mj + 1])  # must search for new minimum sample
            tie = win_max - win_min  # largest error in this window
            if tie > dev:
                dev = tie # store max of all windows as the MTIE result
            ncount += 1

        devs.append(dev)
        deverrs.append(dev / np.sqrt(ncount) )
        ns.append(ncount)
        

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
    var_ab = [x * x for x in dev_ab]
    var_bc = [x * x for x in dev_bc]
    var_ca = [x * x for x in dev_ca]
    assert len(var_ab) == len(var_bc) == len(var_ca)
    var_a = [0.5 * (ab + ca - bc) for (ab, bc, ca) in zip(var_ab, var_bc, var_ca)]
    dev_a = []
    for va in var_a:
        try:
            dev_a.append(np.sqrt(va))
        except:
            dev_a.append(0)
    return tau_ab, dev_a

if __name__ == "__main__":
    print("Nothing to see here.")

