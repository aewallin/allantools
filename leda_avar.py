#!/usr/bin/env python
"""
leda_avar.py
------------

Allan variance code written from the estimate in TMS, to be used on
LEDA data. Written before I realized that this is exactly what the
allantools.adev() function does.

TODO: Integrate this code with main allantools suite. the avar_series
code here seems to run faster by about 60x against allantools.adev(),
and the output is equivalent. Also: need to get sample rates / timestamps
to be handled the same as in allantools.

"""

import numpy as np
import pylab as plt
import allantools
import time

def gen_log_idxs(data, spacing=1.25, minidx=16):
    """ Generate evenly spaced indexes in log space

    data (np.array): numpy array of data
    spacing (float): spacing multiplier between indexes
    mintau (int):    lowest index number

    returns np.array of evenly spaced indexes
    """

    dl = data.shape[0]
    idxs = []
    while dl > minidx:
        idxs.append(dl)
        dl = int(dl / spacing)
    idxs.sort()

    return np.array(idxs)

def avar(data, n=16, m=4):
    """ Compute allan variance of data

    data (np.array): Time series data
    n (int):         Number of samples in tau (int)
    m (int):         Number of samples in estimate

    Notes
    -----
    This is an estimate method from Thompson Moran Swenson (TMS 9.105)
    Interferometry and Synthesis in Radio Astronomy

    AVAR = 1 / 2(m-1) * Sum{k, m-1}(yk+1 - yk)^2

    where yk are averages over a time period. This method computes this
    in a reasonably efficient way by:
    1) Truncate data to be of length n * m
    2) Reshape to (m, n) and average over axis n to get yk values
    3) Difference neighbouring yk values
    4) square data
    5) return allan variance estimate

    TODO: I suspect that for 2^N powers of tau  you could do a recursive
    reshaping to speed up compute times for very long datasets.
    """
    try:
        assert n * m <= data.shape[0]
    except:
        raise RuntimeError("Data length (%i) < n * m (%i * %i)"
                            % (data.shape[0], n, m))

    data  = data[0:n*m]
    data  = data.reshape((m,n))
    avgs  = np.average(data, axis=1)
    #print data.shape, avgs.shape
    diffs = np.diff(avgs)
    avar  = 0.5 * np.average(diffs*diffs)

    return avar

def avar_series(data, spacing=1.25, mintau=16, minest=4, verbose=False):
    """ Compute the allan variance for a series of tau values

    Compute series of allan variances of increasing tau (time), as
    would be needed for making allan variance plots. This calls
    the avar() function several times for tau values that are evenly
    spaced in log space (e.g. tau = 2, 4, 8, 16, ...)

    data (np.array): time series data
    spacing (float): spacing S between time tau(n) =  S * tau(n-1)
    mintau (int): minimum number of time samples to use for yk values
                  as used in the (yk+1 - yk)^2
                  This is the min value of n to be used in avar(data, n , m)
    minest (int): minimum partitioning of time series used in the estimate.
                  This is the m argument in the avar(data, n, m) method.

    returns an array of tau values (time lengths), and prints data to screen
    in tabular format.

    """
    if verbose:
        print "\n----------------------"
        print "  N  |  M   | AVAR "
        print "----------------------"

    avars, taus = [], []
    nidxs = gen_log_idxs(data, spacing=spacing, minidx=mintau)

    for n in nidxs:
        m = int(len(data) / n)
        if m < minest:
            break
        av = avar(data, n=n, m=m)
        if verbose:
            print "%04d | %04d | %2.2e" % (n, m, av)
        avars.append(av)
        taus.append(n)

    if verbose:
        print "----------------------\n"
    return np.array(taus), np.array(avars)

if __name__ == '__main__':

    #########
    # SETUP TEST DATA
    #########

    rate      = 1000.0          # 1000 Hz sample rate
    obs_s     = 60 * 60         # 1 hour
    n_samples = rate * obs_s
    t = np.arange(0, n_samples)

    h  = np.random.random(n_samples) + np.sin(t / n_samples)

    #########
    # COMPUTE & PLOT ALLAN VARIANCES
    #########

    pcol, pcol2  = '#0b566c', '#cc0000'

    t1 = time.time()
    taus, avars = avar_series(h, spacing=1.25, mintau=1, minest=4)
    adevs = np.sqrt(avars)
    t2 = time.time()

    t3 = time.time()
    (tall, ad, ade, adn) = allantools.adev(h, 1.0, taus)
    t4 = time.time()

    assert np.allclose(adevs, ad)

    tstamps = (t[1] - t[0]) * taus / rate
    plt.loglog(np.array(tall)/rate, ad, c=pcol2, label='adev()')
    plt.loglog(tstamps, adevs, c=pcol, label='avar_series()')

    print "Time avar_series():     %2.2e s" % (t2 - t1)
    print "Time allantools.adev(): %2.2e s" % (t4 - t3)
    print "Speedup:                %2.2fx"  % ((t4 - t3) / (t2 - t1))

    plt.legend(loc=1, prop={'size':8})
    plt.xlabel("Time $\\tau$ (s)")
    plt.ylabel("Allan deviation")
    plt.tight_layout()
    plt.show()
