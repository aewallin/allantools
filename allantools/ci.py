"""

this file is part of allantools, https://github.com/aewallin/allantools

- functions for confidence intervals
- functions for noise identification
- functions for computing equivalent degrees of freedom


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

import numpy as np
import scipy.stats  # used in confidence_intervals()
import scipy.signal  # decimation in lag-1 acf
import allantools as at

########################################################################
# Confidence Intervals
ONE_SIGMA_CI = scipy.special.erf(1/np.sqrt(2))
#    = 0.68268949213708585


def confidence_interval(dev, edf, ci=ONE_SIGMA_CI):
    """ returns confidence interval (dev_min, dev_max)
        for a given deviation dev, equivalent degrees of freedom edf,
        and degree of confidence ci.

    Parameters
    ----------
    dev: float
        Mean value (e.g. adev) around which we produce the confidence interval
    edf: float
        Equivalent degrees of freedon
    ci: float, defaults to scipy.special.erf(1/math.sqrt(2))
        for 1-sigma standard error set
        ci = scipy.special.erf(1/math.sqrt(2))
            = 0.68268949213708585

    Returns
    -------
    (dev_min, dev_max): (float, float)
        Confidence interval
    """
    ci_l = min(np.abs(ci), np.abs((ci-1))) / 2
    ci_h = 1 - ci_l

    # function from scipy, works OK, but scipy is large and slow to build
    chi2_l = scipy.stats.chi2.ppf(ci_l, edf)
    chi2_h = scipy.stats.chi2.ppf(ci_h, edf)

    variance = dev*dev
    var_l = float(edf) * variance / chi2_h  # NIST SP1065 eqn (45)
    var_h = float(edf) * variance / chi2_l
    return (np.sqrt(var_l), np.sqrt(var_h))


def confidence_interval_noiseID(x, dev, af, dev_type="adev",
                                data_type="phase", ci=ONE_SIGMA_CI):
    """ returns confidence interval (dev_min, dev_max)
        for a given deviation dev = Xdev( x, tau = af*(1/rate) )

        steps:
        1) identify noise type
        2) compute EDF
        3) compute confidence interval

    Parameters
    ----------
    x: numpy.array
        time-series
    dev: float
        Mean value (e.g. adev) around which we produce the confidence interval
    af: int
        averaging factor
    dev_type: string
        adev, oadev, mdev, tdev, hdev, ohdev
    data_type:
        "phase" or "freq"
    ci: float, defaults to scipy.special.erf(1/math.sqrt(2))
        for 1-sigma standard error set
        ci = scipy.special.erf(1/math.sqrt(2))
            = 0.68268949213708585

    Returns
    -------
    (dev_min, dev_max): (float, float)
        Confidence interval
    """
    # 1) noise ID
    dmax = 2
    if (dev_type == "hdev") or (dev_type == "ohdev"):
        dmax = 3
    alpha_int = autocorr_noise_id(
        x, int(af), data_type=data_type, dmin=0, dmax=dmax)[0]

    # 2) EDF (tes was 'is')
    if dev_type == "adev":
        edf = edf_greenhall(alpha=alpha_int, d=2, m=af, N=len(x),
                            overlapping=False, modified=False)
    elif dev_type == "oadev":
        edf = edf_greenhall(alpha=alpha_int, d=2, m=af, N=len(x),
                            overlapping=True, modified=False)
    elif (dev_type == "mdev") or (dev_type == "tdev"):
        edf = edf_greenhall(alpha=alpha_int, d=2, m=af, N=len(x),
                            overlapping=True, modified=True)
    elif dev_type == "hdev":
        edf = edf_greenhall(alpha=alpha_int, d=3, m=af, N=len(x),
                            overlapping=False, modified=False)
    elif dev_type == "ohdev":
        edf = edf_greenhall(alpha=alpha_int, d=3, m=af, N=len(x),
                            overlapping=True, modified=False)
    else:
        raise NotImplementedError

    # 3) confidence interval
    (low, high) = confidence_interval(dev, edf, ci)
    return (low, high)


########################################################################
# Noise Identification using R(n)


def rn(x, af, rate):
    """ R(n) ratio for noise identification

        ratio of MVAR to AVAR
    """
    (taus, devs, errs, ns) = at.adev(
        x, taus=[af*rate], data_type='phase', rate=rate)
    oadev_x = devs[0]
    (mtaus, mdevs, errs, ns) = at.mdev(
        x, taus=[af*rate], data_type='phase', rate=rate)
    mdev_x = mdevs[0]
    return pow(mdev_x/oadev_x, 2)


def rn_theory(af, b):
    """ R(n) ratio expected from theory for given noise type

        alpha = b + 2
    """
    # From IEEE1139-2008
    #   alpha   beta    ADEV_mu MDEV_mu Rn_mu
    #   -2      -4       1       1       0      Random Walk FM
    #   -1      -3       0       0       0      Flicker FM
    #    0      -2      -1      -1       0      White FM
    #    1      -1      -2      -2       0      Flicker PM
    #    2      0       -2      -3      -1      White PM

    # (a=-3 flicker walk FM)
    # (a=-4 random run FM)
    if b == 0:
        return pow(af, -1)
    elif b == -1:
        # f_h = 0.5/tau0  (assumed!)
        # af = tau/tau0
        # so f_h*tau = 0.5/tau0 * af*tau0 = 0.5*af
        avar = (1.038+3*np.log(2*np.pi*0.5*af)) / (4.0*pow(np.pi, 2))
        mvar = 3*np.log(256.0/27.0)/(8.0*pow(np.pi, 2))
        return mvar/avar
    else:
        return pow(af, 0)


def rn_boundary(af, b_hi):
    """
    R(n) ratio boundary for selecting between [b_hi-1, b_hi]
    alpha = b + 2
    """
    return np.sqrt(rn_theory(af, b_hi)*rn_theory(af, b_hi-1))  # geometric mean

########################################################################
# Noise Identification using B1


def b1(x, af, rate):
    """ B1 ratio for noise identification, [Barnes1974]_ and [Howe2000a]_
        (and bias correction?)

        ratio of Standard Variace to AVAR

    """
    (taus, devs, errs, ns) = at.adev(
        x, taus=[af*rate], data_type="phase", rate=rate)
    oadev_x = devs[0]
    avar = pow(oadev_x, 2.0)

    # variance of y, at given af
    y = np.diff(x)
    y_cut = np.array(y[:len(y)-(len(y) % af)])  # cut to length
    assert len(y_cut) % af == 0
    y_shaped = y_cut.reshape((int(len(y_cut)/af), af))
    y_averaged = np.average(y_shaped, axis=1)  # average
    var = np.var(y_averaged, ddof=1)

    return var/avar


def b1_theory(N, mu):
    """ Expected B1 ratio for given time-series length N and exponent mu

        FIXME: add reference (paper & link)

        The exponents are defined as
        S_y(f) = h_a f^alpha    (power spectrum of y)
        S_x(f) = g_b f^b        (power spectrum of x)
        bias = const * tau^mu

        and (b, alpha, mu) relate to eachother by:
        b    alpha   mu
        0    +2      -2
       -1    +1      -2   resolve between -2 cases with R(n)
       -2     0      -1
       -3    -1       0
       -4    -2      +1
       -5    -3      +2
       -6    -4      +3 for HDEV, by applying B1 to frequency data,
                        and add +2 to resulting mu
    """

    # see Table 3 of Howe 2000
    if mu == 2:
        return float(N)*(float(N)+1.0)/6.0
    elif mu == 1:
        return float(N)/2.0
    elif mu == 0:
        return N*np.log(N)/(2.0*(N-1.0)*np.log(2))
    elif mu == -1:
        return 1
    elif mu == -2:
        return (pow(N, 2)-1.0)/(1.5*N*(N-1.0))
    else:
        up = N*(1.0-pow(N, mu))
        down = 2*(N-1.0)*(1-pow(2.0, mu))
        return up/down

    assert False  # we should never get here


def b1_boundary(b_hi, N):
    """
    B1 ratio boundary for selecting between [b_hi-1, b_hi]
    alpha = b + 2
    """
    b_lo = b_hi-1
    b1_lo = b1_theory(N, b_to_mu(b_lo))
    b1_hi = b1_theory(N, b_to_mu(b_hi))
    if b1_lo >= -4:
        return np.sqrt(b1_lo*b1_hi)  # geometric mean
    else:
        return 0.5*(b1_lo+b1_hi)  # arithemtic mean


def b_to_mu(b):
    """
    return mu, parameter needed for B1 ratio function b1()
    alpha = b + 2
    """
    a = b + 2
    if a == +2:
        return -2
    elif a == +1:
        return -2
    elif a == 0:
        return -1
    elif a == -1:
        return 0
    elif a == -2:
        return 1
    elif a == -3:
        return 2
    elif a == -4:
        return 3
    assert False

########################################################################
# Noise Identification using ACF


def lag1_acf(x, detrend_deg=1):
    """ Lag-1 autocorrelation function
        as defined in [Riley2004]_ eqn (2)
        used by autocorr_noise_id()

        Parameters
        ----------
        x: numpy.array
            time-series
        Returns
        -------
        ACF: float
            Lag-1 autocorrelation for input time-series x

        Notes
        -----
        * a faster algorithm based on FFT might be better!?
        * numpy.corrcoeff() gives similar but not identical results.
            #c = np.corrcoef( np.array(x[:-lag]), np.array(x[lag:]) )
            #r1 = c[0,1] # lag-1 autocorrelation of x
    """
    mu = np.mean(x)
    a = 0
    b = 0
    for n in range(len(x)-1):
        a = a + (x[n]-mu)*(x[n+1]-mu)
    # for n in range(len(x)):
    for xn in x:
        b = b+pow(xn-mu, 2)
    return a/b


def autocorr_noise_id(x, af, data_type="phase", dmin=0, dmax=2):
    """ Lag-1 autocorrelation based noise identification

    Parameters
    ----------
    x: numpy.array
        phase or fractional frequency time-series data
        minimum recommended length is len(x)>30 roughly.
    af: int
        averaging factor
    data_type: string {'phase', 'freq'}
        "phase" for phase data in seconds
        "freq" for fractional frequency data
    dmin: int
        minimum required number of differentiations in the algorithm
    dmax: int
        maximum number of differentiations
        defaults to 2 for ADEV
        set to 3 for HDEV

    Returns
    -------
    alpha_int: int
        noise-slope as integer
    alpha: float
        noise-slope as float
    d: int
        number of differentiations of the time-series performed

    Notes
    -----
        http://www.stable32.com/Auto.pdf
        http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.503.9864&rep=rep1&type=pdf

    Reference [Riley2004]_.

    """
    d = 0  # number of differentiations
    if data_type == "phase":
        if af > 1:
            # x = scipy.signal.decimate(x, af, n=1, ftype='fir')
            x = x[0:len(x):af]  # decimate by averaging factor
            # x = scipy.signal.decimate(x, af, ftype='fir',)
            # resampled_len = int(len(x)/af)
            # x = scipy.signal.resample(x, resampled_len)

        x = detrend(x, deg=2)  # remove quadratic trend (freq offset and drift)
    elif data_type == "freq":
        # average by averaging factor
        y_cut = np.array(x[:len(x)-(len(x) % af)])  # cut to length
        assert len(y_cut) % af == 0
        y_shaped = y_cut.reshape((int(len(y_cut)/af), af))
        x = np.average(y_shaped, axis=1)  # average
        x = detrend(x, deg=1)  # remove frequency drift

    # require minimum length for time-series
    if len(x) < 30:
        print(("autocorr_noise_id() Don't know how to do noise-ID for"
               " time-series length= %d") % len(x))
        raise NotImplementedError

    while True:
        r1 = lag1_acf(x)
        rho = r1/(1.0+r1)
        if d >= dmin and (rho < 0.25 or d >= dmax):
            p = -2*(rho+d)
            # print r1
            # assert r1 < 0
            # assert r1 > -1.0/2.0
            phase_add2 = 0
            if data_type == "phase":
                phase_add2 = 2
            alpha = p+phase_add2
            alpha_int = int(-1.0*np.round(2*rho) - 2.0*d) + phase_add2
            # print "d=",d,"alpha=",p+2
            return alpha_int, alpha, d, rho
        else:
            x = np.diff(x)
            d = d + 1
    assert False  # we should not get here ever.


def detrend(x, deg=1):
    """
    remove polynomial from data.
    used by autocorr_noise_id()

    Parameters
    ----------
    x: numpy.array
        time-series
    deg: int
        degree of polynomial to remove from x

    Returns
    -------
    x_detrended: numpy.array
        detrended time-series
    """
    t = range(len(x))
    p = np.polyfit(t, x, deg)
    residual = x - np.polyval(p, t)
    return residual

########################################################################
# Equivalent Degrees of Freedom


def edf_greenhall_simple(alpha, d, m, S, F, N):
    """ Eqn (13) from [Greenhall2004]_ """
    L = m/F+m*d  # length of filter applied to phase samples
    M = 1 + np.floor(S*(N-L) / m)
    J = min(M, (d+1)*S)
    inv_edf = ((1.0/(pow(greenhall_sz(0, F, alpha, d), 2)*M)) *
               greenhall_BasicSum(J, M, S, F, alpha, d))
    return 1.0/inv_edf


def edf_greenhall(alpha, d, m, N,
                  overlapping=False, modified=False, verbose=False):
    """ returns Equivalent degrees of freedom

        Parameters
        ----------
        alpha: int
            noise type, +2...-4
        d: int
            1 first-difference variance
            2 Allan variance
            3 Hadamard variance
            require alpha+2*d>1
        m: int
            averaging factor
            tau = m*tau0 = m*(1/rate)
        N: int
            number of phase observations (length of time-series)
        overlapping: bool
            True for oadev, ohdev
        modified: bool
            True for mdev, tdev

        Returns
        -------
        edf: float
            Equivalent degrees of freedom

        Notes
        -----
        Reference [Greenhall2004]_.

        Used for the following deviations
        (see [Riley_CI]_ page 8)
        adev()
        oadev()
        mdev()
        tdev()
        hdev()
        ohdev()
    """

    if modified:
        F = 1  # F filter factor, 1 modified variance, m unmodified variance
    else:
        F = int(m)
    if overlapping:  # S stride factor, 1 nonoverlapped estimator,
        S = int(m)  # m overlapped estimator (estimator stride = tau/S )
    else:
        S = 1
    assert(alpha+2*d > 1.0)
    L = m/F+m*d  # length of filter applied to phase samples
    M = 1 + np.floor(S*(N-L) / m)
    J = min(M, (d+1)*S)
    J_max = 100
    r = M/S
    if int(F) == 1 and modified:  # case 1, modified variances, all alpha
        if J <= J_max:
            inv_edf = (1.0/(pow(greenhall_sz(0, 1, alpha, d), 2)*M)) * \
                       greenhall_BasicSum(J, M, S, 1, alpha, d)
            if verbose:
                print("case 1.1 edf= %3f" % float(1.0/inv_edf))
            return 1.0/inv_edf
        elif r > d+1:
            (a0, a1) = greenhall_table1(alpha, d)
            inv_edf = (1.0/r)*(a0-a1/r)
            if verbose:
                print("case 1.2 edf= %3f" % float(1.0/inv_edf))
            return 1.0/inv_edf
        else:
            m_prime = J_max/r
            inv_edf = ((1.0/(pow(greenhall_sz(0, F, alpha, d), 2)*J_max)) *
                       greenhall_BasicSum(J_max, J_max, m_prime, 1, alpha, d))
            if verbose:
                print("case 1.3 edf= %3f" % float(1.0/inv_edf))
            return 1.0/inv_edf
    elif int(F) == int(m) and int(alpha) <= 0 and not modified:
        # case 2, unmodified variances, alpha <= 0
        if J <= J_max:
            if m*(d+1) <= J_max:
                m_prime = m
                variant = "a"
            else:
                m_prime = float('inf')
                variant = "b"

            inv_edf = ((1.0/(pow(greenhall_sz(0, m_prime, alpha, d), 2)*M)) *
                       greenhall_BasicSum(J, M, S, m_prime, alpha, d))
            if verbose:
                print("case 2.1%s edf= %3f" % (variant, float(1.0/inv_edf)))
            return 1.0/inv_edf
        elif r > d+1:
            (a0, a1) = greenhall_table2(alpha, d)
            inv_edf = (1.0/r)*(a0-a1/r)
            if verbose:
                print("case 2.2 edf= %3f" % float(1.0/inv_edf))
            return 1.0/inv_edf
        else:
            m_prime = J_max/r
            inv_edf = (
                (1.0/(pow(greenhall_sz(0, float('inf'), alpha, d), 2)*J_max)) *
                greenhall_BasicSum(
                    J_max, J_max, m_prime, float('inf'), alpha, d))
            if verbose:
                print("case 2.3 edf= %3f" % float(1.0/inv_edf))
            return 1.0/inv_edf
    elif int(F) == int(m) and int(alpha) == 1 and not modified:
        # case 3, unmodified variances, alpha=1
        if J <= J_max:
            # note: m<1e6 to avoid roundoff
            inv_edf = ((1.0/(pow(greenhall_sz(0, m, 1, d), 2)*M)) *
                       greenhall_BasicSum(J, M, S, m, 1, d))
            if verbose:
                print("case 3.1 edf= %3f" % float(1.0/inv_edf))
            return 1.0/inv_edf
        elif r > d+1:
            (a0, a1) = greenhall_table2(alpha, d)
            (b0, b1) = greenhall_table3(alpha, d)
            inv_edf = (1.0/(pow(b0+b1*np.log(m), 2)*r))*(a0-a1/r)
            if verbose:
                print("case 3.2 edf= %3f" % float(1.0/inv_edf))
            return 1.0/inv_edf
        else:
            m_prime = J_max/r
            (b0, b1) = greenhall_table3(alpha, d)
            inv_edf = (
                (1.0/(pow(b0+b1*np.log(m), 2)*J_max)) *
                greenhall_BasicSum(J_max, J_max, m_prime, m_prime, 1, d))
            if verbose:
                print("case 3.3 edf= %3f" % float(1.0/inv_edf))
            return 1.0/inv_edf
    elif int(F) == int(m) and int(alpha) == 2 and not modified:
        # case 4, unmodified variances, alpha=2
        K = np.ceil(r)
        if K <= d:
            # FIXME: add formula from the paper here!
            raise NotImplementedError
        else:
            a0 = (scipy.special.binom(4*d, 2*d) /
                  pow(scipy.special.binom(2*d, d), 2))
            a1 = d/2.0
            inv_edf = (1.0/M)*(a0-a1/r)
            if verbose:
                print("case 4.2 edf= %3f" % float(1.0/inv_edf))
            return 1.0/inv_edf

    print("greenhall_edf() no matching case!")
    raise NotImplementedError


def greenhall_BasicSum(J, M, S, F, alpha, d):
    """ Eqn (10) from [Greenhall2004]_ """
    first = pow(greenhall_sz(0, F, alpha, d), 2)
    second = ((1-float(J)/float(M)) *
              pow(greenhall_sz(float(J)/float(S), F, alpha, d), 2))
    third = 0
    for j in range(1, int(J)):
        third += (2 * (1.0-float(j)/float(M)) *
                  pow(greenhall_sz(float(j) / float(S), F, alpha, d), 2))
    return first+second+third


def greenhall_sz(t, F, alpha, d):
    """ Eqn (9) from [Greenhall2004]_ """
    if d == 1:
        a = 2*greenhall_sx(t, F, alpha)
        b = greenhall_sx(t-1.0, F, alpha)
        c = greenhall_sx(t+1.0, F, alpha)
        return a-b-c
    elif d == 2:
        a = 6*greenhall_sx(t, F, alpha)
        b = 4*greenhall_sx(t-1.0, F, alpha)
        c = 4*greenhall_sx(t+1.0, F, alpha)
        dd = greenhall_sx(t-2.0, F, alpha)
        e = greenhall_sx(t+2.0, F, alpha)
        return a-b-c+dd+e
    elif d == 3:
        a = 20.0*greenhall_sx(t, F, alpha)
        b = 15.0*greenhall_sx(t-1.0, F, alpha)
        c = 15.0*greenhall_sx(t+1.0, F, alpha)
        dd = 6.0*greenhall_sx(t-2.0, F, alpha)
        e = 6.0*greenhall_sx(t+2.0, F, alpha)
        f = greenhall_sx(t-3.0, F, alpha)
        g = greenhall_sx(t+3.0, F, alpha)
        return a-b-c+dd+e-f-g

    assert(0)  # ERROR


def greenhall_sx(t, F, alpha):
    """ Eqn (8) from [Greenhall2004]_
    """
    if F == float('inf'):
        return greenhall_sw(t, alpha+2)
    a = 2*greenhall_sw(t, alpha)
    b = greenhall_sw(t-1.0/float(F), alpha)
    c = greenhall_sw(t+1.0/float(F), alpha)

    return pow(F, 2)*(a-b-c)


def greenhall_sw(t, alpha):
    """ Eqn (7) from [Greenhall2004]_
    """
    alpha = int(alpha)
    if alpha == 2:
        return -np.abs(t)
    elif alpha == 1:
        if t == 0:
            return 0
        else:
            return pow(t, 2)*np.log(np.abs(t))
    elif alpha == 0:
        return np.abs(pow(t, 3))
    elif alpha == -1:
        if t == 0:
            return 0
        else:
            return pow(t, 4)*np.log(np.abs(t))
    elif alpha == -2:
        return np.abs(pow(t, 5))
    elif alpha == -3:
        if t == 0:
            return 0
        else:
            return pow(t, 6)*np.log(np.abs(t))
    elif alpha == -4:
        return np.abs(pow(t, 7))

    assert(0)  # ERROR


def greenhall_table3(alpha, d):
    """ Table 3 from [Greenhall2004]_ """
    assert(alpha == 1)
    idx = d-1
    table3 = [(6.0, 4.0), (15.23, 12.0), (47.8, 40.0)]
    return table3[idx]


def greenhall_table2(alpha, d):
    """ Table 2 from [Greenhall2004]_ """
    row_idx = int(-alpha+2)  # map 2-> row0 and -4-> row6
    assert(row_idx in [0, 1, 2, 3, 4, 5])
    col_idx = int(d-1)
    table2 = [
        # alpha = +2:
        [(3.0/2.0, 1.0/2.0), (35.0/18.0, 1.0), (231.0/100.0, 3.0/2.0)],
        # alpha = +1:
        [(78.6, 25.2), (790.0, 410.0), (9950.0, 6520.0)],
        # alpha = 0:
        [(2.0/3.0, 1.0/6.0), (2.0/3.0, 1.0/3.0), (7.0/9.0, 1.0/2.0)],
        # alpha = -1:
        [(-1, -1), (0.852, 0.375), (0.997, 0.617)],
        # alpha = -2
        [(-1, -1), (1.079, 0.368), (1.033, 0.607)],
        # alpha = -3
        [(-1, -1), (-1, -1), (1.053, 0.553)],
        # alpha = -4
        [(-1, -1), (-1, -1), (1.302, 0.535)]]
    # print("table2 = ", table2[row_idx][col_idx])
    return table2[row_idx][col_idx]


def greenhall_table1(alpha, d):
    """ Table 1 from [Greenhall2004]_ """
    row_idx = int(-alpha+2)  # map 2-> row0 and -4-> row6
    col_idx = int(d-1)
    table1 = [
        # alpha = +2:
        [(2.0/3.0, 1.0/3.0), (7.0/9.0, 1.0/2.0), (22.0/25.0, 2.0/3.0)],
        # alpha = +1:
        [(0.840, 0.345), (0.997, 0.616), (1.141, 0.843)],
        # alpha = 0:
        [(1.079, 0.368), (1.033, 0.607), (1.184, 0.848)],
        # alpha = -1:
        [(-1, -1), (1.048, 0.534), (1.180, 0.816)],
        # alpha = -2
        [(-1, -1), (1.302, 0.535), (1.175, 0.777)],
        # alpha = -3
        [(-1, -1), (-1, -1), (1.194, 0.703)],
        # alpha = -4
        [(-1, -1), (-1, -1), (1.489, 0.702)]]
    # print("table1 = ", table1[row_idx][col_idx])
    return table1[row_idx][col_idx]


def edf_totdev(N, m, alpha):
    """ Equivalent degrees of freedom for Total Deviation
        FIXME: what is the right behavior for alpha outside 0,-1,-2?

        NIST [SP1065]_ page 41, Table 7
    """
    alpha = int(alpha)
    if alpha in [0, -1, -2]:
        # alpha  0 WFM
        # alpha -1 FFM
        # alpha -2 RWFM
        NIST_SP1065_table7 = [(1.50, 0.0), (1.17, 0.22), (0.93, 0.36)]
        (b, c) = NIST_SP1065_table7[int(abs(alpha))]
        return b*(float(N)/float(m))-c
    # alpha outside 0, -1, -2:
    return edf_simple(N, m, alpha)


def edf_mtotdev(N, m, alpha):
    """ Equivalent degrees of freedom for Modified Total Deviation

        NIST [SP1065]_ page 41, Table 8
    """
    assert(alpha in [2, 1, 0, -1, -2])
    NIST_SP1065_table8 = [(1.90, 2.1),
                          (1.20, 1.40),
                          (1.10, 1.2),
                          (0.85, 0.50),
                          (0.75, 0.31)]
    # (b, c) = NIST_SP1065_table8[ abs(alpha-2) ]
    (b, c) = NIST_SP1065_table8[abs(alpha-2)]
    edf = b*(float(N)/float(m))-c
    print("mtotdev b,c= ", (b, c), " edf=", edf)
    return edf


def edf_simple(N, m, alpha):
    """Equivalent degrees of freedom.
    Simple approximate formulae.

    Parameters
    ----------
    N : int
        the number of phase samples
    m : int
        averaging factor, tau = m * tau0
    alpha: int
        exponent of f for the frequency PSD:
        'wp' returns white phase noise.             alpha=+2
        'wf' returns white frequency noise.         alpha= 0
        'fp' returns flicker phase noise.           alpha=+1
        'ff' returns flicker frequency noise.       alpha=-1
        'rf' returns random walk frequency noise.   alpha=-2
        If the input is not recognized, it defaults to idealized, uncorrelated
        noise with (N-1) degrees of freedom.

    Notes
    -----
       See [Stein1985]_

    Returns
    -------
    edf : float
        Equivalent degrees of freedom

    """

    N = float(N)
    m = float(m)
    if alpha in [2, 1, 0, -1, -2]:
        # NIST SP 1065, Table 5
        if alpha == +2:
            edf = (N + 1) * (N - 2*m) / (2 * (N - m))

        if alpha == 0:
            edf = (((3 * (N - 1) / (2 * m)) - (2 * (N - 2) / N)) *
                   ((4*pow(m, 2)) / ((4*pow(m, 2)) + 5)))

        if alpha == 1:
            a = (N - 1)/(2 * m)
            b = (2 * m + 1) * (N - 1) / 4
            edf = np.exp(np.sqrt(np.log(a) * np.log(b)))

        if alpha == -1:
            if m == 1:
                edf = 2 * (N - 2)/(2.3 * N - 4.9)
            if m >= 2:
                edf = 5 * N**2 / (4 * m * (N + (3 * m)))

        if alpha == -2:
            a = (N - 2) / (m * (N - 3)**2)
            b = (N - 1)**2
            c = 3 * m * (N - 1)
            d = 4 * m**2
            edf = a * (b - c + d)

    else:
        edf = (N - 1)
        print("Noise type not recognized."
              " Defaulting to N - 1 degrees of freedom.")

    return edf

########################################################################
# end of ci.py
