import allantools as at
from matplotlib import pyplot as plt
import numpy as np


def b1_noise_id(x, af, rate):
    """ B1 ratio for noise identification

        ratio of Standard Variace to AVAR
    """
    (taus, devs, errs, ns) = at.adev(
        x, taus=[af*rate], data_type="phase", rate=rate)
    oadev_x = devs[0]
    y = np.diff(x)
    y_cut = np.array(y[:len(y)-(len(y) % af)])  # cut to length
    assert len(y_cut) % af == 0
    y_shaped = y_cut.reshape((int(len(y_cut)/af), af))
    y_averaged = np.average(y_shaped, axis=1)  # average

    var = np.var(y_averaged, ddof=1)
    return var/pow(oadev_x, 2.0)


def b1(N, mu):
    """ Expected B1 ratio for given time-series length N and exponent mu

        The exponents are defined as
        S_y(f) = h_a f^alpha
        S_x(f) = g_b f^b
        bias = const * tau^mu


        and relate to eachother by:
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
    if mu == 2:
        return float(N)*(float(N)+1.0)/6.0
        # up = N*(1.0-pow(N, mu))
        # down = 2*(N-1.0)*(1-pow(2.0, mu))
        # return up/down
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


def b_to_mu(b):
    a = b+2
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


ng = at.Noise()
nr = pow(2, 14)
qd = 2e-20
b = -1


def remove_duplicates(l):
    s = []
    for i in l:
        if i not in s:
            s.append(i)
    return s


def b1_boundary(b_hi, N):
    b_lo = b_hi-1
    b1_lo = b1(N, b_to_mu(b_lo))
    b1_hi = b1(N, b_to_mu(b_hi))
    if b1_lo >= -4:
        return np.sqrt(b1_lo*b1_hi)  # geometric mean
    else:
        return 0.5*(b1_lo+b1_hi)  # arithemtic mean


plt.figure()
N_points = pow(2, 6)
# N_points = pow(2,10)
N_points = pow(2, 12)
# N_points = pow(2,15)

nr = pow(2, 14)  # np.random.choice(nrrange)
brange = [0, -1, -2, -3, -4, -5]
# nrrange = [pow(2,10), pow(2,14), pow(2,16)]
pts = []
for n in range(N_points):
    b = np.random.choice(brange)

    taus = np.logspace(0, np.log10(nr/4.0), 200)
    taus = [int(t) for t in taus]
    taus = remove_duplicates(taus)
    af = np.random.choice(taus)
    ng.set_input(nr, qd, b)
    ng.generateNoise()
    x = ng.time_series

    rate = 1.0
    b1_val = b1_noise_id(x, af, rate)
    pts.append((b, af, b1_val))
print("calc done")
colors = {0: 'r', -1: 'g', -2: 'b', -3: 'm', -4: 'k', -5: 'y'}
for p in pts:
    plt.loglog(p[1], p[2], "%s." % colors[p[0]])
print("plot done")
# theory lines
taus = np.logspace(0, np.log10(nr/4.0), 20)
taus = [int(t) for t in taus]
taus = remove_duplicates(taus)

for b in brange:
    x = range(nr)
    Ns = [np.floor(nr/tau) for tau in taus]
    b1s = [b1(N, b_to_mu(b)) for N in Ns]
    plt.loglog(taus, b1s, "%s-" % colors[b], label='B1(b=%d)' % b)

    bound = [b1_boundary(b, N) for N in Ns]
    plt.loglog(taus, bound, "%s--" %
               colors[b], label='b=%d...%d boundary)' % (b, b-1))


plt.xlabel('AF')
plt.ylabel('B1')
plt.title('Power Law Noise Identification with B1(N, mu)\nAW2018-01-13')
plt.legend()
plt.grid()
plt.show()
