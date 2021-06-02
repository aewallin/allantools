import allantools as at
from matplotlib import pyplot as plt
import numpy as np


def rn_noise_id(x, af, rate):
    """ R(n) ratio for noise identification

    """
    # print "len(x) ", len(x), "oadev",af*rate
    (taus, devs, errs, ns) = at.adev(
        x, taus=[af*rate], data_type='phase', rate=rate)
    # print devs
    oadev_x = devs[0]
    (mtaus, mdevs, errs, ns) = at.mdev(
        x, taus=[af*rate], data_type='phase', rate=rate)
    mdev_x = mdevs[0]
    rn = pow(mdev_x/oadev_x, 2)
    return rn


def rn(af, b):
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


plt.figure(figsize=(12, 10))
# N_points = pow(2,6)
# N_points = pow(2,10)
N_points = pow(2, 12)
# N_points = pow(2,15)

nr = pow(2, 14)  # np.random.choice(nrrange)
brange = [0, -1, -2, -3, -4, -5]
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
    rn_val = rn_noise_id(x, af, rate)
    pts.append((b, af, rn_val))

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
    b1s = [rn(af, b) for af in taus]
    plt.loglog(taus, b1s, "%s-" % colors[b], label='R(n) b=%d' % b)

    rnb = [np.sqrt(rn(af, b)*rn(af, b-1)) for af in taus]

    # bound = [ b1_boundary(b, N ) for N in Ns ]
    plt.loglog(taus, rnb, "%s--" %
               colors[b], label='b=%d...%d boundary)' % (b, b-1))


plt.xlabel('AF')
plt.ylabel('R(N)')
plt.title('Power Law Noise Identification with R(n)\nAW2018-01-13')
plt.legend(loc='best')
plt.grid()
plt.show()
