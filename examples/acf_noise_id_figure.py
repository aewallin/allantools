import allantools as at
from matplotlib import pyplot as plt
import numpy as np


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
# N_points = pow(2,12)
N_points = pow(2, 13)
# N_points = pow(2,15)

nr = pow(2, 14)  # np.random.choice(nrrange)
brange = [0, -1, -2, -3, -4, -5]
pts = [[], [], [], [], [], []]  # 6 different b-values
for n in range(N_points):
    b = np.random.choice(brange)

    # note we stop at 30 pts series length
    taus = np.logspace(0, np.log10(nr/30.0), 200)
    taus = [int(t) for t in taus]
    taus = remove_duplicates(taus)
    af = np.random.choice(taus)
    ng.set_input(nr, qd, b)
    ng.generateNoise()
    x = ng.time_series

    alpha_int, alpha, rho, d = at.autocorr_noise_id(
        x, af, data_type="phase", dmin=0, dmax=2)
    pts[-1*b].append((af, alpha, alpha_int))
    print("%d / %d, a=%03d, alpha_id=%03d correct = %d ?" %
          (n, N_points, b+2, alpha_int, b+2 == alpha_int))

print("calc done")

colors = {0: 'r', -1: 'g', -2: 'b', -3: 'm', -4: 'k', -5: 'y'}
for (series, idx) in zip(pts, range(len(pts))):
    b = -idx
    af = [s[0] for s in series]
    alpha = [s[1] for s in series]
    alphai = [s[2] for s in series]
    plt.semilogx(af, alpha, "%s." % colors[b], label='b=%d' % b)
    # plt.semilogx( af, alphai, "%so"%colors[b],label='b=%d'%b) # integer alpha
print("plot done")

plt.ylim((-4, 3))
plt.xlim((.8, 1e3))

plt.xlabel('AF')
plt.ylabel('alpha')  # R(N)
plt.title('Power Law Noise Identification with Lag-1 ACF\nAW2018-01-14')
plt.legend(loc='best')
plt.grid()
plt.show()
