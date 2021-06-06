import allantools as allan
import numpy
# import math
import matplotlib.pyplot as plt

import sys
sys.path.append("..")
sys.path.append("../..")  # hack to import from parent directory
# remove if you have allantools installed in your python path


def read_datafile(filename, column=1):
    p = []
    with open(filename) as f:
        for line in f:

            if line.startswith("#"):  # skip comments
                pass
            else:
                line = line.split()
                p.append(float(line[column]))
    return p


def to_fractional(flist, f0):
    out = []
    for f in flist:
        out.append(f/float(f0) - 1.0)
    return out


fname = "ocxo_frequency.txt"
f10MHz = read_datafile(fname, column=0)

f = to_fractional(f10MHz, 10e6)  # convert to fractional frequency
# log-spaced tau values from 10s and upwards
my_taus = numpy.logspace(1, 5, 40)
rate = 1/float(10)  # data collected with 10s gate time

(oadev_taus, oadev_devs, oadev_errs, ns) = allan.oadev(
    f, data_type='freq', rate=rate, taus=my_taus)
(mdev_taus, mdev_devs, mdev_errs, ns) = allan.mdev(
    f, data_type='freq', rate=rate, taus=my_taus)
(hdev_taus, hdev_devs, hdev_errs, ns) = allan.hdev(
    f, data_type='freq', rate=rate, taus=my_taus)
(ohdev_taus, ohdev_devs, ohdev_errs, ns) = allan.ohdev(
    f, data_type='freq', rate=rate, taus=my_taus)

plt.subplot(111, xscale="log", yscale="log")

plt.errorbar(oadev_taus, oadev_devs, yerr=oadev_errs, label='OADEV')
plt.errorbar(mdev_taus, mdev_devs, yerr=mdev_errs, label='MDEV')
plt.errorbar(hdev_taus, hdev_devs, yerr=hdev_errs, label='HDEV')
plt.errorbar(ohdev_taus, ohdev_devs, yerr=ohdev_errs, label='OHDEV')

plt.xlabel('Taus (s)')
plt.ylabel('ADEV')


plt.legend()
plt.grid()
plt.show()
