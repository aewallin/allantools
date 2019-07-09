"""
    Test for allantools (https://github.com/aewallin/allantools)

    Auto-AF mode: new afs are added automatically as we add more phase points.
"""
import math
import time
import sys
import pytest

import allantools as at
import numpy

n_pts = pow(2,12)
x_series = at.noise.white(n_pts)

def test_oadev_rt_autoAF():
    print("time series length %d"%n_pts)
    oadev_rt = at.realtime.oadev_realtime(tau0=1.0, auto_afs=True)
    for x in x_series:
        oadev_rt.add_phase(x)
    print("number of taus %d"%len(oadev_rt.devs()))
    print("taus: ",oadev_rt.taus())
    print("devs: ",oadev_rt.devs())

def test_ohdev_rt_autoAF():
    print("time series length %d"%n_pts)
    dev_rt = at.realtime.ohdev_realtime(tau0=1.0, auto_afs=True)
    for x in x_series:
        dev_rt.add_phase(x)
    print("number of taus %d"%len(dev_rt.devs()))
    print("taus: ",dev_rt.taus())
    print("devs: ",dev_rt.devs())

def test_tdev_rt_autoAF():
    print("time series length %d"%n_pts)
    dev_rt = at.realtime.tdev_realtime(tau0=1.0, auto_afs=True)
    for x in x_series:
        dev_rt.add_phase(x)
    print("number of taus %d"%len(dev_rt.devs()))
    print("taus: ",dev_rt.taus())
    print("devs: ",dev_rt.devs())

if __name__ == "__main__":
    test_oadev_rt_autoAF()
    test_ohdev_rt_autoAF()
    test_tdev_rt_autoAF()
