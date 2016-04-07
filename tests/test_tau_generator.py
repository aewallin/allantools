#!/usr/bin/python

import sys
sys.path.append("..")

import allantools as at
import numpy as np
import pytest


def test_tau_generator():
    N=1024
    d = np.random.random(N)
    r = 1.0
    t = at.adev(d)
    t = at.adev(d, rate=r, taus="all")
    t = at.adev(d, rate=r, taus="octave")
    t = at.adev(d, rate=r, taus="decade")
    t = at.adev(d, rate=r, taus=[1,2,3,4])

if __name__ == "__main__":
    test_tau_generator()
