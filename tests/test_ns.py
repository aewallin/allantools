#!/usr/bin/python

import sys
sys.path.append("..")

import allantools as allan
from allantools import noise
import numpy

def _test( function, data, rate, taus):

    (taus2,devs2,errs2,ns2) = function(data, rate=rate, taus=taus)
    assert( len(taus2) == len(devs2) )
    assert( len(taus2) == len(errs2) )
    assert( len(taus2) == len(ns2) )
    for n in ns2:
        if n<= 1:
            print("test of ", function, " failed: ", n)
        assert( n > 1 ) # n should be 2 or more for each tau
    print("test_ns of function ",function, " OK.")
    
def test_ns():
    
    # this test asks for results at unreasonable tau-values
    # either zero, not an integer multiple of the data-interval
    # or too large, given the length of the dataset
    N = 500
    rate = 1.0
    phase_white = noise.white(N)
    taus_try = [x for x in numpy.logspace(0,4,4000)] # try insane tau values
    _test( allan.adev, phase_white, rate, taus_try)
    _test( allan.oadev, phase_white, rate, taus_try)
    _test( allan.mdev, phase_white, rate, taus_try)
    _test( allan.tdev, phase_white, rate, taus_try)
    _test( allan.hdev, phase_white, rate, taus_try)
    _test( allan.ohdev, phase_white, rate, taus_try)
    _test( allan.totdev, phase_white, rate, taus_try)
    _test( allan.mtie, phase_white, rate, taus_try)
    _test( allan.tierms, phase_white, rate, taus_try)

if __name__ == "__main__":
    test_ns()
