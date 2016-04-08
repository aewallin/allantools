#!/usr/bin/python

import sys
import allantools as at
import numpy as np

N = 128
d = np.random.random(N)
r = 1.0

expected_all = np.arange(1, 43)
expected_octave = [1.,   2.,   4.,   8.,  16.,  32.]
expected_decade = [1., 2., 4., 10., 20., 40.]

ab = at.noise.white(N)
bc = at.noise.white(N)
ca = at.noise.white(N)

def test_3ch_1():
    (t,d,e,n) = at.three_cornered_hat_phase(ab, bc, ca, rate=r, 
                                       taus='decade', function=at.oadev)
    np.testing.assert_array_equal(t, expected_decade)

if __name__ == "__main__":
    test_3ch_1()
