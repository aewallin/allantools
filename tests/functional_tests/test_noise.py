#!/usr/bin/python

import sys
sys.path.append("..")

from allantools import noise
import numpy
import pytest


def test_noise():

    N = 500
    #rate = 1.0
    w = noise.white(N)
    b = noise.brown(N)
    v = noise.violet(N)
    p = noise.pink(N)
    
    # check output length
    assert len(w) == N
    assert len(b) == N
    assert len(v) == N
    assert len(p) == N
    # check output type
    for x in [w, b, v, p]:
        assert type(x) == numpy.ndarray, "%s is not numpy.ndarray" % (type(x))
        
if __name__ == "__main__":
    test_noise()
