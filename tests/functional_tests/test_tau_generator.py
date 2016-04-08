#!/usr/bin/python

import sys
import allantools as at
import numpy as np
import pytest
#sys.path.append("..")


N = 128
d = np.random.random(N)
r = 1.0

expected_all = np.arange(1, 43)
expected_octave = [1.,   2.,   4.,   8.,  16.,  32.]
expected_decade = [1., 2., 4., 10., 20., 40.]


def test_tau_generator_empty():
    (taus_used, adev, adeverror, adev_n) = at.adev(d)
    np.testing.assert_array_equal(taus_used, expected_octave)

def test_tau_generator_empty_list():
    (taus_used, adev, adeverror, adev_n) = at.adev(d, taus=[])
    np.testing.assert_array_equal(taus_used, expected_octave)

def test_tau_generator_all():
    (taus_used, adev, adeverror, adev_n) = at.adev(d, rate=r, taus="all")
    np.testing.assert_array_equal(taus_used, expected_all)


def test_tau_generator_octave():
    (taus_used, adev, adeverror, adev_n) = at.adev(d, rate=r, taus="octave")
    np.testing.assert_array_equal(taus_used, expected_octave)


def test_tau_generator_decade():
    (taus_used, adev, adeverror, adev_n) = at.adev(d, rate=r, taus="decade")
    np.testing.assert_array_equal(taus_used, expected_decade)


def test_tau_generator_1234():
    wanted_taus = [1, 2, 3, 4]
    (taus_used, adev, adeverror, adev_n) = at.adev(d, rate=r, taus=wanted_taus)
    np.testing.assert_array_equal(taus_used, wanted_taus)

def test_tau_generator_numpy1234():
    wanted_taus = np.array([1, 2, 3, 4])
    (taus_used, adev, adeverror, adev_n) = at.adev(d, rate=r, taus=wanted_taus)
    np.testing.assert_array_equal(taus_used, wanted_taus)

if __name__ == "__main__":
    test_tau_generator_empty()
    test_tau_generator_all()
    test_tau_generator_octave()
    test_tau_generator_decade()
    test_tau_generator_1234()
