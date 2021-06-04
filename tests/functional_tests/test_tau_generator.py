#!/usr/bin/python

import allantools as at
import numpy as np
import pytest


N = 128
d = np.random.random(N)
r = 1.0

expected_all = np.arange(1, 43)
expected_octave = [1.,   2.,   4.,   8.,  16.,  32.]
expected_decade = [1., 2., 4., 10., 20., 40.]
expected_reduced_10 = [1., 2., 3., 4., 5., 7., 8., 11., 14., 17.,
                       22., 28., 35.]
expected_reduced_2 = [1., 5., 17.]


def test_tau_generator_empty():
    (taus_used, adev, adeverror, adev_n) = at.adev(d)
    np.testing.assert_allclose(taus_used, expected_octave)


def test_tau_generator_empty_list():
    (taus_used, adev, adeverror, adev_n) = at.adev(d, taus=[])
    np.testing.assert_allclose(taus_used, expected_octave)


def test_tau_generator_all():
    (taus_used, adev, adeverror, adev_n) = at.adev(d, rate=r, taus="all")
    np.testing.assert_allclose(taus_used, expected_all)


def test_tau_generator_octave():
    (taus_used, adev, adeverror, adev_n) = at.adev(d, rate=r, taus="octave")
    np.testing.assert_allclose(taus_used, expected_octave)


def test_tau_generator_decade():
    (taus_used, adev, adeverror, adev_n) = at.adev(d, rate=r, taus="decade")
    np.testing.assert_allclose(taus_used, expected_decade)


def test_tau_generator_1234():
    wanted_taus = [1, 2, 3, 4]
    (taus_used, adev, adeverror, adev_n) = at.adev(d, rate=r, taus=wanted_taus)
    np.testing.assert_allclose(taus_used, wanted_taus)


def test_tau_generator_numpy1234():
    wanted_taus = np.array([1, 2, 3, 4])
    (taus_used, adev, adeverror, adev_n) = at.adev(d, rate=r, taus=wanted_taus)
    np.testing.assert_allclose(taus_used, wanted_taus)


def test_tau_reduction_10():
    (ms, taus) = at.allantools.tau_reduction(ms=expected_all, rate=r,
                                             n_per_decade=10)
    print(ms, taus)
    np.testing.assert_allclose(expected_reduced_10, ms)


def test_tau_reduction_2():
    (ms, taus) = at.allantools.tau_reduction(ms=expected_all, rate=r,
                                             n_per_decade=2)
    print(ms, taus)
    np.testing.assert_allclose(expected_reduced_2, ms)


def test_zero_rate():
    with pytest.raises(RuntimeError):
        at.adev(d, rate=0.0)


if __name__ == "__main__":
    test_tau_generator_empty()
    test_tau_generator_all()
    test_tau_generator_octave()
    test_tau_generator_decade()
    test_tau_generator_1234()
    test_tau_generator_numpy1234()
    test_tau_reduction_10()
    test_tau_reduction_2()
    test_zero_rate()
