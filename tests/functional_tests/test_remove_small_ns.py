import pytest
import allantools.allantools as at
import numpy as np

minimum_ns = 1

taus = np.arange(1, 499)
devs = np.random.rand(len(taus))
deverrs_l = np.random.rand(len(taus))
deverrs_h = np.random.rand(len(taus))
ns = np.zeros(len(taus))
ns[0] = len(ns)
last_i = 0
for i in range(1, len(ns)):
    ns[i] = max(ns[i-1] // 2, 1.0)
    if ns[i] > minimum_ns:
        last_i = i + 1


def test_4params():
    (o_taus, o_devs, o_deverrs, o_ns) = at.remove_small_ns(
        taus, devs, deverrs_l, ns)
    np.testing.assert_array_equal(o_taus, taus[:last_i])
    np.testing.assert_array_equal(o_devs, devs[:last_i])
    np.testing.assert_array_equal(o_deverrs, deverrs_l[:last_i])
    np.testing.assert_array_equal(o_ns, ns[:last_i])

def test_5params():
    (o_taus, o_devs, [o_deverrs_l, o_deverrs_h], o_ns) = at.remove_small_ns(
        taus, devs, [deverrs_l, deverrs_h], ns)
    np.testing.assert_array_equal(o_taus, taus[:last_i])
    np.testing.assert_array_equal(o_devs, devs[:last_i])
    np.testing.assert_array_equal(o_deverrs_l, deverrs_l[:last_i])
    np.testing.assert_array_equal(o_deverrs_h, deverrs_h[:last_i])
    np.testing.assert_array_equal(o_ns, ns[:last_i])

if __name__ == "__main__":
    pytest.main()
