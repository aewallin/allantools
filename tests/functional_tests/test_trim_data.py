#!/usr/bin/python
import allantools as at
import numpy as np

N = 128
d = np.random.random(N)

array_with_holes = np.array([np.nan, np.nan, 1, 2, 3.4, 2.3,
                             np.nan, 9, 12, 45, np.nan])
expected_trimmed_array = np.array([1, 2, 3.4, 2.3, np.nan, 9, 12, 45])


def test_full_array():
    """ Should be the same """
    output = at.trim_data(d)
    np.testing.assert_array_equal(d, output)


def test_array_with_holes():
    """ Only trim the array, do not remove NaN inside it """
    output = at.trim_data(array_with_holes)
    np.testing.assert_array_equal(expected_trimmed_array, output)

if __name__ == "__main__":
    test_full_array()
    test_array_with_holes()
