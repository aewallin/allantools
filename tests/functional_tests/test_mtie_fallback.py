"""The low-memory path must use the same complete MTIE windows."""

from unittest.mock import patch

import numpy as np
import pytest

import allantools


@pytest.mark.parametrize("phase", [
    [0, 10, 0, 0],
    [0, -10, 0, 0],
    [10, 0, 0, 0],
    [0, 0, 0, 10],
    [0, 1, 2, 3, 4, 5],
    [5, 4, 3, 2, 1, 0],
    [0, 5, 5, -5, -5, 0],
    [3, 3, 3, 3],
])
def test_mtie_fallback_windows(phase):
    phase = np.asarray(phase, dtype=float)
    rate = 2.5
    intervals = np.arange(1, len(phase) - 1)
    taus = intervals / rate
    expected = [
        max(max(phase[i:i + m + 1]) - min(phase[i:i + m + 1])
            for i in range(len(phase) - m))
        for m in intervals
    ]
    normal = allantools.mtie(phase, rate=rate, taus=taus)
    # Trigger the documented NumPy array-size limit without a huge allocation.
    with patch("allantools.allantools.mtie_rolling_window",
               side_effect=ValueError("array is too big")):
        actual = allantools.mtie(phase, rate=rate, taus=taus)
    np.testing.assert_array_equal(actual[0], taus)
    np.testing.assert_array_equal(actual[1], expected)
    np.testing.assert_array_equal(actual[3], len(phase) - intervals)
    np.testing.assert_allclose(actual[2], expected / np.sqrt(actual[3]))
    for result, reference in zip(actual, normal):
        np.testing.assert_array_equal(result, reference)


def test_mtie_fallback_frequency_input():
    # Zero-mean frequency samples integrate to [0, 10, 0, 0] seconds.
    frequency = [20., -20., 0.]
    normal = allantools.mtie(frequency, rate=2, data_type="freq", taus=[0.5])
    with patch("allantools.allantools.mtie_rolling_window",
               side_effect=ValueError("array is too big")):
        actual = allantools.mtie(
            frequency, rate=2, data_type="freq", taus=[0.5])
    np.testing.assert_array_equal(actual[1], [10.])
    for result, reference in zip(actual, normal):
        np.testing.assert_array_equal(result, reference)
