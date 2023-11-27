import allantools
import numpy as np
import pytest  # noqa

import sys
sys.path.append("..")

# Test values have been calculated from the reference documents
# indicated in the source


def test_prc_tdev():
    assert np.isclose(allantools.mask.prc_tdev(50),  3e-9)
    assert np.isclose(allantools.mask.prc_tdev(500),  15e-9)
    assert np.isclose(allantools.mask.prc_tdev(5000),  30e-9)


def test_prc_mtie():
    assert np.isclose(allantools.mask.prc_mtie(50),  0.03875e-6)
    assert np.isclose(allantools.mask.prc_mtie(5000),  0.34e-6)


def test_eprtc_tdev():
    assert np.isclose(allantools.mask.eprtc_tdev(5), 1e-9)
    assert np.isclose(allantools.mask.eprtc_tdev(50), .001666665e-9)
    assert np.isclose(allantools.mask.eprtc_tdev(100e3), 3.33333e-5*1.0e-9*100e3)
    assert np.isclose(allantools.mask.eprtc_tdev(301e3), 10e-9)


def test_eprtc_mtie():
    assert np.isclose(allantools.mask.eprtc_mtie(.5), 4e-9)
    assert np.isclose(allantools.mask.eprtc_mtie(5), 4.4457e-9)
    assert np.isclose(allantools.mask.eprtc_mtie(250), 15.009375e-9)
    assert np.isclose(allantools.mask.eprtc_mtie(500000), 30e-9)


def test_prtcA_tdev():
    assert np.isclose(allantools.mask.prtcA_tdev(5), 3e-9)
    assert np.isclose(allantools.mask.prtcA_tdev(500), 15e-9)
    assert np.isclose(allantools.mask.prtcA_tdev(5000), 30e-9)


def test_prtcB_tdev():
    assert np.isclose(allantools.mask.prtcB_tdev(5), 1e-9)
    assert np.isclose(allantools.mask.prtcB_tdev(250), 2.5e-9)
    assert np.isclose(allantools.mask.prtcB_tdev(5000), 5e-9)


def test_prtcA_mtie():
    assert np.isclose(allantools.mask.prtcA_mtie(5), 0.026375e-6)
    assert np.isclose(allantools.mask.prtcA_mtie(500), .1e-6)


def test_prtcB_mtie():
    assert np.isclose(allantools.mask.prtcB_mtie(5), 0.026375e-6)
    assert np.isclose(allantools.mask.prtcB_mtie(500), .04e-6)


if __name__ == "__main__":
    test_prc_tdev()
    test_prc_mtie()
    test_eprtc_tdev()
    test_eprtc_mtie()
    test_prtcA_tdev()
    test_prtcA_mtie()
    test_prtcB_tdev()
    test_prtcB_mtie()
