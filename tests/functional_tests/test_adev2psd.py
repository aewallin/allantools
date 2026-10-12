"""
Tests for adev2psd_piecewise_approx() / psd_piecewise_to_adev() /
noise.timmer_koenig_from_psd().

Single-slope cases have exact analytic PSD coefficients, which makes them
true independent checks (a round trip alone cannot detect a wrong intensity
coefficient, since both directions evaluate the same kernel integral).
"""

import numpy as np
import pytest

import allantools as at

TAUS = np.geomspace(1.0, 1.0e4, 9)


def test_white_fm_analytic():
    # sigma_y = A tau^-1/2  <->  Sy = h0 = 2 A^2  (alpha = 0)
    A = 1e-13
    f_nodes, _, h, alpha = at.adev2psd_piecewise_approx(A * TAUS ** -0.5, TAUS)
    assert h.size == 1
    np.testing.assert_allclose(alpha[0], 0.0, atol=1e-9)
    np.testing.assert_allclose(h[0], 2.0 * A ** 2, rtol=1e-8)


def test_flicker_fm_analytic():
    # flat ADEV  <->  flicker FM, h = sigma^2 / (2 ln 2)  (alpha = -1)
    sigma = 1e-13
    f_nodes, _, h, alpha = at.adev2psd_piecewise_approx(
        np.full_like(TAUS, sigma), TAUS)
    assert h.size == 1
    np.testing.assert_allclose(alpha[0], -1.0, atol=1e-9)
    np.testing.assert_allclose(h[0], sigma ** 2 / (2.0 * np.log(2.0)),
                               rtol=1e-8)


def test_near_boundary_tail_not_truncated():
    # For mu -> -2 the kernel integral diverges like 3/(8(2+mu)); a naive
    # quad() call plateaus near J ~ 3.3 instead, overestimating h by up to
    # ~100x. Check the implied J is at least the analytic tail term.
    mu = -1.999
    B = 1e-24
    adevs = np.sqrt(B * TAUS ** mu)
    _, _, h, alpha = at.adev2psd_piecewise_approx(adevs, TAUS)
    J_implied = B / (2.0 * h[0] * np.pi ** mu)
    assert J_implied >= 3.0 / (8.0 * (2.0 + mu))   # = 375; naive quad gives ~3.3


def test_round_trip_multi_segment():
    taus = np.array([1, 10, 1e2, 1e3, 5e3, 1e4])
    adevs = np.array([1e-11, 2e-11, 5.5e-11, 1.3e-10, 2.85e-10, 4e-10])
    f_nodes, _, h, alpha = at.adev2psd_piecewise_approx(adevs, taus)
    back = at.psd_piecewise_to_adev(h, alpha, f_nodes, taus)
    # reconstruction error is proportional to local curvature; <10% here
    np.testing.assert_allclose(back, adevs, rtol=0.10)


def test_out_of_range_slopes_rejected():
    t = np.geomspace(1, 1e4, 5)
    with pytest.raises(ValueError):
        at.adev2psd_piecewise_approx(1e-13 * t ** -1.5, t)     # mu = -3
    with pytest.raises(ValueError):
        at.adev2psd_piecewise_approx(1e-13 * t ** 1.5, t)      # mu = +3
    # the Hadamard branch accepts drift-like slopes up to mu < 4
    _, _, h, alpha = at.adev2psd_piecewise_approx(
        1e-13 * t ** 1.5, t, vartype="hdev")
    np.testing.assert_allclose(alpha[0], -4.0, atol=1e-9)


def test_single_segment_f_nodes_accepted():
    # single-slope inputs return two 1/tau endpoint "nodes"; both consumers
    # must accept that form and treat it as "no breaks"
    A = 1e-13
    f_nodes, _, h, alpha = at.adev2psd_piecewise_approx(A * TAUS ** -0.5, TAUS)
    assert f_nodes.size == 2 and h.size == 1
    back = at.psd_piecewise_to_adev(h, alpha, f_nodes, TAUS)
    np.testing.assert_allclose(back, A * TAUS ** -0.5, rtol=1e-6)
    x1 = at.noise.timmer_koenig_from_psd(f_nodes, h, alpha, 4096.0, 1.0,
                                         output="phase", seed=42)
    x2 = at.noise.timmer_koenig_from_psd(np.array([]), h, alpha, 4096.0, 1.0,
                                         output="phase", seed=42)
    np.testing.assert_array_equal(x1, x2)


if __name__ == "__main__":
    test_white_fm_analytic()
    test_flicker_fm_analytic()
    test_near_boundary_tail_not_truncated()
    test_round_trip_multi_segment()
    test_out_of_range_slopes_rejected()
    test_single_segment_f_nodes_accepted()
    print("all adev2psd tests passed")
