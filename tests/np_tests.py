#!/usr/bin/env python
"""
np_tests.py
-----------

Migration from lists to numpy arrays - test routines


1) all functions recreated with _np suffix, e.g.
    tdev_phase -> tdev_phase_np
2) Test all functions against their originals
3) Speed tests
4) Remove originals and rename without _np suffix

"""

import numpy as np
import pylab as plt
import allantools.allantools as alt
import time


if __name__ == "__main__":

    #######################
    # OADEV_PHASE()
    #######################
    print "\ntesting oadev_phase()"
    data = np.random.random(10000)
    taus = [1, 3, 5, 16, 128]
    rates = [1, 20, 10.7]
    strides = [1, 10, 7]

    for rate in rates:
        for stride in strides:
            #print "TAU: %i, RATE: %2.2f, STRIDE: %i" % (tau, rate, stride)
            o_taus, o_dev, o_err, o_n = alt.oadev_phase(data, rate, taus)
            o_taus_, o_dev_, o_err_, o_n_ = alt.oadev_phase_np(data, rate, taus)

            assert np.allclose(o_taus, o_taus_)
            assert np.allclose(o_dev, o_dev_)
            assert np.allclose(o_err, o_err_)

    stride = 1
    tau = 16
    rate = 2.1
    data = np.random.random(100000)
    t1 = time.time()
    o_taus, o_dev, o_err, o_n = alt.oadev_phase(data, rate, taus)
    t2 = time.time()
    t3 = time.time()
    o_taus_, o_dev_, o_err_, o_n_ = alt.oadev_phase_np(data, rate, taus)
    t4 = time.time()

    assert np.allclose(o_taus, o_taus_)
    assert np.allclose(o_dev, o_dev_)
    assert np.allclose(o_err, o_err_)
    print "Original: %2.3fs" % (t2 - t1)
    print "New:      %2.3fs" % (t4 - t3)
    print "Speedup:  %2.2fx" % ((t2 - t1) / (t4 - t3))

    #######################
    # OADEV()
    #######################
    print "\ntesting oadev()"
    data = np.random.random(10000)
    taus = [1, 3, 5, 16, 128]
    rates = [1, 20, 10.7]
    strides = [1, 10, 7]

    for rate in rates:
        for stride in strides:
            #print "TAU: %i, RATE: %2.2f, STRIDE: %i" % (tau, rate, stride)
            o_taus, o_dev, o_err, o_n = alt.oadev(data, rate, taus)
            o_taus_, o_dev_, o_err_, o_n_ = alt.oadev_np(data, rate, taus)

            assert np.allclose(o_taus, o_taus_)
            assert np.allclose(o_dev, o_dev_)
            assert np.allclose(o_err, o_err_)

    stride = 1
    tau = 16
    rate = 2.1
    data = np.random.random(100000)
    t1 = time.time()
    o_taus, o_dev, o_err, o_n = alt.oadev(data, rate, taus)
    t2 = time.time()
    t3 = time.time()
    o_taus_, o_dev_, o_err_, o_n_ = alt.oadev_np(data, rate, taus)
    t4 = time.time()

    assert np.allclose(o_taus, o_taus_)
    assert np.allclose(o_dev, o_dev_)
    assert np.allclose(o_err, o_err_)
    print "Original: %2.3fs" % (t2 - t1)
    print "New:      %2.3fs" % (t4 - t3)
    print "Speedup:  %2.2fx" % ((t2 - t1) / (t4 - t3))

    #######################
    # ADEV_PHASE()
    #######################
    print "\ntesting adev_phase()"
    data = np.random.random(10000)
    taus = [1, 3, 5, 16, 128]
    rates = [1, 20, 10.7]
    strides = [1, 10, 7]

    for rate in rates:
        for stride in strides:
            #print "TAU: %i, RATE: %2.2f, STRIDE: %i" % (tau, rate, stride)
            o_taus, o_dev, o_err, o_n = alt.adev_phase(data, rate, taus)
            o_taus_, o_dev_, o_err_, o_n_ = alt.adev_phase_np(data, rate, taus)

            assert np.allclose(o_taus, o_taus_)
            assert np.allclose(o_dev, o_dev_)
            assert np.allclose(o_err, o_err_)

    stride = 1
    tau = 16
    rate = 2.1
    data = np.random.random(100000)
    t1 = time.time()
    o_taus, o_dev, o_err, o_n = alt.adev_phase(data, rate, taus)
    t2 = time.time()
    t3 = time.time()
    o_taus_, o_dev_, o_err_, o_n_ = alt.adev_phase_np(data, rate, taus)
    t4 = time.time()

    assert np.allclose(o_taus, o_taus_)
    assert np.allclose(o_dev, o_dev_)
    assert np.allclose(o_err, o_err_)
    print "Original: %2.3fs" % (t2 - t1)
    print "New:      %2.3fs" % (t4 - t3)
    print "Speedup:  %2.2fx" % ((t2 - t1) / (t4 - t3))

    #######################
    # ADEV()
    #######################
    print "\ntesting adev()"
    data = np.random.random(10000)
    taus = [1, 3, 5, 16, 128]
    rates = [1, 20, 10.7]
    strides = [1, 10, 7]

    for rate in rates:
        for stride in strides:
            #print "TAU: %i, RATE: %2.2f, STRIDE: %i" % (tau, rate, stride)
            o_taus, o_dev, o_err, o_n = alt.adev(data, rate, taus)
            o_taus_, o_dev_, o_err_, o_n_ = alt.adev_np(data, rate, taus)

            assert np.allclose(o_taus, o_taus_)
            assert np.allclose(o_dev, o_dev_)
            assert np.allclose(o_err, o_err_)

    stride = 1
    tau = 16
    rate = 2.1
    data = np.random.random(100000)
    t1 = time.time()
    o_taus, o_dev, o_err, o_n = alt.adev(data, rate, taus)
    t2 = time.time()
    t3 = time.time()
    o_taus_, o_dev_, o_err_, o_n_ = alt.adev_np(data, rate, taus)
    t4 = time.time()

    assert np.allclose(o_taus, o_taus_)
    assert np.allclose(o_dev, o_dev_)
    assert np.allclose(o_err, o_err_)
    print "Original: %2.3fs" % (t2 - t1)
    print "New:      %2.3fs" % (t4 - t3)
    print "Speedup:  %2.2fx" % ((t2 - t1) / (t4 - t3))


    #######################
    # CALC_ADEV_PHASE()
    #######################
    print "\ntesting calc_adev_phase()"
    data = np.random.random(10000)
    taus = [1, 3, 5, 16, 128]
    rates = [1, 20, 10.7]
    strides = [1, 10, 7]

    for tau in taus:
        for rate in rates:
            for stride in strides:
                #print "TAU: %i, RATE: %2.2f, STRIDE: %i" % (tau, rate, stride)
                mj = tau
                dev, deverr, n = alt.calc_adev_phase(data, rate, mj, stride)
                dev_, deverr_, n_ = alt.calc_adev_phase_np(data, rate, mj, stride)

                assert np.isclose(dev, dev_)
                assert np.isclose(n, n_)
                assert np.isclose(deverr, deverr_)

    stride = 1
    tau = 16
    rate = 2.0
    t1 = time.time()
    dev, deverr, n = alt.calc_adev_phase(data, rate, mj, stride)
    t2 = time.time()
    t3 = time.time()
    dev_, deverr_, n_ = alt.calc_adev_phase_np(data, rate, mj, stride)
    t4 = time.time()

    assert np.isclose(dev, dev_)
    assert np.isclose(n, n_)
    assert np.isclose(deverr, deverr_)
    print "Original: %2.3fs" % (t2 - t1)
    print "New:      %2.3fs" % (t4 - t3)
    print "Speedup:  %2.2fx" % ((t2 - t1) / (t4 - t3))

    #######################
    # TAU_M()
    #######################
    print "\ntesting tau_m()"
    taus = [1, 2, 4, 8, 16, -4, 10000, -3.1, 3.141]
    data = np.random.random(10000)
    rates = [1, 2, 7.1, 123.12]

    for rate in rates:
        m, taus2 = alt.tau_m(data, rate, taus)
        m_, taus2_ = alt.tau_m_np(data, rate, taus)
        assert np.allclose(m, m_)
        assert np.allclose(taus2, taus2_)

    taus = np.random.randint(low=-100, high=10000, size=(10000,))
    rate = 1.234
    t1 = time.time()
    m, taus2 = alt.tau_m(data, rate, taus)
    t2 = time.time()
    t3 = time.time()
    m_, taus2_ = alt.tau_m_np(data, rate, taus)
    t4 = time.time()
    assert np.allclose(m, m_)
    assert np.allclose(taus2, taus2_)
    print "Original: %2.3fs" % (t2 - t1)
    print "New:      %2.3fs" % (t4 - t3)
    print "Speedup:  %2.2fx" % ((t2 - t1) / (t4 - t3))


    #######################
    # FREQUENCY2PHASE()
    #######################
    print "\ntesting frequency2phase()"
    freqdata = np.random.random(10000)

    rates = [1, 2, 7.1, 123.12]

    for rate in rates:
        phase  = alt.frequency2phase(freqdata, rate)
        phase_ = alt.frequency2phase_np(freqdata, rate)
        assert len(phase) == len(phase_)
        assert np.allclose(phase, phase_)

    freqdata = np.random.random(100000)
    t1 = time.time()
    phase  = alt.frequency2phase(freqdata, rate)
    t2 = time.time()
    t3 = time.time()
    phase_ = alt.frequency2phase_np(freqdata, rate)
    t4 = time.time()
    print "Original: %2.3fs" % (t2 - t1)
    print "New:      %2.3fs" % (t4 - t3)
    print "Speedup:  %2.2fx" % ((t2 - t1) / (t4 - t3))

    #######################
    # TDEV_PHASE()
    #######################
    print "\ntesting tdev_phase()"
    rate      = 1.0          # 1000 Hz sample rate
    obs_s     = 10000             # 1 hour
    n_samples = rate * obs_s
    t = np.arange(0, n_samples)
    phase  = np.random.random(n_samples) + np.sin(t / n_samples)

    taus = [4]

    t1 = time.time()
    taus2, td, tde, ns = alt.tdev_phase(phase, rate, taus)
    t2 = time.time()
    t3 = time.time()
    taus2_, td_, tde_, ns_ = alt.tdev_phase_np(phase, rate, taus)
    t4 = time.time()

    assert np.allclose(taus2, taus2_)
    assert np.allclose(td, td_)
    assert np.allclose(tde, tde_)
    assert np.allclose(ns, ns_)

    print "Original: %2.3fs" % (t2 - t1)
    print "New:      %2.3fs" % (t4 - t3)
    print "Speedup:  %2.2fx" % ((t2 - t1) / (t4 - t3))

    #######################
    # TDEV()
    #######################
    print "\ntesting tdev()"
    rate      = 2.0          # 1000 Hz sample rate
    obs_s     = 32768            # 1 hour
    n_samples = rate * obs_s
    t = np.arange(0, n_samples)
    phase  = np.random.random(n_samples) + np.sin(t / n_samples)

    taus = [1, 2, 4]

    t1 = time.time()
    taus2, td, tde, ns = alt.tdev(phase, rate, taus)
    t2 = time.time()

    t3 = time.time()
    taus2_, td_, tde_, ns_ = alt.tdev_np(phase, rate, taus)
    t4 = time.time()

    assert np.allclose(taus2, taus2_)
    assert np.allclose(td, td_)
    assert np.allclose(tde, tde_)
    assert np.allclose(ns, ns_)

    print "Original: %2.3fs" % (t2 - t1)
    print "New:      %2.3fs" % (t4 - t3)
    print "Speedup:  %2.2fx" % ((t2 - t1) / (t4 - t3))

    #######################
    # MDEV_PHASE()
    #######################
    print "\ntesting mdev_phase()"
    rate      = 1.0          # 1000 Hz sample rate
    obs_s     = 10000             # 1 hour
    n_samples = rate * obs_s
    t = np.arange(0, n_samples)
    phase  = np.random.random(n_samples) + np.sin(t / n_samples)

    taus = [4]

    t1 = time.time()
    taus2, td, tde, ns = alt.mdev_phase(phase, rate, taus)
    t2 = time.time()
    t3 = time.time()
    taus2_, td_, tde_, ns_ = alt.mdev_phase_np(phase, rate, taus)
    t4 = time.time()

    assert np.allclose(taus2, taus2_)
    assert np.allclose(td, td_)
    assert np.allclose(tde, tde_)
    assert np.allclose(ns, ns_)

    print "Original: %2.3fs" % (t2 - t1)
    print "New:      %2.3fs" % (t4 - t3)
    print "Speedup:  %2.2fx" % ((t2 - t1) / (t4 - t3))

    #######################
    # mDEV()
    #######################
    print "\ntesting mdev()"
    rate      = 2.0          # 1000 Hz sample rate
    obs_s     = 32768            # 1 hour
    n_samples = rate * obs_s
    t = np.arange(0, n_samples)
    phase  = np.random.random(n_samples) + np.sin(t / n_samples)

    taus = [1, 2, 4]

    t1 = time.time()
    taus2, td, tde, ns = alt.mdev(phase, rate, taus)
    t2 = time.time()

    t3 = time.time()
    taus2_, td_, tde_, ns_ = alt.mdev_np(phase, rate, taus)
    t4 = time.time()

    assert np.allclose(taus2, taus2_)
    assert np.allclose(td, td_)
    assert np.allclose(tde, tde_)
    assert np.allclose(ns, ns_)

    print "Original: %2.3fs" % (t2 - t1)
    print "New:      %2.3fs" % (t4 - t3)
    print "Speedup:  %2.2fx" % ((t2 - t1) / (t4 - t3))

