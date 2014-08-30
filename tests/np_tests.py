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
import allantools.allantools_pure_python as alt
import allantools.allantools as alp
import time


if __name__ == "__main__":


    #######################
    # MTIE()
    #######################
    print "\ntesting mtie()"
    data = np.random.random(1000)
    taus = [1, 3, 5, 16, 128]
    rates = [1, 20, 10.7]
    strides = [1, 10, 7]

    for rate in rates:
        for stride in strides:
            #print "TAU: %i, RATE: %2.2f, STRIDE: %i" % (tau, rate, stride)
            o_taus, o_dev, o_err, o_n = alt.mtie(data, rate, taus)
            o_taus_, o_dev_, o_err_, o_n_ = alp.mtie(data, rate, taus)

            assert np.allclose(o_taus, o_taus_)
            assert np.allclose(o_dev, o_dev_)
            assert np.allclose(o_err, o_err_)

    stride = 1
    tau = 128000
    rate = 2.1
    data = np.random.random(10000)
    t1 = time.time()
    o_taus, o_dev, o_err, o_n = alt.mtie(data, rate, taus)
    t2 = time.time()
    t3 = time.time()
    o_taus_, o_dev_, o_err_, o_n_ = alp.mtie(data, rate, taus)
    t4 = time.time()

    assert np.allclose(o_taus, o_taus_)
    assert np.allclose(o_dev, o_dev_)
    assert np.allclose(o_err, o_err_)
    print "Original: %2.3fs" % (t2 - t1)
    print "New:      %2.3fs" % (t4 - t3)
    print "Speedup:  %2.2fx" % ((t2 - t1) / (t4 - t3))

    #######################
    # MTIE_PHASE()
    #######################
    print "\ntesting mtie_phase()"
    data = np.random.random(1000)
    taus = [1, 3, 5, 16, 128]
    rates = [1, 20, 10.7]
    strides = [1, 10, 7]

    for rate in rates:
        for stride in strides:
            #print "TAU: %i, RATE: %2.2f, STRIDE: %i" % (tau, rate, stride)
            o_taus, o_dev, o_err, o_n = alt.mtie_phase(data, rate, taus)
            o_taus_, o_dev_, o_err_, o_n_ = alp.mtie_phase(data, rate, taus)

            assert np.allclose(o_taus, o_taus_)
            assert np.allclose(o_dev, o_dev_)
            assert np.allclose(o_err, o_err_)

    stride = 1
    tau = 128000
    rate = 2.1
    data = np.random.random(100000)
    t1 = time.time()
    o_taus, o_dev, o_err, o_n = alt.mtie_phase(data, rate, taus)
    t2 = time.time()
    t3 = time.time()
    o_taus_, o_dev_, o_err_, o_n_ = alp.mtie_phase(data, rate, taus)
    t4 = time.time()

    assert np.allclose(o_taus, o_taus_)
    assert np.allclose(o_dev, o_dev_)
    assert np.allclose(o_err, o_err_)
    print "Original: %2.3fs" % (t2 - t1)
    print "New:      %2.3fs" % (t4 - t3)
    print "Speedup:  %2.2fx" % ((t2 - t1) / (t4 - t3))


    #######################
    # THREE_CORNERED_HAT_PHASE()
    #######################
    print "\ntesting three_cornered_hat_phase()"

    stride = 1
    taus = [2, 4, 8, 16]
    rate = 2.1
    pdata_ab = np.random.random(100000)
    pdata_bc = np.random.random(100000)
    pdata_ca = np.random.random(100000)

    t1 = time.time()
    function = alt.adev
    tau, dev_a = alt.three_cornered_hat_phase(pdata_ab, pdata_bc, pdata_ca, rate, taus, function)
    t2 = time.time()
    t3 = time.time()
    function = alp.adev
    tau_, dev_a_ = alp.three_cornered_hat_phase(pdata_ab, pdata_bc, pdata_ca, rate, taus, function)
    t4 = time.time()

    assert np.allclose(tau, tau_)
    assert np.allclose(dev_a, dev_a_)
    print "Original: %2.3fs" % (t2 - t1)
    print "New:      %2.3fs" % (t4 - t3)
    print "Speedup:  %2.2fx" % ((t2 - t1) / (t4 - t3))

    #######################
    # TIERMS_PHASE()
    #######################
    print "\ntesting tierms_phase()"
    data = np.random.random(1000)
    taus = [1, 3, 5, 16, 128]
    rates = [1, 20, 10.7]
    strides = [1, 10, 7]

    for rate in rates:
        for stride in strides:
            #print "TAU: %i, RATE: %2.2f, STRIDE: %i" % (tau, rate, stride)
            o_taus, o_dev, o_err, o_n = alt.tierms_phase(data, rate, taus)
            o_taus_, o_dev_, o_err_, o_n_ = alp.tierms_phase(data, rate, taus)

            assert np.allclose(o_taus, o_taus_)
            assert np.allclose(o_dev, o_dev_)
            assert np.allclose(o_err, o_err_)


    stride = 1
    tau = 16
    rate = 2.1
    data = np.random.random(100000)
    t1 = time.time()
    o_taus, o_dev, o_err, o_n = alt.tierms_phase(data, rate, taus)
    t2 = time.time()
    t3 = time.time()
    o_taus_, o_dev_, o_err_, o_n_ = alp.tierms_phase(data, rate, taus)
    t4 = time.time()

    assert np.allclose(o_taus, o_taus_)
    assert np.allclose(o_dev, o_dev_)
    assert np.allclose(o_err, o_err_)
    print "Original: %2.3fs" % (t2 - t1)
    print "New:      %2.3fs" % (t4 - t3)
    print "Speedup:  %2.2fx" % ((t2 - t1) / (t4 - t3))

    #######################
    # TIERMS()
    #######################
    print "\ntesting tierms()"
    data = np.random.random(1000)
    taus = [1, 3, 5, 16, 128]
    rates = [1, 20, 10.7]
    strides = [1, 10, 7]

    for rate in rates:
        for stride in strides:
            #print "TAU: %i, RATE: %2.2f, STRIDE: %i" % (tau, rate, stride)
            o_taus, o_dev, o_err, o_n = alt.tierms(data, rate, taus)
            o_taus_, o_dev_, o_err_, o_n_ = alp.tierms(data, rate, taus)

            assert np.allclose(o_taus, o_taus_)
            assert np.allclose(o_dev, o_dev_)
            assert np.allclose(o_err, o_err_)


    stride = 1
    tau = 16
    rate = 2.1
    data = np.random.random(100000)
    t1 = time.time()
    o_taus, o_dev, o_err, o_n = alt.tierms(data, rate, taus)
    t2 = time.time()
    t3 = time.time()
    o_taus_, o_dev_, o_err_, o_n_ = alp.tierms(data, rate, taus)
    t4 = time.time()

    assert np.allclose(o_taus, o_taus_)
    assert np.allclose(o_dev, o_dev_)
    assert np.allclose(o_err, o_err_)
    print "Original: %2.3fs" % (t2 - t1)
    print "New:      %2.3fs" % (t4 - t3)
    print "Speedup:  %2.2fx" % ((t2 - t1) / (t4 - t3))

    #######################
    # TOTDEV_PHASE()
    #######################
    print "\ntesting totdev_phase()"
    data = np.random.random(1000)
    taus = [1, 3, 5, 16, 128]
    rates = [1, 20, 10.7]
    strides = [1, 10, 7]

    for rate in rates:
        for stride in strides:
            #print "TAU: %i, RATE: %2.2f, STRIDE: %i" % (tau, rate, stride)
            o_taus, o_dev, o_err, o_n = alt.totdev_phase(data, rate, taus)
            o_taus_, o_dev_, o_err_, o_n_ = alp.totdev_phase(data, rate, taus)

            assert np.allclose(o_taus, o_taus_)
            assert np.allclose(o_dev, o_dev_)
            assert np.allclose(o_err, o_err_)

    stride = 1
    tau = 16
    rate = 2.1
    data = np.random.random(100000)
    t1 = time.time()
    o_taus, o_dev, o_err, o_n = alt.totdev_phase(data, rate, taus)
    t2 = time.time()
    t3 = time.time()
    o_taus_, o_dev_, o_err_, o_n_ = alp.totdev_phase(data, rate, taus)
    t4 = time.time()

    assert np.allclose(o_taus, o_taus_)
    assert np.allclose(o_dev, o_dev_)
    assert np.allclose(o_err, o_err_)
    print "Original: %2.3fs" % (t2 - t1)
    print "New:      %2.3fs" % (t4 - t3)
    print "Speedup:  %2.2fx" % ((t2 - t1) / (t4 - t3))

    #######################
    # TOTDEV()
    #######################
    print "\ntesting totdev()"
    data = np.random.random(1000)
    taus = [1, 3, 5, 16, 128]
    rates = [1, 20, 10.7]
    strides = [1, 10, 7]

    for rate in rates:
        for stride in strides:
            #print "TAU: %i, RATE: %2.2f, STRIDE: %i" % (tau, rate, stride)
            o_taus, o_dev, o_err, o_n = alt.totdev(data, rate, taus)
            o_taus_, o_dev_, o_err_, o_n_ = alp.totdev(data, rate, taus)

            assert np.allclose(o_taus, o_taus_)
            assert np.allclose(o_dev, o_dev_)
            assert np.allclose(o_err, o_err_)

    stride = 1
    tau = 16
    rate = 2.1
    data = np.random.random(100000)
    t1 = time.time()
    o_taus, o_dev, o_err, o_n = alt.totdev(data, rate, taus)
    t2 = time.time()
    t3 = time.time()
    o_taus_, o_dev_, o_err_, o_n_ = alp.totdev(data, rate, taus)
    t4 = time.time()

    assert np.allclose(o_taus, o_taus_)
    assert np.allclose(o_dev, o_dev_)
    assert np.allclose(o_err, o_err_)
    print "Original: %2.3fs" % (t2 - t1)
    print "New:      %2.3fs" % (t4 - t3)
    print "Speedup:  %2.2fx" % ((t2 - t1) / (t4 - t3))

    #######################
    # OHDEV()
    #######################
    print "\ntesting ohdev()"
    data = np.random.random(10000)
    taus = [1, 3, 5, 16, 128]
    rates = [1, 20, 10.7]
    strides = [1, 10, 7]

    for rate in rates:
        for stride in strides:
            #print "TAU: %i, RATE: %2.2f, STRIDE: %i" % (tau, rate, stride)
            o_taus, o_dev, o_err, o_n = alt.ohdev(data, rate, taus)
            o_taus_, o_dev_, o_err_, o_n_ = alp.ohdev(data, rate, taus)

            assert np.allclose(o_taus, o_taus_)
            assert np.allclose(o_dev, o_dev_)
            assert np.allclose(o_err, o_err_)

    stride = 1
    tau = 16
    rate = 2.1
    data = np.random.random(100000)
    t1 = time.time()
    o_taus, o_dev, o_err, o_n = alt.ohdev(data, rate, taus)
    t2 = time.time()
    t3 = time.time()
    o_taus_, o_dev_, o_err_, o_n_ = alp.ohdev(data, rate, taus)
    t4 = time.time()

    assert np.allclose(o_taus, o_taus_)
    assert np.allclose(o_dev, o_dev_)
    assert np.allclose(o_err, o_err_)
    print "Original: %2.3fs" % (t2 - t1)
    print "New:      %2.3fs" % (t4 - t3)
    print "Speedup:  %2.2fx" % ((t2 - t1) / (t4 - t3))

    #######################
    # OHDEV_PHASE()
    #######################
    print "\ntesting ohdev_phase()"
    data = np.random.random(10000)
    taus = [1, 3, 5, 16, 128]
    rates = [1, 20, 10.7]
    strides = [1, 10, 7]

    for rate in rates:
        for stride in strides:
            #print "TAU: %i, RATE: %2.2f, STRIDE: %i" % (tau, rate, stride)
            o_taus, o_dev, o_err, o_n = alt.ohdev_phase(data, rate, taus)
            o_taus_, o_dev_, o_err_, o_n_ = alp.ohdev_phase(data, rate, taus)

            assert np.allclose(o_taus, o_taus_)
            assert np.allclose(o_dev, o_dev_)
            assert np.allclose(o_err, o_err_)

    stride = 1
    tau = 16
    rate = 2.1
    data = np.random.random(100000)
    t1 = time.time()
    o_taus, o_dev, o_err, o_n = alt.ohdev_phase(data, rate, taus)
    t2 = time.time()
    t3 = time.time()
    o_taus_, o_dev_, o_err_, o_n_ = alp.ohdev_phase(data, rate, taus)
    t4 = time.time()

    assert np.allclose(o_taus, o_taus_)
    assert np.allclose(o_dev, o_dev_)
    assert np.allclose(o_err, o_err_)
    print "Original: %2.3fs" % (t2 - t1)
    print "New:      %2.3fs" % (t4 - t3)
    print "Speedup:  %2.2fx" % ((t2 - t1) / (t4 - t3))

    #######################
    # HDEV_PHASE_CALC()
    #######################
    print "\ntesting hdev_phase_calc()"
    data = np.random.random(10000)
    taus = [1, 3, 5, 16, 128]
    rates = [1, 20, 10.7]
    strides = [1, 10, 7]

    for tau in taus:
        for rate in rates:
            for stride in strides:
                #print "TAU: %i, RATE: %2.2f, STRIDE: %i" % (tau, rate, stride)
                mj = tau
                dev, deverr, n = alt.calc_hdev_phase(data, rate, mj, stride)
                dev_, deverr_, n_ = alp.calc_hdev_phase(data, rate, mj, stride)

                assert np.isclose(dev, dev_)
                assert np.isclose(n, n_)
                assert np.isclose(deverr, deverr_)

    stride = 1
    tau = 16
    rate = 2.0
    t1 = time.time()
    dev, deverr, n = alt.calc_hdev_phase(data, rate, mj, stride)
    t2 = time.time()
    t3 = time.time()
    dev_, deverr_, n_ = alp.calc_hdev_phase(data, rate, mj, stride)
    t4 = time.time()

    assert np.isclose(dev, dev_)
    assert np.isclose(n, n_)
    assert np.isclose(deverr, deverr_)
    print "Original: %2.3fs" % (t2 - t1)
    print "New:      %2.3fs" % (t4 - t3)
    print "Speedup:  %2.2fx" % ((t2 - t1) / (t4 - t3))


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
            o_taus_, o_dev_, o_err_, o_n_ = alp.oadev_phase(data, rate, taus)

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
    o_taus_, o_dev_, o_err_, o_n_ = alp.oadev_phase(data, rate, taus)
    t4 = time.time()

    assert np.allclose(o_taus, o_taus_)
    assert np.allclose(o_dev, o_dev_)
    assert np.allclose(o_err, o_err_)
    print "Original: %2.3fs" % (t2 - t1)
    print "New:      %2.3fs" % (t4 - t3)
    print "Speedup:  %2.2fx" % ((t2 - t1) / (t4 - t3))

    #######################
    # HDEV_PHASE()
    #######################
    print "\ntesting hdev_phase()"
    data = np.random.random(10000)
    taus = [1, 3, 5, 16, 128]
    rates = [1, 20, 10.7]
    strides = [1, 10, 7]

    for rate in rates:
        for stride in strides:
            #print "TAU: %i, RATE: %2.2f, STRIDE: %i" % (tau, rate, stride)
            o_taus, o_dev, o_err, o_n = alt.hdev_phase(data, rate, taus)
            o_taus_, o_dev_, o_err_, o_n_ = alp.hdev_phase(data, rate, taus)

            assert np.allclose(o_taus, o_taus_)
            assert np.allclose(o_dev, o_dev_)
            assert np.allclose(o_err, o_err_)

    stride = 1
    tau = 16
    rate = 2.1
    data = np.random.random(100000)
    t1 = time.time()
    o_taus, o_dev, o_err, o_n = alt.hdev_phase(data, rate, taus)
    t2 = time.time()
    t3 = time.time()
    o_taus_, o_dev_, o_err_, o_n_ = alp.hdev_phase(data, rate, taus)
    t4 = time.time()

    assert np.allclose(o_taus, o_taus_)
    assert np.allclose(o_dev, o_dev_)
    assert np.allclose(o_err, o_err_)
    print "Original: %2.3fs" % (t2 - t1)
    print "New:      %2.3fs" % (t4 - t3)
    print "Speedup:  %2.2fx" % ((t2 - t1) / (t4 - t3))

    #######################
    # HDEV()
    #######################
    print "\ntesting hdev()"
    data = np.random.random(10000)
    taus = [1, 3, 5, 16, 128]
    rates = [1, 20, 10.7]
    strides = [1, 10, 7]

    for rate in rates:
        for stride in strides:
            #print "TAU: %i, RATE: %2.2f, STRIDE: %i" % (tau, rate, stride)
            o_taus, o_dev, o_err, o_n = alt.hdev(data, rate, taus)
            o_taus_, o_dev_, o_err_, o_n_ = alp.hdev(data, rate, taus)

            assert np.allclose(o_taus, o_taus_)
            assert np.allclose(o_dev, o_dev_)
            assert np.allclose(o_err, o_err_)

    stride = 1
    tau = 16
    rate = 2.1
    data = np.random.random(100000)
    t1 = time.time()
    o_taus, o_dev, o_err, o_n = alt.hdev(data, rate, taus)
    t2 = time.time()
    t3 = time.time()
    o_taus_, o_dev_, o_err_, o_n_ = alp.hdev(data, rate, taus)
    t4 = time.time()

    assert np.allclose(o_taus, o_taus_)
    assert np.allclose(o_dev, o_dev_)
    assert np.allclose(o_err, o_err_)
    print "Original: %2.3fs" % (t2 - t1)
    print "New:      %2.3fs" % (t4 - t3)
    print "Speedup:  %2.2fx" % ((t2 - t1) / (t4 - t3))


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
            o_taus_, o_dev_, o_err_, o_n_ = alp.oadev_phase(data, rate, taus)

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
    o_taus_, o_dev_, o_err_, o_n_ = alp.oadev_phase(data, rate, taus)
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
            o_taus_, o_dev_, o_err_, o_n_ = alp.oadev(data, rate, taus)

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
    o_taus_, o_dev_, o_err_, o_n_ = alp.oadev(data, rate, taus)
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
            o_taus_, o_dev_, o_err_, o_n_ = alp.adev_phase(data, rate, taus)

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
    o_taus_, o_dev_, o_err_, o_n_ = alp.adev_phase(data, rate, taus)
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
            o_taus_, o_dev_, o_err_, o_n_ = alp.adev(data, rate, taus)

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
    o_taus_, o_dev_, o_err_, o_n_ = alp.adev(data, rate, taus)
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
                dev_, deverr_, n_ = alp.calc_adev_phase(data, rate, mj, stride)

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
    dev_, deverr_, n_ = alp.calc_adev_phase(data, rate, mj, stride)
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
        data_, m_, taus2_ = alp.tau_m(data, rate, taus)
        assert np.allclose(m, m_)
        assert np.allclose(taus2, taus2_)

    taus = np.random.randint(low=-100, high=10000, size=(10000,))
    rate = 1.234
    t1 = time.time()
    m, taus2 = alt.tau_m(data, rate, taus)
    t2 = time.time()
    t3 = time.time()
    data_, m_, taus2_ = alp.tau_m(data, rate, taus)
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
        phase_ = alp.frequency2phase(freqdata, rate)
        assert len(phase) == len(phase_)
        assert np.allclose(phase, phase_)

    freqdata = np.random.random(100000)
    t1 = time.time()
    phase  = alt.frequency2phase(freqdata, rate)
    t2 = time.time()
    t3 = time.time()
    phase_ = alp.frequency2phase(freqdata, rate)
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
    taus2_, td_, tde_, ns_ = alp.tdev_phase(phase, rate, taus)
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
    taus2_, td_, tde_, ns_ = alp.tdev(phase, rate, taus)
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
    taus2_, td_, tde_, ns_ = alp.mdev_phase(phase, rate, taus)
    t4 = time.time()

    assert np.allclose(taus2, taus2_)
    assert np.allclose(td, td_)
    assert np.allclose(tde, tde_)
    assert np.allclose(ns, ns_)

    print "Original: %2.3fs" % (t2 - t1)
    print "New:      %2.3fs" % (t4 - t3)
    print "Speedup:  %2.2fx" % ((t2 - t1) / (t4 - t3))

    #######################
    # MDEV()
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
    taus2_, td_, tde_, ns_ = alp.mdev(phase, rate, taus)
    t4 = time.time()

    assert np.allclose(taus2, taus2_)
    assert np.allclose(td, td_)
    assert np.allclose(tde, tde_)
    assert np.allclose(ns, ns_)

    print "Original: %2.3fs" % (t2 - t1)
    print "New:      %2.3fs" % (t4 - t3)
    print "Speedup:  %2.2fx" % ((t2 - t1) / (t4 - t3))

