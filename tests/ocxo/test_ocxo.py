"""
  Test for allantools (https://github.com/aewallin/allantools)
  Stable32 was used to calculate the deviations we compare against.

  AW2015-06-26
  
  The dataset is from the 10 MHz output at the back of an HP Impedance Analyzer
  measured with Keysight 53230A counter, 1.0s gate, RCON mode, with H-maser 10MHz reference

"""

import os
import testutils
import allantools as allan
import math
import sys
import pytest
import numpy as np

sys.path.append("..")
sys.path.append("../..")  # hack to import from parent directory
# remove if you have allantools installed in your python path


data_file = 'ocxo_frequency.txt'


def change_to_test_dir():
    # hack to run script from its own directory
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)


verbose = 1
tolerance = 1e-4  # relative tolerance
rate = 1/float(1.0)  # stable32 runs were done with this data-interval


class TestOCXO():

    def test_adev_ci(self):
        """ ADEV with confidence intervals, fixed alpha (no noise-ID) """
        change_to_test_dir()
        s32rows = testutils.read_stable32(
            resultfile='adev_octave.txt', datarate=1.0)
        for row in s32rows:
            data = testutils.read_datafile(data_file)
            data = allan.frequency2fractional(data, mean_frequency=1.0e7)
            (taus, devs, errs, ns) = allan.adev(data, rate=rate, data_type="freq",
                                                taus=[row['tau']])
            # NOTE! Here we use alhpa from Stable32-results for the allantools edf computation!
            edf = allan.edf_greenhall(alpha=row['alpha'], d=2, m=row['m'], N=len(
                data), overlapping=False, modified=False, verbose=True)
            print("alpha=", row['alpha'])
            (lo, hi) = allan.confidence_interval(devs[0], edf=edf)
            print("      n check: %d" % testutils.check_equal(ns[0], row['n']))
            print("    dev check: %d" % testutils.check_approx_equal(
                devs[0], row['dev'], tolerance=2e-4))
            print("min dev check: %.4g %.4g %d" % (
                lo, row['dev_min'], testutils.check_approx_equal(lo, row['dev_min'], tolerance=1e-3)))
            print("max dev check: %.4g %.4g %d" % (
                hi, row['dev_max'], testutils.check_approx_equal(hi, row['dev_max'], tolerance=1e-3)))

    def test_adev_ci_and_noiseID(self):
        """ ADEV with confidence intervals, including noise-ID """
        change_to_test_dir()
        s32rows = testutils.read_stable32(
            resultfile='adev_octave.txt', datarate=1.0)
        for row in s32rows:
            data = testutils.read_datafile(data_file)
            data = allan.frequency2fractional(data, mean_frequency=1.0e7)
            (taus, devs, errs, ns) = allan.adev(data, rate=rate, data_type="freq",
                                                taus=[row['tau']])
            dev = devs[0]
            try:
                # CI including noise-ID
                (lo2, hi2) = allan.confidence_interval_noiseID(
                    data, dev, af=int(row['m']), dev_type="adev", data_type="freq")
                assert np.isclose(lo2, row['dev_min'], rtol=1e-2)
                assert np.isclose(hi2, row['dev_max'], rtol=1e-2)
                print(" CI OK! tau= %.0f lo/s32_lo = %f hi/s32_hi = %f " %
                      (row['tau'], lo2/row['dev_min'], hi2/row['dev_max']))
            except NotImplementedError:
                print("can't do CI for tau= %f" % row['tau'])
                pass

    def test_oadev_ci(self):
        """ Overlapping ADEV with confidence intervals """
        s32rows = testutils.read_stable32(
            resultfile='oadev_octave.txt', datarate=1.0)
        for row in s32rows:
            data = testutils.read_datafile(data_file)
            data = allan.frequency2fractional(data, mean_frequency=1.0e7)
            (taus, devs, errs, ns) = allan.oadev(data, rate=rate, data_type="freq",
                                                 taus=[row['tau']])
            # NOTE! Here we use alhpa from Stable32-results for the allantools edf computation!
            edf = allan.edf_greenhall(alpha=row['alpha'], d=2, m=row['m'], N=len(
                data), overlapping=True, modified=False, verbose=True)
            (lo, hi) = allan.confidence_interval(devs[0], edf=edf)
            print("n check: ", testutils.check_equal(ns[0], row['n']))
            print("dev check: ", devs[0], row['dev'], testutils.check_approx_equal(
                devs[0], row['dev'], tolerance=2e-3))
            print("min dev check: ",  lo, row['dev_min'], testutils.check_approx_equal(
                lo, row['dev_min'], tolerance=2e-3))
            print("max dev check: ", hi, row['dev_max'], testutils.check_approx_equal(
                hi, row['dev_max'], tolerance=2e-3))

    def test_mdev_ci(self):
        """ Overlapping ADEV with confidence intervals """
        s32rows = testutils.read_stable32(
            resultfile='mdev_octave.txt', datarate=1.0)
        for row in s32rows:
            data = testutils.read_datafile(data_file)
            data = allan.frequency2fractional(data, mean_frequency=1.0e7)
            (taus, devs, errs, ns) = allan.mdev(data, rate=rate, data_type="freq",
                                                taus=[row['tau']])
            # NOTE! Here we use alhpa from Stable32-results for the allantools edf computation!
            edf = allan.edf_greenhall(alpha=row['alpha'], d=2, m=row['m'], N=len(
                data), overlapping=True, modified=True, verbose=True)
            (lo, hi) = allan.confidence_interval(devs[0], edf=edf)
            print("n check: ", testutils.check_equal(ns[0], row['n']))
            print("dev check: ", devs[0], row['dev'], testutils.check_approx_equal(
                devs[0], row['dev'], tolerance=2e-3))
            print("min dev check: ",  lo, row['dev_min'], testutils.check_approx_equal(
                lo, row['dev_min'], tolerance=2e-3))
            print("max dev check: ", hi, row['dev_max'], testutils.check_approx_equal(
                hi, row['dev_max'], tolerance=2e-3))

    def test_tdev_ci(self):
        """ Time Deviation with confidence intervals """
        s32rows = testutils.read_stable32(
            resultfile='tdev_octave.txt', datarate=1.0)
        for row in s32rows:
            data = testutils.read_datafile(data_file)
            data = allan.frequency2fractional(data, mean_frequency=1.0e7)
            (taus, devs, errs, ns) = allan.tdev(data, rate=rate, data_type="freq",
                                                taus=[row['tau']])
            # NOTE! Here we use alhpa from Stable32-results for the allantools edf computation!
            edf = allan.edf_greenhall(alpha=row['alpha'], d=2, m=row['m'], N=len(
                data), overlapping=True, modified=True, verbose=True)
            (lo, hi) = allan.confidence_interval(devs[0], edf=edf)
            print("n check: ", testutils.check_equal(ns[0], row['n']))
            print("dev check: ", devs[0], row['dev'], testutils.check_approx_equal(
                devs[0], row['dev'], tolerance=2e-3))
            print("min dev check: ",  lo, row['dev_min'], testutils.check_approx_equal(
                lo, row['dev_min'], tolerance=2e-3))
            print("max dev check: ", hi, row['dev_max'], testutils.check_approx_equal(
                hi, row['dev_max'], tolerance=2e-3))

    def test_hdev_ci(self):
        """ Hadamard with confidence intervals """
        s32rows = testutils.read_stable32(
            resultfile='hdev_octave.txt', datarate=1.0)
        for row in s32rows:
            data = testutils.read_datafile(data_file)
            data = allan.frequency2fractional(data, mean_frequency=1.0e7)
            (taus, devs, errs, ns) = allan.hdev(data, rate=rate, data_type="freq",
                                                taus=[row['tau']])
            # NOTE! Here we use alhpa from Stable32-results for the allantools edf computation!
            edf = allan.edf_greenhall(alpha=row['alpha'], d=3, m=row['m'], N=len(
                data), overlapping=False, modified=False, verbose=True)
            (lo, hi) = allan.confidence_interval(devs[0], edf=edf)
            print("n check: ", testutils.check_equal(ns[0], row['n']))
            print("dev check: ", devs[0], row['dev'], testutils.check_approx_equal(
                devs[0], row['dev'], tolerance=2e-3))
            print("min dev check: ",  lo, row['dev_min'], testutils.check_approx_equal(
                lo, row['dev_min'], tolerance=2e-3))
            print("max dev check: ", hi, row['dev_max'], testutils.check_approx_equal(
                hi, row['dev_max'], tolerance=2e-3))

    def test_ohdev_ci(self):
        """ Overlapping Hadamard deviation with confidence intervals  """
        s32rows = testutils.read_stable32(
            resultfile='ohdev_octave.txt', datarate=1.0)
        for row in s32rows:
            data = testutils.read_datafile(data_file)
            data = allan.frequency2fractional(data, mean_frequency=1.0e7)
            (taus, devs, errs, ns) = allan.ohdev(data, rate=rate, data_type="freq",
                                                 taus=[row['tau']])
            # NOTE! Here we use alhpa from Stable32-results for the allantools edf computation!
            edf = allan.edf_greenhall(alpha=row['alpha'], d=3, m=row['m'], N=len(
                data), overlapping=True, modified=False, verbose=True)
            (lo, hi) = allan.confidence_interval(devs[0], edf=edf)
            print("n check: ", testutils.check_equal(ns[0], row['n']))
            print("dev check: ", devs[0], row['dev'], testutils.check_approx_equal(
                devs[0], row['dev'], tolerance=2e-3))
            print("min dev check: ",  lo, row['dev_min'], testutils.check_approx_equal(
                lo, row['dev_min'], tolerance=2e-3))
            print("max dev check: ", hi, row['dev_max'], testutils.check_approx_equal(
                hi, row['dev_max'], tolerance=5e-3))

    # fails
    # totdev() needs bias-correction, depending on alpha(?)
    @pytest.mark.skip(reason="needs bias-correction and noise-ID to work")
    @pytest.mark.xfail
    def test_totdev_ci(self):
        print("totdev()")
        s32rows = testutils.read_stable32(
            resultfile='totdev_octave.txt', datarate=1.0)
        for row in s32rows:
            data = testutils.read_datafile(data_file)
            data = allan.frequency2fractional(data, mean_frequency=1.0e7)
            (taus, devs, errs, ns) = allan.totdev(data, rate=rate, data_type="freq",
                                                  taus=[row['tau']])
            edf = allan.edf_totdev(N=len(data), m=row['m'], alpha=row['alpha'])

            (lo, hi) = allan.confidence_interval(devs[0], edf=edf)
            print("n check: ", testutils.check_equal(ns[0], row['n']))
            print("dev check: ", testutils.check_approx_equal(
                devs[0], row['dev'], tolerance=2e-3))
            print("min dev check: %.4g %.4g %d" % (
                lo, row['dev_min'], testutils.check_approx_equal(lo, row['dev_min'], tolerance=2e-3)))
            print("max dev check: %.4g %.4g %d" % (
                hi, row['dev_max'], testutils.check_approx_equal(hi, row['dev_max'], tolerance=2e-3)))

    def test_noise_id(self):
        """ test for noise-identification """
        s32_rows = testutils.read_stable32('mdev_octave.txt', rate)
        freq = testutils.read_datafile(data_file)
        y_freq = allan.frequency2fractional(freq, mean_frequency=1.0e7)
        phase = allan.frequency2phase(freq, rate)

        for s32 in s32_rows:
            s32_tau, s32_alpha, s32_AF = s32['tau'], s32['alpha'], int(
                s32['m'])

            # noise-ID from frequency
            if len(phase)/s32_AF > 20:
                alpha_int, alpha, d, rho = allan.autocorr_noise_id(
                    freq, data_type='freq', af=s32_AF)
                print("y: ", s32_tau, s32_alpha, alpha_int, alpha, rho, d)
                assert alpha_int == s32_alpha

            # noise-ID from phase
            if len(phase)/s32_AF > 20:
                alpha_int, alpha, d, rho = allan.autocorr_noise_id(
                    phase, data_type='phase', af=s32_AF)
                print("x: ", s32_tau, s32_alpha, alpha_int, alpha, rho, d)
                assert alpha_int == s32_alpha


if __name__ == "__main__":
    # pytest.main()
    t = TestOCXO()
    # t.test_adev_ci_and_noiseID()
    t.test_noise_id()
"""
    t.test_adev_ci()
    t.test_oadev_ci()
    t.test_mdev_ci()
    t.test_tdev_ci()
    t.test_hdev_ci()
    t.test_ohdev_ci()
"""
# t.test_totdev_ci()
