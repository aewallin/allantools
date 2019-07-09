"""
  Test for allantools (https://github.com/aewallin/allantools)
  Stable32 was used to calculate the deviations we compare against.

  GPS tests, AW2016-03-17
"""

import math
import sys
sys.path.append("..")
sys.path.append("../..") # hack to import from parent directory
# remove if you have allantools installed in your python path
import os
import time
import pytest
import numpy as np

import allantools as allan
import testutils

def print_elapsed(start):
    end = time.clock()
    print(" %.2f s"%end-start)
    return time.clock()

def change_to_test_dir():
    # hack to run script from its own directory
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)

data_file = 'gps_1pps_phase_data.txt.gz'
verbose = 1
tolerance = 1e-4 # relative tolerance
rate = 1/float(1.0) # 1 PPS measurements, data interval is 1 s

class TestGPS():
    def test_adev(self):
        self.generic_test(result='stable32_ADEV_decade.txt', fct=allan.adev)
    def test_oadev(self):
        self.generic_test(result='stable32_OADEV_octave.txt', fct=allan.oadev)
    def test_mdev(self):
        self.generic_test(result='stable32_MDEV_octave.txt', fct=allan.mdev)
    def test_tdev(self):
        self.generic_test(result='stable32_TDEV_octave.txt', fct=allan.tdev)
    def test_hdev(self):
        self.generic_test(result='stable32_HDEV_octave.txt', fct=allan.hdev)
    def test_ohdev(self):
        self.generic_test(result='stable32_OHDEV_octave.txt', fct=allan.ohdev)
    def test_totdev(self):
        self.generic_test(result='stable32_TOTDEV_octave.txt', fct=allan.totdev)

    def test_noise_id(self):
        """ test for noise-identification """
        s32_rows = testutils.read_stable32('stable32_ADEV_decade.txt', rate)
        phase = testutils.read_datafile('gps_1pps_phase_data.txt.gz')
        # test noise-ID
        for s32 in s32_rows:
            tau, alpha, AF = s32['tau'], s32['alpha'], int(s32['m'])
            try:
                alpha_int = allan.autocorr_noise_id(phase, af=AF)[0]
                print(tau, alpha, alpha_int)
                assert alpha_int == alpha
            except NotImplementedError:
                print("can't do noise-ID for tau= %f"%s32['tau'])

    def test_adev_ci_and_noiseID(self):
        """ ADEV with confidence intervals, including noise-ID """
        change_to_test_dir()
        phase = testutils.read_datafile('gps_1pps_phase_data.txt.gz')
        
        s32rows = testutils.read_stable32(resultfile='stable32_ADEV_decade.txt', datarate=1.0)
        s32taus = [row['tau'] for row in s32rows]
        (taus, devs, errs, ns) = allan.adev(phase, rate=rate, data_type="phase",
                                                 taus=s32taus)
        for idx, row in enumerate(s32rows):
            
            #(taus, devs, errs, ns) = allan.adev(phase, rate=rate, data_type="phase",
            #                                    taus=[row['tau']])
            dev = devs[idx]
            try:
                # CI including noise-ID
                (lo2, hi2) = allan.confidence_interval_noiseID(phase, dev, af=int(row['m']), dev_type="adev", data_type="phase")
                assert np.isclose(lo2, row['dev_min'], rtol=1e-2, atol=0)
                assert np.isclose(hi2, row['dev_max'], rtol=1e-2, atol=0)
                print(" CI OK! tau= %f  lo/s32_lo = %f hi/s32_hi = %f "% (row['tau'], lo2/row['dev_min'], hi2/row['dev_max']))
            except NotImplementedError:
                print("can't do CI for tau= %f"%row['tau'])
                pass

    def test_oadev_ci_and_noiseID(self):
        """ ADEV with confidence intervals, including noise-ID """
        change_to_test_dir()
        phase = testutils.read_datafile('gps_1pps_phase_data.txt.gz')
        s32rows = testutils.read_stable32(resultfile='stable32_OADEV_octave.txt', datarate=1.0)
        s32taus = [row['tau'] for row in s32rows]
        (taus, devs, errs, ns) = allan.oadev(phase, rate=rate, data_type="phase",
                                                 taus=s32taus)
        for idx, row in enumerate(s32rows):
            dev = devs[idx]
            assert np.isclose(dev, row['dev'], rtol=1e-2, atol=0) # check deviation            
            # now check confidence interval
            try:
                # CI including noise-ID
                (lo2, hi2) = allan.confidence_interval_noiseID(phase, dev, af=int(row['m']), dev_type="oadev", data_type="phase")
                print(" tau= %f  lo/s32_lo = %f hi/s32_hi = %f "% (row['tau'], lo2/row['dev_min'], hi2/row['dev_max']))
                assert np.isclose(lo2, row['dev_min'], rtol=1e-2, atol=0)
                assert np.isclose(hi2, row['dev_max'], rtol=1e-2, atol=0)
            except NotImplementedError:
                print("can't do CI for tau= %f"%row['tau'])
                pass

    def test_mdev_ci_and_noiseID(self):
        """ ADEV with confidence intervals, including noise-ID """
        change_to_test_dir()
        phase = testutils.read_datafile('gps_1pps_phase_data.txt.gz')
        s32rows = testutils.read_stable32(resultfile='stable32_MDEV_octave.txt', datarate=1.0)
        s32taus = [row['tau'] for row in s32rows]
        (taus, devs, errs, ns) = allan.mdev(phase, rate=rate, data_type="phase",
                                                 taus=s32taus)
        for idx, row in enumerate(s32rows):
            
            #(taus, devs, errs, ns) = allan.mdev(phase, rate=rate, data_type="phase",
            #                                    taus=[row['tau']])
            dev = devs[idx]
            #print(dev/row['dev'])
            assert np.isclose(dev, row['dev'], rtol=1e-2, atol=0)
            try:
                # CI including noise-ID
                (lo2, hi2) = allan.confidence_interval_noiseID(phase, dev, af=int(row['m']), dev_type="mdev", data_type="phase")
                print(" tau= %f  lo/s32_lo = %f hi/s32_hi = %f "% (row['tau'], lo2/row['dev_min'], hi2/row['dev_max']))
                assert np.isclose(lo2, row['dev_min'], rtol=1e-2, atol=0)
                assert np.isclose(hi2, row['dev_max'], rtol=1e-2, atol=0)
            except NotImplementedError:
                print("can't do CI for tau= %f"%row['tau'])
                pass

    def test_tdev_ci_and_noiseID(self):
        """ ADEV with confidence intervals, including noise-ID """
        change_to_test_dir()
        phase = testutils.read_datafile('gps_1pps_phase_data.txt.gz')
        s32rows = testutils.read_stable32(resultfile='stable32_TDEV_octave.txt', datarate=1.0)
        s32taus = [row['tau'] for row in s32rows]
        (taus, devs, errs, ns) = allan.tdev(phase, rate=rate, data_type="phase",
                                                 taus=s32taus)
        for idx, row in enumerate(s32rows):
            
            #(taus, devs, errs, ns) = allan.tdev(phase, rate=rate, data_type="phase",
            #                                    taus=[row['tau']])
            dev = devs[idx]
            #print(dev/row['dev'])
            assert np.isclose(dev, row['dev'], rtol=1e-2, atol=0)
            try:
                # CI including noise-ID
                (lo2, hi2) = allan.confidence_interval_noiseID(phase, dev, af=int(row['m']), dev_type="tdev", data_type="phase")
                print(" tau= %f  lo/s32_lo = %f hi/s32_hi = %f "% (row['tau'], lo2/row['dev_min'], hi2/row['dev_max']))
                assert np.isclose(lo2, row['dev_min'], rtol=1e-2, atol=0)
                assert np.isclose(hi2, row['dev_max'], rtol=1e-2, atol=0)
            except NotImplementedError:
                print("can't do CI for tau= %f"%row['tau'])
                pass

    def test_hdev_ci_and_noiseID(self):
        """ ADEV with confidence intervals, including noise-ID """
        change_to_test_dir()
        phase = testutils.read_datafile('gps_1pps_phase_data.txt.gz')
        s32rows = testutils.read_stable32(resultfile='stable32_HDEV_octave.txt', datarate=1.0)
        s32taus = [row['tau'] for row in s32rows]
        (taus, devs, errs, ns) = allan.hdev(phase, rate=rate, data_type="phase",
                                                 taus=s32taus)
        for idx, row in enumerate(s32rows):
            
            #(taus, devs, errs, ns) = allan.hdev(phase, rate=rate, data_type="phase",
            #                                    taus=[row['tau']])
            dev = devs[idx]
            #print(dev/row['dev'])
            assert np.isclose(dev, row['dev'], rtol=1e-2, atol=0)
            try:
                # CI including noise-ID
                (lo2, hi2) = allan.confidence_interval_noiseID(phase, dev, af=int(row['m']), dev_type="hdev", data_type="phase")
                print(" tau= %f  lo/s32_lo = %f hi/s32_hi = %f "% (row['tau'], lo2/row['dev_min'], hi2/row['dev_max']))
                assert np.isclose(lo2, row['dev_min'], rtol=1e-2, atol=0)
                assert np.isclose(hi2, row['dev_max'], rtol=1e-2, atol=0)
            except NotImplementedError:
                print("can't do CI for tau= %f"%row['tau'])
                pass

    def test_ohdev_ci_and_noiseID(self):
        """ ADEV with confidence intervals, including noise-ID """
        change_to_test_dir()
        phase = testutils.read_datafile('gps_1pps_phase_data.txt.gz')
        s32rows = testutils.read_stable32(resultfile='stable32_OHDEV_octave.txt', datarate=1.0)
        s32taus = [row['tau'] for row in s32rows]
        (taus, devs, errs, ns) = allan.ohdev(phase, rate=rate, data_type="phase",
                                                 taus=s32taus)
        for idx, row in enumerate(s32rows):
            
            #(taus, devs, errs, ns) = allan.ohdev(phase, rate=rate, data_type="phase",
            #                                     taus=[row['tau']])
            dev = devs[idx]
            #print(dev/row['dev'])
            assert np.isclose(dev, row['dev'], rtol=1e-2, atol=0)
            try:
                # CI including noise-ID
                (lo2, hi2) = allan.confidence_interval_noiseID(phase, dev, af=int(row['m']), dev_type="ohdev", data_type="phase")
                print(" tau= %f  lo/s32_lo = %f hi/s32_hi = %f "% (row['tau'], lo2/row['dev_min'], hi2/row['dev_max']))
                assert np.isclose(lo2, row['dev_min'], rtol=1e-2, atol=0)
                assert np.isclose(hi2, row['dev_max'], rtol=1e-2, atol=0)
            except NotImplementedError:
                print("can't do CI for tau= %f"%row['tau'])
                pass

    def generic_test(self, datafile=data_file, result="", fct=None):
        change_to_test_dir()
        testutils.test_row_by_row(fct, datafile, 1.0, result, verbose=verbose, tolerance=tolerance)

if __name__ == "__main__":
    t = TestGPS()
    #t.test_ohdev_ci_and_noiseID()
    t.test_oadev_ci_and_noiseID()
    #pytest.main()


