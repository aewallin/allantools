"""
  Test for allantools (https://github.com/aewallin/allantools)
  Stable32 was used to calculate the deviations we compare against.

  AW2015-06-26
  
  The dataset is from the 10 MHz output at the back of an HP Impedance Analyzer
  measured with Keysight 53230A counter, 1.0s gate, RCON mode, with H-maser 10MHz reference

"""

import math
import sys
import pytest

sys.path.append("..")
sys.path.append("../..") # hack to import from parent directory
# remove if you have allantools installed in your python path

import allantools as allan
import testutils

data_file = 'ocxo_frequency.txt'
import os
def change_to_test_dir():
    # hack to run script from its own directory
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)
    
verbose = 1
tolerance = 1e-4 # relative tolerance
rate = 1/float(1.0) # stable32 runs were done with this data-interval
    
class TestOCXO():

    def test_adev_ci(self):
        """ ADEV with confidence intervals """
        change_to_test_dir()
        s32rows = testutils.read_stable32(resultfile='adev_octave.txt', datarate=1.0)
        for row in s32rows:
            data = testutils.read_datafile(data_file)
            data = allan.frequency2fractional(data, mean_frequency=1.0e7)
            (taus, devs, errs, ns) = allan.adev(data, rate=rate, data_type="freq",
                                                  taus=[ row['tau'] ])
            # NOTE! Here we use alhpa from Stable32-results for the allantools edf computation!
            edf = allan.edf_greenhall(alpha=row['alpha'],d=2,m=row['m'],N=len(data),overlapping=False, modified = False, verbose=True)
            print("alpha=",row['alpha'])
            (lo, hi) =allan.confidence_interval(devs[0], edf=edf)
            print("      n check: %d" % testutils.check_equal( ns[0], row['n'] ) )
            print("    dev check: %d" % testutils.check_approx_equal( devs[0], row['dev'], tolerance=2e-4 ) )
            print("min dev check: %.4g %.4g %d" % ( lo, row['dev_min'], testutils.check_approx_equal( lo, row['dev_min'], tolerance=1e-3 ) ) )
            print("max dev check: %.4g %.4g %d" % ( hi, row['dev_max'], testutils.check_approx_equal( hi, row['dev_max'], tolerance=1e-3 ) ) )
    
    def test_oadev_ci(self):
        """ Overlapping ADEV with confidence intervals """
        s32rows = testutils.read_stable32(resultfile='oadev_octave.txt', datarate=1.0)
        for row in s32rows:
            data = testutils.read_datafile(data_file)
            data = allan.frequency2fractional(data, mean_frequency=1.0e7)
            (taus, devs, errs, ns) = allan.oadev(data, rate=rate, data_type="freq",
                                                  taus=[ row['tau'] ])
            # NOTE! Here we use alhpa from Stable32-results for the allantools edf computation!
            edf = allan.edf_greenhall(alpha=row['alpha'],d=2,m=row['m'],N=len(data),overlapping=True, modified = False, verbose=True)
            (lo, hi) =allan.confidence_interval(devs[0], edf=edf)
            print("n check: ", testutils.check_equal( ns[0], row['n'] ) )
            print("dev check: ", devs[0], row['dev'], testutils.check_approx_equal( devs[0], row['dev'], tolerance=2e-3 ) )
            print("min dev check: ",  lo, row['dev_min'], testutils.check_approx_equal( lo, row['dev_min'], tolerance=2e-3 ) )
            print("max dev check: ", hi, row['dev_max'], testutils.check_approx_equal( hi, row['dev_max'], tolerance=2e-3 ) )
    
    def test_mdev_ci(self):
        """ Overlapping ADEV with confidence intervals """
        s32rows = testutils.read_stable32(resultfile='mdev_octave.txt', datarate=1.0)
        for row in s32rows:
            data = testutils.read_datafile(data_file)
            data = allan.frequency2fractional(data, mean_frequency=1.0e7)
            (taus, devs, errs, ns) = allan.mdev(data, rate=rate, data_type="freq",
                                                  taus=[ row['tau'] ])
            # NOTE! Here we use alhpa from Stable32-results for the allantools edf computation!
            edf = allan.edf_greenhall(alpha=row['alpha'],d=2,m=row['m'],N=len(data),overlapping=True, modified = True, verbose=True)
            (lo,hi) =allan.confidence_interval(devs[0], edf=edf)
            print("n check: ", testutils.check_equal( ns[0], row['n'] ) )
            print("dev check: ", devs[0], row['dev'], testutils.check_approx_equal( devs[0], row['dev'], tolerance=2e-3 ) )
            print("min dev check: ",  lo, row['dev_min'], testutils.check_approx_equal( lo, row['dev_min'], tolerance=2e-3 ) )
            print("max dev check: ", hi, row['dev_max'], testutils.check_approx_equal( hi, row['dev_max'], tolerance=2e-3 ) )

    def test_tdev_ci(self):
        """ Time Deviation with confidence intervals """
        s32rows = testutils.read_stable32(resultfile='tdev_octave.txt', datarate=1.0)
        for row in s32rows:
            data = testutils.read_datafile(data_file)
            data = allan.frequency2fractional(data, mean_frequency=1.0e7)
            (taus, devs, errs, ns) = allan.tdev(data, rate=rate, data_type="freq",
                                                  taus=[ row['tau'] ])
            # NOTE! Here we use alhpa from Stable32-results for the allantools edf computation!
            edf = allan.edf_greenhall(alpha=row['alpha'],d=2,m=row['m'],N=len(data),overlapping=True, modified = True, verbose=True)
            (lo,hi) =allan.confidence_interval(devs[0], edf=edf)
            print("n check: ", testutils.check_equal( ns[0], row['n'] ) )
            print("dev check: ", devs[0], row['dev'], testutils.check_approx_equal( devs[0], row['dev'], tolerance=2e-3 ) )
            print("min dev check: ",  lo, row['dev_min'], testutils.check_approx_equal( lo, row['dev_min'], tolerance=2e-3 ) )
            print("max dev check: ", hi, row['dev_max'], testutils.check_approx_equal( hi, row['dev_max'], tolerance=2e-3 ) )

    def test_hdev_ci(self):
        """ Hadamard with confidence intervals """
        s32rows = testutils.read_stable32(resultfile='hdev_octave.txt', datarate=1.0)
        for row in s32rows:
            data = testutils.read_datafile(data_file)
            data = allan.frequency2fractional(data, mean_frequency=1.0e7)
            (taus, devs, errs, ns) = allan.hdev(data, rate=rate, data_type="freq",
                                                  taus=[ row['tau'] ])
            # NOTE! Here we use alhpa from Stable32-results for the allantools edf computation!
            edf = allan.edf_greenhall(alpha=row['alpha'],d=3,m=row['m'],N=len(data),overlapping=False, modified = False, verbose=True)
            (lo,hi) =allan.confidence_interval(devs[0], edf=edf)
            print("n check: ", testutils.check_equal( ns[0], row['n'] ) )
            print("dev check: ", devs[0], row['dev'], testutils.check_approx_equal( devs[0], row['dev'], tolerance=2e-3 ) )
            print("min dev check: ",  lo, row['dev_min'], testutils.check_approx_equal( lo, row['dev_min'], tolerance=2e-3 ) )
            print("max dev check: ", hi, row['dev_max'], testutils.check_approx_equal( hi, row['dev_max'], tolerance=2e-3 ) )
    
    def test_ohdev_ci(self):
        """ Overlapping Hadamard deviation with confidence intervals  """
        s32rows = testutils.read_stable32(resultfile='ohdev_octave.txt', datarate=1.0)
        for row in s32rows:
            data = testutils.read_datafile(data_file)
            data = allan.frequency2fractional(data, mean_frequency=1.0e7)
            (taus, devs, errs, ns) = allan.ohdev(data, rate=rate, data_type="freq",
                                                  taus=[ row['tau'] ])
            # NOTE! Here we use alhpa from Stable32-results for the allantools edf computation!
            edf = allan.edf_greenhall(alpha=row['alpha'],d=3,m=row['m'],N=len(data),overlapping=True, modified = False, verbose=True)
            (lo,hi) =allan.confidence_interval(devs[0], edf=edf)
            print("n check: ", testutils.check_equal( ns[0], row['n'] ) )
            print("dev check: ", devs[0], row['dev'], testutils.check_approx_equal( devs[0], row['dev'], tolerance=2e-3 ) )
            print("min dev check: ",  lo, row['dev_min'], testutils.check_approx_equal( lo, row['dev_min'], tolerance=2e-3 ) )
            print("max dev check: ", hi, row['dev_max'], testutils.check_approx_equal( hi, row['dev_max'], tolerance=5e-3 ) )

    # fails
    # totdev() needs bias-correction, depending on alpha(?)
    @pytest.mark.skip(reason="needs bias-correction and noise-ID to work")
    @pytest.mark.xfail
    def test_totdev_ci(self):
        print("totdev()")
        s32rows = testutils.read_stable32(resultfile='totdev_octave.txt', datarate=1.0)
        for row in s32rows:
            data = testutils.read_datafile(data_file)
            data = allan.frequency2fractional(data, mean_frequency=1.0e7)
            (taus, devs, errs, ns) = allan.totdev(data, rate=rate, data_type="freq",
                                                  taus=[ row['tau'] ])
            edf = allan.edf_totdev(N=len(data),m=row['m'], alpha=row['alpha'])

            (lo,hi) = allan.confidence_interval(devs[0], edf=edf)
            print("n check: ", testutils.check_equal( ns[0], row['n'] ) )
            print("dev check: ", testutils.check_approx_equal( devs[0], row['dev'], tolerance=2e-3 ) )
            print("min dev check: %.4g %.4g %d" %( lo, row['dev_min'], testutils.check_approx_equal( lo, row['dev_min'], tolerance=2e-3 ) ) )
            print("max dev check: %.4g %.4g %d" %( hi, row['dev_max'], testutils.check_approx_equal( hi, row['dev_max'], tolerance=2e-3 ) ) )

if __name__ == "__main__":
    #pytest.main()
    t =TestOCXO()
    t.test_ocxo_adev()
    t.test_adev_ci()
    t.test_oadev_ci()
    t.test_mdev_ci()
    t.test_tdev_ci()
    t.test_hdev_ci()
    t.test_ohdev_ci()

    #t.test_totdev_ci()

