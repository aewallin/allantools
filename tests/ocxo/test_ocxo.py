"""
  Test for allantools (https://github.com/aewallin/allantools)
  Stable32 was used to calculate the deviations we compare against.

  AW2015-06-26
"""

import math
import sys
import pytest

sys.path.append("..")
sys.path.append("../..") # hack to import from parent directory
# remove if you have allantools installed in your python path

import allantools as allan
import testutils

import os
import time

def print_elapsed(start):
    end = time.clock()
    print(" %.2f s"% ( end-start ))
    return time.clock()

def change_to_test_dir():
    # hack to run script from its own directory
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)

data_file = 'ocxo_frequency.txt'

verbose = 1
tolerance = 1e-4 # relative tolerance
rate = 1/float(1.0) # stable32 runs were done with this data-interval
    
class TestOCXO():
    def test_ocxo_adev(self):
        self.generic_test( result= 'stable32_adev_alltau.txt' , fct= allan.adev )
    
    
    def test_adev_ci(self):
        s32rows = testutils.read_stable32(resultfile='adev_octave.txt', datarate=1.0)
        for row in s32rows:
            data = testutils.read_datafile(data_file)
            data = allan.frequency2fractional(data, mean_frequency=1.0e7)
            (taus, devs, errs, ns) = allan.adev(data, rate=rate, data_type="freq",
                                                  taus=[ row['tau'] ])
            edf = allan.edf_greenhall(alpha=row['alpha'],d=2,m=row['m'],N=len(data),overlapping=False, modified = False, verbose=True)
            print("alpha=",row['alpha'])
            (lo,hi) =allan.confidence_intervals(devs[0],ci=0.68268949213708585, edf=edf)
            print("n check: ", testutils.check_equal( ns[0], row['n'] ) )
            print("dev check: ", testutils.check_approx_equal( devs[0], row['dev'], tolerance=2e-4 ) )
            print("min dev check: ",  lo, row['dev_min'], testutils.check_approx_equal( lo, row['dev_min'], tolerance=1e-3 ) )
            print("max dev check: ", hi, row['dev_max'], testutils.check_approx_equal( hi, row['dev_max'], tolerance=1e-3 ) )
    
    def test_ocxo_oadev(self):
        self.generic_test( result= 'stable32_oadev_alltau.txt' , fct= allan.oadev )    
    
    def test_oadev_ci(self):
        s32rows = testutils.read_stable32(resultfile='oadev_octave.txt', datarate=1.0)
        for row in s32rows:
            data = testutils.read_datafile(data_file)
            data = allan.frequency2fractional(data, mean_frequency=1.0e7)
            (taus, devs, errs, ns) = allan.oadev(data, rate=rate, data_type="freq",
                                                  taus=[ row['tau'] ])
            edf = allan.edf_greenhall(alpha=row['alpha'],d=2,m=row['m'],N=len(data),overlapping=True, modified = False, verbose=True)
            (lo,hi) =allan.confidence_intervals(devs[0],ci=0.68268949213708585, edf=edf)
            print("n check: ", testutils.check_equal( ns[0], row['n'] ) )
            print("dev check: ", devs[0], row['dev'], testutils.check_approx_equal( devs[0], row['dev'], tolerance=2e-3 ) )
            print("min dev check: ",  lo, row['dev_min'], testutils.check_approx_equal( lo, row['dev_min'], tolerance=2e-3 ) )
            print("max dev check: ", hi, row['dev_max'], testutils.check_approx_equal( hi, row['dev_max'], tolerance=2e-3 ) )
            
    def test_ocxo_mdev(self):
        self.generic_test( result= 'stable32_mdev_alltau.txt' , fct= allan.mdev )   
    
    def test_mdev_ci(self):
        s32rows = testutils.read_stable32(resultfile='mdev_octave.txt', datarate=1.0)
        for row in s32rows:
            data = testutils.read_datafile(data_file)
            data = allan.frequency2fractional(data, mean_frequency=1.0e7)
            (taus, devs, errs, ns) = allan.mdev(data, rate=rate, data_type="freq",
                                                  taus=[ row['tau'] ])
            edf = allan.edf_greenhall(alpha=row['alpha'],d=2,m=row['m'],N=len(data),overlapping=True, modified = True, verbose=True)
            (lo,hi) =allan.confidence_intervals(devs[0],ci=0.68268949213708585, edf=edf)
            print("n check: ", testutils.check_equal( ns[0], row['n'] ) )
            print("dev check: ", devs[0], row['dev'], testutils.check_approx_equal( devs[0], row['dev'], tolerance=2e-3 ) )
            print("min dev check: ",  lo, row['dev_min'], testutils.check_approx_equal( lo, row['dev_min'], tolerance=2e-3 ) )
            print("max dev check: ", hi, row['dev_max'], testutils.check_approx_equal( hi, row['dev_max'], tolerance=2e-3 ) )
             
    def test_ocxo_tdev(self):
        self.generic_test( result= 'stable32_tdev_alltau.txt' , fct= allan.tdev )  
              
    def test_ocxo_hdev(self):
        self.generic_test( result= 'stable32_hdev_alltau.txt' , fct= allan.hdev )  
    
    def test_hdev_ci(self):
        s32rows = testutils.read_stable32(resultfile='hdev_octave.txt', datarate=1.0)
        for row in s32rows:
            data = testutils.read_datafile(data_file)
            data = allan.frequency2fractional(data, mean_frequency=1.0e7)
            (taus, devs, errs, ns) = allan.hdev(data, rate=rate, data_type="freq",
                                                  taus=[ row['tau'] ])
            edf = allan.edf_greenhall(alpha=row['alpha'],d=3,m=row['m'],N=len(data),overlapping=False, modified = False, verbose=True)
            (lo,hi) =allan.confidence_intervals(devs[0],ci=0.68268949213708585, edf=edf)
            print("n check: ", testutils.check_equal( ns[0], row['n'] ) )
            print("dev check: ", devs[0], row['dev'], testutils.check_approx_equal( devs[0], row['dev'], tolerance=2e-3 ) )
            print("min dev check: ",  lo, row['dev_min'], testutils.check_approx_equal( lo, row['dev_min'], tolerance=2e-3 ) )
            print("max dev check: ", hi, row['dev_max'], testutils.check_approx_equal( hi, row['dev_max'], tolerance=2e-3 ) )
    
    def test_ocxo_ohdev(self):
        self.generic_test( result= 'stable32_ohdev_alltau.txt' , fct= allan.ohdev )
        
    def test_ohdev_ci(self):
        s32rows = testutils.read_stable32(resultfile='ohdev_octave.txt', datarate=1.0)
        for row in s32rows:
            data = testutils.read_datafile(data_file)
            data = allan.frequency2fractional(data, mean_frequency=1.0e7)
            (taus, devs, errs, ns) = allan.ohdev(data, rate=rate, data_type="freq",
                                                  taus=[ row['tau'] ])
            edf = allan.edf_greenhall(alpha=row['alpha'],d=3,m=row['m'],N=len(data),overlapping=True, modified = False, verbose=True)
            (lo,hi) =allan.confidence_intervals(devs[0],ci=0.68268949213708585, edf=edf)
            print("n check: ", testutils.check_equal( ns[0], row['n'] ) )
            print("dev check: ", devs[0], row['dev'], testutils.check_approx_equal( devs[0], row['dev'], tolerance=2e-3 ) )
            print("min dev check: ",  lo, row['dev_min'], testutils.check_approx_equal( lo, row['dev_min'], tolerance=2e-3 ) )
            print("max dev check: ", hi, row['dev_max'], testutils.check_approx_equal( hi, row['dev_max'], tolerance=5e-3 ) )
            
    def test_ocxo_totdev(self):
        self.generic_test( result= 'stable32_totdev_alltau.txt' , fct= allan.totdev )        
    
    """
    def test_totdev_ci(self):
        s32rows = testutils.read_stable32(resultfile='totdev_octave.txt', datarate=1.0)
        for row in s32rows:
            data = testutils.read_datafile(data_file)
            data = allan.frequency2fractional(data, mean_frequency=1.0e7)
            (taus, devs, errs, ns) = allan.totdev(data, rate=rate, data_type="freq",
                                                  taus=[ row['tau'] ])
            edf = allan.edf_totdev(N=len(data),m=row['m'], alpha=row['alpha'])

            (lo,hi) = allan.confidence_intervals(devs[0],ci=0.68268949213708585, edf=edf)
            print("n check: ", testutils.check_equal( ns[0], row['n'] ) )
            print("dev check: ", testutils.check_approx_equal( devs[0], row['dev'] ) )
            print("min dev check: ",  lo, row['dev_min'], testutils.check_approx_equal( lo, row['dev_min'], tolerance=1e-3 ) )
            print("max dev check: ", hi, row['dev_max'], testutils.check_approx_equal( hi, row['dev_max'], tolerance=1e-3 ) )
    """
            
    def generic_test(self, datafile = data_file, result="", fct=None):
        change_to_test_dir()
        testutils.test_row_by_row( fct, datafile, 1.0, result , verbose=verbose, tolerance=tolerance, frequency=True, normalize=True)

if __name__ == "__main__":
    #pytest.main()
    t =TestOCXO()
    #t.test_ocxo_adev()
    t.test_adev_ci()
    t.test_oadev_ci()
    t.test_mdev_ci()
    t.test_hdev_ci()
    t.test_ohdev_ci()
    #t.test_totdev_ci()

