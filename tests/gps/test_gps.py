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

import allantools as allan
import testutils

import os
import time
import pytest

def print_elapsed(start):
    end = time.clock()
    print(" %.2f s"% ( end-start ))
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
        self.generic_test( result='stable32_ADEV_decade.txt' , fct= allan.adev )
    def test_oadev(self):
        self.generic_test( result='stable32_OADEV_octave.txt' , fct= allan.oadev )
    def test_mdev(self):
        self.generic_test( result='stable32_MDEV_octave.txt' , fct= allan.mdev )
    def test_tdev(self):
        self.generic_test( result='stable32_TDEV_octave.txt' , fct= allan.tdev )
    def test_hdev(self):
        self.generic_test( result='stable32_HDEV_octave.txt' , fct= allan.hdev )
    def test_ohdev(self):
        self.generic_test( result='stable32_OHDEV_octave.txt' , fct= allan.ohdev )
    def test_totdev(self):
        self.generic_test( result='stable32_TOTDEV_octave.txt' , fct= allan.totdev )
    
    def generic_test(self, datafile = data_file, result="", fct=None):
        change_to_test_dir()
        testutils.test_row_by_row( fct, datafile, 1.0, result , verbose=verbose, tolerance=tolerance)

if __name__ == "__main__":
    pytest.main()


