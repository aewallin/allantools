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
    def test_ocxo_oadev(self):
        self.generic_test( result= 'stable32_oadev_alltau.txt' , fct= allan.oadev )        
    def test_ocxo_mdev(self):
        self.generic_test( result= 'stable32_mdev_alltau.txt' , fct= allan.mdev )        
    def test_ocxo_tdev(self):
        self.generic_test( result= 'stable32_tdev_alltau.txt' , fct= allan.tdev )        
    def test_ocxo_hdev(self):
        self.generic_test( result= 'stable32_hdev_alltau.txt' , fct= allan.hdev )        
    def test_ocxo_ohdev(self):
        self.generic_test( result= 'stable32_ohdev_alltau.txt' , fct= allan.ohdev )        
    def test_ocxo_totdev(self):
        self.generic_test( result= 'stable32_totdev_alltau.txt' , fct= allan.totdev )        


    def generic_test(self, datafile = data_file, result="", fct=None):
        change_to_test_dir()
        testutils.test_row_by_row( fct, datafile, 1.0, result , verbose=verbose, tolerance=tolerance, frequency=True, normalize=True)

if __name__ == "__main__":
    pytest.main()


