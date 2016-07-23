"""
  Pink frequency noise test for allantools (https://github.com/aewallin/allantools)
  Stable32 was used to calculate the deviations we compare against.

  the pink_frequency.txt was generated with noise.py, for documentation see that file.

"""
import math
import sys
import os
import pytest

sys.path.append("..")
sys.path.append("../..") # hack to import from parent directory
# remove if you have allantools installed in your python path

import allantools as allan
import testutils

def change_to_test_dir():
    # hack to run script from its own directory
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)

data_file = 'pink_frequency.txt'
verbose = 1
tolerance = 1e-4 # relative tolerance
rate = 1/float(42.0) # stable32 runs were done with this data-interval
    
class TestPink():
    def test_adev(self):
        self.generic_test( result= 'adev.txt' , fct= allan.adev )
    def test_oadev(self):
        self.generic_test( result= 'oadev.txt' , fct= allan.oadev )        
    def test_mdev(self):
        self.generic_test( result= 'mdev.txt' , fct= allan.mdev )        
    def test_tdev(self):
        self.generic_test( result= 'tdev.txt' , fct= allan.tdev )        
    def test_hdev(self):
        self.generic_test( result= 'hdev.txt' , fct= allan.hdev )        
    def test_ohdev(self):
        self.generic_test( result= 'ohdev.txt' , fct= allan.ohdev )        
    def test_totdev(self):
        self.generic_test( result= 'totdev_alpha0.txt' , fct= allan.totdev )        

    def generic_test(self, datafile = data_file, result="", fct=None):
        change_to_test_dir()
        testutils.test_row_by_row( fct, datafile, rate, result , verbose=verbose, tolerance=tolerance, frequency=True, normalize=False)
        

if __name__ == "__main__":
    #pytest.main()
    t = TestPink()
    t.test_adev()
