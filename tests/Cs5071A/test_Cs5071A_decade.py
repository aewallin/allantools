"""
  Test for allantools (https://github.com/aewallin/allantools)
  Stable32 was used to calculate the deviations we compare against.

  The 5071A_phase.txt is a dataset collected with a time-interval-counter
  between 1 pulse-per-second outputs from a 5071A Cs clock against a H-maser
 
  first datapoint          1391174210 2014-01-31 13:16:50 UTC 
  last datapoint           1391731199 2014-02-06 23:59:59 UTC
  556990 datapoints in total

  This test uses log-spaced tau-values (i.e. 1, 2, 4, 10, etc.)
  The CS5071A_test_all.py test much more tau-values (1,2,3,4, etc.) but is slower.

  AW2014-02-07
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
    
data_file = '5071A_phase.txt.gz' # input data file
verbose = 1
tolerance = 1e-4
rate = 1/float(1.0) # stable32 runs were done with this data-interval

@pytest.mark.slow
class TestCS():
    def test_adev(self):
        self.generic_test( result= 'adev_decade.txt' , fct= allan.adev )
    def test_oadev(self):
        self.generic_test( result='oadev_decade.txt' , fct= allan.oadev )
    def test_mdev(self):
        self.generic_test( result='mdev_decade.txt' , fct= allan.mdev )
    def test_tdev(self):
        self.generic_test( result='tdev_decade.txt' , fct= allan.tdev )
    def test_hdev(self):
        self.generic_test( result='hdev_decade.txt' , fct= allan.hdev )
    def test_ohdev(self):
        self.generic_test( result='ohdev_decade.txt' , fct= allan.ohdev )
    def test_totdev(self):
        self.generic_test( result='totdev_decade.txt' , fct= allan.totdev )


    #def test_mtie(self):
    #    self.generic_test( result='mtie_fast.txt' , fct= allan.mtie )
    def test_tierms(self):
        self.generic_test( result='tierms_decade.txt' , fct= allan.tierms )
    
    def generic_test(self, datafile = data_file, result="", fct=None):
        change_to_test_dir()
        testutils.test_row_by_row( fct, datafile, 1.0, result , verbose=verbose, tolerance=tolerance)

if __name__ == "__main__":
    pytest.main()


