"""
  Test for allantools (https://github.com/aewallin/allantools)
  Stable32 was used to calculate the deviations we compare against.

  AW2015-03-29
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

def change_to_test_dir():
    # hack to run script from its own directory
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)
    
data_file = 'tic_phase.txt' # input data file
verbose = 1
tolerance = 1e-4
rate = 1/float(1.0) # stable32 runs were done with this data-interval
    
class TestTIC():
    def test_adev(self):
        self.generic_test( result= 'tic_adev.txt' , fct= allan.adev )
    def test_oadev(self):
        self.generic_test( result='tic_oadev.txt' , fct= allan.oadev )
    def test_mdev(self):
        self.generic_test( result='tic_mdev.txt' , fct= allan.mdev )
    def test_tdev(self):
        self.generic_test( result='tic_tdev.txt' , fct= allan.tdev )
    def test_hdev(self):
        self.generic_test( result='tic_hdev.txt' , fct= allan.hdev )
    def test_ohdev(self):
        self.generic_test( result='tic_ohdev.txt' , fct= allan.ohdev )
    def test_totdev(self):
        self.generic_test( result='tic_totdev.txt' , fct= allan.totdev )


    #def test_mtie(self):
    #    self.generic_test( result='mtie_fast.txt' , fct= allan.mtie )
    def test_tierms(self):
        self.generic_test( result='tic_tierms.txt' , fct= allan.tierms )
    
    def generic_test(self, datafile = data_file, result="", fct=None):
        change_to_test_dir()
        testutils.test_row_by_row( fct, datafile, 1.0, result , verbose=verbose, tolerance=tolerance)

    def test_noise_id(self):
        """ test for noise-identification """
        s32_rows = testutils.read_stable32( 'tic_oadev.txt' , 1.0 )
        phase = testutils.read_datafile(data_file)
        for s32 in s32_rows:
            tau, alpha, af = s32['tau'], s32['alpha'], int(s32['m'])
            try:
                alpha_int = allan.autocorr_noise_id( phase , af=af)[0]
                #if len(phase)/af > 30: # noise-id only works for length 30 or longer time-series
                assert alpha_int == alpha
                print( tau, alpha, alpha_int )
            except NotImplementedError:
                print ("no noise-ID: ",tau, alpha, alpha_int )
            
if __name__ == "__main__":
    #pytest.main()
    t=TestTIC()
    t.test_noise_id()
