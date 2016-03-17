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

def print_elapsed(start):
    end = time.clock()
    print(" %.2f s"% ( end-start ))
    return time.clock()
    
def run():
    # hack to run script from its own directory
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)

    data_file = 'gps_1pps_phase_data.txt.gz'
    
    adev_result = 'stable32_ADEV_decade.txt'
    oadev_result = 'stable32_OADEV_octave.txt'
    mdev_result = 'stable32_MDEV_octave.txt'
    tdev_result = 'stable32_TDEV_octave.txt'
    hdev_result = 'stable32_HDEV_octave.txt'
    ohdev_result = 'stable32_OHDEV_octave.txt'
    totdev_result = 'stable32_TOTDEV_octave.txt'
    
    verbose = 1
    tolerance = 1e-4 # relative tolerance
    rate = 1/float(1.0) # 1 PPS measurements, data interval is 1 s

    start0 = time.clock()
    start = time.clock()

    
    testutils.test_row_by_row( allan.adev_phase, data_file, rate, adev_result , verbose, tolerance)
    start = print_elapsed(start)
    
    testutils.test_row_by_row( allan.oadev_phase, data_file, rate, oadev_result, verbose, tolerance )
    start = print_elapsed(start)
    
    testutils.test_row_by_row( allan.mdev_phase, data_file, rate, mdev_result, verbose, tolerance)
    start = print_elapsed(start)
    
    testutils.test_row_by_row( allan.tdev_phase, data_file, rate, tdev_result, verbose, tolerance )
    start = print_elapsed(start)
    
    testutils.test_row_by_row( allan.hdev_phase, data_file, rate, hdev_result, verbose, tolerance )
    start = print_elapsed(start)
    
    testutils.test_row_by_row( allan.ohdev_phase, data_file, rate, ohdev_result, verbose, tolerance )
    start = print_elapsed(start)
    
    testutils.test_row_by_row( allan.totdev_phase, data_file, rate, totdev_result, verbose, tolerance )
    start = print_elapsed(start)
    
    print(" GPS tests took %.2f s" % ( time.clock()-start0 )) 
    
if __name__ == "__main__":
    run()


